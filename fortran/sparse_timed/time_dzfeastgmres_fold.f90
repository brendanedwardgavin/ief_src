subroutine time_dzfeastgmres(UPLO,n,A1,fpm,emin1,emax1,m0)

implicit none

integer :: i,j,k
double precision :: dtmp1,dtmp2
character(len=1) UPLO

!!!!!!!!!!!!!!!!!!!!!!!! Dense matrix:
complex (kind=kind(0.0d0)),dimension(n,*) :: A1

!!!!!!!!!!!!!!!!!!!!!!!!  FEAST:

integer :: n,m0,ijob,loop,info,m
integer, dimension(64)::fpm
double precision :: epsout,emin1,emax1
double precision, dimension(:), allocatable :: res,e
complex (kind=kind(0.0d0)) :: ze
complex (kind=kind(0.0d0)), dimension(:,:), allocatable :: zwork,work,aq,bq,x
complex(kind=(kind(1.0d0))),dimension(:),allocatable :: Zne,Wne
!!!!!!!!!!!!!!!!!!!!!!!!  Linear system solver:

integer :: linjob,linIterations,kdim
! linjob: RCI task returned by linear solver
! linIterations: number of linear system iterations; each iteration requires kdim matvec operations
! kdim: number of Krylov subspace blocks for linear solver
integer, dimension(3) :: linState !internal parameters for linear solver
complex (kind=kind(0.0d0)), dimension(:,:), allocatable :: xsol,V,Av,Ax,zwork2,linwork1,linwork2
! xsol,V,Av,Ax,zwork2: internal storage for linear system solver
! linwork1,linwork2: input and output storage for linear solver
complex (kind=kind(0.0d0)) :: ze2 !FEAST contour point for linear system sovler
!double precision, dimension(:),allocatable :: vaavn
double precision :: lineps,linepsmult ! linear system accuracy
integer :: maxeig !max eigenvector to count in measuring lin sys residual

!!!!!!!!!!!!!!!!!!!!!!!!! Lapack:

double precision, dimension(:), allocatable :: leig
complex (kind=kind(0.0d0)), dimension(:), allocatable :: zhwork
double precision, dimension(:), allocatable :: rwork
integer :: lwork


!!!!!!!!!!!!!!!!!!!!!!! BLAS:
double precision, external :: dznrm2,elapsed_time
character, dimension(6) :: matdescra

!!!!!!!!!!!!!!!!!!!!!!! Timing

integer :: totalc1,totalc2,c1,c2,linc1,linc2
double precision :: times_gm
double precision :: time_ls,time_mm,time_linsys,time_total
double precision, dimension(3) :: times_breakdown_feast
    !1:matrix multiply
    !2:factorize linear system
    !3:solve linear system

!!!!!!!!!!!!!!!!!!!!!!! convergence vs time
double precision, dimension(:),allocatable :: reslist,timelist,linsysreslist
double precision, dimension(:,:),allocatable :: rhsreslistavg
double precision, dimension(:,:,:),allocatable :: rhsreslistcp
double precision, dimension(:),allocatable :: linsysitavg
double precision, dimension(:,:),allocatable :: linsysitcp
integer :: oldloop,cpcount,linitout

!!!!!!!!!!!!!!!!!!!!!!! debug/lin sys accuracy
double precision :: maxres,tempres
double precision, dimension(:),allocatable :: tempreslist
complex (kind=kind(0.0d0)), dimension(:,:), allocatable :: zt1,zt2,zt3 
integer :: maxm !max lin sys index for error calculation
complex (kind=kind(0.0d0)), dimension(:), allocatable :: cplist

!!!!!!!!!!!!!!!!!!!!!!! print flags
integer :: meas_acc,print_times

!!!!!!!!!!!!!!!!!!!!!!! string stuff
character(5) :: m0str,cpstr

!!!!!!!!!!!!!!!!!!!!!!! spectrum folding
double precision :: emin,emax,emid
complex (kind=kind(0.0d0)), dimension(:,:), allocatable :: spectmp1,spectmp2,A,eye

allocate(spectmp1(n,n),spectmp2(n,n),A(n,n),eye(n,n))
emid=(emin1+emax1)/2.0d0
emin=0.0d0
emax=(emax1-emid)**2

do i=1,n
    eye(i,i)=(1.0d0,0.0d0)
end do

spectmp1=A1(1:n,1:n)-emid*eye
call zgemm('N','N',n,n,n,(1.0d0,0.0d0),spectmp1(1,1),n,spectmp1,n,(0.0d0,0.0d0),A,n)

A=A1(1:n,1:n)
emin=emin1
emax=emax1

meas_acc=1
print_times=0

times_gm=0.0d0
time_ls=0.0d0
time_mm=0.0d0
time_linsys=0.0d0
times_breakdown_feast=0.0d0



allocate(zt1(n,m0),zt2(n,m0),zt3(n,m0))

allocate(res(m0))
allocate(e(m0))
allocate(work(n,m0))
allocate(aq(m0,m0))
allocate(bq(m0,m0))
allocate(x(n,m0))
allocate(zwork(n,m0))

allocate(reslist(0:fpm(4)),timelist(0:fpm(4)),linsysreslist(0:fpm(4)),rhsreslistavg(0:fpm(4),1:m0),rhsreslistcp(1:fpm(2),0:fpm(4),1:m0))
allocate(tempreslist(1:m0))
allocate(linsysitavg(0:fpm(4)),linsysitcp(1:fpm(2),0:fpm(4)))
allocate(cplist(1:fpm(2)))
linsysreslist=0.0d0
rhsreslistavg=0.0d0
rhsreslistcp=0.0d0
linsysitavg=0.0d0
linsysitcp=0.0d0
!!!!!!!!!!!!!  Set up linear system solver:
matdescra(1)='H'
matdescra(2)='L'
matdescra(3)='N'
matdescra(4)='F'
linIterations=fpm(50) !number of iterations to use; more linear system iterations allows FEAST to converge faster
kdim=fpm(51) !number of Krylov subspace blocks to use; 1 usually works fine, but more gives faster convergence
linepsmult=10.0d0**(-1.0d0*dble(fpm(53))) !1d-2 !go this far below eigenvalue residual for linear systems
maxeig=4
allocate(xsol(n,m0), V(n,m0*kdim), Av(n,m0*kdim), Ax(n,m0), zwork2(n,m0*kdim), linwork1(n,m0), linwork2(n,m0))

allocate(Zne(fpm(2)),Wne(fpm(2)))
call zfeast_contour(emin,emax,fpm(2),fpm(16),fpm(18),Zne,Wne)
if(fpm(18)==0) then
    do i=1,fpm(2)
       Zne(i)=dble(Zne(i))*(1.0d0,0.0d0)+(0.0d0,1.0d-9)
    end do
end if

if (fpm(11)==0) then


!call zfeast_scsrev_timed(UPLO,n,dsa,isa,jsa,fpm,epsout,loop,emin,emax,m0,e,x,m,res,info,reslist,timelist)
!if (UPLO .ne. 'L') then

!call zfeast_hcsrev(UPLO,n,dsa,isa,jsa,fpm,epsout,loop,emin,emax,m0,e,x,m,res,info)

!call zfeast_hcsrev_timed(UPLO,n,dsa,isa,jsa,fpm,epsout,loop,emin,emax,m0,e,x,m,res,info,reslist,timelist,times_breakdown_feast)

time_total=timelist(loop)
time_mm=times_breakdown_feast(1)
time_linsys=times_breakdown_feast(2)+times_breakdown_feast(3)
!call dfeast_scsrev_timed2(UPLO,n,dsa,isa,jsa,fpm,epsout,loop,emin,emax,m0,e,x,m,res,info,reslist,timelist)
!else
!call dfeast_scsrev_timed2('U',n,ssa,sisa,sjsa,fpm,epsout,loop,emin,emax,m0,e,x,m,res,info,reslist,timelist)
!end if

else


!!!!!!!!!!!!  Run FEAST RCI:
call system_clock(count=totalc1)
ijob=-1
oldloop=0
cpcount=1
do while (ijob .ne. 0)
call zfeast_hrcix(ijob,n,ze,work,zwork,aq,bq,fpm,epsout,loop,emin,emax,m0,e,x,m,res,info,Zne,Wne)

if (oldloop .ne. loop) then
    reslist(oldloop)=epsout
    call system_clock(count=totalc2)
    timelist(oldloop)=elapsed_time(totalc1,totalc2)
    oldloop=loop
    cpcount=1
end if

    select case (ijob)
		case (10) !store complex shift ze
		
            ze2=ze	
		
        case (11) !solve linear system Az*Qz=zwork for Qz, put Qz in zwork 
            

            if(meas_acc) zt1=zwork

            linState(1)=-1 !set linsate(1)=-1 to begin RCI routine dfeast_gmres
            
            
            if(epsout>0.0d0) then
                lineps=epsout*linepsmult
            else
                lineps=linepsmult
            end if
            
            if (m>0) then 
                maxm=min(m,m0)
            else
                maxm=m0
            end if           
            maxeig=maxm

            call system_clock(count=linc1)
            do while (linState(1) .ne. 0)
                     
                call zfeast_gmres(linjob,linState,zwork,xsol,V,Av,Ax,ze2,n,m0,maxeig,lineps,linIterations,kdim,linwork1,linwork2,zwork2,times_gm)
                linitout=linState(2) 

                if(linjob == 30) then
                    
                    !Your matrix multiplication routine goes here! do linwork2=A*linwork1
                    call system_clock(count=c1)
                    call zgemm('N','N',n,m0,n,(1.0d0,0.0d0),A,n,linwork1,n,(0.0d0,0.0d0),linwork2,n)
                    !call mkl_zcsrmm('N', n, m0, n, 1.0d0, matdescra, dsa, jsa, isa, isa(2), linwork1, n, 0.0d0, linwork2, n)
                    !linwork2=matmul(A,linwork1)
                    call system_clock(count=c2)
                    time_mm=time_mm+elapsed_time(c1,c2)
                    if(print_times) print *,'mm gmres 11',elapsed_time(c1,c2)
                end if
            end do
            call system_clock(count=linc2)
            time_ls=time_ls+times_gm
            time_linsys=time_linsys+elapsed_time(linc1,linc2)
            
            if(print_times) print *,'ls gmres 11',times_gm
            
            if(meas_acc>0) then !measure lin sys res!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
            
            !call mkl_zcsrmm('N', n, m0, n, 1.0d0, matdescra, dsa, jsa, isa, isa(2), zwork, n, 0.0d0, zt2, n)
            call zgemm('N','N',n,m0,n,(1.0d0,0.0d0),A,n,zwork,n,(0.0d0,0.0d0),zt2,n)

            zt2=ze2*zwork-zt2
            zt3=zt2-zt1
            maxres=0.0d0
            do i=1,m0
               tempreslist(i)=dznrm2(n,zt3(1:n,i),1)/dznrm2(n,zt1(1:n,i),1)
               rhsreslistavg(loop,i)=rhsreslistavg(loop,i)+tempreslist(i)/dble(fpm(2))
               rhsreslistcp(cpcount,loop,i)=tempreslist(i)
            end do
            
            call quicksort(tempreslist,1,m0) !sort residuals from smallest to biggest
            maxres=tempreslist(maxm)
 
            print *,''
            print *,"lin sys error 1 =",maxres,linState(2)
            !print *,"    tol= ",lineps
            !print *,"    loops=",linitout
            linsysreslist(loop)=linsysreslist(loop)+maxres/dble(fpm(2))
            
            if(maxres>lineps .and. linitout<linIterations) then
                print *,'Error: linear system solver did not converge'
                !stop
            end if

            end if   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            linsysitavg(loop)=linsysitavg(loop)+dble(linitout)/dble(fpm(2))
            linsysitcp(cpcount,loop)=dble(linitout)
            cplist(cpcount)=ze2 
            cpcount=cpcount+1
            
        case(21)

            if(epsout>0.0d0) then
                lineps=epsout*linepsmult
            else
                lineps=linepsmult
            end if

            linState(1)=-1 !set linsate(1)=-1 to begin RCI routine dfeast_gmres
            call system_clock(count=linc1)
            do while (linState(1) .ne. 0)
                if(meas_acc) zt1=zwork
                !call system_clock(count=c1)
                call zfeast_gmres(linjob,linState,zwork,xsol,V,Av,Ax,conjg(ze2),n,m0,maxeig,lineps,linIterations,kdim,linwork1,linwork2,zwork2,times_gm)
                !call system_clock(count=c2)
                !time_gmres=time_gmres+elapsed_time(c1,c2)
                
                if(linjob == 30) then
                    !Your matrix multiplication routine goes here! do linwork2=A*linwork1
                    call system_clock(count=c1)
                    call zgemm('N','N',n,m0,n,(1.0d0,0.0d0),A,n,linwork1,n,(0.0d0,0.0d0),linwork2,n)
                    !call mkl_zcsrmm('N', n, m0, n, 1.0d0, matdescra, dsa, jsa, isa, isa(2), linwork1, n, 0.0d0, linwork2, n)
                    !linwork2=matmul(A,linwork1)
                    call system_clock(count=c2)
                    time_mm=time_mm+elapsed_time(c1,c2)
                    if(print_times) print *,'mm gmres 21',elapsed_time(c1,c2)
                end if
            end do
            call system_clock(count=linc2)
            time_linsys=time_linsys+elapsed_time(linc1,linc2)
            time_ls=time_ls+times_gm
            
            if(print_times) print *,'ls gmres 11',times_gm
            
            if(meas_acc>0) then
            if(m>0) then
                maxm=min(m,m0)
            else
                maxm=m0
            end if

            call zgemm('N','N',n,m0,n,(1.0d0,0.0d0),A,n,zwork,n,(0.0d0,0.0d0),zt2,n)
            !call mkl_zcsrmm('N', n, m0, n, 1.0d0, matdescra, dsa, jsa, isa, isa(2), zwork, n, 0.0d0, zt2, n)
            zt2=conjg(ze2)*zwork-zt2
            zt3=zt2-zt1
            maxres=0.0d0
            do i=1,maxm
               tempres=dznrm2(n,zt3(1:n,i),1)/dznrm2(n,zt1(1:n,i),1)
               if(tempres>maxres) maxres=tempres
            end do
            print *,"lin sys error 2 =",maxres,linState(2)
            print *,''
            end if

        case (30) !A*x	
            
            !Your matrix multiplication routine goes here! Do work=A*x
            call system_clock(count=c1)
            call zgemm('N','N',n,m0,n,(1.0d0,0.0d0),A,n,x,n,(0.0d0,0.0d0),work,n)
            !call mkl_zcsrmm('N', n, m0, n, 1.0d0, matdescra, dsa, jsa, isa, isa(2), x, n, 0.0d0, work, n)
            call system_clock(count=c2)
            time_mm=time_mm+elapsed_time(c1,c2)
            if(print_times) print *,'mm case 30',elapsed_time(c1,c2)
            
		case (40) !B*x
        
            call zlacpy('F',n,m0,x,n,work,n)

	end select
end do

call system_clock(count=totalc2)
time_total=elapsed_time(totalc1,totalc2)

reslist(loop)=epsout
timelist(loop)=elapsed_time(totalc1,totalc2)

end if



!print results:
print *,'FEAST finished; # eigs found = ',m
print *, 'eigenvalues and eigenvector residuals:'
do i=1,m
    print *,i,e(i),res(i)
end do



print *,''
print *,'Total time=',time_total

print *,''
print *,'Matrix multiply time: '!,time_mm,time_mm/time_total
print *,'      Time:',time_mm
print *,'      Percent:',100*time_mm/time_total
print *,''
print *,'Linear system time: '!,time_linsys,time_linsys/time_total
print *,'      Time:',time_linsys
print *,'      Percent:',100*time_linsys/time_total
print *,''
print *,'Least squares time: '
print *,'      Time:',time_ls
print *,'      Percent:',100*time_ls/time_total
print *,''
print *,'Factorize time: '!,time_linsys,time_linsys/time_total
print *,'      Time:',times_breakdown_feast(2)
print *,'      Percent:',100*times_breakdown_feast(2)/time_total
print *,''
print *,'Solve time: '!,time_linsys,time_linsys/time_total
print *,'      Time:',times_breakdown_feast(3)
print *,'      Percent:',100*times_breakdown_feast(3)/time_total
!write residual iterations and times to output file
open(unit=10,file='../output/residualsout.dat',status='REPLACE')
do i=0,loop
    write (10,"(I3, 3ES15.5, F8.2)") i,timelist(i),reslist(i),linsysreslist(i),linsysitavg(i)
    !write(10,*) i,timelist(i),reslist(i),linsysreslist(i)
end do
close(10)

open(unit=10,file='../output/final_vals.dat',status='REPLACE')
    write (10,*) loop
    write (10,*) timelist(loop)
    write (10,*) reslist(loop)
    write (10,*) linsysreslist(loop)
    write (10,*) linsysitavg(loop)
close(10)

write(m0str,"(I5)") m0
open(unit=10,file='../output/rhsresidualsout.dat',status='REPLACE')
open(unit=11,file='../output/linitsout.dat',status='REPLACE')
do i=0,loop
    write (10,"(I3)",advance="no") i
    write (10,"("//m0str//"ES15.5)") (rhsreslistavg(i,j), j=1,m0 )

    write (11,"(I3)",advance="no") i
    write (11,"("//m0str//"F8.2)") (linsysitcp(j,i), j=1,fpm(2) )
    !write(10,*) i,timelist(i),reslist(i),linsysreslist(i)
end do
close(10)
close(11)

open(unit=10,file='../output/contourpoints.dat',status='REPLACE')
do i=1,fpm(2)
    write(10,*) i,cplist(i)
end do
close(10)

deallocate(xsol, V, Av, Ax, zwork2, linwork1, linwork2)
deallocate(res,e,work,aq,bq,x,zwork)

end subroutine time_dzfeastgmres


!double precision function elapsed_time(c1,c2)
!integer :: c1,c2
!integer :: diff
!integer :: maxcount,countrate

!call system_clock(count_rate=countrate,count_max=maxcount)

!if(c2<c1) then
!    diff=maxcount+c2-c1
!else
!    diff=c2-c1
!end if

!elapsed_time= dble(diff)/dble(countrate)

!end function elapsed_time


