subroutine time_szfeastgmres(UPLO,n,dsa,isa,jsa,fpm,emin,emax,m0)

implicit none

integer :: i,j,k
double precision :: dtmp1,dtmp2
character(len=1):: UPLO

!!!!!!!!!!!!!!!!!!!!!!!!  Sparse matrix:
complex (kind=kind(0.0d0)),dimension(*) :: dsa
integer,dimension(*) :: isa,jsa
integer :: nnza

!!!!!!!!!!!!!!!!!!!!!!!!  FEAST:

integer :: n,m0,ijob,loop,info,m
integer, dimension(64)::fpm
double precision :: epsout,emin,emax
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
integer :: dofold
double precision :: eminf,emaxf,emid

dofold=1

emid=(emin+emax)/2.0d0
if(dofold>0) then
    eminf=0.0d0
    emaxf=(emax-emid)**2
else
    eminf=emin
    emaxf=emax
end if

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
if(UPLO=='U' .or. UPLO=='u') then
    matdescra(2)='U'
else
    matdescra(2)='L'
end if
matdescra(3)='N'
matdescra(4)='F'
linIterations=fpm(50) !number of iterations to use; more linear system iterations allows FEAST to converge faster
kdim=fpm(51) !number of Krylov subspace blocks to use; 1 usually works fine, but more gives faster convergence
linepsmult=10.0d0**(-1.0d0*dble(fpm(53))) !1d-2 !go this far below eigenvalue residual for linear systems
maxeig=4
allocate(xsol(n,m0), V(n,m0*kdim), Av(n,m0*kdim), Ax(n,m0), zwork2(n,m0*kdim), linwork1(n,m0), linwork2(n,m0))

allocate(Zne(fpm(2)),Wne(fpm(2)))

call zfeast_contour(eminf,emaxf,fpm(2),fpm(16),fpm(18),Zne,Wne)

if(fpm(18)==0) then
    do i=1,fpm(2)
       Zne(i)=dble(Zne(i))*(1.0d0,0.0d0)+(0.0d0,1.0d-9)
    end do
end if

if (fpm(11)==0) then


!call zfeast_scsrev_timed(UPLO,n,dsa,isa,jsa,fpm,epsout,loop,emin,emax,m0,e,x,m,res,info,reslist,timelist)
!if (UPLO .ne. 'L') then

!call zfeast_hcsrev(UPLO,n,dsa,isa,jsa,fpm,epsout,loop,emin,emax,m0,e,x,m,res,info)

call zfeast_hcsrev_timed(UPLO,n,dsa,isa,jsa,fpm,epsout,loop,emin,emax,m0,e,x,m,res,info,reslist,timelist,times_breakdown_feast)

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

call zfeast_hrcix(ijob,n,ze,work,zwork,aq,bq,fpm,epsout,loop,eminf,emaxf,m0,e,x,m,res,info,Zne,Wne)

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
                    !call zgemm('N','N',n,m0,n,(1.0d0,0.0d0),A,n,linwork1,n,(0.0d0,0.0d0),linwork2,n)
                    
                    !call mkl_zcsrmm('N', n, m0, n, 1.0d0, matdescra, dsa, jsa, isa, isa(2), linwork1, n, 0.0d0, linwork2, n)
                    call specfold_mult(UPLO,n,dsa,isa,jsa,linwork1,linwork2,m0,emid,dofold)
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
            call specfold_mult(UPLO,n,dsa,isa,jsa,zwork,zt2,m0,emid,dofold)
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
                    !call zgemm('N','N',n,m0,n,(1.0d0,0.0d0),A,n,linwork1,n,(0.0d0,0.0d0),linwork2,n)
                    
                    !call mkl_zcsrmm('N', n, m0, n, 1.0d0, matdescra, dsa, jsa, isa, isa(2), linwork1, n, 0.0d0, linwork2, n)
                    
                    call specfold_mult(UPLO,n,dsa,isa,jsa,linwork1,linwork2,m0,emid,dofold)
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
            !call mkl_zcsrmm('N', n, m0, n, 1.0d0, matdescra, dsa, jsa, isa, isa(2), zwork, n, 0.0d0, zt2, n)
            call specfold_mult(UPLO,n,dsa,isa,jsa,zwork,zt2,m0,emid,dofold)
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
            !call mkl_zcsrmm('N', n, m0, n, 1.0d0, matdescra, dsa, jsa, isa, isa(2), x, n, 0.0d0, work, n)
            call specfold_mult(UPLO,n,dsa,isa,jsa,x,work,m0,emid,dofold)
            call system_clock(count=c2)
            time_mm=time_mm+elapsed_time(c1,c2)
            if(print_times) print *,'mm case 30',elapsed_time(c1,c2)
            !call zgemm('N','N',n,m0,n,(1.0d0,0.0d0),A,n,x,n,(0.0d0,0.0d0),work,n)
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

end subroutine time_szfeastgmres




subroutine specfold_mult(UPLO,n,dsa,isa,jsa,x,b,m0,emid,dofold)
implicit none

integer :: dofold
character(len=1) :: UPLO
integer :: n,m0
!!!!!!!!!!!!!!!!!!!!!!!!  Sparse matrix:
complex (kind=kind(0.0d0)),dimension(*) :: dsa
integer,dimension(*) :: isa,jsa

complex (kind=kind(0.0d0)),dimension(n,*) :: x,b
double precision :: emid

character, dimension(6) :: matdescra

complex (kind=kind(0.0d0)),dimension(:,:),allocatable :: temp1,temp2

matdescra(1)='H'
if(UPLO=='U' .or. UPLO=='u') then
    matdescra(2)='U'
else
    matdescra(2)='L'
end if
matdescra(3)='N'
matdescra(4)='F'


if(dofold>0) then
    !print *,'fold'
    allocate(temp1(n,m0),temp2(n,m0))

    !(A-emidI)*(A-emidI)*x=A**2x-2emidAx+emid**2x

    !temp=A*x
    call mkl_zcsrmm('N', n, m0, n, 1.0d0, matdescra, dsa, jsa, isa, isa(2), x(1,1), n, 0.0d0, temp1, n)

    !temp2=A*temp
    call mkl_zcsrmm('N', n, m0, n, 1.0d0, matdescra, dsa, jsa, isa, isa(2), temp1, n, 0.0d0, temp2, n)

    !temp2=temp2-2*emid*A*x
    call mkl_zcsrmm('N', n, m0, n, (-2.0d0,0.0d0)*emid, matdescra, dsa, jsa, isa, isa(2), x(1,1), n, (1.0d0,0.0d0), temp2, n)

    !temp2=temp2-emid**2x
    b(1:n,1:m0)=temp2+emid**2*x(1:n,1:m0)

else
    !print *,'nofold'
    call mkl_zcsrmm('N', n, m0, n, 1.0d0, matdescra, dsa, jsa, isa, isa(2), x(1,1), n, 0.0d0, b(1,1), n)    
end if
end subroutine specfold_mult






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


  subroutine zfeast_hcsrev_timed(UPLO,N,sa,isa,jsa,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,reslist,timelist,times_breakdown_feast)
    !  Purpose 
    !  =======
    !  FEAST SPARSE INTERFACE
    !  Solve the standard Ax=ex eigenvalue problem
    !  
    !  A COMPLEX HERMITIAN :: SPARSE FORMAT 
    ! 
    !  DOUBLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
    !                           triangular part of the matrix(ces) is(are) being supplied.
    !  N          (input)        INTEGER: Size system
    !  sa         (input)        COMPLEX DOUBLE PRECISION (isa(N+1)-1):  Matrix A- CSR format 
    !  isa        (input)        INTEGER(N+1): CSR row array of Matrix A
    !  jsa        (input)        INTEGER(isa(N+1)-1): CSR column array of Matrix A
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  fpm (input/output) INTEGER(*) : FEAST parameters
    !  epsout     (output)       REAL DOUBLE PRECISION : Error on the trace
    !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
    !  Emin,Emax  (input)        REAL DOUBLE PRECISION: search interval
    !  M0         (input/output) INTEGER: Size subspace
    !  E          (output)       REAL DOUBLE PRECISION(M0)   : Eigenvalues -solution
    !  q          (input/output) COMPLEX DOUBLE PRECISION(N,M0) : 
    !                                                       On entry: subspace initial guess if fpm(5)=1 
    !                                                       On exit : Eigenvectors-solution
    !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
    !  res        (output)       REAL DOUBLE PRECISION(M0) : Relative Residual of the solution (1-norm)
    !                                                        if option fpm(6)=1 selected                           
    !  info       (output)       INTEGER: Error handling (0: successful exit)
    !=====================================================================
    ! Eric Polizzi 2009-2015
    ! ====================================================================
    implicit none
    include 'f90_noruntime_interface.fi'
    character(len=1) :: UPLO
    integer :: N
    complex(kind=(kind(1.0d0))),dimension(*) :: sa
    integer,dimension(*) :: isa,jsa
    integer,dimension(*) :: fpm
    double precision :: epsout 
    integer :: loop
    double precision :: Emin,Emax
    integer :: M0
    double precision,dimension(*)  :: E
    complex(kind=(kind(1.0d0))),dimension(N,*):: X
    integer :: mode
    double precision,dimension(*)    :: res
    integer :: info
!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!    Wrapper Routine to expert routine: zfeast_hcsrgvx
!!!    Create Zne, Wne arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    complex(kind=(kind(1.0d0))),dimension(fpm(2)) :: Zne,Wne 
    complex(kind=(kind(1.0d0))),dimension(n) :: sb ! identity
    integer,dimension(n+1) :: isb
    integer,dimension(n) :: jsb
    integer :: i

    !!!!!!!!!!!!!!!!!!!!!!!!timing 
    double precision, dimension(0:fpm(4)) :: reslist,timelist
    double precision, dimension(*) :: times_breakdown_feast

    call zfeast_contour(Emin,Emax,fpm(2),fpm(16),fpm(18),Zne,Wne)

    fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
    ! identity B matrix- option for standard eigenvalue problem
    do i=1,n
       sb(i)=(1.0d0,0.0d0)
       jsb(i)=i
       isb(i)=i
    enddo
    isb(n+1)=n+1

    call zfeast_hcsrgvx_timed(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne,reslist,timelist,times_breakdown_feast)

  end subroutine zfeast_hcsrev_timed





  subroutine zfeast_hcsrgvx_timed(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne,reslist,timelist,times_breakdown_feast)
    !  Purpose 
    !  =======
    !  FEAST SPARSE INTERFACE
    !  Solve the generalized Ax=eBx eigenvalue problem
    !  
    !  A COMPLEX HERMITIAN, B HERMITIAN POSITIVE DEFINITE:: SPARSE FORMAT 
    ! 
    !  DOUBLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
    !                           triangular part of the matrix(ces) is(are) being supplied.
    !  N          (input)        INTEGER: Size system
    !  sa         (input)        COMPLEX DOUBLE PRECISION (isa(N+1)-1):  Matrix A- CSR format 
    !  isa        (input)        INTEGER(N+1): CSR row array of Matrix A
    !  jsa        (input)        INTEGER(isa(N+1)-1): CSR column array of Matrix A
    !  sb         (input)        COMPLEX DOUBLE PRECISION (isb(N+1)-1):  Matrix B- CSR format 
    !  isb        (input)        INTEGER(N+1): CSR row array of Matrix B
    !  jsb        (input)        INTEGER(isb(N+1)-1): CSR column array of Matrix B
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  fpm (input/output) INTEGER(*) : FEAST parameters
    !  epsout     (output)       REAL DOUBLE PRECISION : Error on the trace
    !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
    !  Emin,Emax  (input)        REAL DOUBLE PRECISION: search interval
    !  M0         (input/output) INTEGER: Size subspace
    !  E          (output)       REAL DOUBLE PRECISION(M0)   : Eigenvalues -solution
    !  q          (input/output) COMPLEX DOUBLE PRECISION(N,M0) : 
    !                                                       On entry: subspace initial guess if fpm(5)=1 
    !                                                       On exit : Eigenvectors-solution
    !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
    !  res        (output)       REAL DOUBLE PRECISION(M0) : Relative Residual of the solution (1-norm)
    !                                                        if option fpm(6)=1 selected                           
    !  info       (output)       INTEGER: Error handling (0: successful exit)
    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of ARGUMENTS FOR EXPERT ROUTINE 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  Zne        (input)        COMPLEX DOUBLE PRECISION(fpm(2)): Custom Integration nodes
    !  Wne        (input)        COMPLEX DOUBLE PRECISION(fpm(2)): Custom Integration weights
    !           
    !=====================================================================
    ! Eric Polizzi 2009-2015
    !=====================================================================
    implicit none
    include 'f90_noruntime_interface.fi'
    character(len=1) :: UPLO
    integer :: N
    complex(kind=(kind(1.0d0))),dimension(*) :: sa,sb
    integer,dimension(*) :: isa,jsa,isb,jsb
    integer,dimension(*) :: fpm
    double precision :: epsout 
    integer :: loop
    double precision :: Emin,Emax
    integer :: M0
    double precision,dimension(*)  :: E
    complex(kind=(kind(1.0d0))),dimension(N,*):: X
    integer :: mode
    double precision,dimension(*)    :: res
    integer :: info
    complex(kind=(kind(1.0d0))),dimension(*) :: Zne,Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer :: ijob,infoloc,i,s
    complex(kind=(kind(1.0d0))) :: Ze
    complex(kind=(kind(1.0d0))), dimension(:,:),pointer ::work,workc,zAq,zSq,caux
    complex(kind=(kind(1.0d0))),dimension(:,:),pointer :: saz
    integer,dimension(:),pointer :: isaz,jsaz
    double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
    complex(kind=(kind(1.0d0))),parameter :: ONEC=(DONE,DZERO),ZEROC=(DZERO,DZERO)
    logical :: fact
    integer :: nfact,id
    integer :: rank,code,nb_procs,NEW_COMM_WORLD
!!!!! full csr format
    complex(kind=kind(1.0d0)),dimension(:),pointer :: ssa,ssb
    integer,dimension(:), pointer :: sisa,sjsa,sisb,sjsb
    integer :: opt,nnza,nnzb,nnz
!!!!!for pardiso
    integer(8),dimension(64) :: pt
    integer,dimension(64) :: iparm
    integer :: mtype
    integer :: MAXFCT,MNUM,PHASE,MSGLVL
    integer :: idum

!!!!!!!!!!!!!! for timing
    double precision, dimension(0:fpm(4)) :: reslist,timelist
    double precision, dimension(*) :: times_breakdown_feast
    integer :: totalc1,totalc2,c1,c2
    integer :: oldloop
    double precision, external :: elapsed_time
    character, dimension(6) :: matdescra

    call system_clock(count=totalc1)
    matdescra(1)='S'
    matdescra(2)='U'
    matdescra(3)='N'
    matdescra(4)='F'



    rank=0
    nb_procs=1
    !----------------------------------------------
#ifdef MPI
    NEW_COMM_WORLD=fpm(9)
    call MPI_COMM_RANK(NEW_COMM_WORLD,rank,code)
    call MPI_COMM_SIZE(NEW_COMM_WORLD,nb_procs,code)
#endif
    !----------------------

    INFO = 0
    IF ((UPLO/='F').and.(UPLO/='f').and.(UPLO/='L').and.(UPLO/='l').and.(UPLO/='U').and.(UPLO/='u')) THEN
       INFO=-101
    ELSE IF ( N<=0 ) THEN
       INFO = -102
    END IF
    IF( INFO.NE.0 ) THEN
       CALL XERBLA( 'ZFEAST_HCSRGV', -INFO+100 )
       RETURN
    END IF


    infoloc=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! FORMAT CONVERSION TO full CSR for PARDISO !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Remark: everything should be csr since we work with zS-H which is actually unsymmetric 


    if ((UPLO=='F').or.(UPLO=='f')) then !!! full csr already

       !       nnza=isa(n+1)-1
       !       ssa => sa(1:nnza)
       !       sisa => isa(1:n+1)
       !       sjsa => jsa(1:nnza)

       !       nnzb=isb(n+1)-1
       !       ssb => sb(1:nnzb)
       !       sisb => isb(1:n+1)
       !       sjsb => jsb(1:nnzb)

    else !! upper-csr or lower-csr to full csr

       nnza=2*(isa(n+1)-1) ! may be overestimated
       nnzb=2*(isb(n+1)-1)
       call wallocate_1z(ssa,nnza,infoloc)
       call wallocate_1i(sjsa,nnza,infoloc)
       call wallocate_1z(ssb,nnzb,infoloc)
       call wallocate_1i(sjsb,nnzb,infoloc)
       call wallocate_1i(sisa,n+1,infoloc)
       call wallocate_1i(sisb,n+1,infoloc)
       if (infoloc/=0) then
          info=-1
          return
       end if

       call zhcsr_uplo_to_csr(N,sa,isa,jsa,ssa,sisa,sjsa)
       call zhcsr_uplo_to_csr(N,sb,isb,jsb,ssb,sisb,sjsb)

    endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



    call wallocate_2z(zAq,M0,M0,infoloc)
    call wallocate_2z(zSq,M0,M0,infoloc)
    call wallocate_2z(work,N,M0,infoloc)
    call wallocate_2z(workc,N,M0,infoloc)
    if (infoloc/=0) then
       info=-1
       return
    end if


!!!!!!!!!!!!!!!!!! Factorizations Set-up  
    fact=.true.
    nfact=0
    !! nfact is local (number of total factorization by contour points)
    do i=rank+1,fpm(2),nb_procs
       nfact=nfact+1
    end do

    if (fpm(10)==1) then
       id=0
    else
       id=1
    end if

!!!!!!!!!!!!!!!!! Set up for Az matrix

    call wallocate_1i(isaz,n+1,infoloc)
    call wallocate_2z(saz,1,1,infoloc) ! dummy
    call wallocate_1i(jsaz,1,infoloc)! dummy
    !!>>>
    opt=1
    if ((UPLO=='F').or.(UPLO=='f')) then
       call zaddcsr(N,opt,ONEC,sa,isa,jsa,ONEC,sb,isb,jsb,saz(1,1),isaz,jsaz) !! get isaz
    else
       call zaddcsr(N,opt,ONEC,ssa,sisa,sjsa,ONEC,ssb,sisb,sjsb,saz(1,1),isaz,jsaz) !! get isaz
    end if
    nnz=isaz(n+1)-1
    !!<<<
    call wdeallocate_2z(saz)
    call wdeallocate_1i(jsaz)
    !!>>>
    call wallocate_2z(saz,nnz,nfact,infoloc)
    call wallocate_1i(jsaz,nnz,infoloc)

    opt=2
    if ((UPLO=='F').or.(UPLO=='f')) then
       call zaddcsr(N,opt,ONEC,sa,isa,jsa,ONEC,sb,isb,jsb,saz(1,1),isaz,jsaz) !! get jsaz
    else
       call zaddcsr(N,opt,ONEC,ssa,sisa,sjsa,ONEC,ssb,sisb,sjsb,saz(1,1),isaz,jsaz) !! get jsaz
    endif
!!!!!!!!!!!!!!!
    call wallocate_2z(caux,N,M0,infoloc)

    if (infoloc/=0) then
       info=-1
       return
    end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  pardiso initialization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    MAXFCT=nfact ! Rq: same factorization for (normal+transpose)
    MTYPE=3      ! complex and structurally symmetric 
    call pardisoinit(PT,MTYPE,IPARM)
!!!!!!!!!!!!
    if (fpm(64)==1) then
       do i=1,64
          if (fpm(64+i)/=-9999) iparm(i)=fpm(64+i)
       enddo
    endif
!!!!!!!!!!!!
    IPARM(6)=1 ! solution and rhs are input/output, attention caux is always used
    MSGLVL=0!0 !0- no output, 1- output
    PHASE=11 ! symbolic factorization (do it only once)
    !    call PARDISO(PT,MAXFCT,1,MTYPE,PHASE,n,saz(1,1),isaz,jsaz,idum,M0,IPARM,MSGLVL,workc,caux,infoloc)
!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!
    oldloop=0
    ijob=-1 ! initialization
    do while (ijob/=0) 
       call zfeast_hrcix(ijob,N,Ze,work,workc,zAq,zSq,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)
       if (oldloop .ne. loop) then
            reslist(oldloop)=epsout
            call system_clock(count=totalc2)
            timelist(oldloop)=elapsed_time(totalc1,totalc2)
            print *,'time=',timelist(oldloop)
            oldloop=loop
        end if

       select case(ijob)
       case(10) !! factorize (zeB-A)
          call system_clock(count=c1)
          
          if (fpm(10)==1) then
             id=id+1
             if (id==nfact+1) then
                id=1
                fact=.false.
             endif
          endif

          if (fact) then
             opt=3
             if ((UPLO=='F').or.(UPLO=='f')) then
                call zaddcsr(N,opt,-ONEC,sa,isa,jsa,Ze,sb,isb,jsb,saz(1,id),isaz,jsaz) !! get saz
             else
                call zaddcsr(N,opt,-ONEC,ssa,sisa,sjsa,Ze,ssb,sisb,sjsb,saz(1,id),isaz,jsaz) !! get saz
             end if

             if (PHASE==11) then
                PHASE=12
             else
                PHASE=22
             endif
             !             PHASE=22 ! numerical fact only
             MNUM=id
             call PARDISO(PT,MAXFCT,MNUM,MTYPE,PHASE,n,saz(1,id),isaz,jsaz,idum,fpm(23),IPARM,MSGLVL,workc,caux,infoloc)
             if (infoloc/=0) then
                info=-2
                return
             end if

          end if ! fact true
            call system_clock(count=c2)
            print *,'   factorize time=',elapsed_time(c1,c2)
            times_breakdown_feast(2)=times_breakdown_feast(2)+elapsed_time(c1,c2)
       case(11) !!solve the linear system (ZeB-A)x=workc(1:N,1:fpm(23)) result in to workc
            call system_clock(count=c1)
          IPARM(12)=0 ! normal solve
          PHASE=33 ! solve
          MNUM=id
          call PARDISO(PT,MAXFCT,MNUM,MTYPE,PHASE,n,saz(1,id),isaz,jsaz,idum,fpm(23),IPARM,MSGLVL,workc,caux,infoloc)    

          if (infoloc/=0) then
             info=-2
             return
          end if
            call system_clock(count=c2)
            print *,'   solve time=',elapsed_time(c1,c2)
            times_breakdown_feast(3)=times_breakdown_feast(3)+elapsed_time(c1,c2)
       case(21) !!solve the linear system (ZeB-A)^H x=workc(1:N,1:fpm(23)) result in to workc
            call system_clock(count=c1)
          IPARM(12)=1 ! transpose conjugate solve
          PHASE=33 ! solve
          MNUM=id
          call PARDISO(PT,MAXFCT,MNUM,MTYPE,PHASE,n,saz(1,id),isaz,jsaz,idum,fpm(23),IPARM,MSGLVL,workc,caux,infoloc)    
          if (infoloc/=0) then
             info=-2
             return
          end if
            call system_clock(count=c2)
            print *,'   solve time=',elapsed_time(c1,c2)
            times_breakdown_feast(3)=times_breakdown_feast(3)+elapsed_time(c1,c2)
       case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)
          call system_clock(count=c1)
          if ((UPLO=='F').or.(UPLO=='f')) then
             !call zhcsrmm('F','N',N,N,fpm(25),ONEC,sa,isa,jsa,X(1,fpm(24)),ZEROC,work(1,fpm(24)))
             call mkl_zcsrmm('N', n, fpm(25), n, 1.0d0, matdescra, sa, jsa, isa, isa(2), X(1,fpm(24)), n, 0.0d0, work(1,fpm(24)), n)
          else
             !call zhcsrmm('F','N',N,N,fpm(25),ONEC,ssa,sisa,sjsa,X(1,fpm(24)),ZEROC,work(1,fpm(24)))
             call mkl_zcsrmm('N', n, fpm(25), n, 1.0d0, matdescra, ssa, sjsa, sisa, sisa(2), X(1,fpm(24)), n, 0.0d0, work(1,fpm(24)), n)
          end if
            call system_clock(count=c2)
            print *,'   mm time=',elapsed_time(c1,c2),fpm(24),m0
            times_breakdown_feast(1)=times_breakdown_feast(1)+elapsed_time(c1,c2)
       case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)
            call system_clock(count=c1)
          if ((UPLO=='F').or.(UPLO=='f')) then
             !call zhcsrmm('F','N',N,N,fpm(25),ONEC,sb,isb,jsb,X(1,fpm(24)),ZEROC,work(1,fpm(24)))
             call mkl_zcsrmm('N', n, fpm(25), n, 1.0d0, matdescra, sb, jsb, isb, isb(2), X(1,fpm(24)), n, 0.0d0, work(1,fpm(24)), n)
          else
             !call zhcsrmm('F','N',N,N,fpm(25),ONEC,ssb,sisb,sjsb,X(1,fpm(24)),ZEROC,work(1,fpm(24)))
             call mkl_zcsrmm('N', n, fpm(25), n, 1.0d0, matdescra, ssb, sjsb, sisb, sisb(2), X(1,fpm(24)), n, 0.0d0, work(1,fpm(24)), n)
          endif
            call system_clock(count=c2)
            print *,'   mm time=',elapsed_time(c1,c2),fpm(24),m0
            times_breakdown_feast(1)=times_breakdown_feast(1)+elapsed_time(c1,c2)
       end select
    end do
    reslist(loop)=epsout
    call system_clock(count=totalc2)
    timelist(loop)=elapsed_time(totalc1,totalc2)
!!!!!!!!!!!!!!!!!!!!!!!
!!!! release memory
!!!!!!!!!!!!!!!!!!!!!!!!
    PHASE=-1 

    do MNUM=1,nfact
       call PARDISO(PT,MAXFCT,MNUM,MTYPE,PHASE,n,saz(1,MNUM),isaz,jsaz,idum,fpm(23),IPARM,MSGLVL,workc,caux,infoloc)
       if (infoloc/=0) then
          info=-2
          return
       end if
    end do



    call wdeallocate_2z(zAq)
    call wdeallocate_2z(zSq)
    call wdeallocate_2z(work)
    call wdeallocate_2z(workc)
    call wdeallocate_2z(caux)


    call wdeallocate_2z(saz)
    call wdeallocate_1i(isaz)
    call wdeallocate_1i(jsaz)

    if ((UPLO/='F').and.(UPLO/='f')) then
       call wdeallocate_1z(ssa)
       call wdeallocate_1i(sisa)
       call wdeallocate_1i(sjsa)

       call wdeallocate_1z(ssb)
       call wdeallocate_1i(sisb)
       call wdeallocate_1i(sjsb)
    endif


  end subroutine zfeast_hcsrgvx_timed



