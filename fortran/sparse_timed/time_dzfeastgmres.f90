subroutine time_dzfeastgmres(UPLO,n,A,fpm,emin,emax,m0,outname)

implicit none

character(len=100) :: outname
integer :: i,j,k
double precision :: dtmp1,dtmp2
character(len=1) UPLO

!!!!!!!!!!!!!!!!!!!!!!!! Dense matrix:
complex (kind=kind(0.0d0)),dimension(n,*) :: A

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
double precision, dimension(:,:),allocatable :: linsysitcp,linsysrescp
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


meas_acc=0
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
allocate(linsysitavg(0:fpm(4)),linsysitcp(1:fpm(2),0:fpm(4)),linsysrescp(1:fpm(2),0:fpm(4)))
allocate(cplist(1:fpm(2)))
linsysreslist=0.0d0
rhsreslistavg=0.0d0
rhsreslistcp=0.0d0
linsysitavg=0.0d0
linsysitcp=0.0d0
linsysrescp=0.0d0
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

call zfeast_heev_timed(UPLO,n,A,n,fpm,epsout,loop,emin,emax,m0,e,x,m,res,info,reslist,timelist,times_breakdown_feast)

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
            linsysrescp(cpcount,loop)=maxres 
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
open(unit=10,file='../output/'//trim(outname)//'residualsout.dat',status='REPLACE')
do i=0,loop
    write (10,"(I3, 3ES15.5, F8.2)") i,timelist(i),reslist(i),linsysreslist(i),linsysitavg(i)
    !write(10,*) i,timelist(i),reslist(i),linsysreslist(i)
end do
close(10)

open(unit=10,file='../output/'//trim(outname)//'final_vals.dat',status='REPLACE')
    write (10,*) loop
    write (10,*) timelist(loop)
    write (10,*) reslist(loop)
    write (10,*) linsysreslist(loop)
    write (10,*) linsysitavg(loop)
close(10)

write(m0str,"(I5)") m0
write(cpstr,"(I5)") fpm(2)
open(unit=10,file='../output/'//trim(outname)//'rhsresidualsout.dat',status='REPLACE')
open(unit=11,file='../output/'//trim(outname)//'linitsout.dat',status='REPLACE')
open(unit=12,file='../output/'//trim(outname)//'linrescpout.dat',status='REPLACE')
do i=0,loop
    write (10,"(I3)",advance="no") i
    write (10,"("//m0str//"ES15.5)") (rhsreslistavg(i,j), j=1,m0 )

    write (11,"(I3)",advance="no") i
    write (11,"("//cpstr//"F8.2)") (linsysitcp(j,i), j=1,fpm(2) )
    write (12,"(I3)",advance="no") i
    write (12,"("//cpstr//"ES15.5)") (linsysrescp(j,i), j=1,fpm(2) )
    !write(10,*) i,timelist(i),reslist(i),linsysreslist(i)
end do
close(10)
close(11)
close(12)

open(unit=10,file='../output/'//trim(outname)//'contourpoints.dat',status='REPLACE')
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



subroutine zfeast_heev_timed(UPLO,N,A,LDA,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,reslist,timelist,times_breakdown_feast)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the standard Ax=ex eigenvalue problem
  !  
  !  A COMPLEX HERMITIAN :: DENSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
  !                           triangular part of the matrix(ces) is(are) being supplied.
  !  N          (input)        INTEGER: Size system
  !  A          (input)        COMPLEX DOUBLE PRECISION (LDA,N):  Matrix A 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=N)
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
  integer :: N,LDA
  complex(kind=(kind(1.0d0))),dimension(LDA,*):: A
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
!!!    Wrapper Routine to expert routine: zfeast_hegvx
!!!    Create Zne, Wne arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  complex(kind=(kind(1.0d0))),dimension(fpm(2)) :: Zne,Wne 
  complex(kind=(kind(1.0d0))),dimension(1) :: B ! dummy
  integer :: LDB

    !!!!!!!!!!!!!!!!!!!!!!!!timing 
    double precision, dimension(0:fpm(4)) :: reslist,timelist
    double precision, dimension(*) :: times_breakdown_feast

  call zfeast_contour(Emin,Emax,fpm(2),fpm(16),fpm(18),Zne,Wne)

  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  LDB=-1 ! B is a dummy- option for standard eigenvalue problem
  call zfeast_hegvx_timed(UPLO,N,A,LDA,B,LDB,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne,reslist,timelist,times_breakdown_feast)

end subroutine zfeast_heev_timed





subroutine zfeast_hegvx_timed(UPLO,N,A,LDA,B,LDB,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne,reslist,timelist,times_breakdown_feast)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the generalized Ax=eBx eigenvalue problem
  !  
  !  A COMPLEX HERMITIAN, B HERMITIAN POSITIVE DEFINITE:: DENSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
  !                           triangular part of the matrix(ces) is(are) being supplied.
  !  N          (input)        INTEGER: Size system
  !  A          (input)        COMPLEX DOUBLE PRECISION (LDA,N):  Matrix A 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=N)
  !  B          (input)        COMPLEX DOUBLE PRECISION (LDB,N):  Matrix B 
  !  LDB        (input)        INTEGER: Leading dimension of matrix B (LDB>=N)
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
  integer :: N,LDA,LDB
  complex(kind=(kind(1.0d0))),dimension(LDA,*):: A
  complex(kind=(kind(1.0d0))),dimension(LDB,*):: B
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
  complex(kind=(kind(1.0d0))), dimension(:,:),pointer ::work,workc,zAq,zSq
  complex(kind=(kind(1.0d0))), dimension(:,:,:),pointer :: Az,cAz
  complex(kind=(kind(1.0d0))), dimension(:),pointer ::ztmp
  integer, dimension(:,:),pointer ::ipivloc
  double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
  complex(kind=(kind(1.0d0))),parameter :: ONEC=(DONE,DZERO),ZEROC=(DZERO,DZERO)
  logical :: fact
  integer :: nfact,id
  integer :: rank,code,nb_procs,NEW_COMM_WORLD

!!!!!!!!!!!!!! for timing
    double precision, dimension(0:fpm(4)) :: reslist,timelist
    double precision, dimension(*) :: times_breakdown_feast
    integer :: totalc1,totalc2,c1,c2
    integer :: oldloop
    double precision, external :: elapsed_time


    call system_clock(count=totalc1)


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
  ELSE IF(LDA<N ) THEN
     INFO = -104
  ELSE IF((LDB<N ).and.(LDB/=-1)) THEN
     INFO = -106
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'ZFEAST_HEGV', -INFO+100 )
     RETURN
  END IF
  infoloc=0

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
!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call wallocate_3z(Az,N,N,nfact,infoloc)
  call wallocate_3z(cAz,N,N,nfact,infoloc)

  call wallocate_2i(ipivloc,N,nfact,infoloc)
  call wallocate_1z(ztmp,N,infoloc)
  if (infoloc/=0) then
     info=-1
     return
  end if
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
           if (LDB==-1) then !! standard
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Format CONVERSION

              if ((UPLO=='L').or.(UPLO=='l')) then
                 do i=1,N-1
                    s=N-i
                    Az(i:N,i,id)=-A(i:N,i)
                    Az(i,i,id)=Az(i,i,id)+Ze
                    call ZCOPY(s,Az(i+1,i,id),1,Az(i,i+1,id),N)
                    call ZLACGV(s, Az(i,i+1,id), N )
                 enddo
                 Az(N,N,id)=Ze-A(N,N)*ONEC

              elseif ((UPLO=='U').or.(UPLO=='u')) then
                 do i=1,N-1
                    s=N-i
                    Az(i,i:N,id)=-A(i,i:N)
                    Az(i,i,id)=Az(i,i,id)+Ze
                    call ZCOPY(s,Az(i,i+1,id),N,Az(i+1,i,id),1)
                    call ZLACGV(s, Az(i+1,i,id), 1 )
                 enddo
                 Az(N,N,id)=Ze-A(N,N)*ONEC

              else

                 Az(1:N,1:N,id)=-A(1:N,1:N)
                 do i=1,N
                    Az(i,i,id)=Az(i,i,id)+Ze
                 enddo

              end if


           else 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Format CONVERSION
              if ((UPLO=='L').or.(UPLO=='l')) then
                 do i=1,N-1
                    s=N-i
                    Az(i:N,i,id)=Ze*B(i:N,i)
                    call ZCOPY(s,Az(i+1,i,id),1,Az(i,i+1,id),N)
                    call ZLACGV(s, Az(i,i+1,id), N )
                    call ZSCAL(s,Ze/conjg(Ze),Az(i,i+1,id),N)
                    Az(i:N,i,id)=Az(i:N,i,id)-A(i:N,i)
                 enddo
                 Az(N,N,id)=Ze*B(N,N)-A(N,N)*ONEC
                 do i=1,N-1
                    s=N-i
                    call ZCOPY(s,A(i+1,i),1,ztmp(1),1)
                    call ZLACGV(s,ztmp(1),1)
                    call ZAXPY(s,-ONEC,ztmp(1),1,Az(i,i+1,id),N)
                 enddo


              elseif ((UPLO=='U').or.(UPLO=='u')) then

                 do i=1,N-1
                    s=N-i
                    Az(i,i:N,id)=Ze*B(i,i:N)
                    call ZCOPY(s,Az(i,i+1,id),N,Az(i+1,i,id),1)
                    call ZLACGV(s, Az(i+1,i,id), 1 )
                    call ZSCAL(s,Ze/conjg(Ze),Az(i+1,i,id),1)
                    Az(i,i:N,id)=Az(i,i:N,id)-A(i,i:N)
                 enddo
                 Az(N,N,id)=Ze*B(N,N)-A(N,N)*ONEC
                 do i=1,N-1
                    s=N-i
                    call ZCOPY(s,A(i,i+1),N,ztmp(1),1)
                    call ZLACGV(s,ztmp(1),1)
                    call ZAXPY(s,-ONEC,ztmp(1),1,Az(i+1,i,id),1)
                 enddo
              else ! full 
                 Az(1:N,1:N,id)=Ze*B(1:N,1:N)-A(1:N,1:N)
              end if
           end if

           if(fpm(11)==0) then
           call ZGETRF(N,N,Az(1,1,id),N,IPIVloc(1,id),INFOloc)     
           if (infoloc/=0) then
              info=-2
              return
           end if
           end if
        end if ! fact true
            call system_clock(count=c2)
            print *,'   factorize time=',elapsed_time(c1,c2)
            times_breakdown_feast(2)=times_breakdown_feast(2)+elapsed_time(c1,c2)
     case(11) !!solve the linear system (ZeB-A)x=workc(1:N,1:fpm(23)) result in to workc
            call system_clock(count=c1)

        call ZGETRS( 'N', N, fpm(23), Az(1,1,id), N, IPIVloc(1,id), workc, N, INFOloc )
        if (infoloc/=0) then
           info=-2
           return
        end if

            call system_clock(count=c2)
            print *,'   solve time=',elapsed_time(c1,c2)
            times_breakdown_feast(3)=times_breakdown_feast(3)+elapsed_time(c1,c2)

     case(20) 
        !cAz(1:N,1:N,id)=conjg(Ze)*B(1:N,1:N)-A(1:N,1:N)
        cAz(1:N,1:N,id)=transpose(conjg(Az(1:n,1:n,id)))

     case(21) !!solve the linear system (ZeB-A)^H x=workc(1:N,1:fpm(23)) result in to workc

         call system_clock(count=c1)

        !cAz(1:N,1:N,id)=transpose(conjg(Az(1:n,1:n,id)))

        call ZGETRS( 'C', N, fpm(23), Az(1,1,id), N, IPIVloc(1,id), workc, N, INFOloc )
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
           call ZGEMM('N','N',N,fpm(25),N,ONEC,A,LDA,X(1,fpm(24)),N,ZEROC,work(1,fpm(24)),N)
        else
           call ZHEMM ('L', UPLO, N, fpm(25), ONEC, A, LDA, X(1,fpm(24)), N, ZEROC,work(1,fpm(24)), N)
        endif
            call system_clock(count=c2)
            print *,'   mm time=',elapsed_time(c1,c2),fpm(24),m0
            times_breakdown_feast(1)=times_breakdown_feast(1)+elapsed_time(c1,c2)
     case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)
            call system_clock(count=c1)
        if (LDB==-1) then ! standard
           call ZLACPY( 'F', N, fpm(25),X(1,fpm(24)) , N, work(1,fpm(24)), N )
        else
           if ((UPLO=='F').or.(UPLO=='f')) then
              call ZGEMM('N','N',N,fpm(25),N,ONEC,B,LDB,X(1,fpm(24)),N,ZEROC,work(1,fpm(24)),N)
           else
              call ZHEMM ('L', UPLO, N, fpm(25), ONEC, B, LDB, X(1,fpm(24)), N, ZEROC,work(1,fpm(24)), N)
           endif
        end if
            call system_clock(count=c2)
            print *,'   mm time=',elapsed_time(c1,c2),fpm(24),m0
            times_breakdown_feast(1)=times_breakdown_feast(1)+elapsed_time(c1,c2)

     end select
  end do

    reslist(loop)=epsout
    call system_clock(count=totalc2)
    timelist(loop)=elapsed_time(totalc1,totalc2)

  call wdeallocate_2z(zAq)
  call wdeallocate_2z(zSq)
  call wdeallocate_2z(work)
  call wdeallocate_2z(workc)
  call wdeallocate_3z(Az)
  call wdeallocate_2i(ipivloc)
  call wdeallocate_1z(ztmp)


end subroutine zfeast_hegvx_timed




