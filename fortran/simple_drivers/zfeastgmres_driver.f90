program ziterativeFeast

implicit none

integer :: i,j,k
double precision :: dtmp1,dtmp2
!!!!!!!!!!!!!!!!!!!!!!!!  FEAST:

integer :: n,m0,ijob,loop,info,m
integer, dimension(64)::fpm
double precision :: epsout,emin,emax
double precision, dimension(:), allocatable :: res,e
complex (kind=kind(0.0d0)) :: ze
complex (kind=kind(0.0d0)), dimension(:,:), allocatable :: zwork,work,aq,bq,x

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

!!!!!!!!!!!!!!!!!!!!!!!!! Lapack:

double precision, dimension(:), allocatable :: leig
complex (kind=kind(0.0d0)), dimension(:), allocatable :: zhwork
double precision, dimension(:), allocatable :: rwork
integer :: lwork

!!!!!!!!!!!!!!!!!!!!!!! Artificial eigenvalue problem:
complex (kind=kind(0.0d0)), dimension(:,:), allocatable :: A, ztmpmat1,ztmpmat2,ztmpmat3

!!!!!!!!!!!!!!!!!!!!!!! BLAS:
double precision, external :: dznrm2,elapsed_time

!!!!!!!!!!!!!!!!!!!!!!! Timing

integer :: c1,c2
integer :: totalc1,totalc2
double precision :: time_matvec, time_feast, time_gmres
double precision :: times_gm, time_ls

time_matvec=0.0d0
time_feast=0.0d0
time_gmres=0.0d0
time_ls=0.0d0


!!!!!!!!!!!!! Set up Eigenvalue problem:

n=100
allocate(A(n,n),ztmpmat1(n,n),ztmpmat2(n,n),ztmpmat3(n,n))
call zrandom_subspace(ztmpmat1,n,n)
do i=1,n
    A(i,i)=i*(1.0d0,0.0d0)
end do
call zgemm('N','C',n,n,n,(1.0d0,0.0d0),A,n,ztmpmat1,n,(0.0d0,0.0d0),ztmpmat2,n)
call zgemm('N','N',n,n,n,(1.0d0,0.0d0),ztmpmat1,n,ztmpmat2,n,(0.0d0,0.0d0),A,n)

!!!!!!!!!!!!!! Set up FEAST:

m0=15 !FEAST subspace size
!allocate(vaavn(m0))
allocate(res(m0))
allocate(e(m0))
allocate(work(n,m0))
allocate(aq(m0,m0))
allocate(bq(m0,m0))
allocate(x(n,m0))
allocate(zwork(n,m0))

call feastinit(fpm)
fpm(1)=1 !turn on FEAST output
fpm(4)=50 !maximum FEAST iterations
fpm(2)=4 !number of FEAST contour points, i.e. number of linear systems to solve at each FEAST iteration
fpm(3)=6 !accuracy of FEAST iterations, i.e. stop when residual is less than 10^-fpm(3)

!eigenvalue interval:
emin=0.5d0
emax=10.5d0

!!!!!!!!!!!!!  Set up linear system solver:

linIterations=2 !number of iterations to use; more linear system iterations allows FEAST to converge faster
kdim=1 !number of Krylov subspace blocks to use; 1 usually works fine, but more gives faster convergence
allocate(xsol(n,m0), V(n,m0*kdim), Av(n,m0*kdim), Ax(n,m0), zwork2(n,m0*kdim), linwork1(n,m0), linwork2(n,m0))

!!!!!!!!!!!!  Run FEAST RCI:
call system_clock(count=totalc1)
ijob=-1
do while (ijob .ne. 0)
call zfeast_hrci(ijob,n,ze,work,zwork,aq,bq,fpm,epsout,loop,emin,emax,m0,e,x,m,res,info)
	
    select case (ijob)
		case (10) !store complex shift ze
		
            ze2=ze	
		
        case (11) !solve linear system Az*Qz=zwork for Qz, put Qz in zwork
            ztmpmat3(1:n,1:n)=0.0d0
            do i=1,n
                ztmpmat3(i,i)=ze2
            end do
            ztmpmat3=ztmpmat3-A
            ztmpmat1(1:n,1:m0)=zwork(1:n,1:m0) 
            
            linState(1)=-1 !set linsate(1)=-1 to begin RCI routine dfeast_gmres
            do while (linState(1) .ne. 0)
                call system_clock(count=c1)     
                call zfeast_gmres(linjob,linState,zwork,xsol,V,Av,Ax,ze2,n,m0,linIterations,kdim,linwork1,linwork2,zwork2,times_gm)
                call system_clock(count=c2)
                time_gmres=time_gmres+elapsed_time(c1,c2)
   
                if(linjob == 30) then
                    !Your matrix multiplication routine goes here! do linwork2=A*linwork1
                    call system_clock(count=c1)
                    call zgemm('N','N',n,m0,n,(1.0d0,0.0d0),A,n,linwork1,n,(0.0d0,0.0d0),linwork2,n)
                    !linwork2=matmul(A,linwork1)
                    call system_clock(count=c2)
                    time_matvec=time_matvec+elapsed_time(c1,c2)
                end if
            end do
            time_ls=time_ls+times_gm

            ztmpmat2(1:n,1:m0)=ztmpmat1(1:n,1:m0)-matmul(ztmpmat3,zwork(1:n,1:m0))
            do i=1,m0
                ztmpmat2(1:n,i)=ztmpmat2(1:n,i)/dznrm2(n,ztmpmat1(1:n,i),1)
            end do

            dtmp1=dznrm2(n,ztmpmat2(1:n,1),1)
            do i=1,m
                dtmp2=dznrm2(n,ztmpmat2(1:n,i),1)
                if(dtmp2>dtmp1) then
                    dtmp1=dtmp2
                end if
            end do

            print *,'     Lin sys error=',dtmp1
            

        case(21)

            linState(1)=-1 !set linsate(1)=-1 to begin RCI routine dfeast_gmres
            do while (linState(1) .ne. 0)

                call system_clock(count=c1)
                call zfeast_gmres(linjob,linState,zwork,xsol,V,Av,Ax,conjg(ze2),n,m0,linIterations,kdim,linwork1,linwork2,zwork2,times_gm)
                call system_clock(count=c2)
                time_gmres=time_gmres+elapsed_time(c1,c2)
                
                if(linjob == 30) then
                    !Your matrix multiplication routine goes here! do linwork2=A*linwork1
                    call system_clock(count=c1)
                    call zgemm('N','N',n,m0,n,(1.0d0,0.0d0),A,n,linwork1,n,(0.0d0,0.0d0),linwork2,n)
                    !linwork2=matmul(A,linwork1)
                    call system_clock(count=c2)
                    time_matvec=time_matvec+elapsed_time(c1,c2)
                end if
            end do 
            time_ls=time_ls+times_gm

            
        case (30) !A*x	
            
            !Your matrix multiplication routine goes here! Do work=A*x
            call zgemm('N','N',n,m0,n,(1.0d0,0.0d0),A,n,x,n,(0.0d0,0.0d0),work,n)
		case (40) !B*x
        
            call zlacpy('F',n,m0,x,n,work,n)

	end select
end do

call system_clock(count=totalc2)

!print results:
print *,'FEAST finished; # eigs found = ',m
print *, 'eigenvalues and eigenvector residuals:'
do i=1,m
    print *,i,e(i),res(i)
end do

print *,'total time=',elapsed_time(totalc1,totalc2)
print *,''
print *,'matvec time=',time_matvec
print *,'gmres time=',time_gmres
print *,'ls time=',time_ls
print *,'feast time=',time_feast
print *,'sum =',time_matvec+time_gmres+time_feast
print *,'gmres frac=', time_gmres/elapsed_time(totalc1,totalc2)
print *,'ls frac=',time_ls/time_gmres

deallocate(A,ztmpmat1,ztmpmat2,ztmpmat3)
deallocate(xsol, V, Av, Ax, zwork2, linwork1, linwork2)
deallocate(res,e,work,aq,bq,x,zwork)

end program

subroutine zrandom_subspace(x,n,m)
    integer :: n,m,info
    complex (kind=kind(0.0d0)), dimension(n,m) :: x

    complex (kind=kind(0.0d0)), dimension(m) :: qrtau
    integer :: lwork
    complex (kind=kind(0.0d0)), dimension(:),allocatable :: work

    integer, dimension(4) :: iseed

    lwork=3*n
    allocate(work(lwork))

    !call random_number(x)
    iseed=(/4,3,2,1/)
    call CLARNV(2,iseed,n*m,x)
    !do qr
    call ZGEQRF( n, m, x, n, qrtau, work, lwork, info )
    if (info .ne. 0) then
        print *,'Problem with generating random subspace'
        print *,'ZGEQRF error info = ',info
    end if

    !get q matrix
    !use zungqr for complex case, DORGQR for double precision case
    call ZUNGQR(  n, m, m, x, n, qrtau, work, lwork, info )
    if (info .ne. 0) then
        print *,'Problem with generating random subspace'
        print *,'ZUNGQR error info = ',info
    end if
end subroutine zrandom_subspace

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

