program lapack_dddriver

implicit none

integer :: i,j,k
double precision :: repart,impart

!!!!!!!!!!!!!!!!!!!!!!! read file:
double precision, allocatable :: array(:)

integer :: pc
character(len=100) :: name
character(len=1) UPLO,PRE,SHG,EG
character(len=1) :: cc
complex (kind=kind(0.0d0)),dimension(:),allocatable :: dsa,dca,ssa
integer,dimension(:),allocatable :: isa,jsa,ica,jca,sisa,sjsa
integer :: nnza

complex (kind=kind(0.0d0)),dimension(:,:),allocatable :: A
integer :: row,column
!!!!!!!!!!!!!!!!!!!!!!!!  FEAST:

integer :: n,m0
integer, dimension(64)::fpm
double precision :: emin,emax


!!!!!!!!!!!!!!!!!!! lapack
integer :: info,lwork
double precision, dimension(:),allocatable :: e
complex (kind=kind(0.0d0)),dimension(:),allocatable :: work
double precision, dimension(:),allocatable :: rwork


call feastinit(fpm)
call getarg(1,name)

!!!!!!!!!!!! DRIVER_FEAST_SPARSE input file  
!  open(10,file=trim(name)//'.in',status='old')
!  read(10,*) SHG ! type of eigenvalue problem "General, Hermitian, Symmetric" 
!  read(10,*) EG ! type of eigenvalue probl  g== sparse generalized, e== sparse standard
!  read(10,*) PRE  ! "PRE"==(s,d,c,z) resp. (single real,double real,complex,double complex) 
!  read(10,*) UPLO ! UPLO==(F,L,U) reps. (Full csr, Lower csr, Upper csr)

!  read(10,*) emin
!  read(10,*) emax

!  read(10,*) m0
!  read(10,*) pc ! Some changes from default for fpm
!  do i=1,pc
!     read(10,*) j,fpm(j)
!  enddo

!  close(10)


!goto 100
!!!!!!!!!!!read matrix A
  print *,'reading matrix'
  open(10,file=trim(name)//'_A.mtx',status='old')
  k=0
  cc='%'
  do while(cc=='%')
     k=k+1
     read(10,'(A1)') cc
  end do
  close(10)

  open(10,file=trim(name)//'_A.mtx',status='old')
  do i=1,k-1
     read(10,'(A1)') cc
  enddo
  read(10,*) n,n,nnza

  allocate(A(n,n))

  
  do i=1,nnza
     read(10,*) row,column,repart
     A(row,column)=(1.0d0,0.0d0)*repart
     A(column,row)=A(row,column)
     !print *,ica(i),jca(i),dca(i)
  end do

  close(10)
100 continue

goto 101
n=2315
allocate(array(n*(n+1)/2))
open(24,file=trim('../../matrices/QE/H_matrix_2315'),form='unformatted',status='old',access='sequential')
read(24) array
close(24)
allocate(A(n,n))

do j=1,n
    do i=1,j
        A(j,i)=(1.0d0,0.0d0)*array(j+(i-1)*(2*n-i)/2)
    end do
end do
do j=1,n
    do i=j+1,n
        A(i,j)= A(j,i)
    end do
end do
101 continue




  print *,'done reading.'

!print *,"isa:",isa
!print *,''
!print *,'jsa:',jsa
!print *,''
!print *,'dsa:',dsa
! stop

  lwork=2*n-1
  allocate(work(lwork))
  allocate(rwork(3*n-2))
  allocate(e(n))

    print *,'diagonalizing'
  call zheev('N','U',n,A,n,e,work,lwork,rwork,info)


    !print *,'Eigenvalues:'
    !print *,''
  !do i=1,n
  !  print *,i,e(i)
  !end do
    print *,'writing file'
    open(unit=10,file=trim(name)//'.eigs',status='REPLACE')
   open(unit=11,file=trim(name)//'_plot.eigs',status='REPLACE') 
        write(10,*) n
        do i=1,n
            write(10, *) i,e(i)
            write(11, *) e(i),0.0d0
        end do
    close(10)
    close(11)
end program

