program lapack_dddriver

implicit none

integer :: i,j,k

!!!!!!!!!!!!!!!!!!!!!!! read file:
double precision, allocatable :: array(:)

double precision,dimension(:,:),allocatable :: A

integer :: row,column
!!!!!!!!!!!!!!!!!!!!!!!!  FEAST:

integer :: n,m0
integer, dimension(64)::fpm
double precision :: emin,emax


!!!!!!!!!!!!!!!!!!! lapack
integer :: info,lwork
double precision, dimension(:),allocatable :: e
double precision,dimension(:),allocatable :: work
double precision, dimension(:),allocatable :: rwork


n=2315
allocate(array(n*(n+1)/2))
open(24,file=trim('../../matrices/QE/H_matrix_2315'),form='unformatted',status='old',access='sequential')
read(24) array
close(24)
allocate(A(n,n))

do j=1,n
    do i=1,j
        A(j,i)=array(j+(i-1)*(2*n-i)/2)
    end do
end do
do j=1,n
    do i=j+1,n
        A(j,i)= A(i,j)
    end do
end do

print *,'array'
print *,array(1:3)
print *,''
print *,'Matrix:'
do i=n-4,n
    print '(5ES12.3)',A(i,n-4:n)
end do




  print *,'done reading.'

!print *,"isa:",isa
!print *,''
!print *,'jsa:',jsa
!print *,''
!print *,'dsa:',dsa
! stop

  lwork=5*n-1
  allocate(work(lwork))
  allocate(rwork(3*n-2))
  allocate(e(n))

    print *,'diagonalizing'
  call dsyev('N','U',n,A,n,e,work,lwork,info)
    
    if(info .ne. 1) then
        print *,'lapack error',info
        stop
    end if
    print *,'Eigenvalues:'
    print *,''
  do i=1,n
    print *,i,e(i)
  end do
  
end program

