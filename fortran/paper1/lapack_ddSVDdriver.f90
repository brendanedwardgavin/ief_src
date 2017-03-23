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
double precision,dimension(:),allocatable :: dsa,dca,ssa
integer,dimension(:),allocatable :: isa,jsa,ica,jca,sisa,sjsa
integer :: nnza

double precision,dimension(:,:),allocatable :: A
integer :: row,column
!!!!!!!!!!!!!!!!!!!!!!!!  FEAST:

integer :: n,m,m0
integer, dimension(64)::fpm
double precision :: emin,emax

!!!!!!!!!!!!!! lapack svd
double precision, dimension(:),allocatable :: S
double precision, dimension(:,:), allocatable :: U,VT
integer :: lwork,info
double precision, dimension(:),allocatable :: work

allocate(U(1,1),VT(1,1)) !since we're not calculating U or V


call feastinit(fpm)
call getarg(1,name)


if(iargc()<1) then
    print *,'Too few arguments. Usage:'
    !print *,'time_sddriver <matrix base>  <input file base>  <output file base>'
    print *,'lapack_ddSVDdriver <input/output file base>'
    stop
end if
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
  read(10,*) n,m,nnza

  print *, 'size=',n,m
  allocate(A(n,m))

  
  do i=1,nnza
     read(10,*) row,column,repart
     A(row,column)=repart
     !print *,ica(i),jca(i),dca(i)
  end do

  close(10)
100 continue



  print *,'done reading.'

!print *,"isa:",isa
!print *,''
!print *,'jsa:',jsa
!print *,''
!print *,'dsa:',dsa
! stop

print *,'Doing SVD'
    lwork=max(3*min(n,m)+max(n,m),5*min(n,m))
    allocate(S(min(n,m)),work(lwork))

    call dgesvd('N','N',n,m,A,n,S,U,n,VT,m,work,lwork,info)

    if(info .ne. 0) then
        print *,'DGESVD error: ',info
        stop
    end if

print *,'SVD done'

!    print *,'Singular values:'
!    print *,''
!  do i=1,n
!    print *,i,S(i)
!  end do
 
    print *,'writing file'
    open(unit=10,file=trim(name)//'.svds',status='REPLACE')
   open(unit=11,file=trim(name)//'_plot.svds',status='REPLACE')
        write(10,*) min(n,m)
        do i=1,min(n,m)
            write(10, *) i,S(i)
            write(11, *) S(i),0.0d0
        end do
    close(10)
    close(11)
 
end program

