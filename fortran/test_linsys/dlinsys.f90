program time_sddriver

implicit none
!#ifdef MPI
!include "mpif.h"
!#endif

integer :: i,j,k
double precision :: repart,impart
!!!!!!!!!!!!!!!!!!!!!!! read file:

integer :: pc
character(len=100) :: name
character(len=100) :: outname

character(len=1) UPLO,PRE,SHG,EG
character(len=1) :: cc
complex (kind=kind(0.0d0)),dimension(:),allocatable :: dsa,dca,ssa
integer,dimension(:),allocatable :: isa,jsa,ica,jca,sisa,sjsa
integer :: nnza

!!!!!!!!!!!!!!!!!!!!!!!!  FEAST:

integer :: n,m0
integer, dimension(64)::fpm
double precision :: emin,emax


!!!!!!!!!!!!!!!!!!!!!!! linear system

double precision, dimension(:,:), allocatable :: dB
complex (kind=kind(0.0d0)),dimension(:,:), allocatable :: B,X,R
integer :: nrhs

double precision,external :: dznrm2
double precision :: error,thenorm
character, dimension(6) :: matdescra

call random_seed()

call feastinit(fpm)
call getarg(1,name)
if (iargc()==2) then
    call getarg(2,outname)
else
    outname=""
end if

!!!!!!!!!!!! DRIVER_FEAST_SPARSE input file  
  open(10,file=trim(name)//'.in',status='old')
  read(10,*) SHG ! type of eigenvalue problem "General, Hermitian, Symmetric" 
  read(10,*) EG ! type of eigenvalue probl  g== sparse generalized, e== sparse standard
  read(10,*) PRE  ! "PRE"==(s,d,c,z) resp. (single real,double real,complex,double complex) 
  read(10,*) UPLO ! UPLO==(F,L,U) reps. (Full csr, Lower csr, Upper csr)

  read(10,*) emin
  read(10,*) emax

  read(10,*) m0
  read(10,*) pc ! Some changes from default for fpm
  do i=1,pc
     read(10,*) j,fpm(j)
  enddo

  close(10)



!!!!!!!!!!!read matrix A
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
  allocate(ica(nnza))
  allocate(jca(nnza))

  allocate(dca(nnza))
  do i=1,nnza
     read(10,*) ica(i),jca(i),repart
     dca(i)=(1.0d0,0.0d0)*repart
     !print *,ica(i),jca(i),dca(i)
  end do

  close(10)

  !! create csr format
  allocate(isa(1:n+1))
  allocate(jsa(1:nnza))
  allocate(dsa(1:nnza))
  call zcoo2csr(n,nnza,ica,jca,dca,isa,jsa,dsa)

!print *,"isa:",isa
!print *,''
!print *,'jsa:',jsa
!print *,''
!print *,'dsa:',dsa
! stop

!random rhs
nrhs=m0

allocate(dB(n,m0),B(n,m0),X(n,m0),R(n,m0))

call random_number(dB)
!dB=0.0d0
!dB(1,1)=1.0
!dB(2,1)=2.5

B=(1.0d0,0.0d0)*dB

call zfeast_cgne(UPLO,n,m0,dsa,isa,jsa,(1.0d0,0.0d0),nnza,B,X,100 )


!test convergence
matdescra(1)='H'
matdescra(2)=UPLO
matdescra(3)='N'
matdescra(4)='F'

R=B

call mkl_zcsrmm('N', n, m0, n, (-1.0d0,0.0d0), matdescra, dsa, jsa, isa, isa(2), X, n, (1.0d0,0.0d0), R, n)

error=0.0d0
do i=1,m0
    thenorm=dznrm2(n,R(:,i),1)/dznrm2(n,B(:,i),1)
    if (thenorm>error) error=thenorm
end do

print *,"Lin sys error = ",error

end program

subroutine zcoo2csr(n,nnz,ic,jc,c,isa,jsa,sa)
  implicit none
  integer :: n,nnz
  integer,dimension(*) :: ic,jc,isa,jsa
  complex (kind=kind(0.0d0)),dimension(*) :: c,sa
!!!
  integer :: k,k1,i,j,idum
  integer,dimension(n) :: iloc
  complex (kind=kind(0.0d0)) :: adum


  isa(1:N+1) = 0
  !find how many elements in row
  do  k=1, nnz
     isa(ic(k)) = isa(ic(k))+1
  end do



  !build isa
  k = 1
  do  i=1,N+1
     k1 = isa(i)
     isa(i) = k
     k = k+k1
  end do


  iloc(1:n)=isa(1:n)
  !Build jsa, sa - increment local row counter
  do  k=1, nnz
     sa(iloc(ic(k))) =  c(k)
     jsa(iloc(ic(k))) = jc(k)
     iloc(ic(k)) = iloc(ic(k))+1
  end do
  ! Reorder by increasing column
  do i=1,n
     do k=isa(i),isa(i+1)-1
        do k1=k,isa(i+1)-1
           if (jsa(k1)<jsa(k)) then
              idum=jsa(k)
              jsa(k)=jsa(k1)
              jsa(k1)=idum
              adum=sa(k)
              sa(k)=sa(k1)
              sa(k1)=adum
           endif
        enddo
     enddo
  enddo


end subroutine zcoo2csr
