program time_sddriver
use rundata
implicit none
!#ifdef MPI
!include "mpif.h"
!#endif

integer :: i,j,k
double precision :: repart,impart
!!!!!!!!!!!!!!!!!!!!!!! read file:

integer :: pc
character(len=100) :: name
character(len=100) :: outname1

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
integer :: linloops
double precision, dimension(:,:), allocatable :: dB
complex (kind=kind(0.0d0)),dimension(:,:), allocatable :: B,X,R,temp1,Bs
integer :: nrhs
complex (kind=kind(0.0d0)) :: ze

double precision,external :: dznrm2,zlange
double precision :: error,thenorm
character, dimension(6) :: matdescra
double precision, dimension(:),allocatable :: dwork
integer :: loops

call random_seed()

call feastinit(fpm)
call getarg(1,name)
if (iargc()==2) then
    call getarg(2,outname1)
else
    outname1=""
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

allocate(dB(n,m0),B(n,m0),X(n,m0),R(n,m0),temp1(n,m0),Bs(n,m0))

if(UPLO=='F') then
    matdescra(1)='G'
else
    matdescra(1)='H'
end if
!matdescra(1)='G'
matdescra(2)=UPLO
matdescra(3)='N'
matdescra(4)='F'


!!!!!!!!!!preset rhs:
!Bs=(0.0d0,0.0d0)
!Bs(1,1)=(1.0d0,0.0d0)
!Bs(2,1)=(2.0d0,0.0d0)
!Bs(3,1)=(3.0d0,0.0d0)
!Bs(4,1)=(4.0d0,0.0d0)
!if(m0>1) then
!Bs(5,2)=(1.0d0,0.0d0)
!Bs(6,2)=(2.0d0,0.0d0)
!Bs(7,2)=(3.0d0,0.0d0)
!Bs(8,2)=(4.0d0,0.0d0)
!end if
!call mkl_zcsrmm('N', n, m0, n, (-1.0d0,0.0d0), matdescra, dsa, jsa, isa, isa(2), Bs, n, (0.0d0,0.0d0), B, n)
!!!!!!!!!!!!!!!!!!!!!

!!!!! random rhs:
call random_number(dB)
!dB=0.0d0
!dB(1,1)=1.0
!dB(2,1)=2.5

B=(1.0d0,0.0d0)*dB

!ze=(2.0d0,2.0d0)
ze=(1.0d0,2.0d0)
print *,'ze=',ze
!call zfeast_cgls(UPLO,n,m0,dsa,isa,jsa,ze,nnza,B,X,100 )

!print *,'starting arnoldi'
call initrundata(fpm(2),m0,fpm(4),fpm(51)*fpm(50))
!call blockGMRESarnoldi(UPLO,n,m0,dsa,isa,jsa,(0.0d0,0.0d0),fpm(51),fpm(50),B,X,1.0d-7,loops,1)
!print *,'arnoldi done'

print *,'starting minres'
!call zminres(UPLO,n,dsa,isa,jsa,ze,B,X,1.0d-16,fpm(51),linloops)
call zminresBlock(UPLO,n,m0,dsa,isa,jsa,ze,B,X,1.0d-16,fpm(51),linloops,1)
!call zminres(UPLO,n,   dsa,isa,jsa,ze,B,X,1.0d-16,fpm(51),linloops)
print *,'minres done'

print *,'Sol=',X(1,1)




!test convergence

temp1=X
call mkl_zcsrmm('N', n, m0, n, (-1.0d0,0.0d0), matdescra, dsa, jsa, isa, isa(2), X, n, ze, temp1, n)
R=B-temp1

error=0.0d0
do i=1,m0
    thenorm=dznrm2(n,R(:,i),1)/dznrm2(n,B(:,i),1)
    if (thenorm>error) error=thenorm
end do

allocate(dwork(n))
print *,"Lin sys error = ",zlange('F',n,m0,R,n,dwork)/zlange('F',n,m0,B,n,dwork)

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
