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
complex (kind=kind(0.0d0)),dimension(:),allocatable :: ssa
double precision, dimension(:),allocatable :: dsa,dca
integer,dimension(:),allocatable :: isa,jsa,ica,jca,sisa,sjsa
integer :: nnza

integer :: icain,jcain

!!!!!!!!!!!!!!!!!!!!!!!!  FEAST:

integer :: n,m,m0,loop,mode,info
integer, dimension(64)::fpm
double precision :: emin,emax,epsout
double precision, dimension(:),allocatable :: E,res
double precision, dimension(:,:),allocatable :: X

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

if(UPLO .ne. 'F') then
    print *,'needs full matrix, stopping'
    stop
end if

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
  read(10,*) n,m,nnza

  allocate(ica(2*nnza))
  allocate(jca(2*nnza))
  allocate(dca(2*nnza))
  do i=1,nnza
     read(10,*) icain,jcain,dca(i) 
    ica(i)=icain
    jca(i)=jcain+n
    ica(i+nnza)=jcain+m
    jca(i+nnza)=icain
    dca(i+nnza)=dca(i)
  end do

  close(10)

    nnza=nnza*2
    n=n+m

  !! create csr format
  allocate(isa(1:n+1))
  allocate(jsa(1:nnza))
  allocate(dsa(1:nnza))
  call dcoo2csr(n,nnza,ica,jca,dca,isa,jsa,dsa)

    allocate(E(m0),res(m0),X(n,m0))

    if (abs(emin)==abs(emax)) then
        print *,'singular value range is zero, exiting'
        stop
    end if

    if (abs(emin)<abs(emax)) then
        emin=abs(emin)
        emax=abs(emax)
    else
        emin=abs(emax)
        emax=abs(emin)
    end if

!print *,"isa:",isa
!print *,''
!print *,'jsa:',jsa
!print *,''
!print *,'dsa:',dsa
! stop

  !call time_szfeastgmres(UPLO,n,dsa,isa,jsa,fpm,emin,emax,m0,outname) 

  call dfeast_scsrevit(UPLO,n,dsa,isa,jsa,fpm,epsout,loop,emin,emax,m0,E,X,mode,res,info)

end program



subroutine dcoo2csr(n,nnz,ic,jc,c,isa,jsa,sa)
  implicit none
  integer :: n,nnz
  integer,dimension(*) :: ic,jc,isa,jsa
  double precision,dimension(*) :: c,sa
!!!
  integer :: k,k1,i,j,idum
  integer,dimension(n) :: iloc
  double precision :: adum


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


end subroutine dcoo2csr


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
