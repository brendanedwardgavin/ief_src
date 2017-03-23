program time_sddriver
!use rundata

implicit none
!#ifdef MPI
!include "mpif.h"
!#endif

integer :: i,j,k
double precision :: repart,impart
!!!!!!!!!!!!!!!!!!!!!!! read file:

integer :: pc
character(len=100) :: name,inname
!character(len=100) :: outname

character(len=1) UPLO,PRE,SHG,EG
character(len=1) :: cc
complex (kind=kind(0.0d0)),dimension(:),allocatable :: ssa
double precision, dimension(:),allocatable :: dsa,dca
integer,dimension(:),allocatable :: isa,jsa,ica,jca,sisa,sjsa
integer :: nnza,nnza2

!!!!!!!!!!!!!!!!! SVD:
integer :: na,ma
integer :: ica1,jca1
double precision :: dca1

double precision, dimension(:,:),allocatable :: A

integer, dimension(:),allocatable :: rowcheck,colcheck
!!!!!!!!!!!!!!!!!!!!!!!!  FEAST:

integer :: n,m0,loop,mode,info
integer, dimension(64)::fpm
double precision :: emin,emax,epsout
double precision, dimension(:),allocatable :: E,res
double precision, dimension(:,:),allocatable :: X

call feastinit(fpm)
if(iargc()<1) then
    print *,'Too few arguments. Usage:'
    !print *,'time_sddriver <matrix base>  <input file base>  <output file base>'
    print *,'make_augmat <matrix base>'
    stop
end if
call getarg(1,name)
!call getarg(2,inname)
!call getarg(2,outname)

!if (iargc()==2) then
!    call getarg(2,outname)
!else
!    outname=""
!end if

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
  read(10,*) na,ma,nnza
  
    nnza2=nnza*2+na+ma
      

    allocate(rowcheck(na+ma),colcheck(na+ma))

  allocate(ica(nnza2))
  allocate(jca(nnza2))
  allocate(dca(nnza2))
  k=1
  do i=1,nnza
     read(10,*) ica1,jca1,dca1 
     !ica(i)=ica1
     !jca(i)=na+jca1
     !dca(i)=dca1
     !ica(i)=jca1+na
     !jca(i)=ica1
     !dca(i)=dca1

     ica(k)=ica1
     jca(k)=na+jca1
     dca(k)=dca1
     ica(k+1)=jca1+na
     jca(k+1)=ica1
     dca(k+1)=dca1

    k=k+2

    rowcheck(ica(i))=1
    colcheck(jca(i))=1
    rowcheck(jca(i))=1
    colcheck(ica(i))=1

     if(ica1>na .or. jca1>ma .or. ica(i)>(na+ma) .or. jca(i)>(na+ma) .or. jca1<1 .or. ica1<1) then
        print *,'ica,jca',ica1,jca1
        print *,'na,ma',na,ma
        stop
    end if
    !k=k+2
  end do
  close(10)

    nnza=nnza*2+na+ma

  do i=1,na+ma
     ica(nnza-na-ma+i)=i
     jca(nnza-na-ma+i)=i
     dca(nnza-na-ma+i)=0.0d0 
  end do

  !! create csr format
  n=na+ma

    !ica(nnza+1)=191669
    !jca(nnza+1)=191669
    !dca(nnza+1)=0.0d0


 open(unit=10,file=trim(name)//'aug_A.mtx',status='REPLACE')
        write (10,*) n,n,nnza
        do i=1,nnza
            write (10,*) ica(i),jca(i),dca(i)
        end do
        close(10) 

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
