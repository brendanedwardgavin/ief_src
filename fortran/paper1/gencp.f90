recursive subroutine quicksort(list,lo,hi)
    double precision, dimension(*), intent(inout) :: list
    integer, intent(in) :: lo,hi

    integer :: i,j
    double precision :: pivot,temp

    if(lo<hi) then
    
    pivot=list(hi)
    
    i=lo
    do j=lo,hi-1
        if (list(j)<=pivot) then
            temp=list(i)
            list(i)=list(j)
            list(j)=temp
            i=i+1
        end if
    end do
    temp=list(i)
    list(i)=list(hi)
    list(hi)=temp
    
    call quicksort(list,lo,i-1)
    call quicksort(list,i+1,hi)
    end if

end subroutine quicksort




program lapack_dddriver

implicit none

integer :: i,j,k,l
integer,parameter :: elmin=1,elmax=100,cpmin=1,cpmax=100
double precision :: repart,impart
double precision, dimension(cpmin:cpmax) :: delta_plotcp,trate_plotcp,brate_plotcp
double precision, dimension(elmin:elmax) :: delta_plotel,trate_plotel,brate_plotel
!!!!!!!!!!!!!!!!!!!!!!! read file:
double precision, allocatable :: array(:)

integer :: pc,ncp
character(len=100) :: name, inname
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
complex(kind=(kind(1.0d0))),dimension(:),allocatable :: Zne,Wne

!!!!!!!!!!!!!!!!!!! lapack
integer :: info,lwork
double precision, dimension(:),allocatable :: e,ratfunc
complex (kind=kind(0.0d0)),dimension(:),allocatable :: work
double precision, dimension(:),allocatable :: rwork

!!!!sorting and searching
integer :: ind1,ind2,minind,neigs
double precision:: gammam1,gammai,delta,mrate,alpha,newdist,mindist,brate


call feastinit(fpm)

if(iargc()<1) then
    print *,'Too few arguments. Usage:'
    !print *,'time_sddriver <matrix base>  <input file base>  <output file base>'
    print *,'gencp  <input/output file base>'
    stop
end if
!call getarg(1,inname)
call getarg(1,inname)

!!!!!!!!!!!! DRIVER_FEAST_SPARSE input file  
  open(10,file=trim(inname)//'.in',status='old')
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

    ncp=fpm(2)

!!!!!!!!!!!read eigenvalues
    alpha=1.0d0/(1.0d0*fpm(54))

    allocate(Zne(1:fpm(2)),Wne(1:fpm(2)))

    if(fpm(55)==0) then
        call zfeast_contour(emin,emax,fpm(2),fpm(16),fpm(18),Zne,Wne) 
    elseif (fpm(55)==1) then
        delta=(emax-emin)/fpm(2)
        do i=1,fpm(2)
            call zfeast_contour(emin+(i-1)*delta,emin+i*delta,1,fpm(16),fpm(2)*fpm(18),Zne(i),Wne(i))
        end do
    elseif (fpm(55)==2) then
       call zfeast_contour(Emin,Emax,fpm(2),fpm(16),fpm(18),Zne,Wne)
         do i=1,fpm(2)
            Zne(i)=Zne(i)-(0.0d0,1.0d0)*aimag(Zne(i))
         end do 
    end if


open(unit=10,file=trim(inname)//'_zne.dat',status='REPLACE')
open(unit=11,file=trim(inname)//'_wne.dat',status='REPLACE')
do i=1,fpm(2)
write(10,*) real(Zne(i)),aimag(Zne(i))
write(11,*) Wne(i)
end do
close(11)
close(12)

        open(unit=10,file=trim(inname)//'_cpvals_hplot.dat',status='REPLACE')
        do i=1,ncp
            write (10,*) dble(Zne(i)),aimag(Zne(i))
            write (10,*) dble(Zne(i)),-1.0d0*aimag(Zne(i))
        end do
        close(10)

end program
