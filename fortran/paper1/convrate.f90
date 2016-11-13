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
call getarg(1,name)

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

    alpha=1.0d0/(1.0d0*fpm(54))

    allocate(Zne(1:fpm(2)),Wne(1:fpm(2)))
    if(fpm(55)==0) then
        call zfeast_contour(emin,emax,fpm(2),fpm(16),fpm(18),Zne,Wne)
    else
        delta=(emax-emin)/fpm(2)
        do i=1,fpm(2)
            call zfeast_contour(emin+(i-1)*delta,emin+i*delta,1,fpm(16),fpm(18),Zne(i),Wne(i))
        end do
    end if

!goto 100
!!!!!!!!!!!read eigenvalues
  print *,'reading eigenvalues'
  open(10,file=trim(name)//'.eigs',status='old')
    read(10,*) n
    allocate(e(n),ratfunc(n))
    do i=1,n
        read(10,*) j, e(i)
    end do 
    close(10)


!calculate base convergence rate:
neigs=0
do i=1,n
    if(e(i) >= emin .and. e(i) <= emax) neigs=neigs+1
end do


call dfeast_rationalx(Zne,Wne,fpm(2),e,n,ratfunc)

ratfunc=abs(ratfunc)
call quicksort(ratfunc,1,n)

gammam1=ratfunc(n-(m0+1))
gammai=ratfunc(n-neigs)

print *,'Gammas:',gammam1,gammai

!base convergence rate = gammam1/gammai
brate=gammam1/gammai

print *,"Base rate:",brate

!modified convergence rate:
delta=0.0d0
do i=1,fpm(2)
    minind=1
    mindist=abs(Zne(i)-(1.0d0,0.0d0)*e(1))
    do j=1,n
        newdist=abs(Zne(i)-(1.0d0,0.0d0)*e(j))
        if(newdist<mindist) then
            mindist=newdist
            minind=j
        end if
    end do
    delta=delta+abs(Wne(i)/mindist)
end do

!modified rate: alpha*delta/gammai
mrate=alpha*delta

print *,'Delta=',delta

print *,'Total rate=',mrate/gammai+brate


end program

