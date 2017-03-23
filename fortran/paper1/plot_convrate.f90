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
integer,parameter :: elmin=1,elmax=100,cpmin=1,cpmax=100,nalphas=1000
double precision :: alphamin,alphamax
double precision :: repart,impart
double precision, dimension(cpmin:cpmax) :: delta_plotcp,trate_plotcp,brate_plotcp
double precision, dimension(elmin:elmax) :: delta_plotel,trate_plotel,brate_plotel
double precision, dimension(nalphas) :: trate_plotalpha
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

integer :: n,m0,elpercent
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

alphamin=1.0d-10
alphamax=1.0d0

if(iargc()<2) then
    print *,'Too few arguments. Usage:'
    !print *,'time_sddriver <matrix base>  <input file base>  <output file base>'
    print *,'convrate <matrix base>  <input/output file base>'
    stop
end if
call getarg(1,name)
call getarg(2,inname)

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
    elpercent=fpm(18)

!!!!!!!!!!!read eigenvalues
  print *,'reading eigenvalues'
  open(10,file=trim(name)//'.eigs',status='old')
    read(10,*) n
    allocate(e(n),ratfunc(n))
    do i=1,n
        read(10,*) j, e(i)
    end do 
    close(10)

    alpha=1.0d0/(1.0d0*fpm(54))

    neigs=0
    do i=1,n
        if(e(i) >= emin .and. e(i) <= emax) neigs=neigs+1
    end do

    allocate(Zne(1:cpmax),Wne(1:cpmax))

do l=cpmin,cpmax
    fpm(2)=l
    if(fpm(55)==0) then
        call zfeast_contour(emin,emax,fpm(2),fpm(16),fpm(18),Zne(1:l),Wne(1:l))
    else!if (fpm(55)==1) then
        delta=(emax-emin)/fpm(2)
        do i=1,fpm(2)
            call zfeast_contour(emin+(i-1)*delta,emin+i*delta,1,fpm(16),fpm(2)*fpm(18),Zne(i),Wne(i))
        end do
    end if

!goto 100

!calculate base convergence rate:


call dfeast_rationalx(Zne(1:l),Wne(1:2),fpm(2),e,n,ratfunc)

ratfunc=abs(ratfunc)
call quicksort(ratfunc,1,n)

gammam1=ratfunc(n-(m0+1))
gammai=ratfunc(n-neigs)

!base convergence rate = gammam1/gammai
brate=gammam1/gammai

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

delta_plotcp(l)=delta
trate_plotcp(l)=mrate/gammai+brate
brate_plotcp(l)=brate

end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1! PLOT FOR DIFFERENT ELLIPSE VALUES
do l=elmin,elmax
    fpm(18)=10*l
    fpm(2)=ncp


    if(fpm(55)==0) then
        call zfeast_contour(emin,emax,fpm(2),fpm(16),fpm(18),Zne(1:l),Wne(1:l))
    else!if (fpm(55)==1) then
        delta=(emax-emin)/fpm(2)
        do i=1,fpm(2)
            call zfeast_contour(emin+(i-1)*delta,emin+i*delta,1,fpm(16),fpm(2)*fpm(18),Zne(i),Wne(i))
        end do
    end if

!goto 100

!calculate base convergence rate:


call dfeast_rationalx(Zne(1:l),Wne(1:l),fpm(2),e,n,ratfunc)

ratfunc=abs(ratfunc)
call quicksort(ratfunc,1,n)

gammam1=ratfunc(n-(m0+1))
gammai=ratfunc(n-neigs)

!base convergence rate = gammam1/gammai
brate=gammam1/gammai

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

delta_plotel(l)=delta
trate_plotel(l)=mrate/gammai+brate
brate_plotel(l)=brate

end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!! PLOT FOR DIFFERENT ALPHA VALUES

    fpm(18)=elpercent
    fpm(2)=ncp


    if(fpm(55)==0) then
        call zfeast_contour(emin,emax,fpm(2),fpm(16),fpm(18),Zne(1:l),Wne(1:l))
    else!if (fpm(55)==1) then
        delta=(emax-emin)/fpm(2)
        do i=1,fpm(2)
            call zfeast_contour(emin+(i-1)*delta,emin+i*delta,1,fpm(16),fpm(2)*fpm(18),Zne(i),Wne(i))
        end do
    end if

!goto 100

!calculate base convergence rate:


call dfeast_rationalx(Zne(1:l),Wne(1:l),fpm(2),e,n,ratfunc)

ratfunc=abs(ratfunc)
call quicksort(ratfunc,1,n)

gammam1=ratfunc(n-(m0+1))
gammai=ratfunc(n-neigs)

!base convergence rate = gammam1/gammai
brate=gammam1/gammai

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
do l=1,nalphas
    alpha=exp(log(alphamin)+(l-1)*(log(alphamax)-log(alphamin)/nalphas))
    mrate=alpha*delta
    trate_plotalpha(l)=mrate/gammai+brate
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

print *,'gammai,gammam1',gammai,gammam1
print *,'brate=',brate

open(unit=10,file=trim(inname)//'_plotConvrateEl_delta.dat',status='REPLACE')
open(unit=11,file=trim(inname)//'_plotConvrateEl_trate.dat',status='REPLACE')
open(unit=12,file=trim(inname)//'_plotConvrateEl_brate.dat',status='REPLACE')
do l=1,elmax
write(10,*) l*10,delta_plotel(l)
write(11,*) l*10,trate_plotel(l)
write(12,*) l*10,brate_plotel(l)
end do
close(10)
close(11)
close(12)



open(unit=10,file=trim(inname)//'_plotConvrateCp_delta.dat',status='REPLACE')
open(unit=11,file=trim(inname)//'_plotConvrateCp_trate.dat',status='REPLACE')
open(unit=12,file=trim(inname)//'_plotConvrateCp_brate.dat',status='REPLACE')
do l=1,cpmax
write(10,*) l,delta_plotcp(l)
write(11,*) l,trate_plotcp(l)
write(12,*) l,brate_plotcp(l)
end do
close(10)
close(11)
close(12)

open(unit=12,file=trim(inname)//'_plotConvrateAlpha_trate.dat',status='REPLACE')
do l=1,nalphas
    write(12,*) exp(log(alphamin)+(l-1)*(log(alphamax)-log(alphamin)/nalphas)),trate_plotalpha(l),log(trate_plotalpha(l))/log(10.0d0)
end do
close(12)

end program
