program rational_function
implicit none
    integer :: i
    character(len=32) :: arg
    double precision :: emin,emax
    integer :: samples,ellipse_percent,ncp

    double precision, dimension(:,:),allocatable :: A,X

    integer :: ellipse_samples !number of samples for drawing ellipse in gnuplot

    !!!!!!FEAST
    integer, dimension(1:64) :: fpm
    double precision :: epsout,mid
    double precision, dimension(:),allocatable :: e,res
    integer :: loop,m,m0,info
    complex(kind=(kind(1.0d0))),dimension(:),allocatable :: Zne,Wne,zne2,wne2

    !!!!!!custom contour:
    double precision :: height
    integer :: counter,p
    integer, dimension(:),allocatable :: nedge,tedge
    complex(kind=(kind(1.0d0))),dimension(:),allocatable :: Zedge
    integer :: sides

    if (iargc()<6) then
        print *,"         Missing arguments. Usage: ./custom_contour emin emax samples height #contourpoints #sideints"
        stop
    end if

    call getarg(1,arg)
    read(arg,*) emin

    call getarg(2,arg)
    read(arg,*) emax

    call getarg(3,arg)
    read(arg,*) samples

    call getarg(4,arg)
    read(arg,*) height

    call getarg(5,arg)
    read(arg,*) ncp

    call getarg(6,arg)
    read(arg,*) sides

    m0=samples

    allocate(A(samples,samples),X(samples,samples))
    allocate(e(samples),res(samples))
    !allocate(Zne(ncp),Wne(ncp))

    X=0.0d0
    A=0.0d0
    do i=1,samples
        X(i,i)=1.0d0
        A(i,i)=-1.0d0+dble(i-1)*2/dble(samples)
    end do

    call feastinit(fpm)

    fpm(1)=1 !print feast output
    fpm(5)=1
    fpm(14)=1 !return only subspace
    !fpm(18)=ellipse_percent
    !fpm(2)=ncp

    !call zfeast_contour(emin,emax,ncp,fpm(16),ellipse_percent,Zne,Wne)
    !do i=1,ncp
    !    print *,i,Zne,Wne
    !end do
    !stop
   
    !
    !call zfeast_contour(emin,emax,ncp,fpm(16),ellipse_percent,Zne,Wne)
    !if(ellipse_percent==0) then   
    !    do i=1,ncp
    !        mid=dble(Zne(i))
    !        Zne(i)=mid*(1.0d0,0.0d0)+(0.0d0,1.0d-9)
    !    end do
    !end if

    !custom contour:
    p=4
    allocate(Zedge(p),Nedge(p),Tedge(p))
    Tedge=0
    Zedge(1)=emin*(1.0d0,0.0d0)-height*(0.0d0,1.0d0)
    Zedge(2)=emin*(1.0d0,0.0d0)+height*(0.0d0,1.0d0)
    Zedge(3)=emax*(1.0d0,0.0d0)+height*(0.0d0,1.0d0)
    Zedge(4)=emax*(1.0d0,0.0d0)-height*(0.0d0,1.0d0)

    Nedge(1)=sides
    Nedge(2)=ncp
    Nedge(3)=sides
    Nedge(4)=ncp

    fpm(2)=sum(nedge(1:p))
    allocate(Zne(fpm(2)),Wne(fpm(2)))
    call zfeast_customcontour(fpm(2),p,nedge,tedge,zedge,Zne,Wne)
    allocate(zne2(fpm(2)/2),wne2(fpm(2)/2))

    !get only points above real axis:
    counter=1
    do i=1,fpm(2)
        if (dimag(zne(i))>=0.0d0) then
            zne2(counter)=zne(i)
            wne2(counter)=wne(i)
            counter=counter+1
        end if
    end do

    fpm(2)=fpm(2)/2

    print *,'fpm(2)=',fpm(2)

    call dfeast_syevx('F',samples,A,samples,fpm,epsout,loop,emin,emax,m0,e,X,m,res,info,zne2,wne2)

    open(unit=10,file="rational_function.dat")
    do i=1,samples
        write(10,*) -1.0d0+dble(i-1)*2/dble(samples), abs(X(i,i))
    end do
    close(10)

    open(unit=10,file="contour_points.dat")
    do i=1,fpm(2)
        write(10,*) dble(Zne2(i)),dimag(Zne2(i))
    end do
    do i=1,fpm(2)
        write(10,*) dble(Zne2(i)),-1.0d0*dimag(Zne2(i))
    end do
    close(10)

    goto 100
    ellipse_samples=100
    deallocate(Zne,Wne)
    allocate(Zne(ellipse_samples),Wne(ellipse_samples))

    call zfeast_contour(emin,emax,ellipse_samples,1,ellipse_percent,Zne,Wne) !use trapezoidal contour points to draw actual ellipse

    call zselection_sort(Zne,ellipse_samples)

    open(unit=10,file="ellipse.dat")
    do i=1,ellipse_samples
        write(10,*) dble(Zne(i)),dimag(Zne(i))
    end do
    do i=ellipse_samples,1,-1
        write(10,*) dble(Zne(i)),-1.0d0*dimag(Zne(i))
    end do
    close(10)

    print *,'Rational function generation successful'
    100 continue
end program



subroutine zselection_sort(list,n)
    integer :: n
    complex(kind=(kind(1.0d0))),dimension(*) :: list
    complex(kind=(kind(1.0d0))) :: temp
    integer :: i

    do i=1,n-1
        do j=i+1,n
            if(dble(list(i))>dble(list(j))) then
                temp=list(i)
                list(i)=list(j)
                list(j)=temp
            end if
        end do
    end do
end subroutine zselection_sort
