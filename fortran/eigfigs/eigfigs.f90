program rational_function

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
    complex(kind=(kind(1.0d0))),dimension(:),allocatable :: Zne,Wne

    if (iargc()<5) then
        print *,"         Missing arguments. Usage: ./rational_function emin emax samples ellipse% #contourpoints"
        stop
    end if

    call getarg(1,arg)
    read(arg,*) emin

    call getarg(2,arg)
    read(arg,*) emax

    call getarg(3,arg)
    read(arg,*) samples

    call getarg(4,arg)
    read(arg,*) ellipse_percent

    call getarg(5,arg)
    read(arg,*) ncp

    emin=-0.5d0
    emax=0.5d0
    samples=500
    

    m0=samples

    allocate(A(samples,samples),X(samples,samples))
    allocate(e(samples),res(samples))
    allocate(Zne(ncp),Wne(ncp))

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
    fpm(18)=ellipse_percent
    fpm(2)=ncp

    !call zfeast_contour(emin,emax,ncp,fpm(16),ellipse_percent,Zne,Wne)
    !do i=1,ncp
    !    print *,i,Zne,Wne
    !end do
    !stop
    
    call zfeast_contour(emin,emax,ncp,fpm(16),ellipse_percent,Zne,Wne)
    if(ellipse_percent==0) then   
        do i=1,ncp
            mid=dble(Zne(i))
            Zne(i)=mid*(1.0d0,0.0d0)+(0.0d0,1.0d-9)
        end do
    end if

    call dfeast_syevx('F',samples,A,samples,fpm,epsout,loop,emin,emax,m0,e,X,m,res,info,Zne,Wne)

    open(unit=10,file="rational_function.dat")
    do i=1,samples
        write(10,*) -1.0d0+dble(i-1)*2/dble(samples), abs(X(i,i))
    end do
    close(10)

    open(unit=10,file="contour_points.dat")
    do i=1,ncp
        write(10,*) dble(Zne(i)),dimag(Zne(i))
    end do
    do i=1,ncp
        write(10,*) dble(Zne(i)),-1.0d0*dimag(Zne(i))
    end do
    close(10)

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
