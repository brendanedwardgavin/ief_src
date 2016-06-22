program custom_contour
    
    integer :: i
    character(len=32) :: arg
    double precision :: emin,emax
    integer :: samples,ellipse_percent,ncp

    double precision, dimension(:,:),allocatable :: A,B,X

    integer :: ellipse_samples !number of samples for drawing ellipse in gnuplot

    !!!!!!FEAST
    integer, dimension(1:64) :: fpm
    double precision :: epsout,mid
    double precision, dimension(:),allocatable :: e,res
    integer :: loop,m,m0,info,p,info,m,n,lda,ldb
    complex(kind=(kind(1.0d0))),dimension(:),allocatable :: Zne,Wne,Zedge
    complex(kind=(kind(1.0d0))) :: emid
    integer, dimension(:),allocatable :: nedge,tedge

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

    m0=samples
    n=samples

    allocate(A(n,n),B(n,n),e(m0),X(n,m0),res(2*m0))

    !fill in A


    B=(0.0d0,0.0d0)
    do i=1,n
        B(i,i)=(1.0d0,0.0d0)
    end do

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

    p=4
    allocate(Zedge(p),Nedge(p),Tedge(p))
    Tedge=0
    Zedge(1)=emin*(1.0d0,0.0d0)-h*(0.0d0,1.0d0)
    Zedge(2)=emin*(1.0d0,0.0d0)+h*(0.0d0,1.0d0)
    Zedge(3)=emax*(1.0d0,0.0d0)+h*(0.0d0,1.0d0)
    Zedge(4)=emax*(1.0d0,0.0d0)-h*(0.0d0,1.0d0)

    Nedge(1)=3
    Nedge(2)=ncp
    Nedge(3)=3
    Nedge(4)=ncp

    fpm(2)=sum(nedge(1:p))
    allocate(Zne(fpm(2)),Wne(fpm(2)))
    call zfeast_customcontour(fpm(2),p,nedge,tedge,zedge,Zne,Wne)
    call zfeast_gegvx(n,A,n,B,n,fpm,epsout,loop,emid,r,m0,e,X,m,res,info,Zne,Wne)

end program
