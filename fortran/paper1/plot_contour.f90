program rational_function

    integer :: i,j
    character(len=32) :: arg
    double precision :: emin,emax
    integer :: samples,ellipse_percent,ncp

    !!!!input file stuff:
    character(len=300) :: name,outname,cpstr
    character :: SHG,EG,PRE,UPLO
    integer :: pc
    double precision::realpart,impart
    
    double precision :: delta,demin,demax
    double precision, dimension(:,:),allocatable :: A,X

    integer :: ellipse_samples !number of samples for drawing ellipse in gnuplot

    !!!!!!FEAST
    integer, dimension(1:64) :: fpm
    double precision :: epsout,mid
    double precision, dimension(:),allocatable :: e,res
    integer :: loop,m,m0,info
    complex(kind=(kind(1.0d0))),dimension(:),allocatable :: Zne,Wne

    if (iargc()<2) then
        print *,"         Too few arguments. Usage:"
        print *,"rational_function <input_name> <output_name>"
        stop
    end if
    
    call feastinit(fpm)

    call getarg(1,name)
    call getarg(2,outname)
 



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

       ncp=fpm(2)

    samples=1000
    m0=samples
    allocate(A(samples,samples),X(samples,samples))
    allocate(e(samples),res(samples))
    allocate(Zne(ncp),Wne(ncp))

    !read contour point information from data files
    open(10,file=trim(name)//'_zne.dat',status='old')
    open(11,file=trim(name)//'_wne.dat',status='old')
    do i=1,fpm(2)
       read(10,*) realpart,impart
       Zne(i)=(1.0d0,0.0d0)*realpart+(0.0d0,1.0d0)*impart
       read(11,*) Wne(i)
    end do
    close(10)
    close(11)

    X=0.0d0
    A=0.0d0
    delta=emax-emin
    demin=emin-0.1*delta
    demax=emax+0.1*delta
    delta=(demax-demin)/dble(samples)
    do i=1,samples
        X(i,i)=1.0d0
        A(i,i)=demin+(i-1)*delta
    end do

    fpm(1)=1 !print feast output
    fpm(5)=1
    fpm(14)=1 !return only subspace
    !fpm(18)=ellipse_percent
    !fpm(2)=ncp

    !call zfeast_contour(emin,emax,ncp,fpm(16),fpm(18),Zne,Wne)
    !do i=1,ncp
    !    print *,i,Zne,Wne
    !end do
    !stop
    
    !call zfeast_contour(emin,emax,fpm(2),fpm(16),fpm(18),Zne,Wne)
    !if(ellipse_percent==0) then   
    !    do i=1,ncp
    !        mid=dble(Zne(i))
    !        Zne(i)=mid*(1.0d0,0.0d0)+(0.0d0,1.0d-9)
    !    end do
    !end if

    call dfeast_syevx('F',samples,A,samples,fpm,epsout,loop,emin,emax,m0,e,X,m,res,info,Zne,Wne)
    !call dfeast_syev('F',samples,A,samples,fpm,epsout,loop,emin,emax,m0,e,X,m,res,info)

    print *,trim(outname)//"_ratfunc.dat"
    open(unit=10,file=trim(outname)//"_ratfunc.dat",status='REPLACE')
    do i=1,samples
        write(10,*) demin+(i-1)*delta, abs(X(i,i))
    end do
    close(10)

    !open(unit=10,file="contour_points.dat")
    !do i=1,ncp
    !    write(10,*) dble(Zne(i)),dimag(Zne(i))
    !end do
    !do i=1,ncp
    !    write(10,*) dble(Zne(i)),-1.0d0*dimag(Zne(i))
    !end do
    !close(10)

    ellipse_samples=100
    deallocate(Zne,Wne)
    allocate(Zne(ellipse_samples),Wne(ellipse_samples))
    if(fpm(55)==0) then

        ! 1 : use trapezoidal rule
        call zfeast_contour(emin,emax,ellipse_samples,1,fpm(18),Zne,Wne) !use trapezoidal contour points to draw actual ellipse

        call zselection_sort(Zne,ellipse_samples)
         
        open(unit=10,file=trim(outname)//"_ellipse1.dat",status='REPLACE')
        do i=1,ellipse_samples
            write(10,*) dble(Zne(i)),dimag(Zne(i))
        end do
        do i=ellipse_samples,1,-1
            write(10,*) dble(Zne(i)),-1.0d0*dimag(Zne(i))
        end do
        close(10)
! a change
    elseif(fpm(55)==1) then
        delta=(Emax-Emin)/fpm(2)
        do i=1,fpm(2)
            call zfeast_contour(Emin+(i-1)*delta,Emin+i*delta,ellipse_samples,1,fpm(2)*fpm(18),Zne,Wne)

            call zselection_sort(Zne,ellipse_samples)
        
            write(cpstr,"(I5)") i
            open(unit=10,file=trim(outname)//"_ellipse"//trim(adjustl(cpstr))//".dat",status='REPLACE')
            do j=1,ellipse_samples
                write(10,*) dble(Zne(j)),dimag(Zne(j))
            end do
            do j=ellipse_samples,1,-1
                write(10,*) dble(Zne(j)),-1.0d0*dimag(Zne(j))
            end do
            close(10)
        end do        
    end if

    print *,'Contour and Rational function generation successful'
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
