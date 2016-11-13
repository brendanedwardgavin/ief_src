module rundata
!!!!!!!TODO
!!! save ritz values for each rhs at each feast iteration
!!! save eigenvector error for each rhs as each feast iteration

    !!!!!!!!data by iteration:
    character(len=100) :: outname

    !(feast it, cp#, lin it,rhs#,)
    double precision, dimension(:,:,:,:),allocatable :: linres,lintime
    !(feast it,cp#,rhs #,)
    integer, dimension(:,:,:),allocatable :: linit

    !(feast it,rhs #)
    double precision, dimension(:,:),allocatable :: ritzvals
    double precision, dimension(:), allocatable :: ratfunc

    !(feast it)
    double precision, dimension(:), allocatable :: eigtime,eigres
    double precision, dimension(:,:), allocatable :: eigresall
    
    !(cp #)
    complex (kind=kind(0.0d0)), dimension (:), allocatable :: cpval,znesave
    double precision, dimension(:), allocatable :: wnesave

    !timing data:
    double precision :: qrtime,lstime,gstime,arnolditime,mvtime,totaltime
    integer, dimension(:,:),allocatable :: nmatvec !(feastit,cpnum)

    integer :: startcount,endcount
    integer :: feastit,cpnum !keep track of which feast iteration, contour point we're on 

    integer ::maxfeastit,ncp,m0d,maxlinit

    !!!temporary variables:
    integer :: tc1,tc2

    contains

    double precision function elapsed_time(c1,c2)
        integer :: c1,c2
        integer :: diff
        integer :: maxcount,countrate

        call system_clock(count_rate=countrate,count_max=maxcount)

        if(c2<c1) then
            diff=maxcount+c2-c1
        else
            diff=c2-c1
        end if

        elapsed_time= dble(diff)/dble(countrate)
    end function elapsed_time



    subroutine initrundata(nc,m,mfeastit,mlinit)
        implicit none
        integer :: nc,m,mfeastit,mlinit
        ncp=nc
        m0d=m
        maxfeastit=mfeastit
        maxlinit=mlinit
        
        allocate(linres(0:maxfeastit,ncp,maxlinit,m0d), lintime(0:maxfeastit,ncp,maxlinit,m0d),ratfunc(m0d))
        allocate(linit(0:maxfeastit,ncp,m0d))
        allocate(eigtime(0:maxfeastit),eigres(0:maxfeastit),eigresall(0:maxfeastit,m0d),cpval(ncp),wnesave(ncp),znesave(ncp))
        allocate(nmatvec(0:maxfeastit,ncp))
        allocate(ritzvals(0:maxfeastit,m0d))

        totaltime=0.0d0
        qrtime=0.0d0
        lstime=0.0d0
        mvtime=0.0d0
        gstime=0.0d0
        eigtime=0.0d0
        eigres=0.0d0
        eigresall=0.0d0
        arnolditime=0.0d0
        linres=-1.0d0
        linit=0
        ritzvals=0.0d0
        ratfunc=0.0d0
    end subroutine initrundata

    subroutine printTimes()
        print *,''
        print *,'QR: ',qrtime,qrtime/totaltime
        print *,'Least Squares: ',lstime,lstime/totaltime
        print *,'Matvec: ',mvtime, mvtime/totaltime
        print *,'GS: ',gstime,gstime/totaltime
        print *,'TOTAL TIME: ',totaltime
    end subroutine printTimes

    subroutine printMinCpDist(E)
        implicit none
        double precision, dimension(m0d) :: E
        double precision :: mindist,dist
        integer :: i,j

        print *,''
        do i=1,ncp
            mindist=abs(cpval(i)-dcmplx(E(1)))
            do j=2,m0d
               dist=abs(cpval(i)-dcmplx(E(j))) 
               if(dist>mindist) mindist=dist
            end do
            print *,'cp ',i,mindist
            print *,'      ',cpval(i)
        end do
    end subroutine

    subroutine printCpVal()
        print *,''
        print *,'CP Vals:'
        do i=1,ncp
            print *,i,cpval(i)
        end do
    end subroutine

    subroutine savedata(feastloop)
        implicit none 
        integer :: feastloop
        integer :: i,j,k,temp
        character(len=100) :: cpstr,m0str
  
        double precision, dimension(:),allocatable :: delta

        write(m0str,"(I5)") m0d 

        !save eigenvector residuals:
        open(unit=10,file=trim(outname)//'_eigresiduals.dat',status='REPLACE')
        do i=0,feastloop
            write (10,"(I3, 2ES15.5)") i+1,eigtime(i),eigres(i)
            !write (10,"(I3, 3ES15.5, F8.2)") i,timelist(i),reslist(i),linsysreslist(i),linsysitavg(i)
            !write(10,*) i,timelist(i),reslist(i),linsysreslist(i)
        end do
        close(10) 

        !save all eigenvector residuals
        open(unit=10,file=trim(outname)//'_eigresidualsall.dat',status='REPLACE')
        do i=0,feastloop
            write (10,"(I3,ES15.5)",advance="no") i+1,eigtime(i)
            write (10,"("//m0str//"ES15.5)") (eigresall(i,k), k=1,m0d)
        end do
        close(10)       

        !save eigenvector convergence (from residuals)
        open(unit=10,file=trim(outname)//'_eigconv.dat',status='REPLACE')
        write (10,"(I3, ES15.5)") 1,0.0d0
        do i=1,feastloop
            write (10,"(I3, ES15.5)") i+1,eigres(i)/eigres(i-1)
        end do
        close(10)
   

        !save all eigenvector convergence rates
        open(unit=10,file=trim(outname)//'_eigconvall.dat',status='REPLACE')
        write (10,"(I3)",advance="no") 0
        write (10,"("//m0str//"ES15.5)") (0.0d0, k=1,m0d)
        do i=1,feastloop
            write (10,"(I3)",advance="no") i+1
            write (10,"("//m0str//"ES15.5)") (eigresall(i,k)/eigresall(i-1,k), k=1,m0d)
        end do
        close(10)       


        !save linear system iterations: 
        do i=1,ncp
            write(cpstr,"(I5)") i
            open(unit=10,file=trim(outname)//'_linsysit'//trim(adjustl(cpstr))//'.dat',status='REPLACE')
            do j=0,feastloop
                write (10,"(I3)",advance="no") j+1
                write (10, "("//m0str//"I4)") (linit(j,i,k), k=1,m0d)
            end do
            close(10)
        end do

        !save total and average linear system iterations at each cp:
        do i=1,ncp
            write(cpstr,"(I5)") i
            open(unit=10,file=trim(outname)//'_linsysittotalcp'//trim(adjustl(cpstr))//'.dat',status='REPLACE')
            do j=0,feastloop
                temp=0
                do k=1,m0d
                    temp=temp+linit(j,i,k)
                enddo
                write (10,"(I3)",advance="no") j+1
                write (10, "(I4,F4.1)") temp,(1.0d0*temp)/(1.0d0*m0d)
            end do
            close(10)
        end do


        !save total linear system iterations and total matvecs at each FEAST iteration
        open(unit=10,file=trim(outname)//'_linsystotalmatvec.dat',status='REPLACE')
        do j=0,feastloop     
           temp=0
           do i=1,ncp
                do k=1,m0d
                    temp=temp+linit(j,i,k)
                end do
           end do
           write(10,"(2I4)") j+1,temp
        end do 
        close(10)

        !save ritz values
        open(unit=10,file=trim(outname)//'_ritzvals.dat',status='REPLACE')
        do j=0,feastloop
            write (10,"(I3)",advance="no") j+1
            write (10, "("//m0str//"ES15.5)") (ritzvals(j,k), k=1,m0d)
        end do
        close(10)
        
        !save contour point values
        open(unit=10,file=trim(outname)//'_cpvals.dat',status='REPLACE')
        do i=1,ncp
            write (10,"(2ES15.5)") dble(cpval(i)),aimag(cpval(i))
        end do 
        close(10)

        !measure and save theoretical convergence rate
        !do i=1,ncp
           !find smallest distance of the current cp to an eigenvalue
           !add alpah*abs(wne(i)*1/abs(mindist))
        !end do
        !udpate each entry in delta for each eigenvalue using ratfunc
        !save file
       
    end subroutine savedata
end module
