module rundata

    !!!!!!!!data by iteration:

    !(feast it, cp#, lin it,rhs#,)
    double precision, dimension(:,:,:,:),allocatable :: linres,lintime
    !(feast it,cp#,rhs #,)
    double precision, dimension(:,:,:),allocatable :: linit

    !(feast it)
    double precision, dimension(:), allocatable :: eigtime,eigres

    !(cp #)
    complex (kind=kind(0.0d0)), dimension (:), allocatable :: cpval

    !timing data:
    double precision :: qrtime,lstime,gstime,arnolditime,mvtime
    integer, dimension(:,:),allocatable :: nmatvec !(feastit,cpnum)

    integer :: startCount
    integer :: feastit,cpnum !keep track of which feast iteration, contour point we're on 


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



    subroutine initrundata(ncp,m0,maxfeastit,maxlinit)
        implicit none
        integer :: ncp,m0,maxfeastit,maxlinit
        allocate(linres(0:maxfeastit,ncp,maxlinit,m0), lintime(0:maxfeastit,ncp,maxlinit,m0))
        allocate(linit(0:maxfeastit,ncp,m0))
        allocate(eigtime(0:maxfeastit),eigres(0:maxfeastit),cpval(ncp))
        allocate(nmatvec(0:maxfeastit,ncp))

        qrtime=0.0d0
        lstime=0.0d0
        mvtime=0.0d0
        gstime=0.0d0
        eigtime=0.0d0
        eigres=0.0d0
        arnolditime=0.0d0
        linres=-1.0d0
        
    end subroutine initrundata

end module
