module rundata

    !(cp#,rhs#,feast it, lin it)
    double precision, dimension(:,:,:,:),allocatable :: linres,lintime
    !(cp#,rhs #,feast it)
    double precision, dimension(:,:,:),allocatable :: linit

    !(feast it)
    double precision, dimension(:), allocatable :: eigtime,eigres

    !(cp #)
    complex (kind=kind(0.0d0)), dimension (:), allocatable :: cpval

    integer :: startCount

    contains

    subroutine initrundata(ncp,m0,maxfeastit,maxlinit)
        implicit none
        integer :: ncp,m0,maxfeastit,maxlinit
        allocate(linres(ncp,m0,maxfeastit,maxlinit), lintime(ncp,m0,maxfeastit,maxlinit))
        allocate(eigtime(maxfeastit),eigres(maxfeastit),cpval(ncp))
    end subroutine initrundata

end module
