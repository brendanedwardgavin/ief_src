program testSort

implicit none

double precision, dimension(:),allocatable :: list
integer :: i,n,issorted


n=10

allocate(list(n))
!~call random_seed()
call random_number(list)

print *,'list'

do i=1,n
print *,list(i)
end do

print *,'sorted'

call quicksort(list,1,n)
!call selection_sort(list,n)

do i=1,n
    print *,list(i)
end do

issorted=1
do i=1,n-1
    if (list(i+1)<list(i)) then
        issorted=0
        exit
    end if
enddo

print  *,'issorted=',issorted

end program


