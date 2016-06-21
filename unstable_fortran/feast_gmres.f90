subroutine dfeast_gmres(ijob,stateVars,Brhs,x,V,Av,Ax,ze,n,m,restarts,m0,xwork,workin,Av2)
implicit none
    include 'f90_noruntime_interface.fi'
    
    integer :: ijob,n,m,restarts,m0,i,j,k
    integer, dimension(3)::stateVars

    !ijob: RCI case to call when this routine returns
    !n: dimension of linear sysmte
    !m: number of right hand sides
    !restarts: number of GMRES loops to do
    !m0: number of krylov subspace blocks; total krylov subspace size is m*m0
    !i,j,k: loop variables
    !stateVars(1) = current state of routine
    !stateVars(2) = index for outer loop
    !stateVars(3) = index for inner loop, for building krylov subspace

    double precision, dimension(n,*) :: xwork,workin !workspace variables
    complex (kind=kind(0.0d0)), dimension(n,*) :: x,V,Ax,Av,Brhs,Av2
    !x: solution to linear system
    !V: Krylov subspace in which to solve linear system; has dimension m*m0
    !Ax: A*x
    !Av: A*V; storing this means we only have to do one matrix multiply instead of two
    !Brhs: right hand sides
    complex (kind=kind(0.0d0)) :: ze !complex shift for FEAST

    !!lapack stuff:
    complex (kind=kind(0.0d0)), dimension(:), pointer :: work,qrtau
    integer :: lwork,info
    integer, dimension(:),pointer :: ipiv

    !!least squares stuff:
    complex (kind=kind(0.0d0)), dimension(:,:), pointer :: Rs,Bs
 
    i=stateVars(2)
    j=stateVars(3)
    
    !initialize routine when we start it
    if (stateVars(1)==-1) then
        V(1:n,1:m*m0)=0.0d0
        Av(1:n,1:m*m0)=0.0d0    
        Ax(1:n,1:m)=0.0d0
        x(1:n,1:m)=0.0d0
        stateVars(2)=1
        stateVars(3)=1
        stateVars(1)=1
        xwork(1:n,1:m)=0.0d0
    end if
    

    if (stateVars(1)==3) then !B*re(V(j)) for building V
        call zlacp2('F',n,m,workin(1,1),n,Av2(1,1),n)
        call zscal(n*m,ze,Av2(1,1),1)
        call zlacpy('F',n,m,Av2(1,1),n,V(1,m*j+1),n)
        !V(1:n,m*j+1:m*(j+1))=workin(1:n,1:m)*(1.0d0,0.0d0)*ze
        !xwork(1:n,1:m)=dimag(V(1:n,m*(j-1)+1:m*j)) !if we have B matrix
        workin(1:n,1:m)=dimag(V(1:n,m*(j-1)+1:m*j)) !if we don't have B matrix
        stateVars(1)=31        
        ijob=40
        return
    end if

    if (stateVars(1)==31) then !B*im(V(j)) for building V
        call zlacp2('F',n,m,workin(1,1),n,Av2(1,1),n)
        call zaxpy(n*m,ze*(0.0d0,1.0d0),Av2(1,1),1,V(1,m*j+1),1)
        !V(1:n,m*j+1:m*(j+1))=V(1:n,m*j+1:m*(j+1))+Av2(1:n,1:m)*(0.0d0,1.0d0)*ze!workin(1:n,1:m)*(0.0d0,1.0d0)*ze
        xwork(1:n,1:m)=dble(V(1:n,m*(j-1)+1:m*j))
        stateVars(1)=32        
        ijob=30
        return
    end if

    if (stateVars(1)==32) then !A*re(V(j)) for V
        call zlacp2('F',n,m,workin(1,1),n,Av2(1,1),n)
        call zaxpy(n*m,(-1.0d0,0.0d0),Av2(1,1),1,V(1:n,m*j+1),1)

        !V(1:n,m*j+1:m*(j+1))=V(1:n,m*j+1:m*(j+1))-workin(1:n,1:m)*(1.0d0,0.0d0)
        xwork(1:n,1:m)=dimag(V(1:n,m*(j-1)+1:m*j))
        stateVars(1)=33
        ijob=30
        return
    end if

    if (stateVars(1)==33) then !A*re(V(j)) for V
        call zlacp2('F',n,m,workin(1,1),n,Av2(1,1),n)
        call zaxpy(n*m,(0.0d0,-1.0d0),Av2(1,1),1,V(1,m*j+1),1)

        !V(1:n,m*j+1:m*(j+1))=V(1:n,m*j+1:m*(j+1))-workin(1:n,1:m)*(0.0d0,1.0d0)
    end if


    if (stateVars(1)==4) then !B*re(V(j)) for Av
        call zlacp2('F',n,m,workin(1,1),n,Av2(1,1),n)
        call zscal(n*m,ze,Av2(1,1),1)
        call zlacpy('F',n,m,Av2(1,1),n,Av(1,m*(m0-1)+1),n)

        !V(1:n,m*j+1:m*(j+1))=Av2(1:n,1:m)        
        !Av(1:n,m*(m0-1)+1:m*m0)=Av2(1:n,1:m)
                !Av(1:n,m*(m0-1)+1:m*m0)=workin(1:n,1:m)*(1.0d0,0.0d0)*ze
        !xwork(1:n,1:m)=dimag(V(1:n,m*(m0-1)+1:m*m0)) !if we have B matrix
        workin(1:n,1:m)=dimag(V(1:n,m*(m0-1)+1:m*m0)) !if we don't have B matrix
        stateVars(1)=41
        ijob=40
        return
    end if

    if (stateVars(1)==41) then !B*im(V(j)) for Av
        call zlacp2('F',n,m,workin(1,1),n,Av2(1,1),n)
        call zaxpy(n*m,ze*(0.0d0,1.0d0),Av2(1,1),1,Av(1,(m*(m0-1)+1)),1)
        !Av(1:n,m*(m0-1)+1:m*m0)=Av(1:n,m*(m0-1)+1:m*m0)+Av2(1:n,1:m)*(0.0d0,1.0d0)*ze!workin(1:n,1:m)*(0.0d0,1.0d0)*ze
        xwork(1:n,1:m)=dble(V(1:n,m*(m0-1)+1:m*m0))
        stateVars(1)=42
        ijob=30
        return
    end if

    if (stateVars(1)==42) then !A*re(V(j)) for Av
        call zlacp2('F',n,m,workin(1,1),n,Av2(1,1),n)
        call zaxpy(n*m,(-1.0d0,0.0d0),Av2(1,1),1,Av(1,m*(m0-1)+1),1)
        !Av(1:n,m*(m0-1)+1:m*m0)=Av(1:n,m*(m0-1)+1:m*m0)-workin(1:n,1:m)*(1.0d0,0.0d0)
        xwork(1:n,1:m)=dimag(V(1:n,m*(m0-1)+1:m*m0))
        stateVars(1)=43
        ijob=30
        return
    end if

    if (stateVars(1)==43) then !A*im(V(j)) for Av
        call zlacp2('F',n,m,workin(1,1),n,Av2(1,1),n)
        call zaxpy(n*m,(0.0d0,-1.0d0),Av2(1,1),1,Av(1,m*(m0-1)+1),1)
        !Av(1:n,m*(m0-1)+1:m*m0)=Av(1:n,m*(m0-1)+1:m*m0)-workin(1:n,1:m)*(0.0d0,1.0d0)
    end if


    do i=stateVars(2),restarts
        stateVars(2)=i

        !form right hand side from residual:
        if (stateVars(1)==1) then 

            V(1:n,1:m)=Brhs(1:n,1:m)-Ax(1:n,1:m)
        end if        

        !__________form Krylov subspace V_____________________
        
        do j=stateVars(3),m0-1
            stateVars(3)=j
            if (stateVars(1)==1 .or. stateVars(1)==4) then
                !xwork(1:n,1:m)=dble(V(1:n,m*(j-1)+1:m*j)) !do this when we have a B matrix
                workin(1:n,1:m)=dble(V(1:n,m*(j-1)+1:m*j)) !do this when we don't have B matrix
                stateVars(1)=3     
                ijob=40       
                
                return
            end if

            if(stateVars(1)==33) then
                stateVars(1)=4
            end if
        end do

        if (m0==1) then
            if (stateVars(1)==1) then
                stateVars(1)=4
            end if
        end if

        !____________form reduced system Av=A*V_______________

        if (stateVars(1)==4) then
            if (m0>1) then
                Av(1:n,1:m*(m0-1))=V(1:n,m+1:m0*m)
            end if

            !xwork(1:n,1:m)=dble(V(1:n,m*(m0-1)+1:m*m0)) !if we have B matrix
            workin(1:n,1:m)=dble(V(1:n,m*(m0-1)+1:m*m0)) ! if we don't have B matrix
            ijob=40
            
            return
        end if

        if (stateVars(1)==43) then
            !____________Solve reduced system Av*x=r______________

        
            call wallocate_1i(ipiv,m*m0,info)
            lwork=3*n
            call wallocate_1z(work,lwork,info)
            call wallocate_1z(qrtau,m0*m,info)

            call wallocate_2z(Rs,m0*m,m0*m,info)
            call wallocate_2z(Bs,m0*m,m,info)
                        
            !Av2(1:n,1:m0*m)=Av(1:n,1:m0*m)
            call zlacpy('F',n,m*m0,Av(1,1),n,Av2(1,1),n)

            !get QR factorization
            call ZGEQRF( n, m0*m, Av2, n, qrtau, work, lwork, info )
            if (info .ne. 0) then
                print *,'Problem with least squares solution in GMRES'
                print *,'ZGEQRF error info = ',info
                stop
            end if
    
            !get R matrix
            !Rs(1:m*m0,1:m*m0)=Av2(1:m*m0,1:m*m0)
            call zlacpy('F',m*m0,m*m0,Av2(1,1),n,Rs(1,1),m*m0)

            !get Q matrix
            call ZUNGQR(  n, m0*m, m0*m, Av2, n, qrtau, work, lwork, info )
            if (info .ne. 0) then
                print *,'Problem with least squares solution in GMRES'
                print *,'ZUNGQR error info = ',info
                stop
            end if
            
            !form reduced right hand side matrix:
            !use V(1:n,1:m) since V(1:n,1:m) = r = B-Ax is the right hand side
            call ZGEMM('C','N',m*m0,m,n,(1.0d0,0.0d0),Av2,n,V(1:n,1:m),n,(0.0d0,0.0d0),Bs(1,1),m0*m)
                        
            !solve upper triangular system Rs*x=Q'*Bs
            call ZTRTRS( 'U', 'N', 'N', m*m0, m, Rs, m*m0, Bs, m0*m, info )
            if (info .ne. 0) then
                print *,'Problem with least squares solution in GMRES'
                print *,'ZTRTRS error info = ',info
                stop
            end if
            
            !update Ax
            call ZGEMM('N','N',n,m,m*m0,(1.0d0,0.0d0),Av(1,1),n,Bs(1,1),m0*m,(1.0d0,0.0d0),Ax(1,1),n)
            
            !get full size solution x=V*xr
            call ZGEMM('N','N',n,m,m*m0,(1.0d0,0.0d0),V(1,1),n,Bs(1,1),m0*m,(0.0d0,0.0d0),Av(1,1),n) 
            
            !update solution:
            x(1:n,1:m)=x(1:n,1:m)+Av(1:n,1:m) !reusing Av to save some memory

            call wdeallocate_1i(ipiv)
            call wdeallocate_1z(work)
            call wdeallocate_1z(qrtau)

            call wdeallocate_2z(Rs)
            call wdeallocate_2z(Bs) 
            if  (i<restarts) then
                !xwork(1:n,1:m)=dble(x(1:n,1:m)) !if we have B matrix
                !workin(1:n,1:m)=dble(V(1:n,m*(m0-1)+1:m*m0)) !if we don't have B matrix
                stateVars(1)=1
                
            end if
            
        end if
    end do

stateVars(1)=0!-2
ijob=11
!Brhs(1:n,1:m)=x(1:n,1:m)
call zlacpy('F',n,m,x(1,1),n,Brhs(1,1),n)

end subroutine dfeast_gmres



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


subroutine zprintMat(Mat,n,rowmax,colmax)
integer :: n,colmax,rowmax,j,k
complex (kind=kind(0.0d0)), dimension(n,*) :: Mat
character (len=12) :: fmtstr

write(fmtstr,'(i0)') colmax
print *,''
do j=1,rowmax
!write(*,"(15g15.5)") (real(V(j,k)*conjg(V(j,k))), k=1,m*m0)
!write(*,'('//fmtstr//'g15.5)') (real(Mat(j,k)*conjg(Mat(j,k))), k=1,colmax) !(real(Mat(j,k)*conjg(Mat(j,k))), k=1,rowmax) 
write(*,'('//fmtstr//'g15.5)') (abs(Mat(j,k)), k=1,colmax)
end do
print *,''
end subroutine zprintMat


subroutine selection_sort(list,n) !sort list to be smallest to biggest; not good for really big lists
    double precision, dimension(:) :: list
    integer :: n

    integer :: i,j
    double precision :: temp

    do i=1,n-1

        do j=i+1,n
            if(list(j)<list(i)) then
                temp=list(i)
                list(i)=list(j)
                list(j)=temp
            end if
        end do

    end do

end subroutine selection_sort


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
    end if
end subroutine quicksort



subroutine zfeast_gmres(ijob,stateVars,Brhs,x,V,Av,Ax,ze,n,m,maxm,eps,restarts,m0,xwork,workin,Av2,times)
implicit none
    include 'f90_noruntime_interface.fi'
    
    integer :: ijob,n,m,restarts,m0,i,j,k,l
    integer, dimension(3)::stateVars
    

    !ijob: RCI case to call when this routine returns
    !n: dimension of linear sysmte
    !m: number of right hand sides
    !eps: target linear system accuracy; stop if this is reached
    !restarts: number of GMRES loops to do
    !m0: number of krylov subspace blocks; total krylov subspace size is m*m0
    !i,j,k: loop variables
    !stateVars(1) = current state of routine
    !stateVars(2) = index for outer loop
    !stateVars(3) = index for inner loop, for building krylov subspace

    complex (kind=kind(0.0d0)), dimension(n,*):: xwork,workin !workspace variables
    complex (kind=kind(0.0d0)), dimension(n,*) :: x,V,Ax,Av,Brhs,Av2
    !x: solution to linear system
    !V: Krylov subspace in which to solve linear system; has dimension m*m0
    !Ax: A*x
    !Av: A*V; storing this means we only have to do one matrix multiply instead of two
    !Brhs: right hand sides
    complex (kind=kind(0.0d0)) :: ze !complex shift for FEAST

    !!lapack stuff:
    complex (kind=kind(0.0d0)), dimension(:), pointer :: work,qrtau
    integer :: lwork,info
    integer, dimension(:),pointer :: ipiv

    !!least squares stuff:
    complex (kind=kind(0.0d0)), dimension(:,:), pointer :: Rs,Bs

    !!check residual:
    double precision :: norm 

    !lapack least squares
    complex (kind=kind(0.0d0)), dimension(:,:), allocatable :: B2
    complex (kind=kind(0.0d0)), dimension(:), allocatable :: qwork
    integer :: lwork2
    
    !BLAS:
    complex (kind=kind(0.0d0)),external :: zdotc
    double precision, external :: dznrm2

    !timing:
    double precision :: times
    integer :: c1,c2
    double precision, external :: elapsed_time
  
    !measuring norm:
    integer :: maxm !how many linear systems to look at in determining eps
    double precision :: eps
    double precision :: maxres,tempres
    double precision, dimension(1:m) :: tempreslist    
 
    integer :: debug

    debug=0

    lwork2=n*m*m0+n*m*m0
    allocate(qwork(lwork2),B2(n,m))
    
    i=stateVars(2)
    j=stateVars(3)
    
    !initialize routine when we start it
    if (stateVars(1)==-1) then
        times=0.0d0
        V(1:n,1:m*m0)=0.0d0
        Av(1:n,1:m*m0)=0.0d0    
        Ax(1:n,1:m)=0.0d0
        x(1:n,1:m)=0.0d0
        stateVars(2)=1
        stateVars(3)=1
        stateVars(1)=1
        xwork(1:n,1:m)=0.0d0
    end if
    

    if (stateVars(1)==3) then !B*V(j) for building V
        call zscal(n*m,ze,workin(1,1),1)
        call zlacpy('F',n,m,workin(1,1),n,V(1,m*j+1),n) 
        !V(1:n,m*j+1:m*(j+1))=workin(1:n,1:m)*(1.0d0,0.0d0)*ze
        !xwork(1:n,1:m)=dimag(V(1:n,m*(j-1)+1:m*j)) !if we have B matrix

        call zlacpy('F',n,m,V(1,m*(j-1)+1),n,xwork(1,1),n)
        !xwork(1:n,1:m)=V(1:n,m*(j-1)+1:m*j) !if we don't have B matrix
        stateVars(1)=33
        ijob=30
        return
    end if

    if (stateVars(1)==33) then !A*re(V(j)) for V
        call zaxpy(n*m,(-1.0d0,0.0d0),workin(1:n,1:m),1,V(1:n,m*j+1:m*(j+1)),1)
        !V(1:n,m*j+1:m*(j+1))=V(1:n,m*j+1:m*(j+1))-workin(1:n,1:m)*(1.0d0,0.0d0)
    end if


    if (stateVars(1)==4) then !B*re(V(j)) for Av
        call zscal(n*m,ze,workin(1,1),1)
        call zlacpy('F',n,m,workin(1,1),n,Av(1,m*(m0-1)+1),n)

        !V(1:n,m*j+1:m*(j+1))=Av2(1:n,1:m)        
        !Av(1:n,m*(m0-1)+1:m*m0)=Av2(1:n,1:m)
                !Av(1:n,m*(m0-1)+1:m*m0)=workin(1:n,1:m)*(1.0d0,0.0d0)*ze
        !xwork(1:n,1:m)=dimag(V(1:n,m*(m0-1)+1:m*m0)) !if we have B matrix
        
        !xwork(1:n,1:m)=V(1:n,m*(m0-1)+1:m*m0) !if we don't have B matrix
        call zlacpy('F',n,m,V(1,m*(m0-1)+1),n,xwork(1,1),n)
        stateVars(1)=43
        ijob=30
        return
    end if

    if (stateVars(1)==43) then !A*re(V(j)) for Av
        call zaxpy(n*m,(-1.0d0,0.0d0),workin(1:n,1:m),1,Av(1:n,m*(m0-1)+1:m*m0),1)
        !Av(1:n,m*(m0-1)+1:m*m0)=Av(1:n,m*(m0-1)+1:m*m0)-workin(1:n,1:m)*(1.0d0,0.0d0)
    end if

    do i=stateVars(2),restarts
        stateVars(2)=i
        
        !form right hand side from residual:
        if (stateVars(1)==1) then 
            
            if(debug) print*,'finding residual'
            V(1:n,1:m)=Brhs(1:n,1:m)-Ax(1:n,1:m)
            !print *,'resnorm1=',dznrm2(n*m, V(1:n,1:m), 1)

            !calculate individual vector residuals in order to see if we've converged:
            maxres=0.0d0
            !maxm=m
            do k=1,m
               tempreslist(k)=dznrm2(n,v(1:n,k),1)/dznrm2(n,Brhs(1:n,k),1)
            end do

            call quicksort(tempreslist,1,m)
            
            maxres=tempreslist(maxm)

            if(maxres<eps) then
                !print *,'     lin sys converged:',maxres,eps
                !print *,'     ',i,restarts
                exit
            end if

            !call zlacpy('F',n,m,Brhs(1,1),n,V(1,1),n)
            !call zaxpy(n*m,(-1.0d0,0.0d0),Ax(1,1),1,V(1,1),1)
        end if        

        
        !__________form Krylov subspace V_____________________
        
        do j=stateVars(3),m0-1

            if(debug) print*,'doing krylov',j
            
            stateVars(3)=j
            if (stateVars(1)==1 .or. stateVars(1)==4) then
                !xwork(1:n,1:m)=dble(V(1:n,m*(j-1)+1:m*j)) !do this when we have a B matrix
                
                !workin(1:n,1:m)=V(1:n,m*(j-1)+1:m*j) !do this when we don't have B matrix
                call zlacpy('F',n,m,V(1,m*(j-1)+1),n,workin(1,1),n)
                stateVars(1)=3     
                ijob=40       
                return
            end if

            if(stateVars(1)==33) then
                stateVars(1)=4
            end if
        end do

        if (m0==1) then
            if (stateVars(1)==1) then
                if(debug) print*,'setting state to 4'
                stateVars(1)=4
            end if
        end if

        !____________form reduced system Av=A*V_______________

        if (stateVars(1)==4) then
            if (m0>1) then
                !call zlacpy('F',n,m*(m0-1),V(1,m+1),n,Av(1,1),n)
                Av(1:n,1:m*(m0-1))=V(1:n,m+1:m0*m)
            end if
            if(debug) print*,'forming reduced system, beginning multiply'
            !xwork(1:n,1:m)=dble(V(1:n,m*(m0-1)+1:m*m0)) !if we have B matrix
            !workin(1:n,1:m)=V(1:n,m*(m0-1)+1:m*m0) ! if we don't have B matrix
            call zlacpy('F',n,m,V(1,m*(m0-1)+1),n,workin(1,1),n)
            ijob=40
            return
        end if

        if (stateVars(1)==43) then

            if(debug) print*,'solving reduced system'
            !____________Solve reduced system Av*x=r______________

            call system_clock(count=c1) 
            call wallocate_1i(ipiv,m*m0,info)
            lwork=3*n
            call wallocate_1z(work,lwork,info)
            call wallocate_1z(qrtau,m0*m,info)

            call wallocate_2z(Rs,m0*m,m0*m,info)
            call wallocate_2z(Bs,m0*m,m,info)
           
            !Av2(1:n,1:m0*m)=Av(1:n,1:m0*m) 
            call zlacpy('F',n,m*m0,Av(1,1),n,Av2(1,1),n)

            !goto 111
            !get QR factorization
            call ZGEQRF( n, m0*m, Av2, n, qrtau, work, lwork, info )
            if (info .ne. 0) then
                print *,'Problem with least squares solution in GMRES'
                print *,'ZGEQRF error info = ',info
                stop
            end if
    
            !get R matrix
            !Rs(1:m*m0,1:m*m0)=Av2(1:m*m0,1:m*m0)
            call zlacpy('F',m*m0,m*m0,Av2(1,1),n,Rs(1,1),m*m0) 

            !get Q matrix
            call ZUNGQR(  n, m0*m, m0*m, Av2, n, qrtau, work, lwork, info )
            if (info .ne. 0) then
                print *,'Problem with least squares solution in GMRES'
                print *,'ZUNGQR error info = ',info
                stop
            end if
            
            !form reduced right hand side matrix:
            !use V(1:n,1:m) since V(1:n,1:m) = r = B-Ax is the right hand side
            call ZGEMM('C','N',m*m0,m,n,(1.0d0,0.0d0),Av2,n,V(1:n,1:m),n,(0.0d0,0.0d0),Bs(1,1),m0*m)
         
            !solve upper triangular system Rs*x=Q'*Bs
            call ZTRTRS( 'U', 'N', 'N', m*m0, m, Rs, m*m0, Bs, m0*m, info )
            if (info .ne. 0) then
                print *,'Problem with least squares solution in GMRES'
                print *,'ZTRTRS error info = ',info
                stop
            end if
            111 continue

            !B2(1:n,1:m)=V(1:n,1:m)
            !call ZGELS('N',n,m*m0,m,Av2(1:n,1:m*m0),n,B2(1:n,1:m),n,qwork,lwork2,info)
            !if (info .ne. 0) then
            !    print *,'Problem with least squares solution in GMRES'
            !    print *,'ZGELS error info = ',info
            !    stop
            !end if 
            !Bs(1:m*m0,1:m)=B2(1:m*m0,1:m)

            !update Ax
            call ZGEMM('N','N',n,m,m*m0,(1.0d0,0.0d0),Av(1,1),n,Bs(1,1),m0*m,(1.0d0,0.0d0),Ax(1,1),n)
            
            !get full size solution x=V*xr
            call ZGEMM('N','N',n,m,m*m0,(1.0d0,0.0d0),V(1,1),n,Bs(1,1),m0*m,(0.0d0,0.0d0),Av(1,1),n) 
            
            !update solution:
            x(1:n,1:m)=x(1:n,1:m)+Av(1:n,1:m) !reusing Av to save some memory

            call wdeallocate_1i(ipiv)
            call wdeallocate_1z(work)
            call wdeallocate_1z(qrtau)

            call wdeallocate_2z(Rs)
            call wdeallocate_2z(Bs) 
            call system_clock(count=c2)
            times=times+elapsed_time(c1,c2)            

            if  (i<restarts) then
                !xwork(1:n,1:m)=dble(x(1:n,1:m)) !if we have B matrix
                !workin(1:n,1:m)=V(1:n,m*(m0-1)+1:m*m0) !if we don't have B matrix
                !call zlacpy('F',n,m,V(1,m*(m0-1)+1),n,workin(1,1),n)
                stateVars(1)=1
            end if
        end if
    end do

stateVars(1)=0!-2
ijob=11
!Brhs(1:n,1:m)=x(1:n,1:m)
call zlacpy('F',n,m,x(1,1),n,Brhs(1,1),n)
!stop
end subroutine zfeast_gmres




subroutine zfeast_gmrespre(ijob,stateVars,Brhs,x,V,Av,Ax,ze,n,m,maxm,eps,restarts,m0,xwork,workin,Av2,times)
implicit none
    include 'f90_noruntime_interface.fi'
    
    integer :: ijob,n,m,restarts,m0,i,j,k,l
    integer, dimension(3)::stateVars
    

    !ijob: RCI case to call when this routine returns
    !n: dimension of linear sysmte
    !m: number of right hand sides
    !eps: target linear system accuracy; stop if this is reached
    !restarts: number of GMRES loops to do
    !m0: number of krylov subspace blocks; total krylov subspace size is m*m0
    !i,j,k: loop variables
    !stateVars(1) = current state of routine
    !stateVars(2) = index for outer loop
    !stateVars(3) = index for inner loop, for building krylov subspace

    complex (kind=kind(0.0d0)), dimension(n,*):: xwork,workin !workspace variables
    complex (kind=kind(0.0d0)), dimension(n,*) :: x,V,Ax,Av,Brhs,Av2
    !x: solution to linear system
    !V: Krylov subspace in which to solve linear system; has dimension m*m0
    !Ax: A*x
    !Av: A*V; storing this means we only have to do one matrix multiply instead of two
    !Brhs: right hand sides
    complex (kind=kind(0.0d0)) :: ze !complex shift for FEAST

    !!lapack stuff:
    complex (kind=kind(0.0d0)), dimension(:), pointer :: work,qrtau
    integer :: lwork,info
    integer, dimension(:),pointer :: ipiv

    !!least squares stuff:
    complex (kind=kind(0.0d0)), dimension(:,:), pointer :: Rs,Bs

    !!check residual:
    double precision :: norm 

    !lapack least squares
    complex (kind=kind(0.0d0)), dimension(:,:), allocatable :: B2
    complex (kind=kind(0.0d0)), dimension(:), allocatable :: qwork
    integer :: lwork2
    
    !BLAS:
    complex (kind=kind(0.0d0)),external :: zdotc
    double precision, external :: dznrm2

    !timing:
    double precision :: times
    integer :: c1,c2
    double precision, external :: elapsed_time
  
    !measuring norm:
    integer :: maxm !how many linear systems to look at in determining eps
    double precision :: eps
    double precision :: maxres,tempres
    double precision, dimension(1:m) :: tempreslist    
 
    integer :: debug

    debug=0

    lwork2=n*m*m0+n*m*m0
    allocate(qwork(lwork2),B2(n,m))
    
    i=stateVars(2)
    j=stateVars(3)
    
    !initialize routine when we start it
    if (stateVars(1)==-1) then
        times=0.0d0
        V(1:n,1:m*m0)=0.0d0
        Av(1:n,1:m*m0)=0.0d0    
        Ax(1:n,1:m)=0.0d0
        x(1:n,1:m)=0.0d0
        stateVars(2)=1
        stateVars(3)=1
        stateVars(1)=1
        xwork(1:n,1:m)=0.0d0
    end if
    

    if (stateVars(1)==3) then !B*V(j) for building V
        call zscal(n*m,ze,workin(1,1),1)
        call zlacpy('F',n,m,workin(1,1),n,V(1,m*j+1),n) 
        !V(1:n,m*j+1:m*(j+1))=workin(1:n,1:m)*(1.0d0,0.0d0)*ze
        !xwork(1:n,1:m)=dimag(V(1:n,m*(j-1)+1:m*j)) !if we have B matrix

        call zlacpy('F',n,m,V(1,m*(j-1)+1),n,xwork(1,1),n)
        !xwork(1:n,1:m)=V(1:n,m*(j-1)+1:m*j) !if we don't have B matrix
        stateVars(1)=33
        ijob=30
        return
    end if

    if (stateVars(1)==33) then !A*re(V(j)) for V
        call zaxpy(n*m,(-1.0d0,0.0d0),workin(1:n,1:m),1,V(1:n,m*j+1:m*(j+1)),1)
        !V(1:n,m*j+1:m*(j+1))=V(1:n,m*j+1:m*(j+1))-workin(1:n,1:m)*(1.0d0,0.0d0)
    end if


    if (stateVars(1)==4) then !B*re(V(j)) for Av
        call zscal(n*m,ze,workin(1,1),1)
        call zlacpy('F',n,m,workin(1,1),n,Av(1,m*(m0-1)+1),n)

        !V(1:n,m*j+1:m*(j+1))=Av2(1:n,1:m)        
        !Av(1:n,m*(m0-1)+1:m*m0)=Av2(1:n,1:m)
                !Av(1:n,m*(m0-1)+1:m*m0)=workin(1:n,1:m)*(1.0d0,0.0d0)*ze
        !xwork(1:n,1:m)=dimag(V(1:n,m*(m0-1)+1:m*m0)) !if we have B matrix
        
        !xwork(1:n,1:m)=V(1:n,m*(m0-1)+1:m*m0) !if we don't have B matrix
        call zlacpy('F',n,m,V(1,m*(m0-1)+1),n,xwork(1,1),n)
        stateVars(1)=43
        ijob=30
        return
    end if

    if (stateVars(1)==43) then !A*re(V(j)) for Av
        call zaxpy(n*m,(-1.0d0,0.0d0),workin(1:n,1:m),1,Av(1:n,m*(m0-1)+1:m*m0),1)
        !Av(1:n,m*(m0-1)+1:m*m0)=Av(1:n,m*(m0-1)+1:m*m0)-workin(1:n,1:m)*(1.0d0,0.0d0)
    end if

    do i=stateVars(2),restarts
        stateVars(2)=i
        
        !form right hand side from residual:
        if (stateVars(1)==1) then 
            
            if(debug) print*,'finding residual'
            V(1:n,1:m)=Brhs(1:n,1:m)-Ax(1:n,1:m)
            !print *,'resnorm1=',dznrm2(n*m, V(1:n,1:m), 1)

            !calculate individual vector residuals in order to see if we've converged:
            maxres=0.0d0
            !maxm=m
            do k=1,m
               tempreslist(k)=dznrm2(n,v(1:n,k),1)/dznrm2(n,Brhs(1:n,k),1)
            end do

            call quicksort(tempreslist,1,m)
            
            maxres=tempreslist(maxm)

            if(maxres<eps) then
                !print *,'     lin sys converged:',maxres,eps
                !print *,'     ',i,restarts
                exit
            end if

            !call zlacpy('F',n,m,Brhs(1,1),n,V(1,1),n)
            !call zaxpy(n*m,(-1.0d0,0.0d0),Ax(1,1),1,V(1,1),1)
        end if        

        
        !__________form Krylov subspace V_____________________
        
        do j=stateVars(3),m0-1

            if(debug) print*,'doing krylov',j
            
            stateVars(3)=j
            if (stateVars(1)==1 .or. stateVars(1)==4) then
                !xwork(1:n,1:m)=dble(V(1:n,m*(j-1)+1:m*j)) !do this when we have a B matrix
                call zlacpy('F',n,m,V(1,m*(j-1)+1),n,xwork(1,1),n)

                !workin(1:n,1:m)=V(1:n,m*(j-1)+1:m*j) !do this when we don't have B matrix
                !call zlacpy('F',n,m,V(1,m*(j-1)+1),n,workin(1,1),n)
                
                stateVars(1)=3     
                ijob=40       
                return
            end if

            if(stateVars(1)==33) then
                stateVars(1)=4
            end if
        end do

        if (m0==1) then
            if (stateVars(1)==1) then
                if(debug) print*,'setting state to 4'
                stateVars(1)=4
            end if
        end if

        !____________form reduced system Av=A*V_______________

        if (stateVars(1)==4) then
            if (m0>1) then
                !call zlacpy('F',n,m*(m0-1),V(1,m+1),n,Av(1,1),n)
                Av(1:n,1:m*(m0-1))=V(1:n,m+1:m0*m)
            end if
            if(debug) print*,'forming reduced system, beginning multiply'
            
            !xwork(1:n,1:m)=dble(V(1:n,m*(m0-1)+1:m*m0)) !if we have B matrix
            call zlacpy('F',n,m,V(1,m*(m0-1)+1),n,xwork(1,1),n)
            
            !workin(1:n,1:m)=V(1:n,m*(m0-1)+1:m*m0) ! if we don't have B matrix
            !call zlacpy('F',n,m,V(1,m*(m0-1)+1),n,workin(1,1),n)
            
            ijob=40
            return
        end if

        if (stateVars(1)==43) then

            if(debug) print*,'solving reduced system'
            !____________Solve reduced system Av*x=r______________

            call system_clock(count=c1) 
            call wallocate_1i(ipiv,m*m0,info)
            lwork=3*n
            call wallocate_1z(work,lwork,info)
            call wallocate_1z(qrtau,m0*m,info)

            call wallocate_2z(Rs,m0*m,m0*m,info)
            call wallocate_2z(Bs,m0*m,m,info)
           
            !Av2(1:n,1:m0*m)=Av(1:n,1:m0*m) 
            call zlacpy('F',n,m*m0,Av(1,1),n,Av2(1,1),n)

            !goto 111
            !get QR factorization
            call ZGEQRF( n, m0*m, Av2, n, qrtau, work, lwork, info )
            if (info .ne. 0) then
                print *,'Problem with least squares solution in GMRES'
                print *,'ZGEQRF error info = ',info
                stop
            end if
    
            !get R matrix
            !Rs(1:m*m0,1:m*m0)=Av2(1:m*m0,1:m*m0)
            call zlacpy('F',m*m0,m*m0,Av2(1,1),n,Rs(1,1),m*m0) 

            !get Q matrix
            call ZUNGQR(  n, m0*m, m0*m, Av2, n, qrtau, work, lwork, info )
            if (info .ne. 0) then
                print *,'Problem with least squares solution in GMRES'
                print *,'ZUNGQR error info = ',info
                stop
            end if
            
            !form reduced right hand side matrix:
            !use V(1:n,1:m) since V(1:n,1:m) = r = B-Ax is the right hand side
            call ZGEMM('C','N',m*m0,m,n,(1.0d0,0.0d0),Av2,n,V(1:n,1:m),n,(0.0d0,0.0d0),Bs(1,1),m0*m)
         
            !solve upper triangular system Rs*x=Q'*Bs
            call ZTRTRS( 'U', 'N', 'N', m*m0, m, Rs, m*m0, Bs, m0*m, info )
            if (info .ne. 0) then
                print *,'Problem with least squares solution in GMRES'
                print *,'ZTRTRS error info = ',info
                stop
            end if
            111 continue

            !B2(1:n,1:m)=V(1:n,1:m)
            !call ZGELS('N',n,m*m0,m,Av2(1:n,1:m*m0),n,B2(1:n,1:m),n,qwork,lwork2,info)
            !if (info .ne. 0) then
            !    print *,'Problem with least squares solution in GMRES'
            !    print *,'ZGELS error info = ',info
            !    stop
            !end if 
            !Bs(1:m*m0,1:m)=B2(1:m*m0,1:m)

            !update Ax
            call ZGEMM('N','N',n,m,m*m0,(1.0d0,0.0d0),Av(1,1),n,Bs(1,1),m0*m,(1.0d0,0.0d0),Ax(1,1),n)
            
            !get full size solution x=V*xr
            call ZGEMM('N','N',n,m,m*m0,(1.0d0,0.0d0),V(1,1),n,Bs(1,1),m0*m,(0.0d0,0.0d0),Av(1,1),n) 
            
            !update solution:
            x(1:n,1:m)=x(1:n,1:m)+Av(1:n,1:m) !reusing Av to save some memory

            call wdeallocate_1i(ipiv)
            call wdeallocate_1z(work)
            call wdeallocate_1z(qrtau)

            call wdeallocate_2z(Rs)
            call wdeallocate_2z(Bs) 
            call system_clock(count=c2)
            times=times+elapsed_time(c1,c2)            

            if  (i<restarts) then
                !xwork(1:n,1:m)=dble(x(1:n,1:m)) !if we have B matrix
                !workin(1:n,1:m)=V(1:n,m*(m0-1)+1:m*m0) !if we don't have B matrix
                !call zlacpy('F',n,m,V(1,m*(m0-1)+1),n,workin(1,1),n)
                stateVars(1)=1
            end if
        end if
    end do

stateVars(1)=0!-2
ijob=11
!Brhs(1:n,1:m)=x(1:n,1:m)
call zlacpy('F',n,m,x(1,1),n,Brhs(1,1),n)
!stop
end subroutine zfeast_gmrespre
