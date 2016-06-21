
subroutine linItSolve(ijob,x,V,Av,Ax,ze,n,m,restarts,m0,stateVars,xwork)
implicit none
    include 'f90_noruntime_interface.fi'
    
    integer :: ijob,n,m,restarts,m0,i,j
    integer, dimension(3)::stateVars

    !stateVars(1) = current state of routine
    !stateVars(2) = index for outer loop
    !stateVars(3) = index for inner loop, for building krylov subspace

    double precision, dimension(n,*) :: xwork
    complex (kind=kind(0.0d0)), dimension(n,*) :: x,V,Ax,Av
    complex (kind=kind(0.0d0)) :: ze

    !!lapack stuff:
    complex (kind=kind(0.0d0)), dimension(:), pointer :: work,qrtau
    integer :: lwork,info
    integer, dimension(:),pointer :: ipiv

    !!least squares stuff:
    complex (kind=kind(0.0d0)), dimension(:,:), pointer :: Rs,Bs
 
    i=stateVars(2)
    j=stateVars(3)
 
    if (stateVars(1)==-1) then
        Ax=0.0d0
        x=0.0d0
        stateVars(2)=1
        stateVars(3)=1
        stateVars(1)=1
        xwork=0.0d0
    end if
    
    if (stateVars(1)==3) then !B*re(V(j)) for V
        V(1:n,m*j+1:m*(j+1))=xwork(1:n,1:m)*(1.0d0,0.0d0)*ze
        xwork=aimag(V(1:n,m*(j-1)+1:m*j))
        stateVars(1)=31        
        ijob=40
        return
    end if

    if (stateVars(1)==31) then !B*im(V(j)) for V
        V(1:n,m*j+1:m*(j+1))=V(1:n,m*j+1:m*(j+1))+xwork(1:n,1:m)*(0.0d0,1.0d0)*ze
        xwork=real(V(1:n,m*(j-1)+1:m*j))
        stateVars(1)=32        
        ijob=30
        return
    end if

    if (stateVars(1)==32) then !A*re(V(j)) for V
        V(1:n,m*j+1:m*(j+1))=V(1:n,m*j+1:m*(j+1))-xwork(1:n,1:m)*(1.0d0,0.0d0)
        xwork=aimag(V(1:n,m*(j-1)+1:m*j))
        stateVars(1)=33
        ijob=30
        return
    end if

    if (stateVars(1)==33) then !A*re(V(j)) for V
        V(1:n,m*j+1:m*(j+1))=V(1:n,m*j+1:m*(j+1))-xwork(1:n,1:m)*(0.0d0,1.0d0)
    end if


    if (stateVars(1)==4) then !B*re(V(j)) for Av
        Av(1:n,m*(m0-1)+1:m*m0)=xwork(1:n,1:m)*(1.0d0,0.0d0)*ze
        xwork=aimag(V(1:n,m*(m0-1):m*m0))
        stateVars(1)=41
        ijob=40
        return
    end if

    if (stateVars(1)==41) then !B*im(V(j)) for Av
        Av(1:n,m*(m0-1)+1:m*m0)=Av(1:n,m*(m0-1)+1:m*m0)+xwork(1:n,1:m)*(0.0d0,1.0d0)*ze
        xwork=real(V(1:n,m*(m0-1):m*m0))
        stateVars(1)=42
        ijob=30
        return
    end if

    if (stateVars(1)==42) then !A*re(V(j)) for Av
        Av(1:n,m*(m0-1)+1:m*m0)=Av(1:n,m*(m0-1)+1:m*m0)-xwork(1:n,1:m)*(1.0d0,0.0d0)
        xwork=aimag(V(1:n,m*(m0-1):m*m0))
        stateVars(1)=43
        ijob=30
        return
    end if

    if (stateVars(1)==43) then !A*im(V(j)) for Av
        Av(1:n,m*(m0-1)+1:m*m0)=Av(1:n,m*(m0-1)+1:m*m0)-xwork(1:n,1:m)*(0.0d0,1.0d0)
    end if


    do i=stateVars(2),restarts
        stateVars(2)=i
 
        if (stateVars(1)==1) then 
            V(1:n,1:m)=B(1:n,1:m)-Ax(1:n,1:m)
        end if        

        !__________form Krylov subspace V_____________________
        
        do j=stateVars(3),m0-1
            stateVars(3)=j
            if (stateVars(1)==1 .or. stateVars(1)==4) then
                xwork=real(V(1:n,m*(j-1)+1:m*j))
                stateVars(1)==3     
                ijob=40       
                return
            end if

            if(stateVars(1)==33) then
                stateVars(1)=4
            end if
        end do

        !____________form reduced system Av=A*V_______________

        if (stateVars(1)==4) then
            Av(1:n,1:m*(m0-1))=V(1:n,m+1:m0*m)

            !stateVars(1)=4
            xwork=real(V(1:n,m*(m0-1):m*m0))
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
           
            !get QR factorization
            call ZGEQRF( n, m0*m, Av, n, qrtau, work, lwork, info )
            if (info .ne. 0) then
                print *,'Problem with least squares solution in linItSolve'
                print *,'ZGEQRF error info = ',info
            end if
    
            !get R matrix
            Rs(1:m*m0,1:m*m0)=Av(1:m*m0,1:m*m0)

            !get Q matrix
            call ZUNGQR(  n, m0*m, m0*m, As, n, qrtau, work, lwork, info )
            if (info .ne. 0) then
                print *,'Problem with least squares solution in linItSolve'
                print *,'ZUNGQR error info = ',info
            end if
    
            !form reduced right hand side matrix:
            !use V(1:n,1:m) since V(1:n,1:m) = r = B-Ax is the right hand side
            call ZGEMM('C','N',m*m0,m,n,(1.0d0,0.0d0),Av,n,V(1:n,1:m),n,(0.0d0,0.0d0),Bs(1,1),m0*m)

            !solve upper triangular system Rs*x=Q'*Bs
            call ZTRTRS( 'U', 'N', 'N', m*m0, m, Rs, m*m0, Bs, m0*m, info )
            if (info .ne. 0) then
                print *,'Problem with least squares solution in ssSolve'
                print *,'ZTRTRS error info = ',info
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
                xwork(1:n,1:m)=real(x(1:n,1:m))
                stateVars(1)=1
            end if
        end if
    end do

stateVars(1)=-2


end subroutine linItSolve

