


subroutine refineSolveDense(transA,A,B,n,m,r,m0,usex)
subroutine linItSolve(ijob,B,n,m,r,m0)
implicit none
    include 'f90_noruntime_interface.fi'

    integer :: m0,r,n,m,lda,usex

    character (len=1):: transA
    complex (kind=kind(0.0d0)), dimension(n,*) :: A,B

    integer :: i,infoloc
    complex (kind=kind(0.0d0)), dimension(:,:),pointer :: Bc,rhs,X

    !print *,'entering refinesolve',n,m
    call wallocate_2z(Bc,n,m,infoloc)
    !print *,'allocate again',infoloc
    call wallocate_2z(X,n,m,infoloc)
    call wallocate_2z(rhs,n,m,infoloc)
   !print *,'x='
    !print *,"n,m,r=",n,m,r
    if(usex==1) then
        X(1:n,1:m)=B(1:n,1:m)
    else
        X(1:n,1:m)=(0.0d0,0.0d0)
    end if
    !print *,'bc='
    Bc(1:n,1:m)=(0.0d0,0.0d0)!B(1:n,1:m)

    if (transA=='C') then
        B(1:n,1:m)=conjg(B(1:n,1:m))
    end if

    do i=1,r
        !if (r>1)

            !Bc=B-matmul(A,X)
        !print *,'doing matvec',m
        
        rhs(1:n,1:m)=B(1:n,1:m)
        !print *,'doing zgemm'
        call zgemm('N','N',n,m,n,(1.0d0,0.0d0),A(1,1),n,X(1,1),n,(0.0d0,0.0d0),rhs,n)
        !call zcsrmm(uplo,'N',n,n,m,(1.0d0,0.0d0),sa,isa,jsa,X(1,1),(0.0d0,0.0d0),rhs(1,1))
        !print *,'doing subtraction'
        Bc(1:n,1:m)=B(1:n,1:m)-rhs(1:n,1:m)
        !print *,"first col res",sum(Bc(1:n,1))
        !end if

        !call ssSolve(Atmp,Bc,n,m,m0) 
        !print *,'entering ssSolveDense'
        call ssSolveDense('N',A,Bc,n,m,m0)
        !print *,'done ssSolveDense'
        !rc=Bc-matmul(A,Xc)
        X(1:n,1:m)=X(1:n,1:m)+Bc(1:n,1:m)
    end do

    B(1:n,1:m)=X(1:n,1:m)

    !deallocate(Bc,X)
    call wdeallocate_2z(Bc)
    call wdeallocate_2z(X)
    call wdeallocate_2z(rhs)
end subroutine refineSolveDense

subroutine ssSolveDense(transA,A,B,n,m,m0)
implicit none
    include 'f90_noruntime_interface.fi'

    
    character (len=1) :: transA
    integer :: m0,n,m
    complex (kind=kind(0.0d0)), dimension(n,*), intent(in) :: A
    complex (kind=kind(0.0d0)), dimension(n,*), intent(inout) :: B

    integer :: i,k,j
    complex (kind=kind(0.0d0)), dimension(:,:), pointer :: Q,As
    complex (kind=kind(0.0d0)), dimension(:,:), pointer :: Rs
    complex (kind=kind(0.0d0)), dimension(:,:),pointer :: Bs

    complex (kind=kind(0.0d0)), dimension(:), pointer :: work,qrtau
    complex (kind=kind(0.0d0))::const
    integer :: lwork,info
    integer, dimension(:),pointer :: ipiv

    call wallocate_2z(As,n,m*m0,info)
    call wallocate_2z(Q,n,m*m0,info)
    call wallocate_2z(Rs,m0*m,m0*m,info)
    call wallocate_2z(Bs,m0*m,m,info)

    call wallocate_1i(ipiv,m*m0,info)
    do i=1,m0*m
        ipiv(i)=i
    end do

    lwork=3*n
    !lwork=m*n+m*n
    call wallocate_1z(work,lwork,info)
    call wallocate_1z(qrtau,m0*m,info)

    Q(1:n,1:m)=B(1:n,1:m)

    !call zlacpy('F',n,m,B,n,Q,n)

    do k=1,m0-1
        call zgemm(transA,'N',n,m,n,(1.0d0,0.0d0),A(1,1),n,Q(1,m*(k-1)+1),n,(0.0d0,0.0d0),Q(1,m*k+1),n)
        !call zcsrmm(uplo,'N',n,n,m,(1.0d0,0.0d0),sa,isa,jsa,Q(1,m*(k-1)+1),(0.0d0,0.0d0),Q(1,m*k+1))
        !Q(:,m*k+1:m*(k+1))=matmul(A,Q(:,m*(k-1)+1:m*k)) 
    end do

!least squares problem:

    call zgemm(transA,'N',n,m*m0,n,(1.0d0,0.0d0),A(1,1),n,Q(1,1),n,(0.0d0,0.0d0),As(1,1),n)
    !call zcsrmm(uplo,'N',n,n,m*m0,(1.0d0,0.0d0),sa,isa,jsa,Q(1,1),(0.0d0,0.0d0),As(1,1))
    !As=matmul(A,Q)

    !QR factorization
    !call zgels('N',   n,m0*m,  m, As, n,  B,n, work, lwork, info ) 
    !call ZGEQRF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
    call ZGEQRF( n, m0*m, As, n, qrtau, work, lwork, info )
    if (info .ne. 0) then
        print *,'Problem with least squares solution in ssSolve'
        print *,'ZGEQRF error info = ',info
    end if

    !get R matrix:
    Rs(1:m0*m,1:m0*m)=(0.0d0,0.0d0)
    do i=1,m0*m
        do j=1,i
            Rs(j,i)=As(j,i)
        end do
    end do

    !get Q matrix:
    !call ZUNGQR( M, N, K,       A, LDA, TAU, WORK, LWORK, INFO )
    call ZUNGQR(  n, m0*m, m0*m, As, n, qrtau, work, lwork, info )
    if (info .ne. 0) then
        print *,'Problem with least squares solution in ssSolve'
        print *,'ZUNGQR error info = ',info
    end if

    !Bs(1:m*m0,1:m)=matmul(transpose(conjg(As(1:n,1:m*m0))),B(1:n,1:m))
    call ZGEMM('C','N',m*m0,m,n,(1.0d0,0.0d0),As,n,B(1,1),n,(0.0d0,0.0d0),Bs(1,1),m0*m)


    !solve reduced system
    !call ZGETRS( TRANS, N,    NRHS, A,  LDA, IPIV, B, LDB, INFO )
    !call ZGETRS( 'N',    m0*m,   m,  Rs, m0*m, ipiv, Bs, m0*m, info )
    !call ZTRTRS( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDB, INFO )
    call ZTRTRS( 'U', 'N', 'N', m*m0, m, Rs, m*m0, Bs, m0*m, info )
    !call zgesv(m0*m,m,Rs,m0*m,ipiv,Bs,m0*m,info)

if (info .ne. 0) then
        print *,'Problem with least squares solution in ssSolve'
        print *,'ZTRTRS error info = ',info
    end if


    !B(1:n,1:m)=matmul(Q(1:n,1:m*m0),Bs(1:m*m0,1:m))

    call ZGEMM('N','N',n,m,m*m0,(1.0d0,0.0d0),Q(1,1),n,Bs(1,1),m0*m,(0.0d0,0.0d0),B(1,1),n)

    !deallocate(As,Q,work)
    !deallocate(qrtau)
    call wdeallocate_2z(As)
    call wdeallocate_2z(Q)
    call wdeallocate_2z(Rs)
    call wdeallocate_2z(Bs)
    call wdeallocate_1i(ipiv)
    call wdeallocate_1z(work)
    call wdeallocate_1z(qrtau)
end subroutine ssSolveDense

