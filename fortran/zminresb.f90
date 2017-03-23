subroutine zminresBlock(UPLO,n,m,dsa,isa,jsa,ze,b,x,tol,maxit,linloops,blockstart)

implicit none

integer :: i,j,k

character :: UPLO
integer :: n,m,maxit,linloops,blockstart
integer, dimension(*) :: isa,jsa
complex (kind=kind(0.0d0)) :: ze !complex shift
complex (kind=kind(0.0d0)), dimension(*) :: dsa !unshifted matrix
complex (kind=kind(0.0d0)), dimension(n,m) :: b,x
double precision :: tol

!!!minres stuff
double precision :: normr,normb
complex (kind=kind(0.0d0)), dimension(n,m) :: templ1
complex (kind=kind(0.0d0)), dimension(m,n) :: templ1t
complex (kind=kind(0.0d0)), dimension(n,m) :: tempm1

complex (kind=kind(0.0d0)), dimension(n,m) :: r,Qold,Q,Qnew,P,phi,phiold,phioldold
complex (kind=kind(0.0d0)), dimension(m,m) :: tau,tautilde,Mm,HR,T,Tnew,Rtilde,Rold,Roldold
complex (kind=kind(0.0d0)), dimension(m,m) :: ag,bg,cg,dg,agold,bgold,cgold,dgold,dgoldold,bgoldold
complex (kind=kind(0.0d0)), dimension(2*m,2*m) :: QM,HM,RT 
!!!lapack and blas stuff
complex (kind=kind(0.0d0)), dimension(m) :: qrtau1
complex (kind=kind(0.0d0)), dimension(2*m) :: qrtau2
integer :: lwork,lwork2,info
integer, dimension(m) :: ipiv
complex (kind=kind(0.0d0)), dimension(:),allocatable :: work
double precision, dimension(n) :: dwork
double precision, external :: zlange
character, dimension(6) :: matdescra

    normb=zlange('F',n,m,b,n,dwork)

    lwork=-1
    call zgeqrf(n,m,r,n,qrtau1,work,lwork,info)
    allocate(work(lwork))

    x=(0.0d0,0.0d0)
    r=b !residual, start with x=0.0d0


    !!QR of r!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !(Q,T)=qr(r)
    !tautilde=T
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Q=r
    call zgeqrf(n,m,Q,n,qrtau1,work,lwork,info)
    if(info .ne. 0) then
        print *,'block minres zgeqrf 1 error ',info
        stop
    end if

    T=(0.0d0,0.0d0)
    do i=1,m
        do k=1,i
            T(k,i)=Q(k,i)
        end do
    end do
    tautilde=T

    call zungqr(n,m,m,Q,qrtau1,work,lwork,info)
    if(info .ne. 0) then
        print *,'block minres zungqr 1 error ',info
    end if
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1


    Qold=(0.0d0,0.0d0)
    phioldold=(0.0d0,0.0d0)
    phiold=(0.0d0,0.0d0)
    phi=(0.0d0,0.0d0)

    ag=(0.0d0,0.0d0)
    bg=ag
    dg=ag
    cg=ag
    bgold=ag
    bgoldold=ag
    cgold=ag


    agold=ag
    do i=1,m
        agold(i,i)=(1.0d0,0.0d0)
    end do

    dgold=agold
    dgoldold=agold


    !set up matrix multiplication stuff
    if(UPLO=='F') then
        matdescra(1)='G'
    else
        matdescra(1)='H'
    end if
    !matdescra(1)='G'
    matdescra(2)=UPLO
    matdescra(3)='N'
    matdescra(4)='F'


    do k=1,maxit

        !!!!!!M=Q'A*Q!!!!!!!!!!!!!!!!!!!! 
        !P=A*Q
        call mkl_zcsrmm('N', n, m, n, (-1.0d0,0.0d0), matdescra, dsa, jsa, isa, isa(2), Q, n, (1.0d0,0.0d0)*dble(ze), P, n)
        !M=Q'*P
        call zgemm('C','N',m,m,n,(1.0d0,0.0d0),Q,n,P,n,(0.0d0,0.0d0),Mm,m)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        !!!!!!P=A*Q-Q*M-Qold*T'!!!!!!!!!!       
        !did P=A*Q above
        !P=P-Q*M
        call zgemm('N','N',n,m,m,(-1.0d0,0.0d0),Q,n,Mm,m,(1.0d0,0.0d0),P,n)
        !P=P-Qold*T'
        call zgemm('N','C',n,m,m,(-1.0d0,0.0d0),Qold,n,T,m,(1.0d0,0.0d0),P,n)  
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



        !!!!!(Qnew,Tnew)=qr(P)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        Qnew=P
        call zgeqrf(n,m,Qnew,n,qrtau1,work,lwork,info)
        if(info .ne. 0) then
            print *,'block minres zgeqrf 2 error ',info
            stop
        end if
    
        Tnew=(0.0d0,0.0d0)
        do i=1,m
            do j=1,i
                Tnew(j,i)=Qnew(j,i)
            end do
        end do
   
        call zungqr(n,m,m,Qnew,qrtau1,work,lwork,info)
        if(info .ne. 0) then
            print *,'block minres zungqr 2 error ',info
        end if
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !Roldold=bgoldold*T'
        call zgemm('N','C',m,m,m,(1.0d0,0.0d0),bgoldold,m,T,m,(0.0d0,0.0d0),Roldold,m)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
        !!!!!Rold=agold*dgoldold*T'+bgold*M!!!!!!!!!!!!!!!!!!
        !tempm1=dgoldold*T'
        call zgemm('N','C',m,m,m,(1.0d0,0.0d0),dgoldold,m,T,m,(0.0d0,0.0d0),tempm1,m)
        !Rold=agold*tempm1
        call zgemm('N','N',m,m,m,(1.0d0,0.0d0),agold,m,tempm1,m,(0.0d0,0.0d0),Rold,m)
        !Rold=Rold+bgold*M
        call zgemm('N','N',m,m,m,(1.0d0,0.0d0),bgold,m,Mm,m,(1.0d0,0.0d0),Rold,m)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !Rtilde=cgold*dgoldold*T'+dgold*M
        !Already did tempm1=dgoldold*T' above
        !call zgemm('N','C',m,m,m,(1.0d0,0.0d0),dgoldold,m,T,m,(0.0d0,0.0d0),tempm1,m)
        
        !Rtilde=cgold*tempm1
        call zgemm('N','N',m,m,m,(1.0d0,0.0d0),cgold,m,tempm1,m,(0.0d0,0.0d0),Rtilde,m)
        !Rtilde=Rtilde+dgold*M
        call zgemm('N','N',m,m,m,(1.0d0,0.0d0),dgold,m,Mm,m,(1.0d0,0.0d0),Rtilde,m)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        RT=(0.0d0,0.0d0)
        RT(1:m,1:m)=Rtilde
        RT(m+1:2*m,1:m)=Tnew
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !(QM,HR)=QR(RT)
        !HM=QM'
        !Rkk=HR[1:m,1:m]
        QM=(0.0d0,0.0d0)
        QM(1:2*m,1:m)=RT
        call zgeqrf(2*m,2*m,QM,2*m,qrtau2,work,lwork,info)
        if(info .ne. 0) then
            print *,'block minres zgeqrf 3 error ',info
            stop
        end if

        !form Rkk'
        HR=(0.0d0,0.0d0)
        do i=1,m
            do j=1,i
                HR(i,j)=conjg(QM(j,i))
            end do
        end do

        !call zungqr(2*m,2*m,2*m,Q,qrtau2,work,lwork,info)
        !if(info .ne. 0) then
        !    print *,'block minres zungqr 3 error ',info
        !end if

        HM=(0.0d0,0.0d0)
        do i=1,2*m
            HM(i,i)=(1.0d0,0.0d0)
        end do
        call zunmqr('L','C',2*m,2*m,2*m,QM,2*m,qrtau2,HM,2*m,work,lwork,info)
        if(info .ne. 0) then
            print *,'block minres zunmqr error ',info
            stop
        end if
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1


        ag=HM(1:m,1:m)
        cg=HM(m+1:2*m,1:m)
        bg=HM(1:m,m+1:2*m)
        dg=HM(m+1:2*m,m+1:2*m)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !phi=(Q-phiold*Rold-phioldold*Roldold)*inv(HR[1:m,1:m])
        templ1t=conjg(transpose(Q))
        !templ1t=templ1-Rold'*phiold'
        call zgemm('C','C',m,n,m,(1.0d0,0.0d0),Rold,m,phiold,n,(-1.0d0,0.0d0),templ1t,m)
        !templ1t=templ1t-Roldold'*phioldold'
        call zgemm('C','C',m,n,m,(1.0d0,0.0d0),Roldold,m,phioldold,n,(-1.0d0,0.0d0),templ1t,m)
        !solve system Rkk'*phi'=templ1t
        call zgesv(m,n,HR,m,ipiv,templ1t,m,info)
        if(info .ne. 0) then
            print *,'block minres zgesv error ',info
            stop
        end if
        phi=conjg(transpose(templ1t))
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !tau=ag*tautilde
        call zgemm('N','N',m,m,m,(1.0d0,0.0d0),ag,m,tautilde,m,(0.0d0,0.0d0),tau,m)

        !x=x+phi*tau
        call zgemm('N','N',n,m,m,(1.0d0,0.0d0),phi,n,tau,m,(1.0d0,0.0d0),x,n)

        !Check residual!

        Qold=Q
        Q=Qnew
        T=Tnew

        !tautilde=cg*tautilde
        tempm1=tautilde
        call zgemm('N','N',m,m,m,(1.0d0,0.0d0),cg,m,tempm1,m,(0.0d0,0.0d0),tautilde,m)

        !norm(r)=norm(tautilde)
        normr=zlange('F',m,m,tautilde,m,dwork)  
        print *,'Minres ',k,normr/normb
        if(normr/normb<tol) exit

        phioldold=phiold
        phiold=phi

        agold=ag
        dgoldold=dgold
        dgold=dg
        bgoldold=bgold
        bgold=bg
        cgold=cg

    end do


end subroutine zminresBlock


