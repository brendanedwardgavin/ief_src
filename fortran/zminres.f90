subroutine zminresBlock(UPLO,n,m,dsa,isa,jsa,ze,b,x,tol,maxit,linloops,blockstart)
use rundata
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
complex (kind=kind(0.0d0)), dimension(:,:),allocatable :: templ1,templ2 !(n,m)
complex (kind=kind(0.0d0)), dimension(:,:),allocatable :: templ1t !(m,n)
complex (kind=kind(0.0d0)), dimension(:,:),allocatable :: tempm1 !(m,m)

complex (kind=kind(0.0d0)), dimension(:,:),allocatable :: r,Qold,Q,Qnew,P,phi,phiold,phioldold !(n,m)
complex (kind=kind(0.0d0)), dimension(:,:),allocatable :: tau,tautilde,Mm,HR,T,Tnew,Rtilde,Rold,Roldold !(m,m)

complex (kind=kind(0.0d0)), dimension(:,:),allocatable :: ag,bg,cg,dg,agold,bgold,cgold,dgold,dgoldold,bgoldold !(m,m)
complex (kind=kind(0.0d0)), dimension(:,:),allocatable :: QM,HM,RT !(2*m,2*m)
!!!lapack and blas stuff
complex (kind=kind(0.0d0)), dimension(:),allocatable :: qrtau1 !(m)
complex (kind=kind(0.0d0)), dimension(:),allocatable :: qrtau2 !(2*m)
integer :: lwork,info
integer, dimension(:),allocatable :: ipiv !(m)
complex (kind=kind(0.0d0)), dimension(:),allocatable :: work !(3*n)
double precision, dimension(:),allocatable :: dwork !(n)
double precision, external :: zlange,dznrm2
character, dimension(6) :: matdescra

    allocate(templ1(n,m),templ2(n,m),templ1t(m,n))
    allocate(r(n,m),Qold(n,m),Q(n,m),Qnew(n,m),P(n,m),phi(n,m),phiold(n,m),phioldold(n,m))

    allocate(tempm1(m,m),tau(m,m),tautilde(m,m),Mm(m,m),HR(m,m),T(m,m),Tnew(m,m),Rtilde(m,m),Rold(m,m),Roldold(m,m))

    allocate(ag(m,m),bg(m,m),cg(m,m),dg(m,m),agold(m,m),bgold(m,m),cgold(m,m),dgold(m,m),dgoldold(m,m),bgoldold(m,m))

    allocate(QM(2*m,2*m),HM(2*m,2*m),RT(2*m,2*m))

    allocate(qrtau1(m),qrtau2(2*m),work(3*n),dwork(n),ipiv(m))

    

    normb=zlange('F',n,m,b,n,dwork)
    !call zmatnorm(n,m,b,normb)

    lwork=3*n

    x=(0.0d0,0.0d0)
    r=b !residual, start with x=0.0d0

    call system_clock(count=tc1)
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

    call zungqr(n,m,m,Q,n,qrtau1,work,lwork,info)
    if(info .ne. 0) then
        print *,'block minres zungqr 1 error ',info
    end if
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
    call system_clock(count=tc2)
    qrtime=qrtime+elapsed_time(tc1,tc2)

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

    linloops=k        

        !!!!!!M=Q'A*Q!!!!!!!!!!!!!!!!!!!! 
        !P=zQ-A*Q
        call system_clock(count=tc1)
        P=Q
        call mkl_zcsrmm('N', n, m, n, (-1.0d0,0.0d0), matdescra, dsa, jsa, isa, isa(2), Q, n, ze, P, n)
        call system_clock(count=tc2)
        mvtime=mvtime+elapsed_time(tc1,tc2)
        nmatvec(feastit,cpnum)=nmatvec(feastit,cpnum)+m

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
        call system_clock(count=tc1) 
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
   
        call zungqr(n,m,m,Qnew,n,qrtau1,work,lwork,info)
        if(info .ne. 0) then
            print *,'block minres zungqr 2 error ',info
        end if
        call system_clock(count=tc2)
        qrtime=qrtime+elapsed_time(tc1,tc2)
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
        call system_clock(count=tc1)
        QM=(0.0d0,0.0d0)
        QM=RT

        call zgeqrf(2*m,2*m,QM,2*m,qrtau2,work,lwork,info)
        if(info .ne. 0) then
            print *,'block minres zgeqrf 3 error ',info
            stop
        end if
    
        HR=(0.0d0,0.0d0)
        do i=1,m
            do j=1,i
                HR(i,j)=conjg(QM(j,i))
            end do
        end do
   
        call zungqr(2*m,2*m,2*m,QM,2*m,qrtau2,work,lwork,info)
        if(info .ne. 0) then
            print *,'block minres zungqr 3 error ',info
        end if

        HM=conjg(transpose(QM))
        call system_clock(count=tc2)
        qrtime=qrtime+elapsed_time(tc1,tc2)

        ag=HM(1:m,1:m)
        cg=HM(m+1:2*m,1:m)
        bg=HM(1:m,m+1:2*m)
        dg=HM(m+1:2*m,m+1:2*m)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !phi=(Q-phiold*Rold-phioldold*Roldold)*inv(HR[1:m,1:m])
        templ1t=conjg(transpose(Q))
        !templ1t=templ1-Rold'*phiold'
        call zgemm('C','C',m,n,m,(-1.0d0,0.0d0),Rold,m,phiold,n,(1.0d0,0.0d0),templ1t,m)
        !templ1t=templ1t-Roldold'*phioldold'
        call zgemm('C','C',m,n,m,(-1.0d0,0.0d0),Roldold,m,phioldold,n,(1.0d0,0.0d0),templ1t,m)
        !solve system Rkk'*phi'=templ1t
        
        !call zgemm('N','N',m,n,m,(1.0d0,0.0d0),tempm1,m,templ1t,m,(0.0d0,0.0d0),templ2t,m) 
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
        !call zmatnorm(m,m,tautilde,normr)
        !print *,'Minres ',k,normr/normb
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

    deallocate(templ1,templ2,templ1t)
    deallocate(r,Qold,Q,Qnew,P,phi,phiold,phioldold)

    deallocate(tempm1,tau,tautilde,Mm,HR,T,Tnew,Rtilde,Rold,Roldold)

    deallocate(ag,bg,cg,dg,agold,bgold,cgold,dgold,dgoldold,bgoldold)

    deallocate(QM,HM,RT)

    deallocate(qrtau1,qrtau2,work,dwork,ipiv)


end subroutine zminresBlock



subroutine zminres(UPLO,n,dsa,isa,jsa,ze,b,x,tol,maxit,linloops)
use rundata
implicit none

integer :: i

character :: UPLO
integer :: n,maxit,linloops
integer, dimension(*) :: isa,jsa
complex (kind=kind(0.0d0)) :: ze !complex shift
complex (kind=kind(0.0d0)), dimension(*) :: dsa !unshifted matrix
complex (kind=kind(0.0d0)), dimension(n) :: b,x
double precision :: tol

!!!minres stuff:
complex (kind=kind(0.0d0)), dimension(:),allocatable :: r,v,vold,p,pold,poldold
double precision :: beta,betaold,gtilde,gamma,rnrm,bnrm
complex (kind=kind(0.0d0)) :: c,cm1,s,sm1,expphi,eta1,alpha,epsilon,delta,h,eta
complex (kind=kind(0.0d0)), dimension(:),allocatable :: actr,Tv,Ap,Apold,Apoldold !residual stuff
!complex (kind=kind(0.0d0)), dimension(:),allocatable :: x !solution
!!!BLAS and lapack:
character, dimension(6) :: matdescra
double precision, external :: dznrm2
complex (kind=kind(0.0d0)), external :: zdotc

    allocate(r(n),v(n),vold(n),p(n),pold(n),poldold(n))
    allocate(actr(n),Tv(n),Ap(n),Apold(n),Apoldold(n))
    !allocate(x(n))

    bnrm=dznrm2(n,b,1)

    x=(0.0d0,0.0d0)
    r=b

    v=(0.0d0,0.0d0)
    vold=v
    p=v
    pold=v

    beta=dznrm2(n,r,1)!beta=norm(r)
    betaold=beta
    eta1=beta
    c=(1.0d0,0.0d0)
    cm1=c
    s=(0.0d0,0.0d0)
    sm1=s
    expphi=(1.0d0,0.0d0)

    actr=b
    Tv=(0.0d0,0.0d0)
    Ap=Tv
    Apold=Tv

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

    do i=1,maxit
        linloops=i
        
        vold=v
        v=r/beta

        !Tv=T*v = A*v+re(ze)*v:
        call system_clock(count=tc1)
        Tv=v
        call mkl_zcsrmm('N', n, 1, n, (-1.0d0,0.0d0), matdescra, dsa, jsa, isa, isa(2), v, n, (1.0d0,0.0d0)*dble(ze), Tv, n) 
        call system_clock(count=tc2)
        mvtime=mvtime+elapsed_time(tc1,tc2)
        nmatvec(feastit,cpnum)=nmatvec(feastit,cpnum)+1


        alpha=zdotc(n,Tv,1,v,1)!alpha = Tv'*v
        r=Tv-alpha*v-beta*vold

        betaold=beta
        beta=dznrm2(n,r,1) !beta=norm(r)

        epsilon=sm1*betaold
        delta=s*alpha+c*cm1*betaold*dble(expphi)
        h=(-1.0d0,0.0d0)*s*cm1*betaold*expphi**(-1.0d0,0.0d0)+c*(alpha-(0.0d0,1.0d0)*aimag(ze))
        gtilde=abs(h)
        if(gtilde==0.0d0) then
            expphi=0.0d0
        else
            expphi=h/gtilde
        end if

        gamma=sqrt(gtilde**2+beta**2)

        cm1=c
        sm1=s
        c=gtilde/gamma
        s=beta/gamma

        poldold=pold
        pold=p
        p=(v-delta*pold-epsilon*poldold)/gamma
        eta=c*eta1*expphi

        x=x+eta*p

        !print *,'etap',sum(eta*p)

        !update residual to see if we're done
        Apoldold=Apold
        Apold=Ap
        Ap=(Tv+(0.0d0,1.0d0)*aimag(ze)*v-delta*Apold-epsilon*Apoldold)/gamma
        actr=actr-eta*Ap

        eta1=(-1.0d0,0.0d0)*s*eta1*expphi

        rnrm=dznrm2(n,actr,1)/bnrm
        
        if(rnrm<tol) then
            exit
        end if
    end do

    !b=x
end subroutine zminres



subroutine zminresSVD(UPLO,nrow,mcol,dsa,isa,jsa,ze,b,x,tol,maxit,linloops)
use rundata
implicit none

integer :: i

character :: UPLO
integer :: nrow,mcol
integer :: n,maxit,linloops,m
integer, dimension(*) :: isa,jsa
complex (kind=kind(0.0d0)) :: ze !complex shift
complex (kind=kind(0.0d0)), dimension(*) :: dsa !unshifted matrix
complex (kind=kind(0.0d0)), dimension(min(nrow,mcol)) :: b,x
double precision :: tol

!!!minres stuff:
complex (kind=kind(0.0d0)), dimension(:),allocatable :: r,v,vold,p,pold,poldold
double precision :: beta,betaold,gtilde,gamma,rnrm,bnrm
complex (kind=kind(0.0d0)) :: c,cm1,s,sm1,expphi,eta1,alpha,epsilon,delta,h,eta
complex (kind=kind(0.0d0)), dimension(:),allocatable :: actr,Tv,Ap,Apold,Apoldold !residual stuff

complex (kind=kind(0.0d0)), dimension(:), allocatable :: temp1 
!complex (kind=kind(0.0d0)), dimension(:),allocatable :: x !solution
!!!BLAS and lapack:
character, dimension(6) :: matdescra
double precision, external :: dznrm2
complex (kind=kind(0.0d0)), external :: zdotc

    n=min(nrow,mcol)
    m=max(nrow,mcol)

    allocate(temp1(m))

    allocate(r(n),v(n),vold(n),p(n),pold(n),poldold(n))
    allocate(actr(n),Tv(n),Ap(n),Apold(n),Apoldold(n))
    !allocate(x(n))

    bnrm=dznrm2(n,b,1)

    x=(0.0d0,0.0d0)
    r=b

    v=(0.0d0,0.0d0)
    vold=v
    p=v
    pold=v

    beta=dznrm2(n,r,1)!beta=norm(r)
    betaold=beta
    eta1=beta
    c=(1.0d0,0.0d0)
    cm1=c
    s=(0.0d0,0.0d0)
    sm1=s
    expphi=(1.0d0,0.0d0)

    actr=b
    Tv=(0.0d0,0.0d0)
    Ap=Tv
    Apold=Tv

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

    do i=1,maxit
        linloops=i
        
        vold=v
        v=r/beta

        !Tv=T*v = A*v+re(ze)*v:
        call system_clock(count=tc1)
        
        !Tv=v
        
        if(nrow>mcol) then 
            !temp1=A*v
            !print *,nrow,mcol
            !print *,n
            call mkl_zcsrmm('N', nrow, 1, mcol, (1.0d0,0.0d0), matdescra, dsa, jsa, isa, isa(2), v, n, (0.0d0,0.0d0), temp1, m)
            !Tv=-A'temp1+dble(ze)v
            Tv=v
            call mkl_zcsrmm('C', nrow, 1, mcol, (-1.0d0,0.0d0), matdescra, dsa, jsa, isa, isa(2), temp1, m, (1.0d0,0.0d0)*dble(ze), Tv, n)
        else
            !temp1=A'v
            call mkl_zcsrmm('C', nrow, 1, mcol, (1.0d0,0.0d0), matdescra, dsa, jsa, isa, isa(2), v, n, (0.0d0,0.0d0), temp1, m)
            !Tv=-Atemp1+dble(ze)v
            Tv=v
            call mkl_zcsrmm('N', nrow, 1, mcol, (-1.0d0,0.0d0), matdescra, dsa, jsa, isa, isa(2), temp1, m, (1.0d0,0.0d0)*dble(ze), Tv, n) 
        end if

        !call mkl_zcsrmm('N', n, 1, n, (-1.0d0,0.0d0), matdescra, dsa, jsa, isa, isa(2), v, n, (1.0d0,0.0d0)*dble(ze), Tv, n) 
        call system_clock(count=tc2)
        mvtime=mvtime+elapsed_time(tc1,tc2)
        nmatvec(feastit,cpnum)=nmatvec(feastit,cpnum)+2


        alpha=zdotc(n,Tv,1,v,1)!alpha = Tv'*v
        r=Tv-alpha*v-beta*vold

        betaold=beta
        beta=dznrm2(n,r,1) !beta=norm(r)

        epsilon=sm1*betaold
        delta=s*alpha+c*cm1*betaold*dble(expphi)
        h=(-1.0d0,0.0d0)*s*cm1*betaold*expphi**(-1.0d0,0.0d0)+c*(alpha-(0.0d0,1.0d0)*aimag(ze))
        gtilde=abs(h)
        if(gtilde==0.0d0) then
            expphi=0.0d0
        else
            expphi=h/gtilde
        end if

        gamma=sqrt(gtilde**2+beta**2)

        cm1=c
        sm1=s
        c=gtilde/gamma
        s=beta/gamma

        poldold=pold
        pold=p
        p=(v-delta*pold-epsilon*poldold)/gamma
        eta=c*eta1*expphi

        x=x+eta*p

        !print *,'etap',sum(eta*p)

        !update residual to see if we're done
        Apoldold=Apold
        Apold=Ap
        Ap=(Tv+(0.0d0,1.0d0)*aimag(ze)*v-delta*Apold-epsilon*Apoldold)/gamma
        actr=actr-eta*Ap

        eta1=(-1.0d0,0.0d0)*s*eta1*expphi

        rnrm=dznrm2(n,actr,1)/bnrm
        
        if(rnrm<tol) then
            exit
        end if
    end do

    deallocate(temp1)

    deallocate(r,v,vold,p,pold,poldold)
    deallocate(actr,Tv,Ap,Apold,Apoldold)

    !b=x
end subroutine zminresSVD


