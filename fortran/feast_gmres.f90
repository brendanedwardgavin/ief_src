





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
    double precision, dimension(*),intent(inout) :: list
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
    
    call quicksort(list,lo,i-1)
    call quicksort(list,i+1,hi)
    end if

end subroutine quicksort



subroutine blockGivens(n,V,M) 
implicit none

!V: 2n x n    input
!M: 2n x 2n   output
integer :: n
complex (kind=kind(0.0d0)), dimension(2*n,n) :: V
complex (kind=kind(0.0d0)), dimension(2*n,2*n) :: M

!X=conjg(transpose(V(1:n,1:n)))
!Y=conjg(transpose(V(n+1:2*n,1:n)))

!B=\(X,Y)

!call zgetrf(n,n,X,)

!A=F, F'F=BB'
!C=F, F'F=B'B

!M(1:n,1:n)=A

!M(1:n,n+1:2*n)=A*B'
!M(n+1:2*n,1:n)=-1.0d0*C*B
!M(n+1:2*n,n+1:2*n)=C
end subroutine blockGivens





subroutine blockGMRESarnoldi(UPLO,n,m,dsa,isa,jsa,ze,kmax,restarts,Brhs,Xlhs,eps,loops,blockstart)
use rundata
implicit none

character :: UPLO
integer :: n,m,kmax,restarts,loops,blockstart
integer, dimension(*) :: isa,jsa
complex (kind=kind(0.0d0)) :: ze
complex (kind=kind(0.0d0)), dimension(*) :: dsa
complex (kind=kind(0.0d0)), dimension(n,m) :: Brhs,Xlhs
double precision :: eps !target residual error

!!!!Arnoldi stuff
complex (kind=kind(0.0d0)), dimension(:,:),allocatable :: V,H,Bsm,ym,Htmp,R,R2,Xtmp

!!!BLAS and lapack:
character, dimension(6) :: matdescra
integer :: lwork,info
complex (kind=kind(0.0d0)), dimension(:),allocatable :: work
double precision, external :: dznrm2

double precision :: error,error2

integer :: i,j,jmax,l

!print *,kmax,restarts

allocate(V(n,(kmax+1)*m),H(1:(kmax+1)*m,1:kmax*m), Htmp((kmax+1)*m,kmax*m)  ,Bsm((kmax+1)*m,m), ym((kmax+1)*m,m), R(1:n,1:m), R2(1:n,1:m), Xtmp(1:n,1:m) )



lwork=m*(kmax+1)*m*2
allocate(work(lwork))

if(UPLO=='F') then
    matdescra(1)='G'
else
    matdescra(1)='H'
end if
!matdescra(1)='G'
matdescra(2)=UPLO
matdescra(3)='N'
matdescra(4)='F'

R=Brhs
Xlhs=(0.0d0,0.0d0)

loops=0
do i=1,restarts

    if(i>1) then
        !R=Brhs-A*Xlhs
        !R=Brhs(1:n,1:m)
        R=R2
        !call mkl_zcsrmm('N', n, m, n, (-1.0d0,0.0d0), matdescra, dsa, jsa, isa, isa(2), Xlhs(1,1), n, (1.0d0,0.0d0), R(1,1), n)
    end if

    V(1:n,1:m)=R(1:n,1:m)

    do j=1,kmax
        loops=loops+1
        !next arnoldi step
        call system_clock(count=tc1)
        call blockArnoldiIt(UPLO,n,m,dsa,isa,jsa,ze,kmax,j,V,H,Bsm)
        call system_clock(count=tc2)
        arnolditime=arnolditime+elapsed_time(tc1,tc2)

        !solve system H(1:(j+1)*m,1:j*m)ym=Bsm(1:(j+1)*m,1:m)
        ym(1:(j+1)*m,1:m)=Bsm(1:(j+1)*m,1:m)
        Htmp(1:(j+1)*m,1:j*m)=H(1:(j+1)*m,1:j*m)

        !print *,'H=',Htmp(1,2)
        !print *,'Bsm=',sum(Bsm)
        !print *,'V=',sum(V)

        call system_clock(count=tc1)      
        call zgels('N',(j+1)*m,j*m,m,Htmp(1,1),(kmax+1)*m,ym(1,1),(kmax+1)*m,work,lwork,info)
        if(info .ne. 0) then
            print *,'Error in blockGMRESarnoldi'
            print *,'ZGELS error ',info
            stop
        end if
        call system_clock(count=tc2)
        lstime=lstime+elapsed_time(tc1,tc2)

        !print *,'Ym=',ym(1,1)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !use QR explicitly
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            goto 116
           !get QR factorization
            !call ZGEQRF( (j+1)*m, m, Htmp(1:(j+1)*m,1:j*m),(j+1)*m , qrtau, work, lwork, info )
            if (info .ne. 0) then
                print *,'Problem with least squares solution in GMRES'
                print *,'ZGEQRF error info = ',info
                stop
            end if
    
            !get R matrix
            !Rs(1:m*m0,1:m*m0)=Av2(1:m*m0,1:m*m0)
            !call zlacpy('F',m*m0,m*m0,Av2(1,1),n,Rs(1,1),m*m0) 

            !get Q matrix
            !call ZUNGQR(  n, m0*m, m0*m, Av2, n, qrtau, work, lwork, info )
            if (info .ne. 0) then
                print *,'Problem with least squares solution in GMRES'
                print *,'ZUNGQR error info = ',info
                stop
            end if
            
            !form reduced right hand side matrix:
            !use V(1:n,1:m) since V(1:n,1:m) = r = B-Ax is the right hand side
            !call ZGEMM('C','N',m*m0,m,n,(1.0d0,0.0d0),Av2,n,V(1:n,1:m),n,(0.0d0,0.0d0),Bs(1,1),m0*m)
         
            !solve upper triangular system Rs*x=Q'*Bs
            !call ZTRTRS( 'U', 'N', 'N', m*m0, m, Rs, m*m0, Bs, m0*m, info )
            if (info .ne. 0) then
                print *,'Problem with least squares solution in GMRES'
                print *,'ZTRTRS error info = ',info
                stop
            end if

        116 continue
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



        !measure residual:
        !error=norm(R-V[1:n,1:(j+1)*m]*H[1:(j+1)*m,1:j*m]*ym)/norm(Brhs)
        R2=R
        call zgemm('N','N',(j+1)*m,m,j*m,(1.0d0,0.0d0),H(1,1),(kmax+1)*m,ym(1,1),(kmax+1)*m,(0.0d0,0.0d0),Htmp(1,1),(kmax+1)*m)
        call zgemm('N','N',n,m,(j+1)*m,(-1.0d0,0.0d0),V(1,1),n,Htmp(1,1),(kmax+1)*m,(1.0d0,0.0d0),R2(1,1),n)

        !Xtmp=Xlhs(1:n,1:m)
        !call zgemm('N','N',n,m,j*m,(1.0d0,0.0d0),V(1:n,1:j*m),n,ym(1:j*m,1:m),j*m,(1.0d0,0.0d0),Xtmp(1:n,1:m),n)
        !R2=Brhs(1:n,1:m)
        !call mkl_zcsrmm('N', n, m, n, (-1.0d0,0.0d0), matdescra, dsa, jsa, isa, isa(2), Xtmp(1,1), n, (1.0d0,0.0d0), R2(1,1), n)

        !print *,i,j
        error=0.0d0
        do l=1,m
            error2=dznrm2(n,R2(:,l),1)/dznrm2(n,Brhs(:,l),1)
            if (error2>error) error=error2
            linres(feastit,cpnum,loops,blockstart+l-1)=error2
        !    print *,'    ',error2
        end do

        !print *,i,j,error

        jmax=j
        if(error<eps) then 
            !Xlhs(1:n,1:m)=Xtmp
            !return
            !print *,'   error=',error,eps
            exit
        end if
    end do

    R=R2

    !update solution
    !Xlhs=Xlhs+V[1:n,1:j*m]*ym
    j=jmax
    call zgemm('N','N',n,m,j*m,(1.0d0,0.0d0),V(1,1),n,ym(1,1),(kmax+1)*m,(1.0d0,0.0d0),Xlhs(1,1),n)

    if(error<eps) then 
        exit
    end if
end do


    !Xtmp=Xlhs(1:n,1:m)
    !call zgemm('N','N',n,m,j*m,(1.0d0,0.0d0),V(1:n,1:j*m),n,ym(1:j*m,1:m),j*m,(1.0d0,0.0d0),Xtmp(1:n,1:m),n)
    !R2=Brhs(1:n,1:m)
    !call mkl_zcsrmm('N', n, m, n, (-1.0d0,0.0d0), matdescra, dsa, jsa, isa, isa(2), Xtmp(1,1), n, (1.0d0,0.0d0), R2(1,1), n)

    !error=0.0d0
    !do l=1,m
    !    error2=dznrm2(n,R2(:,l),1)/dznrm2(n,Brhs(:,l),1)
    !    if (error2>error) error=error2
    !end do
    
    !print *,'its ',i,j
    !print *,'     Final error=',error

deallocate(V,H, Htmp  ,Bsm, ym, R, R2 )

end subroutine blockGMRESarnoldi



subroutine blockArnoldiIt(UPLO,n,m,dsa,isa,jsa,ze,k,k0,V,H,Bsm)
use rundata
implicit none

character :: UPLO
integer :: n,m,k,k0
integer, dimension(*) :: isa,jsa
complex (kind=kind(0.0d0)) :: ze
complex (kind=kind(0.0d0)), dimension(*) :: dsa
complex (kind=kind(0.0d0)), dimension(n,*) :: V
complex (kind=kind(0.0d0)), dimension((k+1)*m,*) ::H,Bsm

integer :: nnza

!complex (kind=kind(0.0d0)), dimension(:,:), allocatable :: Hnew,Vnew
complex (kind=kind(0.0d0)), dimension(m,m) :: Hnew
complex (kind=kind(0.0d0)), dimension(n,m) :: Vnew

!!!BLAS and lapack:
character, dimension(6) :: matdescra
!!lapack stuff:
!complex (kind=kind(0.0d0)), dimension(:), allocatable :: work,qrtau
complex (kind=kind(0.0d0)), dimension(3*n) :: work
complex (kind=kind(0.0d0)), dimension(m) :: qrtau
integer :: lwork,info

integer :: j,i

nnza=isa(n+1)+1

!allocate(Vnew(n,m),Hnew(m,m))
lwork=3*n
!allocate(work(lwork),qrtau(m))

if(UPLO=='F') then
    matdescra(1)='G'
else
    matdescra(1)='H'
end if
!matdescra(1)='G'
matdescra(2)=UPLO
matdescra(3)='N'
matdescra(4)='F'

if(k0==1) then !initialize everything
    H(1:(k+1)*m,1:k*m)=(0.0d0,0.0d0)    
    Bsm(1:(k+1)*m,1:m)=(0.0d0,0.0d0)
    
    !QR=V
    !Bsm(1:m,1:m)=R(1:m,1:m)
    !V(:,1:m)=Q
    
    call system_clock(count=tc1)
    !get QR factorization
    call ZGEQRF( n, m, V(1,1), n, qrtau, work, lwork, info )
    if (info .ne. 0) then
        print *,'Problem with least squares solution in blockArnoldiIt'
        print *,'ZGEQRF error info = ',info
        stop
    end if

    !put R matrix into H:
    do i=1,m
        do j=i,m
            Bsm(i,j)=V(i,j)
        end do
    end do

    !put Q matrix into V:
    call ZUNGQR(  n, m, m, V(1,1), n, qrtau, work, lwork, info )
    if (info .ne. 0) then
        print *,'Problem with least squares solution in ArnoldiIt'
        print *,'ZUNGQR error info = ',info
        stop
    end if
    call system_clock(count=tc2)
    qrtime=qrtime+elapsed_time(tc1,tc2)

end if


!Do matrix multiply:

!Vnew=A*V0(:,(i-1)*m+1:i*m)
call system_clock(count=tc1)
Vnew=V(:,(k0-1)*m+1:k0*m)
call mkl_zcsrmm('N', n, m, n, (-1.0d0,0.0d0), matdescra, dsa, jsa, isa, isa(2), V(1,(k0-1)*m+1), n, ze, Vnew, n)
!call mkl_zcsrmm('N', n, m, n, (1.0d0,0.0d0), matdescra, dsa, jsa, isa, isa(2), V(1,(k0-1)*m+1), n, (0.0d0,0.0d0), Vnew, n)
call system_clock(count=tc2)
mvtime=mvtime+elapsed_time(tc1,tc2)
nmatvec(feastit,cpnum)=nmatvec(feastit,cpnum)+m

call system_clock(count=tc1)
do j=1,k0 !Orthogonalize with respect to previous basis vectors:
    !Hnew=V(:,(j-1)*m+1:j*m)'*Vnew
    call zgemm('C','N',m,m,n,(1.0d0,0.0d0),V(1,(j-1)*m+1),n,Vnew,n,(0.0d0,0.0d0),Hnew,m)

    !H((j-1)*m+1:j*m,(k0-1)*m+1:k0*m)=Hnew
    H((j-1)*m+1:j*m,(k0-1)*m+1:k0*m)=Hnew
    !Vnew=Vnew-V(:,(j-1)*m+1:j*m)*Hnew

    call zgemm('N','N',n,m,m,(-1.0d0,0.0d0),V(:,(j-1)*m+1),n,Hnew,m,(1.0d0,0.0d0),Vnew,n)

end do
call system_clock(count=tc2)
gstime=gstime+elapsed_time(tc1,tc2)


!Use QR to orthonormalize new vectors:

call system_clock(count=tc1)
!get QR factorization
call ZGEQRF( n, m, Vnew, n, qrtau, work, lwork, info )
if (info .ne. 0) then
    print *,'Problem with least squares solution in blockArnoldiIt'
    print *,'ZGEQRF error info = ',info
    stop
end if

!put R matrix into H:
!H(k0*m+1:(k0+1)*m,(k0-1)*m+1:k0*m)=Vnew(1:m,1:m)
do i=1,m
    do j=i,m
        H(k0*m+i,(k0-1)*m+j)=Vnew(i,j)
    end do
end do


!put Q matrix into V:
call ZUNGQR(  n, m, m, Vnew, n, qrtau, work, lwork, info )
if (info .ne. 0) then
    print *,'Problem with least squares solution in ArnoldiIt'
    print *,'ZUNGQR error info = ',info
    stop
end if

V(:,k0*m+1:(k0+1)*m)=Vnew
call system_clock(count=tc2)
qrtime=qrtime+elapsed_time(tc1,tc2)
!V(:,k0*m+1:(k0+1)*m)=Q
!H(k0*m+1:(k0+1)*m,(k0-1)*m+1:k0*m)=R(1:m,1:m)
!deallocate(Vnew,Hnew)
!deallocate(work,qrtau)

end subroutine blockArnoldiIt









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




subroutine zfeast_gmres(ijob,stateVars,Brhs,x,V,Av,Ax,ze,n,m,maxm,eps,restarts,m0,xwork,workin,Av2,times)
    use rundata
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
    !double precision, external :: elapsed_time
  
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




subroutine zfeast_gmres_norm(ijob,stateVars,Brhs,x,V,Av,Ax,ze,n,m,maxm,eps,restarts,m0,xwork,workin,Av2,times)
use rundata
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
    !double precision, external :: elapsed_time
  
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
            !call ZGEQRF( n, m0*m, Av2, n, qrtau, work, lwork, info )
            !if (info .ne. 0) then
            !    print *,'Problem with least squares solution in GMRES'
            !    print *,'ZGEQRF error info = ',info
            !    stop
            !end if
    
            !get R matrix
            !Rs(1:m*m0,1:m*m0)=Av2(1:m*m0,1:m*m0)
            !call zlacpy('F',m*m0,m*m0,Av2(1,1),n,Rs(1,1),m*m0) 

            !get Q matrix
            !call ZUNGQR(  n, m0*m, m0*m, Av2, n, qrtau, work, lwork, info )
            !if (info .ne. 0) then
            !    print *,'Problem with least squares solution in GMRES'
            !    print *,'ZUNGQR error info = ',info
            !    stop
            !end if
            
            !form normal equations:
            call ZGEMM('C','N',m*m0,m*m0,n,(1.0d0,0.0d0),Av2,n,Av2,n,(0.0d0,0.0d0),Rs(1,1),m0*m)
         
            !form reduced right hand side matrix:
            call ZGEMM('C','N',m*m0,m,n,(1.0d0,0.0d0),Av2,n,V,n,(0.0d0,0.0d0),Bs(1,1),m0*m)

            !solve upper triangular system Rs*x=Q'*Bs
            !call ZTRTRS( 'U', 'N', 'N', m*m0, m, Rs, m*m0, Bs, m0*m, info )
            call ZPOSV( 'U', m*m0, m, Rs, m*m0, Bs, m*m0, info )
            if (info .ne. 0) then
                print *,'Problem with least squares solution in GMRES'
                print *,'ZPOSV error info = ',info
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
end subroutine zfeast_gmres_norm




subroutine zfeast_BiCGSTABRes(UPLO,n,dsa,isa,jsa,ze,nnza,B,X,maxit,eps,neigs,error)
implicit none
!A=Az=(ze*I-A) in this routine

    integer :: n,maxit,neigs
    complex (kind=kind(0.0d0)) :: ze
    double precision :: eps,error  !target error, output error
    character :: UPLO

    !!!!!!!!!!!!!!!!!!!!!!!!  Sparse matrix:
    complex (kind=kind(0.0d0)),dimension(*) :: dsa
    integer,dimension(*) :: isa,jsa
    integer :: nnza

    !!! RHS, solution
    complex (kind=kind(0.0d0)), dimension(n) :: B,X

    !!! CG stuff
    integer :: m
    complex (kind=kind(0.0d0)), dimension(:), allocatable :: R,Rnew,P,Rstar,S,tv1,tv2
    complex (kind=kind(0.0d0)) :: tz1,tz2,alpha,beta,omega
    
    !!!BLAS and lapack:
    character, dimension(6) :: matdescra
    
    integer :: i,j,debug,its
    double precision :: dtemp
    double precision, external :: dznrm2
    complex (kind=kind(0.0d0)) :: zdotc
    double precision :: errorprint

    m=1

    !if (m<neigs) then
    !    print *,'zfeast_bicgstabRes error: neigs > m',neigs,m
    !    stop
    !endif

    
    debug=0

    if(UPLO=='F') then
        matdescra(1)='G'
    else
        matdescra(1)='H'
    end if
    matdescra(2)=UPLO
    matdescra(3)='N'
    matdescra(4)='F'

    !all this allocating is probably slow; maybe have user allocate once and for all?
    allocate(R(n),Rnew(n),Rstar(n),P(n),tv1(n),tv2(n),S(n))

    !X=0.0
    X(1:n)=0.0d0

    R=B(1:n)
    !call mkl_zcsrmm('C', n, m, n, (-1.0d0,0.0d0), matdescra, dsa, jsa, isa, isa(2), B, n, (1.0d0,0.0d0)*conjg(ze), R, n)

    P=R
    Rstar=R

    its=0

    do i=1,maxit
        its=its+1
        
        !alpha=r'*rstar/(rstar'*A*p)

        !----tv1=A*p
        tv1=P
        call mkl_zcsrmm('N', n, m, n, (-1.0d0,0.0d0), matdescra, dsa, jsa, isa, isa(2), P, n, ze, tv1, n)
        !----tz1=rstar'*tv1
        tz1=zdotc(n,tv1,1,Rstar,1)
        !----tz2=r'*rstar
        tz2=zdotc(n,R,1,Rstar,1)
        !----alpha=tz2/tz1
        alpha=tz2/tz1

        !S=R-alpha*A*P
        S=R-alpha*tv1

        !omega=(A*S)'*S/([A*S]'*[A*S])
        
        !----tv2=A*S
        tv2=S
        call mkl_zcsrmm('N', n, m, n, (-1.0d0,0.0d0), matdescra, dsa, jsa, isa, isa(2), S, n, ze, tv2, n)
        !----tz1=tv2'*S
        tz1=zdotc(n,tv2,1,S,1)
        !----tz2=tv2'*tv2
        tz2=zdotc(n,tv2,1,tv2,1)
        !----omega=tz1/tz2
        omega=tz1/tz2

        X=X+alpha*P+omega*S

        !Rnew=S-omega*A*S
        Rnew=S-omega*tv2

        !beta=(rnew'*rstar)/(r'*rstar) * alpha/omega
        !----tz1=rnew'*rstar
        tz1=zdotc(n,Rnew,1,Rstar,1)
        !----tz2=R'*Rstar
        tz2=zdotc(n,R,1,Rstar,1)
        !----beta=(tz1/tz2)*(alpha/omega)
        beta=(tz1/tz2)*(alpha/omega)

        !P=Rnew+beta*(p-omega*A*P)
        P=Rnew+beta*(p-omega*tv1)

        R=Rnew

        !measure residual
        error=dznrm2(n,R(:),1)/dznrm2(n,B(:),1)

        if(error<eps) exit !if error is low enough, end loop

    end do  
    
    print *,'   linits=',its
    print *,'      errors=',error
    print *,'beta = ',beta
    deallocate(R,Rnew,Rstar,P,tv1,tv2,S)

end subroutine zfeast_BiCGSTABRes



  subroutine  zbicgstab(UPLO,N,saz,isa,jsa,M0,fj,xj,res,ares,nbit_out,epso,comd,info) 
    implicit none
    character(len=1) :: UPLO
    integer :: N,M0
    complex(kind=kind(1.0d0)),dimension(*) :: saz
    integer,dimension(*) :: isa,jsa

    complex(kind=kind(1.0d0)), dimension(N,*) :: fj
    complex(kind=kind(1.0d0)), dimension(N,*):: xj
    double precision, dimension(M0) :: res
    integer :: info
    integer :: nbit_out
    double precision ::epso
    double precision ::ares
    logical :: comd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

    integer :: rank,nbit_out0
    integer :: p,k
    integer :: s,sizeAj,i,etiquette
    integer :: tca,tcb,tc1,tc2,tim
    logical :: fail,loop,half
    character(len=1) :: ASTRU,AFAC
    complex(kind=kind(1.0d0)) :: zdotc
    complex(kind=kind(1.0d0)), dimension(:,:),allocatable :: rb,pp,pb,v,ss,sb,t,rj
    complex(kind=kind(1.0d0)),dimension(:),allocatable ::rho_1,prho_1,pbeta,palpha,pomega,prho_2,aux0,paux0,aux1,paux1
    double precision,dimension(:),allocatable :: resf,nss

    complex(kind=kind(1.0d0)) :: ONEC,ZEROC

    character, dimension(6) :: matdescra


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
    ONEC=(1.0d0,0.0d0)
    ZEROC=(0.0d0,0.0d0)

if(UPLO=='F') then
        matdescra(1)='G'
    else
        matdescra(1)='H'
    end if
    matdescra(2)=UPLO
    matdescra(3)='N'
    matdescra(4)='F'


rank=0


    ! Remove Fortran runtime dependency
    allocate(rj(N,M0))
    allocate(rb(N,M0))
    allocate(pp(N,M0))
    allocate(pb(N,M0))
    allocate(v(N,M0))
    allocate(ss(N,M0))
    allocate(sb(N,M0))
    allocate(t(N,M0))
    allocate(rho_1(M0))
    allocate(prho_1(M0))
    allocate(prho_2(M0))
    allocate(pbeta(M0))
    allocate(palpha(M0))
    allocate(pomega(M0))
    allocate(aux0(M0))
    allocate(paux0(M0))
    allocate(aux1(M0))
    allocate(paux1(M0))
    allocate(nss(M0))
    allocate(resf(M0))
!!!!!!!!!!!!!


    s=M0
    sizeAj=N

    info=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Initial rj,xj,res and resf !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!! SOLVE M*xj=fj ==> xj
!    call ZCOPY(sizeAj*s,fj(1,1),1,xj(1,1),1)   
 !!   call zSOLVER_for_BICGSTAB(pre,xj,info)


      
!!!!!!!!!!!!!!! MATVEC for the initial residual r=f-Ax
    call ZCOPY(sizeAj*s,fj(1,1),1,rj(1,1),1)

!    call zMATVEC_for_BICGSTAB(-ONEC,mat,xj,ONEC,rj)
! call zcsrmm(UPLO,'N',N,N,M0,-ONEC,saz,isa,jsa,xj,(1.0d0,0.0d0),rj)
call mkl_zcsrmm('N', n, M0, n, -ONEC, matdescra, saz, jsa, isa, isa(2), xj, n, (1.0d0,0.0d0), rj, n)

!!!
    do i=1,s
       res(i)=maxval(abs(rj(1:N,i)))
       resf(i)=maxval(abs(fj(1:N,i))) ! cte all along
    enddo
    ares=maxval(res/resf)
    if (comd) print *,'rel. residual before iteration',ares

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!! BIPCG-STAB 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if ((comd).and.(rank==0)) then
       tca=0
       tcb=0
    end if

!!!!!!!!!! choose rb
    call ZCOPY(sizeAj*s,rj(1,1),1,rb(1,1),1)

    nbit_out0=0
    k=0
    fail=.false.
    loop=.true.

    !loop=.false.


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111
    if (.not.((nbit_out0<=nbit_out).and.(ares.gt.epso).and.(fail.eqv..false.))) loop=.false. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111

    Do While (loop) 
       k=k+1
       nbit_out0=k
       if (k>1) prho_2=prho_1 

!!!!!!!!!!! SCALAR PRODUCT on different RHS
       do i=1,s
          rho_1(i)=zdotc(sizeAj,rb(1,i),1,rj(1,i),1)
       end do
       prho_1=rho_1
!!!!!!!!!!!! TEST
       do i=1,s
          if (prho_1(i)==ZEROC) then
             info= -201
             if ((comd).and.(rank==0)) print*,'ATTENTION ----- BICG-STAB OUT FAILED, prho_1=0 !!'
             return 
             fail=.true.
          end if
       end do

!!!!!!!!!!!! CONDITION
       IF (k==1) THEN
          call ZCOPY(sizeAj*s,rj(1,1),1,pp(1,1),1)
       ELSE
          pbeta=(prho_1/prho_2)*(palpha/pomega)
          do i=1,s
             call ZAXPY(sizeAj,-pomega(i),v(1,i),1,pp(1,i),1)
             call ZSCAL(sizeAj,pbeta(i),pp(1,i),1)
             call ZAXPY(sizeAj,ONEC,rj(1,i),1,pp(1,i),1)
             !      pp(:,i)=r(:,i)+pbeta(i)*(pp(:,i)-pomega(i)*v(:,i)) ! by processor
          end do
       END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!! SOLVE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !! solve the spike M*pb=pp(k) ---> pb
       call ZCOPY(sizeAj*s,pp(1,1),1,pb(1,1),1)
       if (comd) call system_clock(tc1,tim) 
       !call zSOLVER_for_BICGSTAB(pre,pb,info)

       if (info/=0) return 
       if (comd) then
          call system_clock(tc2,tim)
          tcb=tcb+tc2-tc1
       end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! MAT-VEC A by pb, results in v(k)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !v=ZERO
       if (comd) call system_clock(tc1,tim)
       !call zMATVEC_for_BICGSTAB(ONEC,mat,pb,ZEROC,v)
       !call zcsrmm(UPLO,'N',N,N,M0,ONEC,saz,isa,jsa,pb,(0.0d0,0.0d0),v)
call mkl_zcsrmm('N', n, M0, n, ONEC, matdescra, saz, jsa, isa, isa(2), pb, n, (0.0d0,0.0d0), v, n)


       if (comd) then
          call system_clock(tc2,tim)
          tca=tca+tc2-tc1
       end if

!!!!!! SCALAR PRODUCT on different RHS
       do i=1,s!! different RHS
          aux0(i)=zdotc(sizeAj,rb(1,i),1,v(1,i),1)
       end do
       paux0=aux0
!!!!!!
       palpha=prho_1/paux0
       do i=1,s
          !ss(:,i)=r(:,i)-palpha(i)*v(:,i) ! by processors
          call ZCOPY(sizeAj,rj(1,i),1,ss(1,i),1)
          call ZAXPY(sizeAj,-palpha(i),v(1,i),1,ss(1,i),1) 
       end do

!!!! CHECK THE  NORM Loo of ss
!!! call pvnorm_Loo(nss,ss,p,nbpart,pm%new_comm_world)
       do i=1,s
          nss(i)=maxval(abs(ss(1:N,i)))
       enddo

       half=.false.
       ares=maxval(nss/resf)
       if (ares<epso) half=.true. !

       IF (half) then
          do i=1,s
             !      xj(:,i)=xj(:,i)+palpha(i)*pb(:,i) 
             call ZAXPY(sizeAj,palpha(i),pb(1,i),1,xj(1,i),1) 
          end do
          res=nss
          if ((comd).and.(rank==0)) print *,(k-1)*1.0d0+0.5d0,ares


       ELSE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
!!!!!!!!!!!! SOLVE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! solve the spike M*sb=ss ---> sb
          call ZCOPY(sizeAj*s,ss(1,1),1,sb(1,1),1)
          if (comd) call system_clock(tc1,tim)
          !call zSOLVER_for_BICGSTAB(pre,sb,info)

          if (info/=0) return 

          if (comd) then
             call system_clock(tc2,tim)
             tcb=tcb+tc2-tc1
          end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! MAT-VEC A by sb, results in t
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !t=ZERO
          if (comd) call system_clock(tc1,tim)
          !call zMATVEC_for_BICGSTAB(ONEC,mat,sb,ZEROC,t)
          !call zcsrmm(UPLO,'N',N,N,M0,ONEC,saz,isa,jsa,sb,(0.0d0,0.0d0),t)
          call mkl_zcsrmm('N', n, M0, n, ONEC, matdescra, saz, jsa, isa, isa(2), sb, n, (0.0d0,0.0d0), rj, n)


          if (comd) then
             call system_clock(tc2,tim)
             tca=tca+tc2-tc1
          end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! SCALAR PRODUCTS on different RHS
          do i=1,s!! different RHS
             aux0(i)=zdotc(sizeAj,t(1,i),1,ss(1,i),1)
             aux1(i)=zdotc(sizeAj,t(1,i),1,t(1,i),1)
          end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          paux0=aux0
          paux1=aux1
!!!!!!
          pomega=paux0/paux1
          do i=1,s
             if (pomega(i)==ZEROC) then
		info= -202 
                if ((comd).and.(rank==0)) print *,'ATTENTION ----- BICG-STAB OUT FAILED, pomega=0'
		return 
                fail=.true.
             end if
          end do

          do i=1,s
             ! xj(:,i)=xj(:,i)+palpha(i)*pb(:,i)+pomega(i)*sb(:,i)
             call ZAXPY(sizeAj,palpha(i),pb(1,i),1,xj(1,i),1) 
             call ZAXPY(sizeAj,pomega(i),sb(1,i),1,xj(1,i),1) 
          end do

          do i=1,s
             !   r(:,i)=ss(:,i)-pomega(i)*t(:,i)
             call ZCOPY(sizeAj,ss(1,i),1,rj(1,i),1)
             call ZAXPY(sizeAj,-pomega(i),t(1,i),1,rj(1,i),1) 
          end do

!!!! NORM Loo RESIDUAL (||r||oo)
!!! call pvnorm_Loo(res,rj,p,nbpart,pm%new_comm_world)
          do i=1,s
             res(i)=maxval(abs(rj(1:N,i)))
          enddo

          ares=maxval(res/resf)
          if ((comd).and.(rank==0)) print *,k,ares

       END IF



!!!! test for the other loop
       if (.not.((nbit_out0<=nbit_out).and.(ares.gt.epso).and.(fail.eqv..false.))) loop=.false. 


    end Do



    if ((comd).and.(rank==0)) then
       print *,'TIME postprocess MATMUL ',tca*1.0d0/tim
       print *,'TIME postprocess SOLVE ',tcb*1.0d0/tim
       print *,''
    end if


    nbit_out=nbit_out0


    deallocate(rj,rb,pp,pb,v,ss,sb,t)
    deallocate(nss,rho_1,prho_1,pbeta,palpha,pomega,prho_2,aux0,paux0,aux1,paux1,resf)

  end subroutine zbicgstab


subroutine zfeast_cglsRes(UPLO,n,m,dsa,isa,jsa,ze,nnza,B,X,maxit,eps,neigs,error,its)
use rundata
implicit none
!A=Az=(ze*I-A) in this routine

    integer :: n,m,maxit,neigs
    complex (kind=kind(0.0d0)) :: ze
    double precision :: eps,error  !target error, output error
    character :: UPLO

    !!!!!!!!!!!!!!!!!!!!!!!!  Sparse matrix:
    complex (kind=kind(0.0d0)),dimension(*) :: dsa
    integer,dimension(*) :: isa,jsa
    integer :: nnza

    !!! RHS, solution
    complex (kind=kind(0.0d0)), dimension(n,*) :: B,X

    !!! CG stuff
    complex (kind=kind(0.0d0)), dimension(:,:), allocatable :: R,Rnew,P,lambda,psi,T,D
    complex (kind=kind(0.0d0)), dimension(:,:), allocatable :: temp1,temp2,sqtemp1,sqtemp2
    
    !!!BLAS and lapack:
    character, dimension(6) :: matdescra
    integer :: info
    integer, dimension(m) :: ipiv

    integer :: i,j,debug,its
    double precision :: dtemp
    double precision, external :: dznrm2
    double precision, dimension(:),allocatable :: errorlist
    double precision :: errorprint

    if (m<neigs) then
        print *,'zfeast_cglsRes error: neigs > m',neigs,m
        stop
    endif

    allocate(errorlist(m))

    debug=0

    if(UPLO=='F') then
        matdescra(1)='G'
    else
        matdescra(1)='H'
    end if
    matdescra(2)=UPLO
    matdescra(3)='N'
    matdescra(4)='F'

    !all this allocating is probably slow; maybe have user allocate once and for all?
    allocate(R(n,m),Rnew(n,m),P(n,m),lambda(m,m),psi(m,m),temp1(n,m),sqtemp1(m,m),sqtemp2(m,m),temp2(n,m),D(n,m),T(n,m))

    !X=0.0
    X(1:n,1:m)=0.0d0

    D=B(1:n,1:m)

    !R=A'*B
    call system_clock(count=tc1)
    R=B(1:n,1:m)
    call mkl_zcsrmm('C', n, m, n, (-1.0d0,0.0d0), matdescra, dsa, jsa, isa, isa(2), B, n, (1.0d0,0.0d0)*conjg(ze), R, n)
    !call mkl_zcsrmm('C', n, m, n, (1.0d0,0.0d0), matdescra, dsa, jsa, isa, isa(2), B, n, (0.0d0,0.0d0), R, n)
    !R(1:n,1:m)=-1.0*B(1:n,1:m)
    call system_clock(count=tc2)
    nmatvec(feastit,cpnum)=nmatvec(feastit,cpnum)+m
    mvtime=mvtime+elapsed_time(tc1,tc2)
    !call mkl_zcsrmm('N', n, m, n, 1.0d0, matdescra, dsa, jsa, isa, isa(2), X, n, 0.0d0, linwork2, n)
    
    !P=-R
    P=R

    if(debug>0) then
        print *,'X=',X(1:n,1)
        print *,'P=',P(1:n,1)
        print *,'R=',R(1:n,1)
    end if

    its=0

    do i=1,maxit
        its=its+1
        
        !lambda=inv(P'*A'*A*P)*R'*R
        !-----T=A*P
        
        call system_clock(count=tc1)
        T=P
        call mkl_zcsrmm('N', n, m, n, (-1.0d0,0.0d0), matdescra, dsa, jsa, isa, isa(2), P, n, ze, T, n) 
        !call mkl_zcsrmm('N', n, m, n, (1.0d0,0.0d0), matdescra, dsa, jsa, isa, isa(2), P, n, (0.0d0,0.0d0), T, n)
        call system_clock(count=tc2)
        nmatvec(feastit,cpnum)=nmatvec(feastit,cpnum)+m
        mvtime=mvtime+elapsed_time(tc1,tc2)

        call system_clock(count=tc1)
        !-----sqtemp=T'*T
        call zgemm('C','N',m,m,n,(1.0d0,0.0d0),T,n,T,n,(0.0d0,0.0d0),sqtemp1,m)
        !-----lambda=R'*R    !might be better to do (inv(P'A'AP)R')R, not sure...
        call zgemm('C','N',m,m,n,(1.0d0,0.0d0),R,n,R,n,(0.0d0,0.0d0),lambda,m)
        !-----lambda=\(sqtemp1,lambda)
        !call zposv('U',m,m,sqtemp1,m,lambda,m,info)
        call zgesv(m,m,sqtemp1,m,ipiv,lambda,m,info)
        if(info .ne. 0) then
            print *,'CGLS error: ZGESV info ',info
            stop
        end if

        if(debug>0) print *,'lambda = ', lambda(1,1)

        !X=X+P*lambda
        call zgemm('N','N',n,m,m,(1.0d0,0.0d0),P,n,lambda,m,(1.0d0,0.0d0),X,n)
        
        if(debug>0) print *,'Xnew = ',X(1:n,1)

        !D=D-T*lambda
        call zgemm('N','N',n,m,m,(-1.0d0,0.0d0),T,n,lambda,m,(1.0d0,0.0d0),D,n)
        call system_clock(count=tc2)
        gstime=gstime+elapsed_time(tc1,tc2)

        !Rnew=A'*D 
        
        call system_clock(count=tc1)
        Rnew=D
        call mkl_zcsrmm('C', n, m, n, (-1.0d0,0.0d0), matdescra, dsa, jsa, isa, isa(2), D, n, conjg(ze), Rnew, n)
        !call mkl_zcsrmm('C', n, m, n, (1.0d0,0.0d0), matdescra, dsa, jsa, isa, isa(2), D, n, (0.0d0,0.0d0), Rnew, n)
        call system_clock(count=tc2)
        nmatvec(feastit,cpnum)=nmatvec(feastit,cpnum)+m
        mvtime=mvtime+elapsed_time(tc1,tc2)
        !measure residual
        error=0.0d0
        do j=1,m
            errorlist(j)=dznrm2(n,D(:,j),1)/dznrm2(n,B(:,j),1)
        end do

        call quicksort(errorlist,1,m)
        !call selection_sort(errorlist,m)
        !error=errorlist(neigs)
        error=errorlist(m)!errorlist(neigs)!(m)
        errorprint=errorlist(neigs)

        !print *,error

        if(error<eps) exit !if error is low enough, end loop

        if(debug>0) print *,'Rnew = ', Rnew(1:n,1)

        call system_clock(count=tc1)
        !psi=inv(R'*R)*Rnew'*Rnew
        !-----sqtemp1=R'*R
        call zgemm('C','N',m,m,n,(1.0d0,0.0d0),R,n,R,n,(0.0d0,0.0d0),sqtemp1,m)
        !-----sqtemp2=Rnew'*Rnew
        call zgemm('C','N',m,m,n,(1.0d0,0.0d0),Rnew,n,Rnew,n,(0.0d0,0.0d0),psi,m)
        !-----psi=\(sqtemp1,psi)
        !call zposv('U',m,m,sqtemp1,m,psi,m,info)
        call zgesv(m,m,sqtemp1,m,ipiv,psi,m,info)
        if(info .ne. 0) then
            print *,'CGLS error: second ZGESV info ',info
            stop
        end if

        if(debug>0) print *,'psi = ', psi(1,1)

        !P=Rnew+P*psi
        temp1=P
        P=Rnew
        call zgemm('N','N',n,m,m,(1.0d0,0.0d0),temp1,n,psi,m,(1.0d0,0.0d0),P,n)
        call system_clock(count=tc2)
        gstime=gstime+elapsed_time(tc1,tc2)

        if(debug>0) print *,'Pnew = ', P(1:n,1)

        R=Rnew
       
        !temp1=X(1:n,1:m)
        !call mkl_zcsrmm('N', n, m, n, (-1.0d0,0.0d0), matdescra, dsa, jsa, isa, isa(2), X, n, ze, temp1, n)
        !temp2=B(1:n,1:m)-temp1
        !error=0.0d0
        !do j=1,m
        !    dtemp=dznrm2(n,temp2(1:n,j),1)/dznrm2(n,B(1:n,j),1)
        !   !dtemp=dznrm2(n,R(1:n,j),1)/dznrm2(n,B(1:n,j),1)
        !    if (dtemp>error) error=dtemp
        !end do
        !print *,i,error

        if(debug>0 .and. i>1) stop
    end do  
    
    !print *,'   linits=',its
    !print *,'      errors=',error,errorprint
    if (error<errorprint) then
        !do i=1,m
        !print *,'   ',errorlist(i)
        !end do
    end if

    !measure actual error:
    !temp1=X(1:n,1:m)
    !call mkl_zcsrmm('N', n, m, n, (-1.0d0,0.0d0), matdescra, dsa, jsa, isa, isa(2), X, n, ze, temp1, n)
    !temp2=B(1:n,1:m)-temp1
    !do j=1,m
    !    errorlist(j)=dznrm2(n,temp2(1:n,j),1)/dznrm2(n,B(1:n,j),1)
    !end do
    !call quicksort(errorlist,1,m)
    !errorprint=errorlist(neigs)
    !error=errorlist(m)
    !print *,'      acterr=',error,errorprint

    !stop

    deallocate(errorlist,R,Rnew,P,lambda,psi,temp1,sqtemp1,sqtemp2,temp2,T,D)
end subroutine zfeast_cglsRes




subroutine zfeast_cgls(UPLO,n,m,dsa,isa,jsa,ze,nnza,B,X,maxit)
implicit none
!A=Az=(ze*I-A) in this routine

    integer :: n,m,maxit
    complex (kind=kind(0.0d0)) :: ze
    character :: UPLO

    !!!!!!!!!!!!!!!!!!!!!!!!  Sparse matrix:
    complex (kind=kind(0.0d0)),dimension(*) :: dsa
    integer,dimension(*) :: isa,jsa
    integer :: nnza

    !!! RHS, solution
    complex (kind=kind(0.0d0)), dimension(n,*) :: B,X

    !!! CG stuff
    complex (kind=kind(0.0d0)), dimension(:,:), allocatable :: R,Rnew,P,lambda,psi,T,D
    complex (kind=kind(0.0d0)), dimension(:,:), allocatable :: temp1,temp2,sqtemp1,sqtemp2
    
    !!!BLAS and lapack:
    character, dimension(6) :: matdescra
    integer :: info
    integer, dimension(m) :: ipiv

    integer :: i,j,debug
    double precision :: error,dtemp
    double precision, external :: dznrm2

    debug=0

    if(UPLO=='F') then
        matdescra(1)='G'
    else
        matdescra(1)='H'
    end if
    matdescra(2)=UPLO
    matdescra(3)='N'
    matdescra(4)='F'

    !all this allocating is probably slow; maybe have user allocate once and for all?
    allocate(R(n,m),Rnew(n,m),P(n,m),lambda(m,m),psi(m,m),temp1(n,m),sqtemp1(m,m),sqtemp2(m,m),temp2(n,m),D(n,m),T(n,m))

    !X=0.0
    X(1:n,1:m)=0.0d0

    D=B(1:n,1:m)

    !R=A'*B
    R=B(1:n,1:m)
    call mkl_zcsrmm('C', n, m, n, (-1.0d0,0.0d0), matdescra, dsa, jsa, isa, isa(2), B, n, (1.0d0,0.0d0)*conjg(ze), R, n)
    !R(1:n,1:m)=-1.0*B(1:n,1:m)
    
    !call mkl_zcsrmm('N', n, m, n, 1.0d0, matdescra, dsa, jsa, isa, isa(2), X, n, 0.0d0, linwork2, n)
    
    !P=-R
    P=R

    if(debug>0) then
        print *,'X=',X(1:n,1)
        print *,'P=',P(1:n,1)
        print *,'R=',R(1:n,1)
    end if

    do i=1,maxit

        !lambda=inv(P'*A'*A*P)*R'*R
        !-----T=A*P
        T=P
        call mkl_zcsrmm('N', n, m, n, (-1.0d0,0.0d0), matdescra, dsa, jsa, isa, isa(2), P, n, ze, T, n)
        !-----sqtemp=T'*T
        call zgemm('C','N',m,m,n,(1.0d0,0.0d0),T,n,T,n,(0.0d0,0.0d0),sqtemp1,m)
        !-----lambda=R'*R    !might be better to do (inv(P'A'AP)R')R, not sure...
        call zgemm('C','N',m,m,n,(1.0d0,0.0d0),R,n,R,n,(0.0d0,0.0d0),lambda,m)
        !-----lambda=\(sqtemp1,lambda)
        !call zposv('U',m,m,sqtemp1,m,lambda,m,info)
        call zgesv(m,m,sqtemp1,m,ipiv,lambda,m,info)
        if(info .ne. 0) then
            print *,'CGLS error: ZGESV info ',info
            stop
        end if

        if(debug>0) print *,'lambda = ', lambda(1,1)

        !X=X+P*lambda
        call zgemm('N','N',n,m,m,(1.0d0,0.0d0),P,n,lambda,m,(1.0d0,0.0d0),X,n)
        
        if(debug>0) print *,'Xnew = ',X(1:n,1)

        !D=D-T*lambda
        call zgemm('N','N',n,m,m,(-1.0d0,0.0d0),T,n,lambda,m,(1.0d0,0.0d0),D,n)

        !Rnew=A'*D 
        Rnew=D
        call mkl_zcsrmm('C', n, m, n, (-1.0d0,0.0d0), matdescra, dsa, jsa, isa, isa(2), D, n, conjg(ze), Rnew, n)

        if(debug>0) print *,'Rnew = ', Rnew(1:n,1)

        !psi=inv(R'*R)*Rnew'*Rnew
        !-----sqtemp1=R'*R
        call zgemm('C','N',m,m,n,(1.0d0,0.0d0),R,n,R,n,(0.0d0,0.0d0),sqtemp1,m)
        !-----sqtemp2=Rnew'*Rnew
        call zgemm('C','N',m,m,n,(1.0d0,0.0d0),Rnew,n,Rnew,n,(0.0d0,0.0d0),psi,m)
        !-----psi=\(sqtemp1,psi)
        !call zposv('U',m,m,sqtemp1,m,psi,m,info)
        call zgesv(m,m,sqtemp1,m,ipiv,psi,m,info)
        if(info .ne. 0) then
            print *,'CGLS error: second ZGESV info ',info
            stop
        end if

        if(debug>0) print *,'psi = ', psi(1,1)

        !P=Rnew+P*psi
        temp1=P
        P=Rnew
        call zgemm('N','N',n,m,m,(1.0d0,0.0d0),temp1,n,psi,m,(1.0d0,0.0d0),P,n)

        if(debug>0) print *,'Pnew = ', P(1:n,1)

        R=Rnew
       
        temp1=X(1:n,1:m)
        !call mkl_zcsrmm('N', n, m, n, (-1.0d0,0.0d0), matdescra, dsa, jsa, isa, isa(2), X, n, ze, temp1, n)
        temp2=B(1:n,1:m)-temp1

        error=0.0d0
        do j=1,m
            dtemp=dznrm2(n,temp2(1:n,j),1)/dznrm2(n,B(1:n,j),1)
            !dtemp=dznrm2(n,R(1:n,j),1)/dznrm2(n,B(1:n,j),1)
            if (dtemp>error) error=dtemp
        end do
        !print *,i,error

        if(debug>0 .and. i>1) stop
    end do

    !stop

    deallocate(R,Rnew,P,lambda,psi,temp1,sqtemp1,sqtemp2,temp2,T,D)
end subroutine zfeast_cgls



subroutine zfeast_cgne(UPLO,n,m,dsa,isa,jsa,ze,nnza,B,X,maxit)
implicit none
!A=Az=(ze*I-A) in this routine

    integer :: n,m,maxit
    complex (kind=kind(0.0d0)) :: ze
    character :: UPLO

    !!!!!!!!!!!!!!!!!!!!!!!!  Sparse matrix:
    complex (kind=kind(0.0d0)),dimension(*) :: dsa
    integer,dimension(*) :: isa,jsa
    integer :: nnza

    !!! RHS, solution
    complex (kind=kind(0.0d0)), dimension(n,*) :: B,X

    !!! CG stuff
    complex (kind=kind(0.0d0)), dimension(:,:), allocatable :: R,Rnew,P,lambda,psi
    complex (kind=kind(0.0d0)), dimension(:,:), allocatable :: temp1,temp2,sqtemp1,sqtemp2,AP
    
    !!!BLAS and lapack:
    character, dimension(6) :: matdescra
    integer :: info
    integer, dimension(m) :: ipiv

    integer :: i,j,debug
    double precision :: error,dtemp
    double precision, external :: dznrm2

    debug=0

    if(UPLO=='F') then
        matdescra(1)='G'
    else
        matdescra(1)='H'
    end if
    matdescra(2)=UPLO
    matdescra(3)='N'
    matdescra(4)='F'

    !all this allocating is probably slow; maybe have user allocate once and for all?
    allocate(R(n,m),Rnew(n,m),P(n,m),lambda(m,m),psi(m,m),temp1(n,m),sqtemp1(m,m),sqtemp2(m,m),temp2(n,m),AP(n,m))

    !X=0.0
    X(1:n,1:m)=0.0d0

    !R=A'*A*X-B
    R=B(1:n,1:m)
    call mkl_zcsrmm('C', n, m, n, (1.0d0,0.0d0), matdescra, dsa, jsa, isa, isa(2), B, n, (-1.0d0,0.0d0)*conjg(ze), R, n)
    !R(1:n,1:m)=-1.0*B(1:n,1:m)
    
    !call mkl_zcsrmm('N', n, m, n, 1.0d0, matdescra, dsa, jsa, isa, isa(2), X, n, 0.0d0, linwork2, n)
    
    !P=-R
    P=-1.0*R

    if(debug>0) then
        print *,'X=',X(1:n,1)
        print *,'P=',P(1:n,1)
        print *,'R=',R(1:n,1)
    end if

    do i=1,maxit

        !lambda=inv(P'*A'*A*P)*R'*R
        !(P'*A'*A*P)x=R'*R
        
        !-----temp1=A*P
        AP=P
        call mkl_zcsrmm('N', n, m, n, (-1.0d0,0.0d0), matdescra, dsa, jsa, isa, isa(2), P, n, (1.0d0,0.0d0)*ze, AP, n)
        !-----sqtemp=AP'*AP
        call zgemm('C','N',m,m,n,(1.0d0,0.0d0),AP,n,AP,n,(0.0d0,0.0d0),sqtemp1,m)
        !-----lambda=R'*R    !might be better to do (inv(P'A'AP)R')R, not sure...
        call zgemm('C','N',m,m,n,(1.0d0,0.0d0),R,n,R,n,(0.0d0,0.0d0),lambda,m)
        !-----lambda=\(sqtemp1,lambda)
        !call zposv('U',m,m,sqtemp1,m,lambda,m,info)
        call zgesv(m,m,sqtemp1,m,ipiv,lambda,m,info)
        if(info .ne. 0) then
            print *,'CGNE error: ZPOSV info ',info
            stop
        end if

        if(debug>0) print *,'lambda = ', lambda(1,1)

        !X=X+P*lambda
        call zgemm('N','N',n,m,m,(1.0d0,0.0d0),P,n,lambda,m,(1.0d0,0.0d0),X,n)
        
        if(debug>0) print *,'Xnew = ',X(1:n,1)

        !Rnew=R+A'*A*P*lambda
        !-----temp1=AP*lambda    !maybe do this first and then add it to X?
        call zgemm('N','N',n,m,m,(1.0d0,0.0d0),AP,n,lambda,m,(0.0d0,0.0d0),temp1,n)
        !-----temp2=A*temp1
        !call mkl_zcsrmm('N', n, m, n, (1.0d0,0.0d0), matdescra, dsa, jsa, isa, isa(2), temp1, n, (0.0d0,0.0d0), temp2, n)
        !-----Rnew=R+A'*temp1
        temp2=temp1
        call mkl_zcsrmm('C', n, m, n, (-1.0d0,0.0d0), matdescra, dsa, jsa, isa, isa(2), temp1, n, conjg(ze), temp2, n)
        Rnew=R+temp2

        if(debug>0) print *,'Rnew = ', Rnew(1:n,1)

        !psi=inv(R'*R)*Rnew'*Rnew
        !-----sqtemp1=R'*R
        call zgemm('C','N',m,m,n,(1.0d0,0.0d0),R,n,R,n,(0.0d0,0.0d0),sqtemp1,m)
        !-----sqtemp2=Rnew'*Rnew
        call zgemm('C','N',m,m,n,(1.0d0,0.0d0),Rnew,n,Rnew,n,(0.0d0,0.0d0),psi,m)
        !-----psi=\(sqtemp1,psi)
        !call zposv('U',m,m,sqtemp1,m,psi,m,info)
        call zgesv(m,m,sqtemp1,m,ipiv,psi,m,info)
        if(info .ne. 0) then
            print *,'CGNE error: second ZPOSV info ',info
            stop
        end if

        if(debug>0) print *,'psi = ', psi(1,1)

        !P=-Rnew+P*psi
        temp1=P
        P=Rnew
        call zgemm('N','N',n,m,m,(1.0d0,0.0d0),temp1,n,psi,m,(-1.0d0,0.0d0),P,n)

        if(debug>0) print *,'Pnew = ', P(1:n,1)

        R=Rnew
       
        temp1=X(1:n,1:m)
        call mkl_zcsrmm('N', n, m, n, (-1.0d0,0.0d0), matdescra, dsa, jsa, isa, isa(2), X, n, ze, temp1, n)
        temp2=B(1:n,1:m)-temp1

        error=0.0d0
        do j=1,m
            dtemp=dznrm2(n,temp2(1:n,j),1)/dznrm2(n,B(1:n,j),1)
            !dtemp=dznrm2(n,R(1:n,j),1)/dznrm2(n,B(1:n,j),1)
            if (dtemp>error) error=dtemp
        end do
        !print *,i,error

        if(debug>0 .and. i>1) stop
    end do

    !stop

    deallocate(R,Rnew,P,lambda,psi,temp1,sqtemp1,sqtemp2,temp2,AP)
end subroutine zfeast_cgne


subroutine zfeast_gmrespre(ijob,stateVars,Brhs,x,V,Av,Ax,ze,n,m,maxm,eps,restarts,m0,xwork,workin,Av2,times)
use rundata
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
    !double precision, external :: elapsed_time
  
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


subroutine pre_dmult1(n,m,A,X,Y,Mout,ze)
    implicit none
    integer :: n,m
    complex (kind=kind(0.0d0)) :: ze
    complex (kind=kind(0.0d0)), dimension(n,n), intent(in) :: A
    complex (kind=kind(0.0d0)), dimension(n,m), intent(in) :: X,Y
    complex (kind=kind(0.0d0)), dimension(n,m), intent(out) ::Mout

    complex (kind=kind(0.0d0)), dimension(:,:),allocatable :: xta1,xta,xax,ixax,invixaxxta,B

    integer :: i,j,k

    !!!!lapack:
    integer :: info
    integer, dimension(:),allocatable :: ipiv

    integer :: debug

    debug=1

    allocate(xta1(n,m),xta(m,n),xax(m,m),invixaxxta(m,n),B(m,m))

    allocate(ipiv(m))

    call zgemm('N','N',n,m,n,(1.0d0,0.0d0),A,n,X,n,(0.0d0,0.0d0),xta1,n)
    xta=conjg(transpose(xta1))

    xax=(0.0d0,0.0d0)
    do i=1,m
        xax(i,i)=(1.0d0,0.0d0)
    end do

    call zgemm('N','N',m,m,n,(1.0d0,0.0d0),xta,m,X,n,(1.0d0,0.0d0),xax,m)
    !xax=matmul(xta,X)
    !do i=1,m
    !    xax(i,i)=xax(i,i)+(1.0d0,0.0d0)
    !end do

    call zgetrf(m,m,xax,m,ipiv,info)

    if (info .ne. 0) then
        print *,'pre_dmult zgetrf error',info
        stop
    end if

    invixaxxta=xta
    call zgetrs('N',m,n,xax,m,ipiv,invixaxxta,m,info)

    if (info .ne. 0) then
        print *,'pre_dmult zgetrs error',info
        stop
    end if

    call zgemm('N','N',m,m,n,(1.0d0,0.0d0),invixaxxta,m,Y,n,(0.0d0,0.0d0),B,m) 

    Mout=Y

    call zgemm('N','N',n,m,m,(-1.0d0,0.0d0),X,n,B,m,(1.0d0,0.0d0),Mout,n)

end subroutine pre_dmult1




subroutine pre_dmult(n,m,A,X,Y,Mout,ze)
    implicit none
    integer :: n,m
    complex (kind=kind(0.0d0))::ze
    complex (kind=kind(0.0d0)), dimension(n,n), intent(in) :: A
    complex (kind=kind(0.0d0)), dimension(n,m), intent(in) :: X,Y
    complex (kind=kind(0.0d0)), dimension(n,m), intent(out) ::Mout

    complex (kind=kind(0.0d0)), dimension(:,:),allocatable :: xta1,xta,xax,ixax,invixaxxta,B

    integer :: i,j,k

    !!!!lapack:
    integer :: info
    integer, dimension(:),allocatable :: ipiv

    integer :: debug

    debug=1

    allocate(xta1(n,m),xta(m,n),xax(m,m),invixaxxta(m,n),B(m,m))

    allocate(ipiv(m))

    !xta=transpose(A*X)
    call zgemm('N','N',n,m,n,(1.0d0,0.0d0),A,n,X,n,(0.0d0,0.0d0),xta1,n)
    xta=conjg(transpose(xta1))

    !xax=I-transpose(A*X)*X/ze
    xax=(0.0d0,0.0d0)
    do i=1,m
        xax(i,i)=(1.0d0,0.0d0)
    end do

    call zgemm('N','N',m,m,n,(-1.0d0,0.0d0)/ze,xta,m,X,n,(1.0d0,0.0d0),xax,m)
    !xax=matmul(xta,X)
    !do i=1,m
    !    xax(i,i)=xax(i,i)+(1.0d0,0.0d0)
    !end do

    !invixaxxta=inv(I-X'A'X/ze)*X'A'
    call zgetrf(m,m,xax,m,ipiv,info)

    if (info .ne. 0) then
        print *,'pre_dmult zgetrf error',info
        stop
    end if

    invixaxxta=xta
    call zgetrs('N',m,n,xax,m,ipiv,invixaxxta,m,info)

    if (info .ne. 0) then
        print *,'pre_dmult zgetrs error',info
        stop
    end if

    !B=inv(I-X'A'X/ze)*X'A'Y/ze^2    
    call zgemm('N','N',m,m,n,(1.0d0,0.0d0),invixaxxta,m,Y,n,(0.0d0,0.0d0),B,m) 
    B=B/ze**2

    Mout=Y/ze

    call zgemm('N','N',n,m,m,(1.0d0,0.0d0),X,n,B,m,(1.0d0,0.0d0),Mout,n)

    !Mout=Y

end subroutine pre_dmult





subroutine dfeast_gmres_norm(ijob,stateVars,Brhs,x,V,Av,Ax,ze,n,m,restarts,m0,xwork,workin,Av2)
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
            !call ZGEQRF( n, m0*m, Av2, n, qrtau, work, lwork, info )
            !if (info .ne. 0) then
            !    print *,'Problem with least squares solution in GMRES'
            !    print *,'ZGEQRF error info = ',info
            !    stop
            !end if
    
            !get R matrix
            !Rs(1:m*m0,1:m*m0)=Av2(1:m*m0,1:m*m0)
            !call zlacpy('F',m*m0,m*m0,Av2(1,1),n,Rs(1,1),m*m0)

            !get Q matrix
            !call ZUNGQR(  n, m0*m, m0*m, Av2, n, qrtau, work, lwork, info )
            !if (info .ne. 0) then
            !    print *,'Problem with least squares solution in GMRES'
            !    print *,'ZUNGQR error info = ',info
            !    stop
            !end if
            
            !!!!!normal equations:
            call ZGEMM('C','N',m*m0,m,n,(1.0d0,0.0d0),Av2,n,Av2(1:n,1:m),n,(0.0d0,0.0d0),Rs(1,1),m0*m)

            !form reduced right hand side matrix:
            !use V(1:n,1:m) since V(1:n,1:m) = r = B-Ax is the right hand side
            call ZGEMM('C','N',m*m0,m,n,(1.0d0,0.0d0),Av2,n,V(1:n,1:m),n,(0.0d0,0.0d0),Bs(1,1),m0*m)
                        
            !solve upper triangular system Rs*x=Q'*Bs
            !call ZTRTRS( 'U', 'N', 'N', m*m0, m, Rs, m*m0, Bs, m0*m, info )
            call ZPOSV( 'U', m*m0, m*m0, Rs, m*m0, Bs, m*m0, info )
            if (info .ne. 0) then
                print *,'Problem with least squares solution in GMRES'
                print *,'ZPOSV error info = ',info
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

end subroutine dfeast_gmres_norm
