program iterativeFeast

implicit none

integer :: i,j,k
double precision :: dtmp1,dtmp2

!!!!!!!!!!!!!!!!!!!!!!!!  FEAST:

integer :: n,m0,ijob,loop,info,m
integer, dimension(64)::fpm
double precision :: epsout,emin,emax
double precision, dimension(:), allocatable :: res,e
double precision, dimension(:,:),allocatable :: work,aq,bq,x
complex (kind=kind(0.0d0)) :: ze
complex (kind=kind(0.0d0)), dimension(:,:), allocatable :: zwork
complex (kind=kind(0.0d0)), dimension(:), allocatable :: zne,wne

!!!!!!!!!!!!!!!!!!!!!!!!  Linear system solver:

integer :: linjob,linIterations,kdim
! linjob: RCI task returned by linear solver
! linIterations: number of linear system iterations; each iteration requires kdim matvec operations
! kdim: number of Krylov subspace blocks for linear solver
integer, dimension(3) :: linState !internal parameters for linear solver
complex (kind=kind(0.0d0)), dimension(:,:), allocatable :: xsol,V,Av,Ax,zwork2
! xsol,V,Av,Ax,zwork2: internal storage for linear system solver
complex (kind=kind(0.0d0)) :: ze2 !FEAST contour point for linear system sovler
double precision, dimension (:,:), allocatable :: linwork1,linwork2 !more internal storage for linear system solver

!!!!!!!!!!!!!!!!!!!!!!! Artificial eigenvalue problem:
complex (kind=kind(0.0d0)), dimension(:,:), allocatable :: ztmpmat1,ztmpmat2,ztmpmat3

!!!!!!!!!!!!!!!!!!!!!!! External BLAS:
double precision, external :: dznrm2
double precision, external :: elapsed_time
character, dimension(6) :: matdescra

!!!!!!!!!!!!!!!!!!!!!!! Timing

integer :: totalc1,totalc2,c1,c2

!!!!!!!!!!!!!!!!!!!!!!! convergence vs time
double precision, dimension(:),allocatable :: reslist,timelist
integer :: oldloop

!!!!!!!!!!!!!!!!!!!!!!! read file:

integer :: pc
character(len=100) :: name
character(len=1) UPLO,PRE,SHG,EG
character(len=1) :: cc
double precision,dimension(:),allocatable :: dsa,dca,ssa
integer,dimension(:),allocatable :: isa,jsa,ica,jca,sisa,sjsa
integer :: nnza


!allocate(ztmpmat1(n,n),ztmpmat2(n,n),ztmpmat3(n,n))

call feastinit(fpm)
call getarg(1,name)

!!!!!!!!!!!! DRIVER_FEAST_SPARSE input file  
  open(10,file=trim(name)//'.in',status='old')
  read(10,*) SHG ! type of eigenvalue problem "General, Hermitian, Symmetric" 
  read(10,*) EG ! type of eigenvalue probl  g== sparse generalized, e== sparse standard
  read(10,*) PRE  ! "PRE"==(s,d,c,z) resp. (single real,double real,complex,double complex) 
  read(10,*) UPLO ! UPLO==(F,L,U) reps. (Full csr, Lower csr, Upper csr)

  read(10,*) emin
  read(10,*) emax

  read(10,*) m0
  read(10,*) pc ! Some changes from default for fpm
  do i=1,pc
     read(10,*) j,fpm(j)
  enddo

  close(10)



!!!!!!!!!!!read matrix A
  open(10,file=trim(name)//'_A.mtx',status='old')
  k=0
  cc='%'
  do while(cc=='%')
     k=k+1
     read(10,'(A1)') cc
  end do
  close(10)

  open(10,file=trim(name)//'_A.mtx',status='old')
  do i=1,k-1
     read(10,'(A1)') cc
  enddo
  read(10,*) n,n,nnza
  allocate(ica(nnza))
  allocate(jca(nnza))

  allocate(dca(nnza))
  do i=1,nnza
     read(10,*) ica(i),jca(i),dca(i)
  end do

  close(10)

  !! create csr format
  allocate(isa(1:n+1))
  allocate(jsa(1:nnza))
  allocate(dsa(1:nnza))
  call dcoo2csr(n,nnza,ica,jca,dca,isa,jsa,dsa)

!!!!!!!!!!!!!! Set up FEAST:

allocate(res(m0))
allocate(e(m0))
allocate(work(n,m0))
allocate(aq(m0,m0))
allocate(bq(m0,m0))
allocate(x(n,m0))
allocate(zwork(n,m0))

!fpm(1)=1 !turn on FEAST output
!fpm(4)=10 !maximum FEAST iterations
!fpm(2)=4 !number of FEAST contour points, i.e. number of linear systems to solve at each FEAST iteration
!fpm(3)=6 !accuracy of FEAST iterations, i.e. stop when residual is less than 10^-fpm(3)

allocate(reslist(0:fpm(4)),timelist(0:fpm(4)))

!eigenvalue interval:
!emin=0.5d0
!emax=10.5d0

!!!!!!!!!!!!!  Set up linear system solver:
matdescra(1)='S'
matdescra(2)='L'
matdescra(3)='N'
matdescra(4)='F'
linIterations=3 !number of iterations to use; more linear system iterations allows FEAST to converge faster
kdim=2 !number of Krylov subspace blocks to use; 1 usually works fine, but more gives faster convergence
allocate(xsol(n,m0), V(n,m0*kdim), Av(n,m0*kdim), Ax(n,m0), zwork2(n,m0*kdim), linwork1(n,m0), linwork2(n,m0))


if (fpm(11)==0) then


call dfeast_scsrev_timed(UPLO,n,dsa,isa,jsa,fpm,epsout,loop,emin,emax,m0,e,x,m,res,info,reslist,timelist)
!if (UPLO .ne. 'L') then

!call dfeast_scsrev_timed2(UPLO,n,dsa,isa,jsa,fpm,epsout,loop,emin,emax,m0,e,x,m,res,info,reslist,timelist)
!else
!call dfeast_scsrev_timed2('U',n,ssa,sisa,sjsa,fpm,epsout,loop,emin,emax,m0,e,x,m,res,info,reslist,timelist)
!end if

else

!!!!!!!!!!!!  Run FEAST RCI:
call system_clock(count=totalc1)
ijob=-1
oldloop=0
do while (ijob .ne. 0)
call dfeast_srci(ijob,n,ze,work,zwork,aq,bq,fpm,epsout,loop,emin,emax,m0,e,x,m,res,info)

if (oldloop .ne. loop) then
    reslist(oldloop)=epsout
    call system_clock(count=totalc2)
    timelist(oldloop)=elapsed_time(totalc1,totalc2)
    oldloop=loop
end if


select case (ijob)
		case (10) !store complex shift ze
		
            ze2=ze	
		
        case (11) !solve linear system Az*Qz=zwork for Qz, put Qz in zwork
            !ztmpmat3(1:n,1:n)=0.0d0
            !do i=1,n
            !    ztmpmat3(i,i)=ze2
            !end do
            !ztmpmat3=ztmpmat3-A
            !ztmpmat1(1:n,1:m0)=zwork(1:n,1:m0)
			
            linState(1)=-1 !set linsate(1)=-1 to begin RCI routine dfeast_gmres
            do while (linState(1) .ne. 0)
                !print *,'linstate',linstate
                !call dfeast_gmres(linjob,linState,zwork,xsol,V,Av,Ax,ze2,n,m0,linIterations,kdim,linwork1,linwork2,zwork2)                
                !print *,'linstate',linstate
                call dfeast_gmres(linjob,linState,zwork,xsol,V,Av,Ax,ze2,n,m0,linIterations,kdim,linwork1,linwork2,zwork2)
                if(linjob == 30) then
                    !Your matrix multiplication routine goes here! do linwork2=A*linwork1
                    !call dgemm('N','N',n,m0,n,1.0d0,A,n,linwork1,n,0.0d0,linwork2,n)
                    !print *,'doing matmul'
                    call system_clock(count=c1)
                    !call dcsrmm(UPLO,'N',N,N,m0,1.0d0,dsa,isa,jsa,linwork1,0.0d0,linwork2)
                    call mkl_dcsrmm('N', n, m0, n, 1.0d0, matdescra, dsa, jsa, isa, isa(2), linwork1, n, 0.0d0, linwork2, n)

                    !if (UPLO .ne. 'L') then
                    !    call mkl_dcsrmm('N', n, m0, n, 1.0d0, matdescra, dsa, jsa, isa, isa(2), linwork1, n, 0.0d0, linwork2, n)
                    !else
                    !    call mkl_dcsrmm('N', n, m0, n, 1.0d0, matdescra, ssa, sjsa, sisa, sisa(2), linwork1, n, 0.0d0, linwork2, n)
                    !end if
                    call system_clock(count=c2)
                    !print *,'done matmul', elapsed_time(c1,c2)
                end if

                if(linjob==50) then
                    linwork2(1:n,1:m0)=linwork1(1:n,1:m0)
                end if

            end do
            !print *, 'done lin sys'	
            !ztmpmat2(1:n,1:m0)=ztmpmat1(1:n,1:m0)-matmul(ztmpmat3,zwork(1:n,1:m0))
            !do i=1,m0
            !    ztmpmat2(1:n,i)=ztmpmat2(1:n,i)/dznrm2(n,ztmpmat1(1:n,i),1)
            !end do

            !dtmp1=dznrm2(n,ztmpmat2(1:n,1),1)
            !do i=1,m
            !    dtmp2=dznrm2(n,ztmpmat2(1:n,i),1)
            !    if(dtmp2>dtmp1) then
            !        dtmp1=dtmp2
            !    end if
            !end do

            !print *,'     Lin sys error=',dtmp1

            !stop

        case (30) !A*x	
           !print *,'case30' 
            !Your matrix multiplication routine goes here! Do work=A*x
            !work(:,:)=matmul(A,x(:,:))
            !call dgemm('N','N',n,m0,n,1.0d0,A,n,x,n,0.0d0,work,n)
            !call dcsrmm(UPLO,'N',N,N,m0,1.0d0,dsa,isa,jsa,x,0.0d0,work)
            call system_clock(count=c1)
            call mkl_dcsrmm('N', n, m0, n, 1.0d0, matdescra, dsa, jsa, isa, isa(2), x, n, 0.0d0, work, n)
            call system_clock(count=c2)
            !print *,'done matmul2', elapsed_time(c1,c2)
		case (40) !B*x
        
            call dlacpy('F',n,m0,x,n,work,n)

	end select
end do
call system_clock(count=totalc2)

reslist(oldloop)=epsout
timelist(oldloop)=elapsed_time(totalc1,totalc2)
end if


open(unit=1,file="dfeast_out.dat")
do i=0,loop
    write(1,*) i,timelist(i),reslist(i)
end do
close(1)

!print results:
print *,'FEAST finished; # eigs found = ',m
print *, 'eigenvalues and eigenvector residuals:'
if (fpm(14)==0) then
do i=1,m
    print *,i,e(i),res(i)
end do
end if

if (fpm(14)==2) then
    do i=1,m0
        print *,res(i)
    end do
end if

print *,'total time=',elapsed_time(totalc1,totalc2)

!deallocate(A)
!deallocate(ztmpmat1,ztmpmat2,ztmpmat3)
deallocate(res,e,work,aq,bq,x,zwork)
deallocate(xsol, V, Av, Ax, zwork2, linwork1, linwork2)

deallocate(reslist,timelist)

deallocate(dsa,dca,isa,jsa,ica,jca)

!if (UPLO=='L') then
!    deallocate(ssa,sjsa,sisa)
!end if
end program



subroutine dcoo2csr(n,nnz,ic,jc,c,isa,jsa,sa)
  implicit none
  integer :: n,nnz
  integer,dimension(*) :: ic,jc,isa,jsa
  double precision,dimension(*) :: c,sa
!!!
  integer :: k,k1,i,j,idum
  integer,dimension(n) :: iloc
  double precision :: adum


  isa(1:N+1) = 0
  !find how many elements in row
  do  k=1, nnz
     isa(ic(k)) = isa(ic(k))+1
  end do



  !build isa
  k = 1
  do  i=1,N+1
     k1 = isa(i)
     isa(i) = k
     k = k+k1
  end do


  iloc(1:n)=isa(1:n)
  !Build jsa, sa - increment local row counter
  do  k=1, nnz
     sa(iloc(ic(k))) =  c(k)
     jsa(iloc(ic(k))) = jc(k)
     iloc(ic(k)) = iloc(ic(k))+1
  end do
  ! Reorder by increasing column
  do i=1,n
     do k=isa(i),isa(i+1)-1
        do k1=k,isa(i+1)-1
           if (jsa(k1)<jsa(k)) then
              idum=jsa(k)
              jsa(k)=jsa(k1)
              jsa(k1)=idum
              adum=sa(k)
              sa(k)=sa(k1)
              sa(k1)=adum
           endif
        enddo
     enddo
  enddo


end subroutine dcoo2csr

subroutine dfeast_scsrev_timed(UPLO,N,sa,isa,jsa,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,reslist,timelist)
    implicit none
    character(len=1) :: UPLO
    integer :: N
    double precision,dimension(*) :: sa
    integer,dimension(*) :: isa,jsa
    integer,dimension(*) :: fpm
    double precision :: epsout 
    integer :: loop
    double precision :: Emin,Emax
    integer :: M0
    double precision,dimension(*)  :: E
    double precision,dimension(N,*):: X
    integer :: mode
    double precision,dimension(*)    :: res
    integer :: info

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!    Wrapper Routine to expert routine: dfeast_scsrgvx
!!!    Create Zne, Wne arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    complex(kind=(kind(1.0d0))),dimension(fpm(2)) :: Zne,Wne 
    double precision,dimension(n) :: sb ! identity
    integer,dimension(n+1) :: isb
    integer,dimension(n) :: jsb
    integer :: i

!!!!!!!!!!!!!!!!!!!!!!!! timing
    double precision, dimension(0:fpm(4)) :: reslist,timelist

    call zfeast_contour(Emin,Emax,fpm(2),fpm(16),fpm(18),Zne,Wne)

    fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
    ! identity B matrix- option for standard eigenvalue problem
    do i=1,n
       sb(i)=1.0d0
       jsb(i)=i
       isb(i)=i
    enddo
    isb(n+1)=n+1

    call dfeast_scsrgvx_timed(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne,reslist,timelist)

  end subroutine dfeast_scsrev_timed



  subroutine dfeast_scsrgvx_timed(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne,reslist,timelist)
    !  Purpose 
    !  =======
    !  FEAST SPARSE INTERFACE
    !  Solve the generalized Ax=eBx eigenvalue problem
    !  
    !  A REAL SYMMETRIC, B SYMMETRIC POSITIVE DEFINITE:: SPARSE FORMAT 
    ! 
    !  DOUBLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
    !                           triangular part of the matrix(ces) is(are) being supplied.
    !  N          (input)        INTEGER: Size system
    !  sa         (input)        REAL DOUBLE PRECISION (isa(N+1)-1):  Matrix A- CSR format 
    !  isa        (input)        INTEGER(N+1): CSR row array of Matrix A
    !  jsa        (input)        INTEGER(isa(N+1)-1): CSR column array of Matrix A
    !  sb         (input)        REAL DOUBLE PRECISION (isb(N+1)-1):  Matrix B- CSR format 
    !  isb        (input)        INTEGER(N+1): CSR row array of Matrix B
    !  jsb        (input)        INTEGER(isb(N+1)-1): CSR column array of Matrix B
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  fpm (input/output) INTEGER(*) : FEAST parameters
    !  epsout     (output)       REAL DOUBLE PRECISION : Error on the trace
    !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
    !  Emin,Emax  (input)        REAL DOUBLE PRECISION: search interval
    !  M0         (input/output) INTEGER: Size subspace
    !  E          (output)       REAL DOUBLE PRECISION(M0)   : Eigenvalues -solution
    !  q          (input/output) REAL DOUBLE PRECISION(N,M0) : 
    !                                                       On entry: subspace initial guess if fpm(5)=1 
    !                                                       On exit : Eigenvectors-solution
    !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
    !  res        (output)       REAL DOUBLE PRECISION(M0) : Relative Residual of the solution (1-norm)
    !                                                        if option fpm(6)=1 selected                           
    !  info       (output)       INTEGER: Error handling (0: successful exit)
    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of ARGUMENTS FOR EXPERT ROUTINE 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  Zne        (input)        COMPLEX DOUBLE PRECISION(fpm(2)): Custom Integration nodes
    !  Wne        (input)        COMPLEX DOUBLE PRECISION(fpm(2)): Custom Integration weights
    !           
    !=====================================================================
    ! Eric Polizzi 2009-2015
    !=====================================================================
    implicit none
    include 'f90_noruntime_interface.fi'
    character(len=1) :: UPLO
    integer :: N
    double precision,dimension(*) :: sa,sb
    integer,dimension(*) :: isa,jsa,isb,jsb
    integer,dimension(*) :: fpm
    double precision :: epsout 
    integer :: loop
    double precision :: Emin,Emax
    integer :: M0
    double precision,dimension(*)  :: E
    double precision,dimension(N,*):: X
    integer :: mode
    double precision,dimension(*)    :: res
    integer :: info
    complex(kind=(kind(1.0d0))),dimension(*) :: Zne,Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer :: ijob,infoloc,i,s
    complex(kind=(kind(1.0d0))) :: Ze
    complex(kind=(kind(1.0d0))), dimension(:,:),pointer ::workc,caux
    double precision, dimension(:,:),pointer ::work,Aq,Sq
    complex(kind=(kind(1.0d0))),dimension(:,:),pointer :: saz
    integer,dimension(:),pointer :: isaz,jsaz
    double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
    complex(kind=(kind(1.0d0))),parameter :: ONEC=(DONE,DZERO)
    logical :: fact
    integer :: nfact,id
    integer :: rank,code,nb_procs,NEW_COMM_WORLD
!!!!! csr-upper format
    double precision,dimension(:),pointer :: ssa,ssb
    integer,dimension(:), pointer :: sisa,sjsa,sisb,sjsb
    integer :: opt,nnza,nnzb,nnz
!!!!!for pardiso
    integer(8),dimension(64) :: pt
    integer,dimension(64) :: iparm
    integer :: mtype
    integer :: MAXFCT,MNUM,PHASE,MSGLVL
    integer :: idum

!!!!!!!!!!!!!! for timing
    double precision, dimension(0:fpm(4)) :: reslist,timelist
    integer :: totalc1,totalc2,c1,c2
    integer :: oldloop
    double precision, external :: elapsed_time

    character, dimension(6) :: matdescra

    call system_clock(count=totalc1)

    matdescra(1)='S'
    matdescra(2)='U'
    matdescra(3)='N'
    matdescra(4)='F'

    rank=0
    nb_procs=1
    !----------------------------------------------
#ifdef MPI
    NEW_COMM_WORLD=fpm(9)
    call MPI_COMM_RANK(NEW_COMM_WORLD,rank,code)
    call MPI_COMM_SIZE(NEW_COMM_WORLD,nb_procs,code)
#endif
    !----------------------

    INFO = 0
    IF ((UPLO/='F').and.(UPLO/='f').and.(UPLO/='L').and.(UPLO/='l').and.(UPLO/='U').and.(UPLO/='u')) THEN
       INFO=-101
    ELSE IF ( N<=0 ) THEN
       INFO = -102
    END IF
    IF( INFO.NE.0 ) THEN
       CALL XERBLA( 'DFEAST_SCSRGV', -INFO+100 )
       RETURN
    END IF


    infoloc=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! FORMAT CONVERSION TO CSR-UPPER for PARDISO !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if ((UPLO=='F').or.(UPLO=='f')) then !!! full csr to upper-csr
       call wallocate_1i(sisa,n+1,infoloc)
       call wallocate_1i(sisb,n+1,infoloc)

       !!<<<
       call wallocate_1d(ssa,1,infoloc) ! dummy
       call wallocate_1i(sjsa,1,infoloc) !dummy
       call wallocate_1d(ssb,1,infoloc) !dummy
       call wallocate_1i(sjsb,1,infoloc) !dummy
       if (infoloc/=0) then
          info=-1
          return
       end if
       !!>>>
       opt=1
       call dcsr2csr_up(opt,N,sa,isa,jsa,ssa,sisa,sjsa)
       nnza=sisa(n+1)-1 
       call dcsr2csr_up(opt,N,sb,isb,jsb,ssb,sisb,sjsb)
       nnzb=sisb(n+1)-1 
       !!<<<
       call wdeallocate_1d(ssa)
       call wdeallocate_1i(sjsa)
       call wdeallocate_1d(ssb)
       call wdeallocate_1i(sjsb)
       !!>>>

       call wallocate_1d(ssa,nnza,infoloc)
       call wallocate_1i(sjsa,nnza,infoloc)
       call wallocate_1d(ssb,nnzb,infoloc)
       call wallocate_1i(sjsb,nnzb,infoloc)
       if (infoloc/=0) then
          info=-1
          return
       end if

       opt=2
       call dcsr2csr_up(opt,N,sa,isa,jsa,ssa,sisa,sjsa)
       call dcsr2csr_up(opt,N,sb,isb,jsb,ssb,sisb,sjsb)

!!!!!!!!!!!!!!!!!!!!!!!!!
    elseif ((UPLO=='U').or.(UPLO=='u')) then !! upper-csr already 

       !       nnza=isa(n+1)-1
       !       ssa => sa(1:nnza)
       !       sisa => isa(1:n+1)
       !       sjsa => jsa(1:nnza)

       !       nnzb=isb(n+1)-1
       !       ssb =>  sb(1:nnzb)
       !       sisb => isb(1:n+1)
       !       sjsb => jsb(1:nnzb)


    elseif ((UPLO=='L').or.(UPLO=='l')) then !!! lower-csr to upper-csr

       nnza=isa(n+1)-1
       nnzb=isb(n+1)-1
       call wallocate_1d(ssa,nnza,infoloc)
       call wallocate_1i(sjsa,nnza,infoloc)
       call wallocate_1d(ssb,nnzb,infoloc)
       call wallocate_1i(sjsb,nnzb,infoloc)
       call wallocate_1i(sisa,n+1,infoloc)
       call wallocate_1i(sisb,n+1,infoloc)
       if (infoloc/=0) then
          info=-1
          return
       end if

       call dcsr_transpose(N,sa,isa,jsa,ssa,sisa,sjsa)
       call dcsr_transpose(N,sb,isb,jsb,ssb,sisb,sjsb)

    endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



    call wallocate_2d(Aq,M0,M0,infoloc)
    call wallocate_2d(Sq,M0,M0,infoloc)
    call wallocate_2d(work,N,M0,infoloc)
    call wallocate_2z(workc,N,M0,infoloc)
    if (infoloc/=0) then
       info=-1
       return
    end if


!!!!!!!!!!!!!!!!!! Factorizations Set-up  
    fact=.true.
    nfact=0
    !! nfact is local (number of total factorization by contour points)
    do i=rank+1,fpm(2),nb_procs
       nfact=nfact+1
    end do

    if (fpm(10)==1) then
       id=0
    else
       id=1
    end if

!!!!!!!!!!!!!!!!! Set up for Az matrix

    call wallocate_1i(isaz,n+1,infoloc)
    call wallocate_2z(saz,1,1,infoloc) ! dummy
    call wallocate_1i(jsaz,1,infoloc)! dummy
    !!>>>
    opt=1
    if ((UPLO=='U').or.(UPLO=='u')) then
       call zdaddcsr(N,opt,ONEC,sa,isa,jsa,ONEC,sb,isb,jsb,saz,isaz,jsaz) !! get isaz
    else
       call zdaddcsr(N,opt,ONEC,ssa,sisa,sjsa,ONEC,ssb,sisb,sjsb,saz,isaz,jsaz) !! get isaz
    end if

    nnz=isaz(n+1)-1
    !!<<<
    call wdeallocate_2z(saz)
    call wdeallocate_1i(jsaz)
    !!>>>
    call wallocate_2z(saz,nnz,nfact,infoloc)
    call wallocate_1i(jsaz,nnz,infoloc)

    opt=2
    if ((UPLO=='U').or.(UPLO=='u')) then
       call zdaddcsr(N,opt,ONEC,sa,isa,jsa,ONEC,sb,isb,jsb,saz,isaz,jsaz) !! get jsaz
    else
       call zdaddcsr(N,opt,ONEC,ssa,sisa,sjsa,ONEC,ssb,sisb,sjsb,saz,isaz,jsaz) !! get jsaz
    end if
!!!!!!!!!!!!!!!
    call wallocate_2z(caux,N,M0,infoloc)


    if (infoloc/=0) then
       info=-1
       return
    end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  pardiso initialization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    MAXFCT=nfact ! Rq: same factorization for (normal+transpose)
    MTYPE=6      ! complex and symmetric 
    call pardisoinit(PT,MTYPE,IPARM)
!!!!!!!!!!!!
    if (fpm(64)==1) then
       do i=1,64
          if (fpm(64+i)/=-9999) iparm(i)=fpm(64+i)
       enddo
    endif
!!!!!!!!!!!!
    IPARM(6)=1 ! solution and rhs are input/output, attention caux is always used
    MSGLVL=0!0 !0- no output, 1- output
    PHASE=11 ! symbolic factorization (do it only once)
    !    call PARDISO(PT,MAXFCT,1,MTYPE,PHASE,n,saz(1,1),isaz,jsaz,idum,M0,IPARM,MSGLVL,workc,caux,infoloc)
!!!!!!!!!!


    ijob=-1 ! initialization
    oldloop=0
    do while (ijob/=0)
       call dfeast_srcix(ijob,N,Ze,work,workc,Aq,Sq,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)    

       if (oldloop .ne. loop) then
            reslist(oldloop)=epsout
            call system_clock(count=totalc2)
            timelist(oldloop)=elapsed_time(totalc1,totalc2)
            print *,'time=',timelist(oldloop)
            oldloop=loop
        end if
     

       select case(ijob)
       case(10) !! factorize (zeB-A)
            call system_clock(count=c1)
            

          if (fpm(10)==1) then
             id=id+1
             if (id==nfact+1) then
                id=1
                fact=.false.
             endif
          endif

          if (fact) then
             opt=3
             if ((UPLO=='U').or.(UPLO=='u')) then
                call zdaddcsr(N,opt,-ONEC,sa,isa,jsa,Ze,sb,isb,jsb,saz(1,id),isaz,jsaz) !! get saz
             else
                call zdaddcsr(N,opt,-ONEC,ssa,sisa,sjsa,Ze,ssb,sisb,sjsb,saz(1,id),isaz,jsaz) !! get saz
             end if
             if (PHASE==11) then
                PHASE=12
             else
                PHASE=22
             endif
             !             PHASE=22 ! numerical fact only
             MNUM=id
             call PARDISO(PT,MAXFCT,MNUM,MTYPE,PHASE,n,saz(1,id),isaz,jsaz,idum,fpm(23),IPARM,MSGLVL,workc,caux,infoloc)
             if (infoloc/=0) then
                info=-2
                return
             end if

          endif ! fact true
            call system_clock(count=c2)
            print *,'   factorize time=',elapsed_time(c1,c2)
       case(11) !!solve the linear system (ZeB-A)x=workc(1:N,1:fpm(23)) result in to workc
          call system_clock(count=c1)
          
          PHASE=33 ! solve
          MNUM=id
          call PARDISO(PT,MAXFCT,MNUM,MTYPE,PHASE,n,saz(1,id),isaz,jsaz,idum,fpm(23),IPARM,MSGLVL,workc,caux,infoloc)        
          if (infoloc/=0) then
             info=-2
             return
          end if
            call system_clock(count=c2)
            print *,'   solve time=',elapsed_time(c1,c2)
       case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)
            call system_clock(count=c1)
          print *,'begining matmul'  
          if ((UPLO=='U').or.(UPLO=='u')) then
             !call dcsrmm('U','N',N,N,fpm(25),DONE,sa,isa,jsa,X(1,fpm(24)),DZERO,work(1,fpm(24)))
             call mkl_dcsrmm('N', n, fpm(25), n, 1.0d0, matdescra, sa, jsa, isa, isa(2), X(1,fpm(24)), n, 0.0d0, work(1,fpm(24)), n)
          else
             !call dcsrmm('L','N',N,N,fpm(25),DONE,sa,isa,jsa,X(1,fpm(24)),DZERO,work(1,fpm(24)))
             !call dcsrmm('U','N',N,N,fpm(25),DONE,ssa,sisa,sjsa,X(1,fpm(24)),DZERO,work(1,fpm(24)))
             call mkl_dcsrmm('N', n, fpm(25), n, 1.0d0, matdescra, ssa, sjsa, sisa, sisa(2), X(1,fpm(24)), n, 0.0d0, work(1,fpm(24)), n)
          end if
            call system_clock(count=c2)
            print *,'   mm time=',elapsed_time(c1,c2),fpm(24),m0
       case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

          if ((UPLO=='U').or.(UPLO=='u')) then
             call dcsrmm('U','N',N,N,fpm(25),DONE,sb,isb,jsb,X(1,fpm(24)),DZERO,work(1,fpm(24)))
          else
             !call dcsrmm('L','N',N,N,fpm(25),DONE,sb,isb,jsb,X(1,fpm(24)),DZERO,work(1,fpm(24)))
             call dcsrmm('U','N',N,N,fpm(25),DONE,ssb,sisb,sjsb,X(1,fpm(24)),DZERO,work(1,fpm(24)))
          end if


       end select
    end do

!!!!!!!!!!!!!!!!!!!!!!!
!!!! release memory
!!!!!!!!!!!!!!!!!!!!!!!!
    PHASE=-1 

    do MNUM=1,nfact
       call PARDISO(PT,MAXFCT,MNUM,MTYPE,PHASE,n,saz(1,MNUM),isaz,jsaz,idum,fpm(23),IPARM,MSGLVL,workc,caux,infoloc)
       if (infoloc/=0) then
          info=-2
          return
       end if
    end do

    call system_clock(count=totalc2)

    reslist(oldloop)=epsout
    timelist(oldloop)=elapsed_time(totalc1,totalc2)

    call wdeallocate_2d(Aq)
    call wdeallocate_2d(Sq)
    call wdeallocate_2d(work)
    call wdeallocate_2z(workc)
    call wdeallocate_2z(caux)


    call wdeallocate_2z(saz)
    call wdeallocate_1i(isaz)
    call wdeallocate_1i(jsaz)


    if ((UPLO/='U').and.(UPLO/='u')) then
       call wdeallocate_1d(ssa)
       call wdeallocate_1i(sisa)
       call wdeallocate_1i(sjsa)

       call wdeallocate_1d(ssb)
       call wdeallocate_1i(sisb)
       call wdeallocate_1i(sjsb)
    endif

  end subroutine dfeast_scsrgvx_timed



