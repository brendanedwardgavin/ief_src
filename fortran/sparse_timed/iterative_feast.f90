  subroutine dfeast_scsrevit(UPLO,N,sa,isa,jsa,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
    !  Purpose 
    !  =======
    !  FEAST SPARSE INTERFACE
    !  Solve the standard Ax=ex eigenvalue problem
    !  
    !  A REAL SYMMETRIC:: SPARSE FORMAT 
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
    !=====================================================================
    ! Eric Polizzi 2009-2015
    ! ====================================================================
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

    call zfeast_contour(Emin,Emax,fpm(2),fpm(16),fpm(18),Zne,Wne)

    fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
    ! identity B matrix- option for standard eigenvalue problem
    do i=1,n
       sb(i)=1.0d0
       jsb(i)=i
       isb(i)=i
    enddo
    isb(n+1)=n+1

    call dfeast_scsrgvxit(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)

  end subroutine dfeast_scsrevit



subroutine dfeast_scsrgvxit(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)
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
    use rundata
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
double precision,dimension(:),pointer :: ssa,ssb,nres
double precision :: ares
logical :: comb
integer :: infob,itmax

    integer,dimension(:), pointer :: sisa,sjsa,sisb,sjsb
    integer :: opt,nnza,nnzb,nnz
!!!!!for pardiso
    integer(8),dimension(64) :: pt
    integer,dimension(64) :: iparm
    integer :: mtype
    integer :: MAXFCT,MNUM,PHASE,MSGLVL
    integer :: idum

    !!!!!!!!!! iterative solver
    complex(kind=(kind(1.0d0))) :: ze2 
    complex(kind=(kind(1.0d0))), dimension(:,:), pointer :: ztempmat
    complex(kind=(kind(1.0d0))),dimension(:),pointer :: zsa
    integer :: linresindex !number of RHS to use in measuring linear system error
    integer :: blockits,blockremainder
    double precision :: lintargeterror,linsyserror !goal error for linear system, return linear system error
    integer :: linloops
    double precision :: linepsout
    integer :: k

    integer :: oldloop

    !!!BLAS and lapack:
    character, dimension(6) :: matdescra

    if(UPLO=='F') then
        matdescra(1)='G'
    else
        matdescra(1)='S'
    end if
    !matdescra(1)='G'
    matdescra(2)=UPLO
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

    if (fpm(11)>0) then
        call wallocate_1z(zsa,nnza,infoloc)
        call wallocate_2z(ztempmat,n,m0,infoloc)
        zsa=(1.0d0,0.0d0)*sa(1:nnza)

call wallocate_1d(nres,M0,infoloc) ! dummy
    end if

    ijob=-1 ! initialization

    !!!!!!!!!! start keeping track of stuff
    cpnum=0
    oldloop=0
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do while (ijob/=0)
       call dfeast_srcix(ijob,N,Ze,work,workc,Aq,Sq,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)    
       
        feastit=loop
        if(oldloop .ne. loop) then !new iteration
            call system_clock(count=tc1)
            eigtime(oldloop)=elapsed_time(startcount,tc1)
            eigres(oldloop)=epsout
            !call quicksort(E,1,m0)
            ritzvals(oldloop,1:m0)=E(1:m0)
            eigresall(oldloop,1:m0)=res(1:m0)
            oldloop=loop
        end if
        
        !print *,epsout
       select case(ijob)
       case(10) !! factorize (zeB-A)
          
          !!!!!! keep track of contour points
          if (cpnum<fpm(2)) then
              cpnum=cpnum+1
          else
              cpnum=1
          end if
          cpval(cpnum)=ze
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          if (fpm(11)==0) then

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

          else !user iterative solver, don't need factorization

             ze2=Ze

          end if

       case(11) !!solve the linear system (ZeB-A)x=workc(1:N,1:fpm(23)) result in to workc

        if(fpm(11)==0) then
              PHASE=33 ! solve
              MNUM=id
              call PARDISO(PT,MAXFCT,MNUM,MTYPE,PHASE,n,saz(1,id),isaz,jsaz,idum,fpm(23),IPARM,MSGLVL,workc,caux,infoloc)        
              if (infoloc/=0) then
                 info=-2
                 return
              end if
        else
              !call zfeast_cgls(UPLO,n,m0,zsa,isa,jsa,ze2,nnza,workc,ztempmat,fpm(50)) 
              if (mode>0 .and. mode<=m0) then
                  linresindex=mode
              else
                  linresindex=m0
              end if

              if(epsout <=1.0d-1 .and. loop>0) then
                  lintargeterror=1d-2*epsout
              else
                  lintargeterror=1d-1
              endif          

              !if (loop<3) linresindex=m0
              !if(fpm(11)==2) then !block CGLS
              !  call zfeast_cglsRes(UPLO,n,m0,zsa,isa,jsa,ze2,nnza,workc,ztempmat,fpm(50),lintargeterror,linresindex,linsyserror,linloops) 
                !linit(feastit,cpnum,1+(i-1)*fpm(53):m0)=linloops
              !end if

              !if(fpm(11)==3) then !single vector BICGSTAB
                !zsa(1:nnz)=-(1.0d0,0.0d0)*sa(1:nnz)
                !do  i=1,n
                !do k=isa(i),isa(i+1)-1
                !if (jsa(k)==i) zsa(k)=zsa(k)+Ze2
                !enddo
                !enddo

                !do i=1,m0
                !call zfeast_BiCGSTABRes(UPLO,n,zsa,isa,jsa,ze2,nnza,workc(:,i),ztempmat(:,i),fpm(50),lintargeterror,linresindex,linsyserror) 
               ! end do

              !  ztempmat(1:N,1:M0)=(0.0d0,0.0d0)
              !   comb=.true.
              !   itmax=fpm(50)
              !   call zbicgstab(UPLO,N,zsa,isa,jsa,M0,workc,ztempmat,nres,ares,itmax,lintargeterror,comb,infob) 
              !      print *,'lin sys error=',ares
              !     print *,'lin sys it = ',itmax
              !end if

              !if(fpm(11)==1) then
              !  zsa(1:nnz)=-(1.0d0,0.0d0)*sa(1:nnz)
              !  do  i=1,n
              !  do k=isa(i),isa(i+1)-1
              !  if (jsa(k)==i) zsa(k)=zsa(k)+Ze2
              !  enddo
              !  enddo
              !end if

                !call blockGMRESarnoldi(UPLO,n,m0,zsa,isa,jsa,fpm(51),fpm(50),workc,ztempmat,lintargeterror) 
              if(fpm(11)>0) then
                !use fpm(53) as block size
                blockits=m0/fpm(53)
                blockremainder=m0-blockits*fpm(53)
                
                do i=1,blockits

                    linresindex=fpm(53)
                    if(fpm(11)==1) call blockGMRESarnoldi(UPLO,n,fpm(53),zsa,isa,jsa,ze2,fpm(51),fpm(50),workc(1,1+(i-1)*fpm(53)),ztempmat(1,1+(i-1)*fpm(53)),lintargeterror,linloops,1+(i-1)*fpm(53)) 
                    if(fpm(11)==2) call zfeast_cglsRes(UPLO,n,fpm(53),zsa,isa,jsa,ze2,nnza,workc(1,1+(i-1)*fpm(53)),ztempmat(1,1+(i-1)*fpm(53)),fpm(51)*fpm(50),lintargeterror,linresindex,linsyserror,linloops)
                    !workc(1+(i-1)*fpm(53):i*fpm(53),1:n)=ztempmat(1:fpm(53),1:n)
                    linit(feastit,cpnum,1+(i-1)*fpm(53):i*fpm(53))=linloops
                    !print *,'block ',i,linloops

                end do
                !print *,'i ',i,blockits
                if(blockremainder .ne. 0) then
                    linresindex=blockremainder
                    if(fpm(11)==1) call blockGMRESarnoldi(UPLO,n,blockremainder,zsa,isa,jsa,ze2,fpm(51),fpm(50),workc(1,1+(i-1)*fpm(53)),ztempmat(1,1+(i-1)*fpm(53)),lintargeterror,linloops,1+(i-1)*fpm(53))
                    if(fpm(11)==2) call zfeast_cglsRes(UPLO,n,blockremainder,zsa,isa,jsa,ze2,nnza,workc(1,1+(i-1)*fpm(53)),ztempmat(1,1+(i-1)*fpm(53)),fpm(51)*fpm(50),lintargeterror,linresindex,linsyserror,linloops)
                    linit(feastit,cpnum,1+(i-1)*fpm(53):m0)=linloops
                end if
              end if

              print *, 'lin sys',cpnum,sum(linit(feastit,cpnum,1:m0))/m0
              workc=ztempmat
        end if

       case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)
          !matdescra(2)='U'
          if ((UPLO=='U').or.(UPLO=='u')) then
          !   call dcsrmm('U','N',N,N,fpm(25),DONE,sa,isa,jsa,X(1,fpm(24)),DZERO,work(1,fpm(24)))
          !    call mkl_dcsrmm('N', n, m0, n, 0.0d0, matdescra, sa, jsa, isa, isa(2), X, n, 0.0d0, work, n)
          else
          !   call dcsrmm('U','N',N,N,fpm(25),DONE,ssa,sisa,sjsa,X(1,fpm(24)),DZERO,work(1,fpm(24)))
          !    call mkl_dcsrmm('N', n, m0, n, 0.0d0, matdescra, ssa, sjsa, sisa, sisa(2), X, n, 0.0d0, work, n)
          end if

            call mkl_dcsrmm('N', n, m0, n, 1.0d0, matdescra, sa, jsa, isa, isa(2), X, n, 0.0d0, work, n)

        case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

          !if ((UPLO=='U').or.(UPLO=='u')) then
          !   call dcsrmm('U','N',N,N,fpm(25),DONE,sb,isb,jsb,X(1,fpm(24)),DZERO,work(1,fpm(24)))
          !else
          !   call dcsrmm('U','N',N,N,fpm(25),DONE,ssb,sisb,sjsb,X(1,fpm(24)),DZERO,work(1,fpm(24)))
          !end if
          work(1:n,1:m0)=X(1:n,1:m0)
          !call mkl_zcsrmm('N', n, fpm(25), n, (1.0d0,0.0d0), matdescra, sb, jsb, isb, isb(2), X(1,fpm(24)), n, (0.0d0,0.0d0), work(1,fpm(24)), n)


       end select
    end do

    call system_clock(count=tc1)
    eigtime(loop)=elapsed_time(startcount,tc1)
    eigres(loop)=epsout
    !call quicksort(E,1,m0)
    ritzvals(loop,1:m0)=E(1:m0)
    eigresall(loop,1:m0)=res(1:m0)

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

  end subroutine dfeast_scsrgvxit


