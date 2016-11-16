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
    double precision :: delta

    if(fpm(55)==0) then
        call zfeast_contour(Emin,Emax,fpm(2),fpm(16),fpm(18),Zne,Wne)
    else
        delta=(Emax-Emin)/fpm(2)
        do i=1,fpm(2)
            call zfeast_contour(Emin+(i-1)*delta,Emin+i*delta,1,fpm(16),fpm(18),Zne(i),Wne(i))
        end do
    end if

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

    if(fpm(54)==0) print *,'WARNING: using default 1d-2*epsout for linear system accuracy. Set feastparam(54)/=0 to change this.'

    do while (ijob/=0)
       call dfeast_myrcix(ijob,N,Ze,work,workc,Aq,Sq,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)    
       !call dfeast_srcix(ijob,N,Ze,work,workc,Aq,Sq,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)    

        feastit=loop
        if(oldloop .ne. loop) then !new iteration
            call system_clock(count=tc1)
            eigtime(oldloop)=elapsed_time(startcount,tc1)
            print *,''
            print *,'Time elapsed ',eigtime(oldloop)
            eigres(oldloop)=epsout
            !call quicksort(E,1,m0)
            !ritzvals(oldloop,1:m0)=E(1:m0)
            !eigresall(oldloop,1:m0)=res(1:m0)
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
          znesave=Zne(1:ncp)
          wnesave=Wne(1:ncp)
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
                  if(fpm(54)==0) then
                      lintargeterror=1.0d-2*epsout
                  else
                      !lintargeterror=epsout*10.0d0**(-1.0d0*fpm(54))
                      lintargeterror=epsout/(1.0d0*fpm(54))

                  end if
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
              if(cpnum==1) print *,''
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





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine dfeast_myrcix(ijob,N,Ze,work,workc,Aq,Sq,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info,Zne,Wne)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use rundata
!  Purpose 
  !  =======
  !  FEAST RCI (Reverse Communication Interfaces) - Includes option for custom integration nodes/weight
  !  Solve generalized Aq=lambda Bq eigenvalue problems
  !  
  !  A REAL SYMMETRIC, B SYMMETRIC POSITIVE DEFINITE (or B Identity) 
  !  DOUBLE PRECISION version 
  ! 
  !  Arguments
  !  =========
  !
  !  ijob       (input/output) INTEGER :: ID of the RCI
  !                            INPUT on first entry: ijob=-1 
  !                            OUTPUT Return values (0,10,20,21,30,40)-- see FEAST documentation
  !  N          (input)        INTEGER: Size system
  !  Ze         (output)       COMPLEX DOUBLE PRECISION        : Integration node
  !  work       (input/output) REAL DOUBLE PRECISION (N,M0)   :  Workspace 
  !  workc      (input/output) COMPLEX DOUBLE PRECISION (N,M0):  Workspace 
  !  Aq,Sq      (input/output) REAL DOUBLE PRECISION (M0,M0)  :  Workspace for Reduced Eigenvalue System
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  fpm     (input/output)    INTEGER(*) : FEAST parameters (see FEAST documentation)
  !  epsout     (output)       REAL DOUBLE PRECISION : Error on the trace
  !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
  !  Emin,Emax  (input)        REAL DOUBLE PRECISION: search interval
  !  M0         (input/output) INTEGER: Size subspace
  !  lambda     (output)       REAL DOUBLE PRECISION(M0)   : Eigenvalues -solution
  !  q          (input/output) REAL DOUBLE PRECISION(N,M0) : 
  !                                                       On entry: subspace initial guess if fpm(5)=1 
  !                                                       On exit : Eigenvectors-solution
  !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
  !  res        (output)       REAL DOUBLE PRECISION(M0) : Relative Residual of the solution (1-norm)
  !                                                               
  !  info       (output)       INTEGER: Error handling (0: successful exit) -- see FEAST documentation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of ARGUMENTS FOR EXPERT ROUTINE 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Zne        (input)        COMPLEX DOUBLE PRECISION(fpm(2)): Custom Integration nodes
  !  Wne        (input)        COMPLEX DOUBLE PRECISION(fpm(2)): Custom Integration weights
  !                
  !                            Expert comment: if fpm(29)=1 Zne, Wne have already been generated using default contour   
  !===========================================================================================================
  ! Eric Polizzi 2009-2015
  ! ====================================================================

  implicit none
  !-------------------------------------
#ifdef MPI
  include 'mpif.h'
#endif
  !-------------------------------------
  include "f90_noruntime_interface.fi"
  integer :: ijob,N,M0
  complex(kind=(kind(1.0d0))) :: Ze
  double precision, dimension(N,*) ::work
  complex(kind=(kind(1.0d0))),dimension(N,*):: workc
  integer,dimension(*) :: fpm
  double precision,dimension(M0,*):: Aq,Sq
  double precision :: epsout 
  integer :: loop
  double precision :: Emin,Emax
  double precision,dimension(*)  :: lambda
  double precision,dimension(N,*):: q
  integer :: mode
  double precision,dimension(*) :: res
  integer :: info
  complex(kind=(kind(1.0d0))),dimension(*) :: Zne,Wne
  !! parameters
  double precision, Parameter :: pi=3.1415926535897932d0
  double precision, Parameter :: DONE=1.0d0, DZERO=0.0d0
  complex(kind=(kind(1.0d0))),parameter :: ONEC=(DONE,DZERO)
  complex(kind=(kind(1.0d0))),parameter :: ZEROC=(DZERO,DZERO)
  double precision, parameter :: ba=-pi/2.0d0, ab=pi/2.0d0
  integer(8),parameter :: fout =6
  !! variable for FEAST
  integer :: i,e,j,k
  integer,dimension(4) :: iseed
  double precision :: theta,r,Emid
  complex(kind=(kind(1.0d0))) :: zxe,zwe,aux
  double precision, dimension(:,:),pointer :: Sqo
  integer, dimension(:),pointer :: fpm_default
  logical :: testconv
  double precision :: trace
  !! Lapack variable (reduced system)
  character(len=1) :: JOBZ,UPLO
  double precision, dimension(:),pointer :: work_loc,work_loc2
  integer :: lwork_loc,info_lap,infoloc
  !! MPI compatibility variables
  integer :: rank,code,nb_procs,NEW_COMM_WORLD

  rank=0
  nb_procs=1
  !----------------------------------------------
#ifdef MPI
  NEW_COMM_WORLD=fpm(9)
  call MPI_COMM_RANK(NEW_COMM_WORLD,rank,code)
  call MPI_COMM_SIZE(NEW_COMM_WORLD,nb_procs,code)
  if (rank/=0) fpm(1)=0 ! comment only in rank 0 if any
#endif
  !---------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! Initialization!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (ijob==-1) then 
     info=0 ! default value
     if (fpm(1)==1) then
        call wwrite_n(fout)
        call wwrite_s(fout, '***********************************************')  
        call wwrite_n(fout) 
        call wwrite_s(fout, '*********** FEAST- BEGIN **********************')
        call wwrite_n(fout) 
        call wwrite_s(fout, '***********************************************')  
        call wwrite_n(fout)
        call wwrite_s(fout, 'Routine DFEAST_S{}{}')
        if (fpm(29)==0) call wwrite_s(fout,'X') 
        call wwrite_n(fout)
!!!!!!!!!!!! Print the FEAST parameters which has been changed from default
        call wwrite_s(fout, 'List of input parameters fpm(1:64)-- if different from default')
        call wwrite_n(fout)

        call wallocate_1i(fpm_default,64,infoloc)
        call feastinit(fpm_default)
        do i=1,19
           if (fpm(i)/=fpm_default(i)) then
              call wwrite_s(fout, '   fpm(')
              call wwrite_i(fout, i)
              call wwrite_s(fout, ')=')
              call wwrite_i(fout, fpm(i))
              call wwrite_n(fout)
           endif
        enddo
        call wdeallocate_1i(fpm_default)

        call wwrite_s(fout, 'Search interval [')
        call wwrite_d(fout,Emin)
        call wwrite_s(fout, '; ')
        call wwrite_d(fout,Emax)
        call wwrite_s(fout, ']')
        call wwrite_n(fout)
     end if
     call check_feast_fpm_input(fpm,info)
     call dcheck_feast_srci_input(Emin,Emax,M0,N,info)

     if (info/=0) fpm(21)=100 ! The End
  
!!!!!!!!!!!!!!!!
     IF (info==0) then
        fpm(22)=fpm(2) ! only one half contour necessary
        loop=0
        if (fpm(1)==1) then
           call wwrite_s(fout, 'Size system    ')  
           call wwrite_t(fout) 
           call wwrite_i(fout,N)
           call wwrite_n(fout)
           call wwrite_s(fout, 'Size subspace  ')  
           call wwrite_t(fout) 
           call wwrite_i(fout,M0)
           call wwrite_n(fout)
           call wwrite_s(fout, '#Linear systems')  
           call wwrite_t(fout) 
           call wwrite_i(fout,fpm(22))
           call wwrite_n(fout)
           call wwrite_s(fout, '-----------------------------------------------------------------------------------')
           call wwrite_n(fout)
           call wwrite_s(fout, '#Loop | #Eig  |       Trace           |     Error-Trace       |     Max-Residual')  
           call wwrite_n(fout)
           call wwrite_s(fout, '-----------------------------------------------------------------------------------')
           call wwrite_n(fout)
        endif
        fpm(23)=min(M0,N) ! 'current M0' size (global value)
        fpm(25)=fpm(23) !! 'current M0' size (by default)
        fpm(24)=1  !! origin first column number of vector q for parallel mat-vec (default)
        !----------------------------------------------
#ifdef MPI
        if (fpm(23)/nb_procs>=1) then ! criteria for parallelism of mat-vec
           fpm(25)=fpm(23)/nb_procs ! local size of 'current M0'
           if (rank==nb_procs-1) fpm(25)=fpm(25)+fpm(23)-(fpm(23)/nb_procs)*nb_procs 
           fpm(24)=1+rank*(fpm(23)/nb_procs) ! local origin of first column for vector q for parallel mat-vec 
        end if
#endif
        !-----------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        fpm(21)=1 ! prepare reentry
        if (fpm(5)==0) then !!! random vectors (option 2 in DLARNV)
           iseed=(/56,890,3456,2333/)
           ! copy work into q for multiplication by B matrix, if stochastic estimate is "on" 
           if (fpm(14)==2) then
              call DLARNV(2,iseed,N*fpm(23),q(1,1)) 
           else
              call DLARNV(2,iseed,N*fpm(23),work(1,1)) 
           endif
        end if

        if ((fpm(5)==1).or.(fpm(14)==2)) then !!!!!! q is the initial guess
           !----------------------------------------
#ifdef MPI
           work(1:N,1:fpm(23))=DZERO 
#endif
           !------------------------------------------
           ijob=40 !! B*q=>work
           return
        end if
     end IF ! info=0
!!!!!!!!!!!!!!
  end if   !ijob=-1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! CONTOUR INTEGRATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  IF (fpm(21)==1) then !! we initialize a new contour integration
     !------------------------------------------------------------------------
#ifdef MPI
     if ((loop>0).or.(fpm(5)==1).or.(fpm(14)==2)) then        
        if (fpm(23)/nb_procs>=1) call MPI_ALLREDUCE(MPI_IN_PLACE,work,N*fpm(23),MPI_DOUBLE_PRECISION,MPI_SUM,NEW_COMM_WORLD,code)
     endif

     if ((fpm(29)==1).and.(fpm(16)==2)) then ! Zolotarev 
        if ((loop>0).and.(fpm(23)/nb_procs>=1)) call MPI_ALLREDUCE(MPI_IN_PLACE,q,N*fpm(23),MPI_DOUBLE_PRECISION,MPI_SUM,NEW_COMM_WORLD,code)
     end if
#endif
     !------------------------------------------------------------------------

     if ((fpm(29)==1).and.(fpm(16)==2)) then !! Zolotarev
        call dset_feast_zolotarev(fpm(2),0,zxe,zwe)
        !   r=(Emax-Emin)/2.0d0
        if ((loop==0).and.(fpm(5)==0)) then 
        !   !q(1:N,fpm(23))=-dble(zwe)*r*work(1:N,1:fpm(23)) 
            q(1:N,1:fpm(23))=DZERO
        else
            if (rank==0) then
             q(1:N,1:fpm(23))=-dble(zwe)*q(1:N,1:fpm(23))
            else 
             q(1:N,1:fpm(23))=DZERO
            endif
        end if
     else
        q(1:N,1:fpm(23))=DZERO
     end if

     fpm(20)=1
     fpm(21)=2
     ijob=-2 ! just initialization 
  end IF


!!!!!!!!!!!!
  IF (fpm(21)==2) then !! we start or pursue the contour integration

     IF (info==0) then !! will end up checking info errors returned by FEAST drivers
        do e=fpm(20)+rank,fpm(22),nb_procs !!!! loop over the contour 

           if (ijob==-2) then !!Factorize the linear system (complex) (zS-A)
              Ze=Zne(e)
              fpm(20)=e-rank
              ijob=10 ! for fact
              if ((loop==0).or.(fpm(22)>nb_procs)) return ! no need to factorize again if one linear system per processor
           endif

           if (ijob==10) then !!Solve the linear system (complex) (zS-A)q=v 
              call ZLACP2( 'F', N, fpm(23),work , N, workc, N )
              ijob=11 ! for solve
              return
           endif

           if (ijob==11) then 
              !! summation              
              aux=2.0d0*Wne(e)            
              !! aux has been multiplied by 2 to account for the 2 half-contours
              q(1:N,1:fpm(23))=q(1:N,1:fpm(23))-dble(aux*workc(1:N,1:fpm(23))) 

              ijob=-2 ! just for identification
           end if
        end do
     end IF  !! info=0

     !------------------------------------------------
#ifdef MPI
     call MPI_BCAST(info,1,MPI_INTEGER,0,NEW_COMM_WORLD,code)
#endif
     !-----------------------------------------------  
     if (info/=0) fpm(21)=100 ! the end

     if (info==0) then
        fpm(21)=4 
        !------------------------------------------------
#ifdef MPI
        call MPI_ALLREDUCE(MPI_IN_PLACE,q(1,1),N*fpm(23),MPI_DOUBLE_PRECISION,MPI_SUM,NEW_COMM_WORLD,code)
#endif
        !-----------------------------------------------  
     end if
  end IF     ! fpm(21)==1
!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if ((fpm(21)==4).and.(fpm(14)==1)) then !! only q vectors has been computed and is returned
     info=4
     if (info/=0) fpm(21)=100 ! The End
  end if


!!!!!!!!!!!!!!!!!!!!!!!!!!! Stochastic estimates
  if ((fpm(21)==4).and.(fpm(14)==2)) then !! only q is returned with stochastic estimation in res
     iseed=(/56,890,3456,2333/)
     call DLARNV(2,iseed,N*fpm(23),work(1,1)) 
     call DGEMM('T','N',fpm(23),fpm(23),N,-DONE,work,N,q,N,DZERO,Aq,M0) ! projection
     call DGEMM('T','N',fpm(23),fpm(23),N,DONE,work,N,work,N,DZERO,Sq,M0) ! normalization
     theta=DZERO
     do i=1,fpm(23)
        theta=theta+Aq(i,i)/Sq(i,i)
        res(i)=abs(theta*(N/(DONE*i)))
     enddo
     mode=int(real(res(fpm(23))))+1
     info=5
     if (info/=0) fpm(21)=100 ! The End
  end if

!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Form the reduced eigenvalue problem
!!!!!!! Aq xq=eq Bq xq     with Aq=Q^TAQ Bq=Q^TBQ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!form  Bq=> Bq=Q^T B Q
  if (fpm(21)==4) then 
#ifdef MPI 
     work(1:N,1:fpm(23))=DZERO
#endif
     fpm(21)=5 ! preparing reenty
     ijob=40 
     return! mat-vec B*q => work
  end if


  if (fpm(21)==5) then
     !----------------------------------------------
#ifdef MPI
     Sq(1:M0,1:fpm(23))=DZERO
#endif
     !-------------------------------------------------

     call DGEMM('T','N',fpm(23),fpm(25),N,DONE,q(1,1),N,work(1,fpm(24)),N,DZERO,Sq(1,fpm(24)),M0) 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! Spurious test 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     mode=0
     if (loop>0) then
        do i=fpm(24),fpm(24)+fpm(25)-1
           if (res(i)==DONE) then ! indicator for being inside the search interval
              if (abs(sqrt(abs(Sq(i,i)))-lambda(i))/sqrt(abs(Sq(i,i)))<0.1d0) mode=mode+1
           endif
        enddo
     endif
!!!!!!!!!!!!!!!!!!
#ifdef MPI
     call MPI_ALLREDUCE(MPI_IN_PLACE,mode,1,MPI_INTEGER,MPI_SUM,NEW_COMM_WORLD,code)
#endif
!!!!!!!!!!!!!!!!!!
     fpm(21)=6
  end if


!!!!!!!!!form Aq=> Aq=Q^T A Q 
  if (fpm(21)==6) then
     !------------------------------------------------
#ifdef MPI 
     work(1:N,1:fpm(23))=DZERO
#endif
     fpm(21)=7 ! preparing reentry
     ijob=30 
     return  ! mat-vec A*q => work
  endif


  if (fpm(21)==7) then 
     !------------------------------------------------
#ifdef MPI 
     Aq(1:M0,1:fpm(23))=DZERO
#endif
     !-------------------------------------------------
     call DGEMM('T','N',fpm(23),fpm(25),N,DONE,q(1,1),N,work(1,fpm(24)),N,DZERO,Aq(1,fpm(24)),M0) !q^tAq
     fpm(21)=8
  endif


  if (fpm(21)==8) then

     !---------------------------------------- !(Aq,Sq known to all processors) 
#ifdef MPI
     if (fpm(23)/nb_procs>=1) then
        call MPI_ALLREDUCE(MPI_IN_PLACE,Aq(1,1),M0*fpm(23),MPI_DOUBLE_PRECISION,MPI_SUM,NEW_COMM_WORLD,code)
        call MPI_ALLREDUCE(MPI_IN_PLACE,Sq(1,1),M0*fpm(23),MPI_DOUBLE_PRECISION,MPI_SUM,NEW_COMM_WORLD,code)
     end if
#endif
     !---------------------------------------

     !---------------------------------------
     if (fpm(12)==1) then ! customize eigenvalue solver
        fpm(21)=9 ! preparing reentry - could return new value of M0 in fpm(23) if reduced subspace is needed
        ijob=50
        return
     endif
  endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Solve the reduced eigenvalue problem
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (fpm(21)==8) then
     ! if using FEAST-MPI ==> solve on a single proc
     if (rank==0) then
        JOBZ='V'
        UPLO='L'
        info_lap=1 ! initialization
        i=1
        LWORK_LOC=3*fpm(23)-1 !! for lapack eig reduced system
        call wallocate_1d(WORK_LOC,LWORK_LOC,infoloc)
        if (infoloc/=0) info=-1


        do while ((info_lap/=0).and.(info==0))
           i=i+1
           if (i==10) info=-3 ! arbitrary maximum
           call wallocate_2d(Sqo,fpm(23),fpm(23),infoloc)
           if (infoloc/=0) info=-1
           call DLACPY( 'F', fpm(23), fpm(23),Sq, M0, Sqo, fpm(23) )
           call DSYGV(1, JOBZ, UPLO, fpm(23),Aq,M0,Sqo,fpm(23),lambda,work_loc,Lwork_loc,INFO_lap)
            ritzvals(feastit,1:m0)=lambda(1:m0)
           if ((info_lap<=fpm(23)).and.(info_lap/=0)) info=-3
           if (info_lap>fpm(23)) then !! Sqo is not spd (a posteriori resize subspace)
              fpm(23)=info_lap-fpm(23)-1 
              if (fpm(1)==1) then
                 call wwrite_s(fout, 'Resize subspace')  
                 call wwrite_t(fout) 
                 call wwrite_i(fout,fpm(23))
                 call wwrite_n(fout)
              end if
           end if
           call wdeallocate_2d(Sqo)
        end do
        call wdeallocate_1d(work_loc) 
     end if !(rank 0)
     !-------------------------------- !(info common to all processors)
#ifdef MPI
     call MPI_BCAST(info,1,MPI_INTEGER,0,NEW_COMM_WORLD,code)
#endif
     !--------------------------------
     if (info/=0) fpm(21)=100 ! the end

     if (info==0) then
        fpm(25)=fpm(23)!! current M0 size (by default) -- global
        fpm(24)=1  !! origin first column number of vector q for parallel mat-vec (default)
        !----------------------------------------!(Aq==> vectors, lambda and fpm(23), known to all processors) 
#ifdef MPI 
        call MPI_BCAST(fpm(23),1,MPI_INTEGER,0,NEW_COMM_WORLD,code)
        call MPI_BCAST(Aq,fpm(23)*M0,MPI_DOUBLE_PRECISION,0,NEW_COMM_WORLD,code)
        call MPI_BCAST(lambda,fpm(23),MPI_DOUBLE_PRECISION,0,NEW_COMM_WORLD,code)
        if (fpm(23)/nb_procs>=1) then ! criteria for parallelism of mat-vec
           fpm(25)=fpm(23)/nb_procs ! local size of current M0
           if (rank==nb_procs-1) fpm(25)=fpm(25)+fpm(23)-(fpm(23)/nb_procs)*nb_procs 
           fpm(24)=1+rank*(fpm(23)/nb_procs) ! local origin of first column for vector q for parallel mat-vec 
        end if
#endif
        !-----------------------------
        fpm(21)=9
     end if
  end if !! fpm(21)=8



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Compute Ritz vectors X=Qxq  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (fpm(21)==9) then
     call DLACPY( 'F', N, fpm(23),q , N, work, N ) 
     !! option - non shifted
#ifdef MPI 
     q(1:N,1:fpm(23))=DZERO
#endif
     call DGEMM('N','N',N,fpm(25),fpm(23),DONE,work(1,1),N,Aq(1,fpm(24)),M0,DZERO,q(1,fpm(24)),N)
  endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Compute Residual ||AX-lambda*BX||
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (fpm(21)==9) then
     fpm(21)=10 ! preparing reentry
     ijob=30
     return  ! mat-vec A*q => work
  endif

  if (fpm(21)==10) then
     call ZLACP2( 'F', N, fpm(25),work(1,fpm(24)), N, workc(1,fpm(24)), N )
     fpm(21)=11 ! preparing reentry
     ijob=40
     !----------------------------------------
#ifdef MPI
     work(1:N,1:fpm(23))=DZERO  !! work is also needed for outer-loop if any 
#endif
     !----------------------------------------
     return  ! mat-vec S*q => work 
  endif

  if (fpm(21)==11) then
     !----------------------------------------
#ifdef MPI
     res(1:fpm(23))=DZERO
#endif
     !-------- Absolute residual
     do i=fpm(24),fpm(24)+fpm(25)-1
           res(i)=sum(abs(dble(workc(1:N,i))-lambda(i)*work(1:N,i)))/sum(abs(max(abs(Emin),abs(Emax))*(work(1:N,i)))) 
     end do
     eigresall(feastit,fpm(24):fpm(24)+fpm(25)-1)=res(fpm(24):fpm(24)+fpm(25)-1)
     !----------------------------------------
#ifdef MPI 
     if (fpm(23)/nb_procs>=1) call MPI_ALLREDUCE(MPI_IN_PLACE,res,fpm(23),MPI_DOUBLE_PRECISION,MPI_SUM,NEW_COMM_WORLD,code)
#endif
     !-----------------------------------------
  end if



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Compute Trace / count eigenvalues (remove spurious if any)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  if (fpm(21)==11) then
     k=0
     trace=DZERO
     theta=DZERO
!!! count how many eigenvalues have converged + trace + max residual
     do i=1,fpm(23)
        if ((lambda(i)>Emin).and.(lambda(i)<Emax)) then ! inside the search interval
           k=k+1 ! number of eigenvalues (could include spurious)
           trace=trace+lambda(i)
           if (res(i)>theta) theta=res(i) ! max residual 
        endif
     enddo

!!!!!!!! remove spurious if any       
     if ((mode==0).and.(k>0)) mode=k ! wait before looking into mode 
     ! Rq: if all eigenvalues k are spurious...FEAST will not converge
     if (mode<k) then
        do j=1,k-mode
           theta=DZERO
           e=1
           do i=1,fpm(23)
              if ((lambda(i)>Emin).and.(lambda(i)<Emax)) then ! inside the search interval
                 if (res(i)>theta) then
                    e=i
                    theta=res(i) !max
                 endif
              end if
           enddo
           trace=trace-lambda(e) ! update trace
           res(e)=-DONE !spurious
        end do
!!! max residual
        theta=DZERO
        do i=1,fpm(23)
           if ((lambda(i)>Emin).and.(lambda(i)<Emax)) then ! inside the search interval
              if (res(i)>theta) theta=res(i)
           end if
        enddo
     end if

     !if (mode==0) info=1  ! no eigenvalue detected in the interval
     if (loop>1) then ! wait second iteration (spurious related)
        if ((mode==M0).and.(mode/=N)) info=3 ! size subspace too small
     endif
     if (info/=0) fpm(21)=100 ! The End
  end if



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Check FEAST Convergence
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       

  if (fpm(21)==11) then
     testconv=.false. ! initialization

     !! trace convergence
     if (loop==0) epsout=DONE
     if (loop>0) then
        epsout=(abs(trace-epsout))/max(abs(Emin),abs(Emax))           
        if ((fpm(6)==0).and.(log10(epsout)<(-fpm(3)))) testconv=.true.
     end if

     !! residual convergence
     if ((fpm(6)/=0).and.(log10(theta)<(-fpm(3)))) testconv=.true.

     !! Two loops minimum if spurious are found
     if ((loop<=1).and.(k>mode)) testconv=.false.

     if (fpm(1)==1) then
        call wwrite_i(fout,loop)
        call wwrite_t(fout) 
        call wwrite_i(fout,mode)
        call wwrite_t(fout)
        call wwrite_d(fout,trace)
        call wwrite_t(fout) 
        call wwrite_d(fout,epsout)
        call wwrite_t(fout) 
        call wwrite_d(fout,theta)
        call wwrite_n(fout) 
     end if

     if (.not.testconv) then
        epsout=trace
        if (loop==fpm(4)) then
           info=2 ! FEAST did not converge (#loop reaches maximum)
           testconv=.true. ! return final eigenvector anyway
        endif
     endif

    if (mode==0) testconv=.false. !assume we'll find something, don't stop until we do

    !make epsout be the eigenvector residual rather than the trace residual`    
     if (fpm(6)/=0) then 
         epsout=theta
     end if
 
  end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! FEAST exit IF Convergence - FEAST iteration IF NOT 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!            

  if (fpm(21)==11) then

     if (.not.testconv) then !!! need FEAST refinement loop


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! rational function for the inner eigenvalues
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        res(1:fpm(23))=DZERO ! temp indicator for being inside the contour
        call wallocate_1d(work_loc,fpm(23),infoloc)
        call dfeast_rationalx(Zne,Wne,fpm(2),lambda,fpm(23),work_loc)
        if ((fpm(29)==1).and.(fpm(16)==2)) then ! zolotarev case
        call dset_feast_zolotarev(fpm(2),0,zxe,zwe)
        work_loc(1:fpm(23))=work_loc(1:fpm(23))+zwe
        endif
        do i=1,fpm(23)
           if ((lambda(i)>Emin).and.(lambda(i)<Emax))  res(i)=1.0d0
           lambda(i)=work_loc(i)
        enddo
        call wdeallocate_1d(work_loc)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!! Remark: q and work=S*q is known already, using FEAST-MPI work and q are already distributed


        fpm(21)=1   ! prepare reentry- reloop (with contour integration)
        !fpm(21)=4 ! reloop (without contour integration) -in this case work=q (actually does not need "work")
        loop=loop+1
        ijob=-2 ! do nothing
        return  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     else    !!!!!!! final eigenvectors/eigenvalues
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef MPI       
        if (fpm(23)/nb_procs>=1) call MPI_ALLREDUCE(MPI_IN_PLACE,q,N*fpm(23),MPI_DOUBLE_PRECISION,MPI_SUM,NEW_COMM_WORLD,code)
#endif

!!! reorder (shift) lambda, eigenvector and residual
        call DLACPY( 'F', N, fpm(23),q , N, work, N ) 
        call wallocate_1d(work_loc,fpm(23),infoloc)
        call wallocate_1d(work_loc2,fpm(23),infoloc)
        call DCOPY(fpm(23),lambda , 1, work_loc, 1 )
        call DCOPY(fpm(23),res , 1, work_loc2, 1 )
        q(1:N,1:fpm(23))=DZERO
        e=0
        j=0
        k=0
        do i=1,fpm(23)
           if ((work_loc(i)>Emin).and.(work_loc(i)<Emax)) then ! inside the search interval
              if (work_loc2(i)/=-DONE) then ! spurious 
                 e=e+1
                 q(1:N,e)=work(1:N,i)
                 lambda(e)=work_loc(i)
                 res(e)=work_loc2(i)
              endif
           else
              j=j+1
              q(1:N,mode+j)=work(1:N,i)
              lambda(mode+j)=work_loc(i)
              res(mode+j)=work_loc2(i)
           end if
           if (work_loc2(i)==-DONE) then ! spurious at the end
              k=k+1
              q(1:N,fpm(23)-k+1)=work(1:N,i)
              lambda(fpm(23)-k+1)=work_loc(i)
              res(fpm(23)-k+1)=work_loc2(i)
           endif
        enddo
        call wdeallocate_1d(work_loc)
        call wdeallocate_1d(work_loc2)

!!!!!!!!!!!!!!!!!!!!!!!!!!
        M0=fpm(23)  ! update value of M0 (new subspace)
        fpm(21)=100 ! The End

     end if ! test convergence
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (fpm(21)==100) then !! THE END (ijob=0) 
     ijob=0 !! exit FEAST

     if (fpm(1)==1) then !! Print  Information

        if (info>=200) then
           call wwrite_s(fout, 'PROBLEM with input parameters')
           call wwrite_n(fout) 
        end if

        if ((info>100).and.(info<200)) then
           call wwrite_s(fout, 'PROBLEM with FEAST array parameters')
           call wwrite_n(fout) 
        end if

        if (info==-3) then
           call wwrite_s(fout, 'ERROR with reduced system')  
           call wwrite_n(fout) 
        end if

        if (info==-2) then
           call wwrite_s(fout, 'ERROR from Inner Linear System Solver in FEAST driver')  
           call wwrite_n(fout) 
        end if

        if (info==-1) then
           call wwrite_s(fout, 'ERROR with Internal memory allocation')  
           call wwrite_n(fout) 
        end if

        if (info==1) then
           call wwrite_s(fout, '==>WARNING: No eigenvalue has been found in the proposed search interval')
           call wwrite_n(fout)
        endif

        if (info==3) then
           call wwrite_s(fout, '==>WARNING: Size subspace M0 too small')  
           call wwrite_n(fout)
        end if

        if (info==4) then
           call wwrite_s(fout, '==>WARNING: Only the subspace has been returned')  
           call wwrite_n(fout)
        end if

        if (info==5) then
           call wwrite_s(fout, '==>WARNING: Only stochastic estimation of #eigenvalues returned')  
           call wwrite_n(fout)
        end if

        if (info==2) then
           call wwrite_s(fout, '==>WARNING: FEAST did not converge "yet" (#loop reaches maximum allowed)')  
           call wwrite_n(fout)
        end if

        if (info==0) then
           call wwrite_s(fout, '==>FEAST has successfully converged (to desired tolerance)')  
           call wwrite_n(fout) 
        else
           call wwrite_s(fout, '==>INFO code = ') 
           call wwrite_i(fout,info)
           call wwrite_n(fout)
        end if


        call wwrite_s(fout, '***********************************************')  
        call wwrite_n(fout) 
        call wwrite_s(fout, '*********** FEAST- END*************************')
        call wwrite_n(fout) 
        call wwrite_s(fout, '***********************************************')  
        call wwrite_n(fout) 
        call wwrite_n(fout)
     endif
  end if

end subroutine dfeast_myrcix





