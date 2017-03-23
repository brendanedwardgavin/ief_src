  subroutine dfeast_scsrevit(UPLO,nrow,ncol,sa,isa,jsa,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
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
    integer :: Nrow,Ncol
    double precision,dimension(*) :: sa
    integer,dimension(*) :: isa,jsa
    integer,dimension(*) :: fpm
    double precision :: epsout 
    integer :: loop
    double precision :: Emin,Emax
    integer :: M0
    double precision,dimension(*)  :: E
    double precision,dimension(min(nrow,ncol),*):: X
    integer :: mode
    double precision,dimension(*)    :: res
    integer :: info

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!    Wrapper Routine to expert routine: dfeast_scsrgvx
!!!    Create Zne, Wne arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    complex(kind=(kind(1.0d0))),dimension(fpm(2)) :: Zne,Wne 
    double precision,dimension(nrow) :: sb ! identity !n 
    integer,dimension(nrow+1) :: isb !n+1
    integer,dimension(nrow) :: jsb !n
    integer :: i
    double precision :: delta


    
    if(fpm(55)<=0) then ! include the case -1 where the nodes are forced to real axis
        call zfeast_contour(Emin,Emax,fpm(2),fpm(16),fpm(18),Zne,Wne)
    elseif(fpm(55)==1) then
        delta=(Emax-Emin)/fpm(2)
        do i=1,fpm(2)
            call zfeast_contour(Emin+(i-1)*delta,Emin+i*delta,1,fpm(16),fpm(2)*fpm(18),Zne(i),Wne(i))
        end do
    elseif(fpm(55)==2) then
         !set contour points, then flatten them on to real axis
         call zfeast_contour(Emin,Emax,fpm(2),fpm(16),fpm(18),Zne,Wne)
         do i=1,fpm(2)
            Zne(i)=Zne(i)-(0.0d0,1.0d0)*aimag(Zne(i))
         end do
!elseif(fpm(55)==3) then ! real axis - one point
!        call zfeast_contour(Emin,Emax,fpm(2),fpm(16),fpm(18),Zne,Wne)
!!Wne=aimag(Wne)*(1.0d0,0.0d0)
!Zne=dble(Zne)*(1.0d0,0.0d0)


    end if

    fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
    ! identity B matrix- option for standard eigenvalue problem
    do i=1,nrow
       sb(i)=1.0d0
       jsb(i)=i
       isb(i)=i
    enddo
    isb(nrow+1)=nrow+1

    call dfeast_scsrgvxit(UPLO,Nrow,Ncol,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)

  end subroutine dfeast_scsrevit



  subroutine dfeast_scsrgvit(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
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
    !=====================================================================
    ! Eric Polizzi 2009-2015
    !===================================================================== 
    implicit none
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!    Wrapper Routine to expert routine: dfeast_scsrgvx
!!!    Create Zne, Wne arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    complex(kind=(kind(1.0d0))),dimension(fpm(2)) :: Zne,Wne 

    call zfeast_contour(Emin,Emax,fpm(2),fpm(16),fpm(18),Zne,Wne)

    fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
     call dfeast_scsrgvxit(UPLO,N,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)
   
  end subroutine dfeast_scsrgvit


  
  
  
  
  
  
  
  
  subroutine dfeast_scsrgvxit(UPLO,Nrow,Ncol,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)
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
    integer :: Nrow,Ncol
    double precision,dimension(*) :: sa,sb
    integer,dimension(*) :: isa,jsa,isb,jsb
    integer,dimension(*) :: fpm
    double precision :: epsout 
    integer :: loop
    double precision :: Emin,Emax
    integer :: M0
    double precision,dimension(*)  :: E
    double precision,dimension(min(nrow,ncol),*):: X
    integer :: mode
    double precision,dimension(*)    :: res
    integer :: info
    complex(kind=(kind(1.0d0))),dimension(*) :: Zne,Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer :: ijob,infoloc,i,s,nnz
    complex(kind=(kind(1.0d0))) :: Ze
    complex(kind=(kind(1.0d0))), dimension(:,:),pointer ::workc,caux
    double precision, dimension(:,:),pointer ::work,Aq,Sq,aux
    complex(kind=(kind(1.0d0))),dimension(:,:),pointer :: saz  
    double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
    complex(kind=(kind(1.0d0))),parameter :: ONEC=(DONE,DZERO)
  

    logical :: fact
    integer :: nfact,id
    integer :: rank,code,nb_procs,NEW_COMM_WORLD
!!!!! csr-upper format
    double precision,dimension(:),pointer :: nres
    double precision :: ares
    logical :: comb
    integer :: infob,itmax

!!!!!!!!!! iterative solver
    complex(kind=(kind(1.0d0))) :: ze2 
    complex(kind=(kind(1.0d0))), dimension(:,:), pointer :: ztempmat,zX,Rvec
    integer :: linresindex !number of RHS to use in measuring linear system error
    integer :: blockits,blockremainder
    double precision :: lintargeterror,linsyserror !goal error for linear system, return linear system error
    integer :: linloops
    double precision :: linepsout
    integer :: k,ii,j

    integer :: nmin,nmax

!!!MKL mat-vec:
    character, dimension(6) :: matdescra


    if(UPLO=='F') then
       matdescra(1)='G'
    else
       matdescra(1)='S'
    end if
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
    ELSE IF ( Nrow<=0 ) THEN
       INFO = -102
    END IF
    IF( INFO.NE.0 ) THEN
       CALL XERBLA( 'DFEAST_SCSRIT', -INFO+100 )
       RETURN
    END IF


    infoloc=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    nmin=min(nrow,ncol)
    nmax=max(nrow,ncol)


    call wallocate_2d(Aq,M0,M0,infoloc)
    call wallocate_2d(Sq,M0,M0,infoloc)
    call wallocate_2d(work,nmin,M0,infoloc)
    call wallocate_2z(workc,nmin,M0,infoloc)    
    call wallocate_2z(ztempmat,nmin,m0,infoloc)



    call wallocate_2d(aux,nmax,M0,infoloc)


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
    nnz=isa(nrow+1)-1
    call wallocate_2z(saz,nnz,nfact,infoloc)
    do i=1,nfact
       saz(1:nnz,i)=sa(1:nnz)*(1.0d0,0.0d0)   
    enddo


    if (infoloc/=0) then
       info=-1
       return
    end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!! solution saved at ijob=-2
    call wallocate_2z(zX,nmin,m0,infoloc)
    call wallocate_2z(Rvec,nmin,M0,infoloc)   
    zX(1:nrow,1:M0)=(0.0d0,0.0d0)



    ijob=-1 ! initialization


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do while (ijob/=0)
       call dfeast_myrcix(ijob,nmin,Ze,work,workc,Aq,Sq,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)    

       !print *,ijob

       select case(ijob)

       case(-2)  !! retrieve current solution, residual after first loop (eigenvalues are lambda)
          zX(1:nrow,1:M0)=X(1:nrow,1:M0)*(1.0d0,0.0d0)

          !! get the  residual
          do i=1,M0
             call DSCAL(nmin,E(i),work(1,i),1)
          enddo
          call ZLACP2( 'F', nmin, M0,work , nmin, Rvec, nmin )
          call ZAXPY(nmin*M0,(-1.0d0,0.0d0),workc,1,Rvec,1)


       case(10) !! factorize (zeB-A)
          !
          !  print *,Ze     
          !user iterative solver, don't need factorization
          !            ze2=dble(Ze)*(1.0d0,0.0d0)
          if (fpm(55)==-1) then !!! real axis
             ze2=dble(Ze)*(1.0d0,0.0d0)
          else
             ze2=Ze
          endif
          !print *,ze2!,Zne(1),Wne(1)

       case(11) !!solve the linear system (ZeB-A)x=workc(1:N,1:fpm(23)) result in to workc

          if (loop>0) then   
             if (fpm(11)<0) workc(1:nrow,1:M0)=-Rvec(1:nrow,1:M0) 
          end if

!!!!!! Set convergence criteria
             if (fpm(54)==0) then
                lintargeterror=1d-2
             else
                lintargeterror=1.0d0/fpm(54)
             end if

             if ((loop>0).and.(fpm(11)>0)) lintargeterror=lintargeterror*epsout !depends on feast convergence
             
!!$          
!!$          if(loop>0 .and. (epsout .ne. 0.0d0)) then
!!$             if(fpm(54)==0) then
!!$                lintargeterror=1.0d-2*epsout
!!$                if (fpm(11)<0)  lintargeterror=1.0d-2 !! error indep of FEAST since act on the residual
!!$             else
!!$                lintargeterror=epsout/(1.0d0*fpm(54))
!!$                if (fpm(11)<0) lintargeterror=1.0d0/(1.0d0*fpm(54))!! error indep of FEAST since act on the residual
!!$
!!$             end if
!!$          else
!!$             if (fpm(54)==0) then
!!$                lintargeterror=1d-2
!!$             else
!!$                lintargeterror=1.0d0/fpm(54)
!!$             end if
!!$          endif


          
          !!!!!!!use fpm(53) as block size
          !!!!!! Attention: M0 should a multiple of fpm(53)
          blockits=m0/fpm(53)
          blockremainder=m0-blockits*fpm(53)

          if (blockremainder/=0) then
             print *,'Input error: M0 should be a multiple of fpm(53)'
             stop
          end if


!!!!!!!!!!! main loop over the block rhs
          do i=1,blockits

!!$if (fpm(62)>0) then
!!$if (loop>0) then
!!$if (E(i)<fpm(62)*0.1d0) then ! outside interval E--> res-200
!!$fpm(51)=2 
!!$else 
!!$fpm(51)=fpm51
!!$endif
!!$endif
!!$end if

             if ((fpm(11)==1).or.(fpm(11)==-1))   call blockBGMRESarnoldi(UPLO,'N',nrow,ncol,fpm(53),saz(1,id),isa,jsa,sb,isb,jsb,ze2,fpm(51),fpm(50),workc(1,1+(i-1)*fpm(53)),ztempmat(1,1+(i-1)*fpm(53)),lintargeterror,linloops,1+(i-1)*fpm(53)) 



             if(fpm(11)==2) then    
                linresindex=fpm(53)
                print *,'to check'!!call zfeast_cglsRes(UPLO,n,fpm(53),zsa,isa,jsa,ze2,nnza,workc(1,1+(i-1)*fpm(53)),ztempmat(1,1+(i-1)*fpm(53)),fpm(51)*fpm(50),lintargeterror,linresindex,linsyserror,linloops)
                !workc(1+(i-1)*fpm(53):i*fpm(53),1:n)=ztempmat(1:fpm(53),1:n)
                !!linit(feastit,cpnum,1+(i-1)*fpm(53):i*fpm(53))=linloops
                !print *,'block ',i,linloops
             end if

             
             if(fpm(11)==3) call zminres(UPLO,nrow,ncol,saz(1,id),isa,jsa,ze2,workc(1,1+(i-1)*fpm(53)),ztempmat(1,1+(i-1)*fpm(53)),lintargeterror,fpm(51),linloops)

             
             if ((fpm(11)==4).or.(fpm(11)==-4)) call zminresBlock(UPLO,'N',nrow,ncol,fpm(53),saz(1,id),isa,jsa,sb,isb,jsb,ze2,workc(1,1+(i-1)*fpm(53)),ztempmat(1,1+(i-1)*fpm(53)),lintargeterror,fpm(51),linloops,1+(i-1)*fpm(53))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


             print *,i,linloops

          end do


          call ZCOPY( nmin*M0, ztempmat, 1, workc, 1 ) 
          !workc=ztempmat
          if (loop>0) then
             if (fpm(11)<0) then
                ! workc<=(ztempmat+initial_guess)*(1/(ze-Ei))
                call ZAXPY(nmin*M0,(1.0d0,0.0d0),zX,1,workc,1)
                do i=1,M0
                   call ZSCAL(nmin,(1.0d0,0.0d0)/(ze2-E(i)*(1.0d0,0.0d0)), workc(1,i), 1)
                enddo
             endif
          end if



       case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

          if (nrow==ncol) then
             call mkl_dcsrmm('N', nrow, m0, ncol, 1.0d0, matdescra, sa, jsa, isa, isa(2), X, ncol, 0.0d0, work, ncol)
          elseif (nrow>ncol) then
             call mkl_dcsrmm('N', nrow, m0, ncol, 1.0d0, matdescra, sa, jsa, isa, isa(2), X, ncol, 0.0d0, aux, nrow)
             call mkl_dcsrmm('T', nrow, m0, ncol, 1.0d0, matdescra, sa, jsa, isa, isa(2), aux, nrow, 0.0d0, work, ncol)
          elseif (nrow<ncol) then
             call mkl_dcsrmm('T', nrow, m0, ncol, 1.0d0, matdescra, sa, jsa, isa, isa(2), X, nrow, 0.0d0, aux, ncol)
             call mkl_dcsrmm('N', nrow, m0, ncol, 1.0d0, matdescra, sa, jsa, isa, isa(2), aux, ncol, 0.0d0, work, nrow)
          endif



       case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)


          call mkl_dcsrmm('N', nrow, m0, ncol, 1.0d0, matdescra, sb, jsb, isb, isb(2), X, ncol, 0.0d0, work, ncol)
          ! work(1:nmin,1:m0)=X(1:nmin,1:m0)
          !            print *,'why??? X',X(1:nrow,1:m0)




       end select
    end do



    call wdeallocate_2d(Aq)
    call wdeallocate_2d(Sq)
    call wdeallocate_2d(work)
    call wdeallocate_2z(workc)
    call wdeallocate_2d(aux)

    call wdeallocate_2z(saz)


  end subroutine dfeast_scsrgvxit





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine dfeast_myrcix(ijob,N,Ze,work,workc,Aq,Sq,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info,Zne,Wne)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
    use rundata

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
           if (res(i)>100d0) then !(res(i)==DONE) then ! indicator for being inside the search interval
              !if (abs(sqrt(abs(Sq(i,i)))-lambda(i))/sqrt(abs(Sq(i,i)))<0.1d0) mode=mode+1
              if (abs(sqrt(abs(Sq(i,i)))-(res(i)-200d0))/sqrt(abs(Sq(i,i)))<0.1d0) mode=mode+1
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
            !ritzvals(feastit,1:m0)=lambda(1:m0)
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
     !eigresall(feastit,fpm(24):fpm(24)+fpm(25)-1)=res(fpm(24):fpm(24)+fpm(25)-1)
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

    if (mode==0) testconv=.false. !assume we'll find something, don't stop until we do !<<<<<<<<<<<

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
           if ((lambda(i)>Emin).and.(lambda(i)<Emax))  res(i)=200d0+work_loc(i)!res(i)=1.0d0 !!!!<<<< trick not to modify lambda
           !lambda(i)=work_loc(i)
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





subroutine blockGMRESarnoldi(UPLO,TRANS,nrow,ncol,m,dsa,isa,jsa,ze,kmax,restarts,Brhs,Xlhs,eps,loops,blockstart)
!n is ncol
use rundata

implicit none

character :: UPLO,TRANS
integer :: nrow,ncol,m,kmax,restarts,loops,blockstart
integer, dimension(*) :: isa,jsa
complex (kind=kind(0.0d0)) :: ze
complex (kind=kind(0.0d0)), dimension(*) :: dsa
complex (kind=kind(0.0d0)), dimension(min(nrow,ncol),*) :: Brhs,Xlhs
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

!double precision, external::elapsed_time
integer :: nmin,nmax

nmin=min(nrow,ncol)
nmax=max(nrow,ncol)


!print *,kmax,restarts

allocate(V(nmin,(kmax+1)*m),H(1:(kmax+1)*m,1:kmax*m), Htmp((kmax+1)*m,kmax*m)  ,Bsm((kmax+1)*m,m), ym((kmax+1)*m,m), R(1:nmin,1:m), R2(1:nmin,1:m), Xtmp(1:nmin,1:m) )



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

R=Brhs(1:nmin,1:m)
Xlhs(1:nmin,1:m)=(0.0d0,0.0d0)

loops=0
do i=1,restarts

    if(i>1) then
        !R=Brhs-A*Xlhs
        !R=Brhs(1:n,1:m)
        R=R2
        !call mkl_zcsrmm('N', n, m, n, (-1.0d0,0.0d0), matdescra, dsa, jsa, isa, isa(2), Xlhs(1,1), n, (1.0d0,0.0d0), R(1,1), n)
    end if

    V(1:nmin,1:m)=R(1:nmin,1:m)

    do j=1,kmax
        loops=loops+1
        !next arnoldi step
        call system_clock(count=tc1)
        call blockArnoldiIt(UPLO,TRANS,nrow,ncol,m,dsa,isa,jsa,ze,kmax,j,V,H,Bsm)
        !call blockArnoldiIt(UPLO,'N',nrow,ncol,m,dsa,isa,jsa,ze,kmax,j,V,H,Bsm)
        
        call system_clock(count=tc2)
        !arnolditime=arnolditime+elapsed_time(tc1,tc2)

        !solve system H(1:(j+1)*m,1:j*m)ym=Bsm(1:(j+1)*m,1:m)
        ym(1:(j+1)*m,1:m)=Bsm(1:(j+1)*m,1:m)
        Htmp(1:(j+1)*m,1:j*m)=H(1:(j+1)*m,1:j*m)

        !print *,'here'  
        !print *,'H=',Htmp(1,2)
        !print *,'Bsm=',sum(Bsm)
        !print *,'V=',sum(V)

        call system_clock(count=tc1)      
        call zgels('N',(j+1)*m,j*m,m,Htmp(1,1),(kmax+1)*m,ym(1,1),(kmax+1)*m,work,lwork,info)
        !call zgels(TRANS,(j+1)*m,j*m,m,Htmp(1,1),(kmax+1)*m,ym(1,1),(kmax+1)*m,work,lwork,info)
        
        if(info .ne. 0) then
            print *,'Error in blockGMRESarnoldi'
            print *,'ZGELS error ',info
            stop
        end if
        call system_clock(count=tc2)
        !lstime=lstime+elapsed_time(tc1,tc2)

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
        call zgemm('N','N',nmin,m,(j+1)*m,(-1.0d0,0.0d0),V(1,1),nmin,Htmp(1,1),(kmax+1)*m,(1.0d0,0.0d0),R2(1,1),nmin)

        !Xtmp=Xlhs(1:n,1:m)
        !call zgemm('N','N',n,m,j*m,(1.0d0,0.0d0),V(1:n,1:j*m),n,ym(1:j*m,1:m),j*m,(1.0d0,0.0d0),Xtmp(1:n,1:m),n)
        !R2=Brhs(1:n,1:m)
        !call mkl_zcsrmm('N', n, m, n, (-1.0d0,0.0d0), matdescra, dsa, jsa, isa, isa(2), Xtmp(1,1), n, (1.0d0,0.0d0), R2(1,1), n)

        !print *,i,j
        error=0.0d0
        do l=1,m
            error2=dznrm2(nmin,R2(:,l),1)/dznrm2(nmin,Brhs(:,l),1)
            if (error2>error) error=error2
            !!linres(feastit,cpnum,loops,blockstart+l-1)=error2
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
    call zgemm('N','N',nmin,m,j*m,(1.0d0,0.0d0),V(1,1),nmin,ym(1,1),(kmax+1)*m,(1.0d0,0.0d0),Xlhs(1,1),nmin)

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



subroutine blockArnoldiIt(UPLO,TRANS,nrow,ncol,m,dsa,isa,jsa,ze,k,k0,V,H,Bsm)
use rundata

implicit none


character :: UPLO,TRANS
integer :: nrow,ncol,m,k,k0
integer, dimension(*) :: isa,jsa
complex (kind=kind(0.0d0)) :: ze
complex (kind=kind(0.0d0)), dimension(*) :: dsa
complex (kind=kind(0.0d0)), dimension(min(nrow,ncol),*) :: V
complex (kind=kind(0.0d0)), dimension((k+1)*m,*) ::H,Bsm

integer :: nnza

complex (kind=kind(0.0d0)), dimension(:,:), allocatable :: Hnew,Vnew,aux
!complex (kind=kind(0.0d0)), dimension(m,m) :: Hnew
!complex (kind=kind(0.0d0)), dimension(ncol,m) :: Vnew
!complex (kind=kind(0.0d0)), dimension(max(nrow,ncol),m) :: aux


!!!BLAS and lapack:
character, dimension(6) :: matdescra
!!lapack stuff:
complex (kind=kind(0.0d0)), dimension(:), allocatable :: work,qrtau
!complex (kind=kind(0.0d0)), dimension(3*ncol) :: work
!complex (kind=kind(0.0d0)), dimension(m) :: qrtau
integer :: lwork,info

integer :: j,i

!double precision, external::elapsed_time
integer :: nmin,nmax

nmin=min(nrow,ncol)
nmax=max(nrow,ncol)


nnza=isa(nrow+1)+1

allocate(Vnew(nmin,m),Hnew(m,m),aux(nmax,m))
lwork=3*nmin
allocate(work(lwork),qrtau(m))

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
    
    !call system_clock(count=tc1)
    !get QR factorization
    call ZGEQRF( nmin, m, V(1,1), nmin, qrtau, work, lwork, info )
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
    call ZUNGQR(  nmin, m, m, V(1,1), nmin, qrtau, work, lwork, info )
    if (info .ne. 0) then
        print *,'Problem with least squares solution in ArnoldiIt'
        print *,'ZUNGQR error info = ',info
        stop
    end if
    !call system_clock(count=tc2)
    !qrtime=qrtime+elapsed_time(tc1,tc2)

end if

!print *,'here',k,k0
!Do matrix multiply:

!Vnew=A*V0(:,(i-1)*m+1:i*m)
call system_clock(count=tc1)
Vnew=V(1:nmin,(k0-1)*m+1:k0*m)
!print *,'here',k,k0

if (ncol==nrow) then
call mkl_zcsrmm(TRANS, nrow, m, ncol, (-1.0d0,0.0d0), matdescra, dsa, jsa, isa, isa(2), V(1,(k0-1)*m+1), ncol, ze, Vnew, nrow)

elseif (nrow>ncol) then
call mkl_zcsrmm('N', nrow, m, ncol, (1.0d0,0.0d0), matdescra, dsa, jsa, isa, isa(2), V(1,(k0-1)*m+1), ncol, (0.0d0,0.0d0), aux, nrow)
call mkl_zcsrmm('T', nrow, m, ncol, -(1.0d0,0.0d0), matdescra, dsa, jsa, isa, isa(2), aux, nrow, ze,Vnew, ncol)

elseif (nrow<ncol) then
call mkl_zcsrmm('T', nrow, m, ncol, (1.0d0,0.0d0), matdescra, dsa, jsa, isa, isa(2), V(1,(k0-1)*m+1), nrow, (0.0d0,0.0d0), aux, ncol)
call mkl_zcsrmm('N', nrow, m, ncol, -(1.0d0,0.0d0), matdescra, dsa, jsa, isa, isa(2), aux, ncol, ze,Vnew, nrow)
   
endif   


!print *,'here',k,k0
!call mkl_zcsrmm('N', n, m, n, (1.0d0,0.0d0), matdescra, dsa, jsa, isa, isa(2), V(1,(k0-1)*m+1), n, (0.0d0,0.0d0), Vnew, n)
!call system_clock(count=tc2)
!mvtime=mvtime+elapsed_time(tc1,tc2)
!!nmatvec(feastit,cpnum)=nmatvec(feastit,cpnum)+m


!call system_clock(count=tc1)
do j=1,k0 !Orthogonalize with respect to previous basis vectors:
!print *,j,k0
   !Hnew=V(:,(j-1)*m+1:j*m)'*Vnew
    call zgemm('C','N',m,m,nmin,(1.0d0,0.0d0),V(1,(j-1)*m+1),nmin,Vnew,nmin,(0.0d0,0.0d0),Hnew,m)
!print *,j,k0
    !H((j-1)*m+1:j*m,(k0-1)*m+1:k0*m)=Hnew
    H((j-1)*m+1:j*m,(k0-1)*m+1:k0*m)=Hnew
    !Vnew=Vnew-V(:,(j-1)*m+1:j*m)*Hnew

    call zgemm('N','N',nmin,m,m,(-1.0d0,0.0d0),V(1,(j-1)*m+1),nmin,Hnew,m,(1.0d0,0.0d0),Vnew,nmin)
!print *,j,k0
 end do
 
! print *,'hhere',k,k0
!call system_clock(count=tc2)
!gstime=gstime+elapsed_time(tc1,tc2)


!Use QR to orthonormalize new vectors:

!call system_clock(count=tc1)
!get QR factorization
call ZGEQRF( nmin, m, Vnew, nmin, qrtau, work, lwork, info )
if (info .ne. 0) then
    print *,'Problem with least squares solution in blockArnoldiIt'
    print *,'ZGEQRF error info = ',info
    stop
end if

!print *,'hereop',k,k0

!put R matrix into H:
!H(k0*m+1:(k0+1)*m,(k0-1)*m+1:k0*m)=Vnew(1:m,1:m)
do i=1,m
    do j=i,m
        H(k0*m+i,(k0-1)*m+j)=Vnew(i,j)
    end do
end do

!print *,'herer',k,k0
!put Q matrix into V:
call ZUNGQR(  nmin, m, m, Vnew, nmin, qrtau, work, lwork, info )
if (info .ne. 0) then
    print *,'Problem with least squares solution in ArnoldiIt'
    print *,'ZUNGQR error info = ',info
    stop
end if

!print *,'herej',k,k0


V(1:nmin,k0*m+1:(k0+1)*m)=Vnew
!call system_clock(count=tc2)
!qrtime=qrtime+elapsed_time(tc1,tc2)
!V(:,k0*m+1:(k0+1)*m)=Q
!H(k0*m+1:(k0+1)*m,(k0-1)*m+1:k0*m)=R(1:m,1:m)
deallocate(Vnew,Hnew,aux)
deallocate(work,qrtau)

end subroutine blockArnoldiIt



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



subroutine blockBGMRESarnoldi(UPLO,TRANS,nrow,ncol,m,dsa,isa,jsa,sb,isb,jsb,ze,kmax,restarts,Brhs,Xlhs,eps,loops,blockstart)
!n is ncol
use rundata

implicit none

character :: UPLO,TRANS
integer :: nrow,ncol,m,kmax,restarts,loops,blockstart
integer, dimension(*) :: isa,jsa,isb,jsb
complex (kind=kind(0.0d0)) :: ze
complex (kind=kind(0.0d0)), dimension(*) :: dsa
double precision, dimension(*) :: sb
complex (kind=kind(0.0d0)), dimension(min(nrow,ncol),*) :: Brhs,Xlhs
double precision :: eps !target residual error

!!!!Arnoldi stuff
complex (kind=kind(0.0d0)), dimension(:,:),allocatable :: V,H,Bsm,ym,Htmp,R,R2,Xtmp

!!!BLAS and lapack:
integer :: lwork,info
complex (kind=kind(0.0d0)), dimension(:),allocatable :: work
double precision, external :: dznrm2

double precision :: error,error2

integer :: i,j,jmax,l

!double precision, external::elapsed_time
integer :: nmin,nmax

nmin=min(nrow,ncol)
nmax=max(nrow,ncol)


!print *,kmax,restarts

allocate(V(nmin,(kmax+1)*m),H(1:(kmax+1)*m,1:kmax*m), Htmp((kmax+1)*m,kmax*m)  ,Bsm((kmax+1)*m,m), ym((kmax+1)*m,m), R(1:nmin,1:m), R2(1:nmin,1:m), Xtmp(1:nmin,1:m) )



lwork=m*(kmax+1)*m*2
allocate(work(lwork))


R=Brhs(1:nmin,1:m)
Xlhs(1:nmin,1:m)=(0.0d0,0.0d0)

loops=0
do i=1,restarts

    if(i>1) R=R2
   

    V(1:nmin,1:m)=R(1:nmin,1:m)

    do j=1,kmax
       !print *,j
        loops=loops+1
        !next arnoldi step
        call system_clock(count=tc1)
        call blockBArnoldiIt(UPLO,TRANS,nrow,ncol,m,dsa,isa,jsa,sb,isb,jsb,ze,kmax,j,V,H,Bsm)
        !call blockArnoldiIt(UPLO,'N',nrow,ncol,m,dsa,isa,jsa,ze,kmax,j,V,H,Bsm)
        
        call system_clock(count=tc2)
        !arnolditime=arnolditime+elapsed_time(tc1,tc2)

        !solve system H(1:(j+1)*m,1:j*m)ym=Bsm(1:(j+1)*m,1:m)
        ym(1:(j+1)*m,1:m)=Bsm(1:(j+1)*m,1:m)
        Htmp(1:(j+1)*m,1:j*m)=H(1:(j+1)*m,1:j*m)


        call system_clock(count=tc1)      
        call zgels('N',(j+1)*m,j*m,m,Htmp(1,1),(kmax+1)*m,ym(1,1),(kmax+1)*m,work,lwork,info)
        if(info .ne. 0) then
            print *,'Error in blockGMRESarnoldi'
            print *,'ZGELS error ',info
            stop
        end if
        call system_clock(count=tc2)
        !lstime=lstime+elapsed_time(tc1,tc2)

       
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        !measure residual:
        !error=norm(R-V[1:n,1:(j+1)*m]*H[1:(j+1)*m,1:j*m]*ym)/norm(Brhs)
        R2=R
        call zgemm('N','N',(j+1)*m,m,j*m,(1.0d0,0.0d0),H(1,1),(kmax+1)*m,ym(1,1),(kmax+1)*m,(0.0d0,0.0d0),Htmp(1,1),(kmax+1)*m)
        call zgemm('N','N',nmin,m,(j+1)*m,(-1.0d0,0.0d0),V(1,1),nmin,Htmp(1,1),(kmax+1)*m,(1.0d0,0.0d0),R2(1,1),nmin)

       
        error=0.0d0
        do l=1,m
            error2=dznrm2(nmin,R2(:,l),1)/dznrm2(nmin,Brhs(:,l),1)
            if (error2>error) error=error2
            !!linres(feastit,cpnum,loops,blockstart+l-1)=error2
        !    print *,'    ',error2
        end do

        print *,i,error

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
    call zgemm('N','N',nmin,m,j*m,(1.0d0,0.0d0),V(1,1),nmin,ym(1,1),(kmax+1)*m,(1.0d0,0.0d0),Xlhs(1,1),nmin)

    if(error<eps) then 
        exit
    end if
end do


deallocate(V,H, Htmp  ,Bsm, ym, R, R2 )

end subroutine blockBGMRESarnoldi



subroutine blockBArnoldiIt(UPLO,TRANS,nrow,ncol,m,dsa,isa,jsa,sb,isb,jsb,ze,k,k0,V,H,Bsm)
use rundata

implicit none


character :: UPLO,TRANS
integer :: nrow,ncol,m,k,k0
integer, dimension(*) :: isa,jsa,isb,jsb
complex (kind=kind(0.0d0)) :: ze
complex (kind=kind(0.0d0)), dimension(*) :: dsa
double precision,dimension(*) :: sb
complex (kind=kind(0.0d0)), dimension(min(nrow,ncol),*) :: V
complex (kind=kind(0.0d0)), dimension((k+1)*m,*) ::H,Bsm

integer :: nnza

complex (kind=kind(0.0d0)), dimension(:,:), allocatable :: Hnew,Vnew,aux,Vold,zt,Vold0
!complex (kind=kind(0.0d0)), dimension(m,m) :: Hnew
!complex (kind=kind(0.0d0)), dimension(ncol,m) :: Vnew
!complex (kind=kind(0.0d0)), dimension(max(nrow,ncol),m) :: aux


!!!BLAS and lapack:
character, dimension(6) :: matdescra,matdescrb
!!lapack stuff:
complex (kind=kind(0.0d0)), dimension(:), allocatable :: work,qrtau,zsb
!complex (kind=kind(0.0d0)), dimension(3*ncol) :: work
!complex (kind=kind(0.0d0)), dimension(m) :: qrtau
integer :: lwork,info

double precision :: zmax
integer :: imax

integer :: j,i

!double precision, external::elapsed_time
integer :: nmin,nmax,nnz,n

nmin=min(nrow,ncol)
nmax=max(nrow,ncol)

if (nrow==ncol) n=nrow

nnza=isa(nrow+1)+1

allocate(Vnew(nmin,m),Hnew(m,m),aux(nmax,m))
lwork=3*nmin
allocate(work(lwork),qrtau(m))




if(UPLO=='F') then
    matdescra(1)='G'
else
    matdescra(1)='H'
end if
!matdescra(1)='G'
matdescra(2)=UPLO
matdescra(3)='N'
matdescra(4)='F'

!!!!!!!!!! Bmatrix
if (nrow==ncol) then !! only square matrix
n=nrow
nnz=isb(n+1)-1
allocate(zsb(nnz))
zsb(1:nnz)=sb(1:nnz)*(1.0d0,0.0d0)
allocate(Vold(nmin,m))



if(UPLO=='F') then
    matdescrb(1)='G'
else
    matdescrb(1)='H'
end if
!matdescra(1)='G'
matdescrb(2)=UPLO
matdescrb(3)='N'
matdescrb(4)='F'

end if


if(k0==1) then !initialize everything
    H(1:(k+1)*m,1:k*m)=(0.0d0,0.0d0)    
    Bsm(1:(k+1)*m,1:m)=(0.0d0,0.0d0)

    !get QR factorization
  call ZGEQRF( nmin, m, V(1,1), nmin, qrtau, work, lwork, info )
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
    call ZUNGQR(  nmin, m, m, V(1,1), nmin, qrtau, work, lwork, info )
    if (info .ne. 0) then
        print *,'Problem with least squares solution in ArnoldiIt'
        print *,'ZUNGQR error info = ',info
        stop
    end if

end if

!print *,'here',k,k0
!Do matrix multiply:

!Vnew=A*V0(:,(i-1)*m+1:i*m)
call system_clock(count=tc1)
Vnew=V(1:nmin,(k0-1)*m+1:k0*m)
!print *,'here',k,k0

if (ncol==nrow) then

   
      Vold(1:N,1:m)=Vnew(1:N,1:m)
call mkl_zcsrmm(TRANS, nrow, m, ncol, (1.0d0,0.0d0), matdescrb, zsb, jsb, isb, isb(2), Vold, ncol, (0.0d0,0.0d0), Vnew, nrow)

call mkl_zcsrmm(TRANS, nrow, m, ncol, (-1.0d0,0.0d0), matdescra, dsa, jsa, isa, isa(2), V(1,(k0-1)*m+1), ncol, ze, Vnew, nrow)



elseif (nrow>ncol) then
call mkl_zcsrmm('N', nrow, m, ncol, (1.0d0,0.0d0), matdescra, dsa, jsa, isa, isa(2), V(1,(k0-1)*m+1), ncol, (0.0d0,0.0d0), aux, nrow)
call mkl_zcsrmm('T', nrow, m, ncol, -(1.0d0,0.0d0), matdescra, dsa, jsa, isa, isa(2), aux, nrow, ze,Vnew, ncol)

elseif (nrow<ncol) then
call mkl_zcsrmm('T', nrow, m, ncol, (1.0d0,0.0d0), matdescra, dsa, jsa, isa, isa(2), V(1,(k0-1)*m+1), nrow, (0.0d0,0.0d0), aux, ncol)
call mkl_zcsrmm('N', nrow, m, ncol, -(1.0d0,0.0d0), matdescra, dsa, jsa, isa, isa(2), aux, ncol, ze,Vnew, nrow)
   
endif   


!print *,'here',k,k0
!call mkl_zcsrmm('N', n, m, n, (1.0d0,0.0d0), matdescra, dsa, jsa, isa, isa(2), V(1,(k0-1)*m+1), n, (0.0d0,0.0d0), Vnew, n)
!call system_clock(count=tc2)
!mvtime=mvtime+elapsed_time(tc1,tc2)
!!nmatvec(feastit,cpnum)=nmatvec(feastit,cpnum)+m


!call system_clock(count=tc1)
do j=1,k0 !Orthogonalize with respect to previous basis vectors:
!print *,j,k0
   !Hnew=V(:,(j-1)*m+1:j*m)'*Vnew
    call zgemm('C','N',m,m,nmin,(1.0d0,0.0d0),V(1,(j-1)*m+1),nmin,Vnew,nmin,(0.0d0,0.0d0),Hnew,m)
!print *,j,k0
    !H((j-1)*m+1:j*m,(k0-1)*m+1:k0*m)=Hnew
    H((j-1)*m+1:j*m,(k0-1)*m+1:k0*m)=Hnew
    !Vnew=Vnew-V(:,(j-1)*m+1:j*m)*Hnew

    call zgemm('N','N',nmin,m,m,(-1.0d0,0.0d0),V(1,(j-1)*m+1),nmin,Hnew,m,(1.0d0,0.0d0),Vnew,nmin)
!print *,j,k0
 end do
 
! print *,'hhere',k,k0
!call system_clock(count=tc2)
!gstime=gstime+elapsed_time(tc1,tc2)


!Use QR to orthonormalize new vectors:

!call system_clock(count=tc1)
!get QR factorization
call ZGEQRF( nmin, m, Vnew, nmin, qrtau, work, lwork, info )
if (info .ne. 0) then
    print *,'Problem with least squares solution in blockArnoldiIt'
    print *,'ZGEQRF error info = ',info
    stop
end if

!print *,'hereop',k,k0

!put R matrix into H:
!H(k0*m+1:(k0+1)*m,(k0-1)*m+1:k0*m)=Vnew(1:m,1:m)
do i=1,m
    do j=i,m
        H(k0*m+i,(k0-1)*m+j)=Vnew(i,j)
    end do
end do

!print *,'herer',k,k0
!put Q matrix into V:
call ZUNGQR(  nmin, m, m, Vnew, nmin, qrtau, work, lwork, info )
if (info .ne. 0) then
    print *,'Problem with least squares solution in ArnoldiIt'
    print *,'ZUNGQR error info = ',info
    stop
end if



V(1:nmin,k0*m+1:(k0+1)*m)=Vnew
!call system_clock(count=tc2)
!qrtime=qrtime+elapsed_time(tc1,tc2)
!V(:,k0*m+1:(k0+1)*m)=Q
!H(k0*m+1:(k0+1)*m,(k0-1)*m+1:k0*m)=R(1:m,1:m)
deallocate(Vnew,Hnew,aux)
deallocate(work,qrtau)

if (nrow==ncol) then
deallocate(zsb)
deallocate(Vold)
endif

end subroutine blockBArnoldiIt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!














!!$
!!$subroutine zfeast_cglsRes(UPLO,n,m,dsa,isa,jsa,ze,nnza,B,X,maxit,eps,neigs,error,its)
!!$use rundata
!!$implicit none
!!$!A=Az=(ze*I-A) in this routine
!!$
!!$    integer :: n,m,maxit,neigs
!!$    complex (kind=kind(0.0d0)) :: ze
!!$    double precision :: eps,error  !target error, output error
!!$    character :: UPLO
!!$
!!$    !!!!!!!!!!!!!!!!!!!!!!!!  Sparse matrix:
!!$    complex (kind=kind(0.0d0)),dimension(*) :: dsa
!!$    integer,dimension(*) :: isa,jsa
!!$    integer :: nnza
!!$
!!$    !!! RHS, solution
!!$    complex (kind=kind(0.0d0)), dimension(n,*) :: B,X
!!$
!!$    !!! CG stuff
!!$    complex (kind=kind(0.0d0)), dimension(:,:), allocatable :: R,Rnew,P,lambda,psi,T,D
!!$    complex (kind=kind(0.0d0)), dimension(:,:), allocatable :: temp1,temp2,sqtemp1,sqtemp2
!!$    
!!$    !!!BLAS and lapack:
!!$    character, dimension(6) :: matdescra
!!$    integer :: info
!!$    integer, dimension(m) :: ipiv
!!$
!!$    integer :: i,j,debug,its
!!$    double precision :: dtemp
!!$    double precision, external :: dznrm2
!!$    double precision, dimension(:),allocatable :: errorlist
!!$    double precision :: errorprint
!!$
!!$    if (m<neigs) then
!!$        print *,'zfeast_cglsRes error: neigs > m',neigs,m
!!$        stop
!!$    endif
!!$
!!$    allocate(errorlist(m))
!!$
!!$    debug=0
!!$
!!$    if(UPLO=='F') then
!!$        matdescra(1)='G'
!!$    else
!!$        matdescra(1)='H'
!!$    end if
!!$    matdescra(2)=UPLO
!!$    matdescra(3)='N'
!!$    matdescra(4)='F'
!!$
!!$    !all this allocating is probably slow; maybe have user allocate once and for all?
!!$    allocate(R(n,m),Rnew(n,m),P(n,m),lambda(m,m),psi(m,m),temp1(n,m),sqtemp1(m,m),sqtemp2(m,m),temp2(n,m),D(n,m),T(n,m))
!!$
!!$    !X=0.0
!!$    X(1:n,1:m)=0.0d0
!!$
!!$    D=B(1:n,1:m)
!!$
!!$    !R=A'*B
!!$    call system_clock(count=tc1)
!!$    R=B(1:n,1:m)
!!$    call mkl_zcsrmm('C', n, m, n, (-1.0d0,0.0d0), matdescra, dsa, jsa, isa, isa(2), B, n, (1.0d0,0.0d0)*conjg(ze), R, n)
!!$    !call mkl_zcsrmm('C', n, m, n, (1.0d0,0.0d0), matdescra, dsa, jsa, isa, isa(2), B, n, (0.0d0,0.0d0), R, n)
!!$    !R(1:n,1:m)=-1.0*B(1:n,1:m)
!!$    call system_clock(count=tc2)
!!$    nmatvec(feastit,cpnum)=nmatvec(feastit,cpnum)+m
!!$    mvtime=mvtime+elapsed_time(tc1,tc2)
!!$    !call mkl_zcsrmm('N', n, m, n, 1.0d0, matdescra, dsa, jsa, isa, isa(2), X, n, 0.0d0, linwork2, n)
!!$    
!!$    !P=-R
!!$    P=R
!!$
!!$    if(debug>0) then
!!$        print *,'X=',X(1:n,1)
!!$        print *,'P=',P(1:n,1)
!!$        print *,'R=',R(1:n,1)
!!$    end if
!!$
!!$    its=0
!!$
!!$    do i=1,maxit
!!$        its=its+1
!!$        
!!$        !lambda=inv(P'*A'*A*P)*R'*R
!!$        !-----T=A*P
!!$        
!!$        call system_clock(count=tc1)
!!$        T=P
!!$        call mkl_zcsrmm('N', n, m, n, (-1.0d0,0.0d0), matdescra, dsa, jsa, isa, isa(2), P, n, ze, T, n) 
!!$        !call mkl_zcsrmm('N', n, m, n, (1.0d0,0.0d0), matdescra, dsa, jsa, isa, isa(2), P, n, (0.0d0,0.0d0), T, n)
!!$        call system_clock(count=tc2)
!!$        nmatvec(feastit,cpnum)=nmatvec(feastit,cpnum)+m
!!$        mvtime=mvtime+elapsed_time(tc1,tc2)
!!$
!!$        call system_clock(count=tc1)
!!$        !-----sqtemp=T'*T
!!$        call zgemm('C','N',m,m,n,(1.0d0,0.0d0),T,n,T,n,(0.0d0,0.0d0),sqtemp1,m)
!!$        !-----lambda=R'*R    !might be better to do (inv(P'A'AP)R')R, not sure...
!!$        call zgemm('C','N',m,m,n,(1.0d0,0.0d0),R,n,R,n,(0.0d0,0.0d0),lambda,m)
!!$        !-----lambda=\(sqtemp1,lambda)
!!$        !call zposv('U',m,m,sqtemp1,m,lambda,m,info)
!!$        call zgesv(m,m,sqtemp1,m,ipiv,lambda,m,info)
!!$        if(info .ne. 0) then
!!$            print *,'CGLS error: ZGESV info ',info
!!$            stop
!!$        end if
!!$
!!$        if(debug>0) print *,'lambda = ', lambda(1,1)
!!$
!!$        !X=X+P*lambda
!!$        call zgemm('N','N',n,m,m,(1.0d0,0.0d0),P,n,lambda,m,(1.0d0,0.0d0),X,n)
!!$        
!!$        if(debug>0) print *,'Xnew = ',X(1:n,1)
!!$
!!$        !D=D-T*lambda
!!$        call zgemm('N','N',n,m,m,(-1.0d0,0.0d0),T,n,lambda,m,(1.0d0,0.0d0),D,n)
!!$        call system_clock(count=tc2)
!!$        gstime=gstime+elapsed_time(tc1,tc2)
!!$
!!$        !Rnew=A'*D 
!!$        
!!$        call system_clock(count=tc1)
!!$        Rnew=D
!!$        call mkl_zcsrmm('C', n, m, n, (-1.0d0,0.0d0), matdescra, dsa, jsa, isa, isa(2), D, n, conjg(ze), Rnew, n)
!!$        !call mkl_zcsrmm('C', n, m, n, (1.0d0,0.0d0), matdescra, dsa, jsa, isa, isa(2), D, n, (0.0d0,0.0d0), Rnew, n)
!!$        call system_clock(count=tc2)
!!$        nmatvec(feastit,cpnum)=nmatvec(feastit,cpnum)+m
!!$        mvtime=mvtime+elapsed_time(tc1,tc2)
!!$        !measure residual
!!$        error=0.0d0
!!$        do j=1,m
!!$            errorlist(j)=dznrm2(n,D(:,j),1)/dznrm2(n,B(:,j),1)
!!$        end do
!!$
!!$        call quicksort(errorlist,1,m)
!!$        !call selection_sort(errorlist,m)
!!$        !error=errorlist(neigs)
!!$        error=errorlist(m)!errorlist(neigs)!(m)
!!$        errorprint=errorlist(neigs)
!!$
!!$        !print *,error
!!$
!!$        if(error<eps) exit !if error is low enough, end loop
!!$
!!$        if(debug>0) print *,'Rnew = ', Rnew(1:n,1)
!!$
!!$        call system_clock(count=tc1)
!!$        !psi=inv(R'*R)*Rnew'*Rnew
!!$        !-----sqtemp1=R'*R
!!$        call zgemm('C','N',m,m,n,(1.0d0,0.0d0),R,n,R,n,(0.0d0,0.0d0),sqtemp1,m)
!!$        !-----sqtemp2=Rnew'*Rnew
!!$        call zgemm('C','N',m,m,n,(1.0d0,0.0d0),Rnew,n,Rnew,n,(0.0d0,0.0d0),psi,m)
!!$        !-----psi=\(sqtemp1,psi)
!!$        !call zposv('U',m,m,sqtemp1,m,psi,m,info)
!!$        call zgesv(m,m,sqtemp1,m,ipiv,psi,m,info)
!!$        if(info .ne. 0) then
!!$            print *,'CGLS error: second ZGESV info ',info
!!$            stop
!!$        end if
!!$
!!$        if(debug>0) print *,'psi = ', psi(1,1)
!!$
!!$        !P=Rnew+P*psi
!!$        temp1=P
!!$        P=Rnew
!!$        call zgemm('N','N',n,m,m,(1.0d0,0.0d0),temp1,n,psi,m,(1.0d0,0.0d0),P,n)
!!$        call system_clock(count=tc2)
!!$        gstime=gstime+elapsed_time(tc1,tc2)
!!$
!!$        if(debug>0) print *,'Pnew = ', P(1:n,1)
!!$
!!$        R=Rnew
!!$       
!!$        !temp1=X(1:n,1:m)
!!$        !call mkl_zcsrmm('N', n, m, n, (-1.0d0,0.0d0), matdescra, dsa, jsa, isa, isa(2), X, n, ze, temp1, n)
!!$        !temp2=B(1:n,1:m)-temp1
!!$        !error=0.0d0
!!$        !do j=1,m
!!$        !    dtemp=dznrm2(n,temp2(1:n,j),1)/dznrm2(n,B(1:n,j),1)
!!$        !   !dtemp=dznrm2(n,R(1:n,j),1)/dznrm2(n,B(1:n,j),1)
!!$        !    if (dtemp>error) error=dtemp
!!$        !end do
!!$        !print *,i,error
!!$
!!$        if(debug>0 .and. i>1) stop
!!$    end do  
!!$    
!!$    !print *,'   linits=',its
!!$    !print *,'      errors=',error,errorprint
!!$    if (error<errorprint) then
!!$        !do i=1,m
!!$        !print *,'   ',errorlist(i)
!!$        !end do
!!$    end if
!!$
!!$    !measure actual error:
!!$    !temp1=X(1:n,1:m)
!!$    !call mkl_zcsrmm('N', n, m, n, (-1.0d0,0.0d0), matdescra, dsa, jsa, isa, isa(2), X, n, ze, temp1, n)
!!$    !temp2=B(1:n,1:m)-temp1
!!$    !do j=1,m
!!$    !    errorlist(j)=dznrm2(n,temp2(1:n,j),1)/dznrm2(n,B(1:n,j),1)
!!$    !end do
!!$    !call quicksort(errorlist,1,m)
!!$    !errorprint=errorlist(neigs)
!!$    !error=errorlist(m)
!!$    !print *,'      acterr=',error,errorprint
!!$
!!$    !stop
!!$
!!$    deallocate(errorlist,R,Rnew,P,lambda,psi,temp1,sqtemp1,sqtemp2,temp2,T,D)
!!$end subroutine zfeast_cglsRes
!!$
!!$
!!$
!!$recursive subroutine quicksort(list,lo,hi)
!!$    double precision, dimension(*), intent(inout) :: list
!!$    integer, intent(in) :: lo,hi
!!$
!!$    integer :: i,j
!!$    double precision :: pivot,temp
!!$
!!$    if(lo<hi) then
!!$    
!!$    pivot=list(hi)
!!$    
!!$    i=lo
!!$    do j=lo,hi-1
!!$        if (list(j)<=pivot) then
!!$            temp=list(i)
!!$            list(i)=list(j)
!!$            list(j)=temp
!!$            i=i+1
!!$        end if
!!$    end do
!!$    temp=list(i)
!!$    list(i)=list(hi)
!!$    list(hi)=temp
!!$    
!!$    call quicksort(list,lo,i-1)
!!$    call quicksort(list,i+1,hi)
!!$    end if
!!$
!!$end subroutine quicksort
!!$



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
 


subroutine zminres(UPLO,nrow,ncol,dsa,isa,jsa,ze,b,x,tol,maxit,linloops)

use rundata

implicit none

integer :: i

character :: UPLO
integer :: nrow,ncol,maxit,linloops
integer, dimension(*) :: isa,jsa
complex (kind=kind(0.0d0)) :: ze !complex shift
complex (kind=kind(0.0d0)), dimension(*) :: dsa !unshifted matrix
complex (kind=kind(0.0d0)), dimension(min(nrow,ncol)) :: b,x
double precision :: tol

!!!minres stuff:
complex (kind=kind(0.0d0)), dimension(:),allocatable :: r,v,vold,p,pold,poldold,aux
double precision :: beta,betaold,gtilde,gamma,rnrm,bnrm
complex (kind=kind(0.0d0)) :: c,cm1,s,sm1,expphi,eta1,alpha,epsilon,delta,h,eta
complex (kind=kind(0.0d0)), dimension(:),allocatable :: actr,Tv,Ap,Apold,Apoldold !residual stuff
!complex (kind=kind(0.0d0)), dimension(:),allocatable :: x !solution
!!!BLAS and lapack:
character, dimension(6) :: matdescra
double precision, external :: dznrm2
complex (kind=kind(0.0d0)), external :: zdotc
integer :: nmin,nmax

nmin=min(nrow,ncol)
nmax=max(nrow,ncol)


    allocate(r(nmin),v(nmin),vold(nmin),p(nmin),pold(nmin),poldold(nmin),aux(nmax))
    allocate(actr(nmin),Tv(nmin),Ap(nmin),Apold(nmin),Apoldold(nmin))
    !allocate(x(n))

    bnrm=dznrm2(nmin,b,1)

    x=(0.0d0,0.0d0)
    r=b

    v=(0.0d0,0.0d0)
    vold=v
    p=v
    pold=v

    beta=dznrm2(nmin,r,1)!beta=norm(r)
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
    matdescra(2)=UPLO
    matdescra(3)='N'
    matdescra(4)='F'

    do i=1,maxit
        linloops=i
        vold=v
        v=r/beta

        !Tv=T*v = A*v+re(ze)*v:
        Tv=v
if (nrow==ncol) then
        call mkl_zcsrmm('N', nrow, 1, ncol, (-1.0d0,0.0d0), matdescra, dsa, jsa, isa, isa(2), v, nrow, (1.0d0,0.0d0)*dble(ze), Tv, nrow)
elseif (nrow>ncol) then
   call mkl_zcsrmm('N', nrow, 1, ncol, (1.0d0,0.0d0), matdescra, dsa, jsa, isa, isa(2), v, ncol, (0.0d0,0.0d0), aux, nrow)
   call mkl_zcsrmm('T', nrow, 1, ncol, -(1.0d0,0.0d0), matdescra, dsa, jsa, isa, isa(2), aux, nrow, (1.0d0,0.0d0)*dble(ze),Tv, ncol)
elseif (nrow<ncol) then
   call mkl_zcsrmm('T', nrow, 1, ncol, (1.0d0,0.0d0), matdescra, dsa, jsa, isa, isa(2), v, nrow, (0.0d0,0.0d0), aux, ncol)
   call mkl_zcsrmm('N', nrow, 1, ncol, -(1.0d0,0.0d0), matdescra, dsa, jsa, isa, isa(2), aux, ncol, (1.0d0,0.0d0)*dble(ze),Tv, nrow)
   endif
        alpha=zdotc(nmin,Tv,1,v,1)!alpha = Tv'*v
        r=Tv-alpha*v-beta*vold

        betaold=beta
        beta=dznrm2(nmin,r,1) !beta=norm(r)

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

        !update residual to see if we're done
        Apoldold=Apold
        Apold=Ap
        Ap=(Tv+(0.0d0,1.0d0)*aimag(ze)*v-delta*Apold-epsilon*Apoldold)/gamma
        actr=actr-eta*Ap

        eta1=(-1.0d0,0.0d0)*s*eta1*expphi

        rnrm=dznrm2(nmin,actr,1)
        if(rnrm<tol) then
            exit
        end if
    end do

    !b=x
   ! print *,'residual',rnrm
end subroutine zminres




subroutine zminresBlock(UPLO,TRANS,nrow,ncol,m,dsa,isa,jsa,sb,isb,jsb,ze,b,x,tol,maxit,linloops,blockstart)

implicit none

integer :: i,j,k

character :: UPLO,TRANS
integer :: nrow,ncol,m,maxit,linloops,blockstart
integer, dimension(*) :: isa,jsa,isb,jsb
complex (kind=kind(0.0d0)) :: ze !complex shift
complex (kind=kind(0.0d0)), dimension(*) :: dsa !unshifted matrix
double precision,dimension(*) :: sb
complex (kind=kind(0.0d0)), dimension(min(nrow,ncol),*) :: b,x
double precision :: tol

!!!minres stuff
double precision :: normr,normb
complex (kind=kind(0.0d0)), dimension(:,:),allocatable :: templ1,templ2 !(n,m)
complex (kind=kind(0.0d0)), dimension(:,:),allocatable :: templ1t !(m,n)
complex (kind=kind(0.0d0)), dimension(:,:),allocatable :: tempm1 !(m,m)

complex (kind=kind(0.0d0)), dimension(:,:),allocatable :: r,Qold,Q,Qnew,P,phi,phiold,phioldold,aux !(n,m)
complex (kind=kind(0.0d0)), dimension(:,:),allocatable :: tau,tautilde,Mm,HR,T,Tnew,Rtilde,Rold,Roldold !(m,m)

complex (kind=kind(0.0d0)), dimension(:,:),allocatable :: ag,bg,cg,dg,agold,bgold,cgold,dgold,dgoldold,bgoldold !(m,m)
complex (kind=kind(0.0d0)), dimension(:,:),allocatable :: QM,HM,RT !(2*m,2*m)
!!!lapack and blas stuff
complex (kind=kind(0.0d0)), dimension(:),allocatable :: qrtau1,zsb !(m)
complex (kind=kind(0.0d0)), dimension(:),allocatable :: qrtau2 !(2*m)
integer :: lwork,info
integer, dimension(:),allocatable :: ipiv !(m)
complex (kind=kind(0.0d0)), dimension(:),allocatable :: work !(3*n)
double precision, dimension(:),allocatable :: dwork !(n)
double precision, external :: zlange,dznrm2
character, dimension(6) :: matdescra,matdescrb
integer :: nmin,nmax,nnz,n

nmin=min(nrow,ncol)
nmax=max(nrow,ncol)


!!!!!!!!!!!!!!!!!!!!!!!! only square matrix
n=nrow
nnz=isb(n+1)-1
allocate(zsb(nnz))
zsb(1:nnz)=sb(1:nnz)*(1.0d0,0.0d0)

if(UPLO=='F') then
    matdescrb(1)='G'
else
    matdescrb(1)='H'
end if
!matdescra(1)='G'
matdescrb(2)=UPLO
matdescrb(3)='N'
matdescrb(4)='F'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    allocate(templ1(nmin,m),templ2(nmin,m),templ1t(m,nmin))
    allocate(r(nmin,m),Qold(nmin,m),Q(nmin,m),Qnew(nmin,m),P(nmin,m),phi(nmin,m),phiold(nmin,m),phioldold(nmin,m),aux(nmax,m))

    allocate(tempm1(m,m),tau(m,m),tautilde(m,m),Mm(m,m),HR(m,m),T(m,m),Tnew(m,m),Rtilde(m,m),Rold(m,m),Roldold(m,m))

    allocate(ag(m,m),bg(m,m),cg(m,m),dg(m,m),agold(m,m),bgold(m,m),cgold(m,m),dgold(m,m),dgoldold(m,m),bgoldold(m,m))

    allocate(QM(2*m,2*m),HM(2*m,2*m),RT(2*m,2*m))

    allocate(qrtau1(m),qrtau2(2*m),work(3*nmin),dwork(nmin),ipiv(m))

    

    normb=zlange('F',nmin,m,b,nmin,dwork)
    !call zmatnorm(n,m,b,normb)

    lwork=3*nmin

    x(1:nmin,1:m)=(0.0d0,0.0d0)
    r=b(1:nmin,1:m) !residual, start with x=0.0d0


    !!QR of r!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !(Q,T)=qr(r)
    !tautilde=T
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Q=r
    call zgeqrf(nmin,m,Q,nmin,qrtau1,work,lwork,info)
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

    call zungqr(nmin,m,m,Q,nmin,qrtau1,work,lwork,info)
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

    linloops=k        

        !!!!!!M=Q'A*Q!!!!!!!!!!!!!!!!!!!! 
        !P=zQ-A*Q
    !P=Q
    P=(0.0d0,0.0d0)
    if (nrow==ncol) then
      
        call mkl_zcsrmm(TRANS, nmin, m, nmin, (1.0d0,0.0d0), matdescrb, zsb, jsb, isb, isb(2), Q, nmin, (0.0d0,0.0d0), P, nmin)

        call mkl_zcsrmm(TRANS, nmin, m, nmin, (-1.0d0,0.0d0), matdescra, dsa, jsa, isa, isa(2), Q, nmin, ze, P, nmin)



     elseif (nrow>ncol) then
   call mkl_zcsrmm('N', nrow, m, ncol, (1.0d0,0.0d0), matdescra, dsa, jsa, isa, isa(2), Q, ncol, (0.0d0,0.0d0), aux, nrow)
   call mkl_zcsrmm('T', nrow, m, ncol, -(1.0d0,0.0d0), matdescra, dsa, jsa, isa, isa(2), aux, nrow,ze,P, ncol)
elseif (nrow<ncol) then
   call mkl_zcsrmm('T', nrow, m, ncol, (1.0d0,0.0d0), matdescra, dsa, jsa, isa, isa(2), Q, nrow, (0.0d0,0.0d0), aux, ncol)
   call mkl_zcsrmm('N', nrow, m, ncol, -(1.0d0,0.0d0), matdescra, dsa, jsa, isa, isa(2), aux, ncol,ze,P, nrow)
   endif


        !M=Q'*P
        call zgemm('C','N',m,m,nmin,(1.0d0,0.0d0),Q,nmin,P,nmin,(0.0d0,0.0d0),Mm,m)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        !!!!!!P=A*Q-Q*M-Qold*T'!!!!!!!!!!       
        !did P=A*Q above
        !P=P-Q*M
        call zgemm('N','N',nmin,m,m,(-1.0d0,0.0d0),Q,nmin,Mm,m,(1.0d0,0.0d0),P,nmin)
        !P=P-Qold*T'
        call zgemm('N','C',nmin,m,m,(-1.0d0,0.0d0),Qold,nmin,T,m,(1.0d0,0.0d0),P,nmin)  
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



        !!!!!(Qnew,Tnew)=qr(P)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        Qnew=P
        call zgeqrf(nmin,m,Qnew,nmin,qrtau1,work,lwork,info)
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
   
        call zungqr(nmin,m,m,Qnew,nmin,qrtau1,work,lwork,info)
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


        ag=HM(1:m,1:m)
        cg=HM(m+1:2*m,1:m)
        bg=HM(1:m,m+1:2*m)
        dg=HM(m+1:2*m,m+1:2*m)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !phi=(Q-phiold*Rold-phioldold*Roldold)*inv(HR[1:m,1:m])
        templ1t=conjg(transpose(Q))
        !templ1t=templ1-Rold'*phiold'
        call zgemm('C','C',m,nmin,m,(-1.0d0,0.0d0),Rold,m,phiold,nmin,(1.0d0,0.0d0),templ1t,m)
        !templ1t=templ1t-Roldold'*phioldold'
        call zgemm('C','C',m,nmin,m,(-1.0d0,0.0d0),Roldold,m,phioldold,nmin,(1.0d0,0.0d0),templ1t,m)
        !solve system Rkk'*phi'=templ1t
        
        !call zgemm('N','N',m,n,m,(1.0d0,0.0d0),tempm1,m,templ1t,m,(0.0d0,0.0d0),templ2t,m) 
        call zgesv(m,nmin,HR,m,ipiv,templ1t,m,info)
        if(info .ne. 0) then
            print *,'block minres zgesv error ',info
            stop
        end if
        phi=conjg(transpose(templ1t))
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !tau=ag*tautilde
        call zgemm('N','N',m,m,m,(1.0d0,0.0d0),ag,m,tautilde,m,(0.0d0,0.0d0),tau,m)

        !x=x+phi*tau
        call zgemm('N','N',nmin,m,m,(1.0d0,0.0d0),phi,nmin,tau,m,(1.0d0,0.0d0),x,nmin)
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
    deallocate(r,Qold,Q,Qnew,P,phi,phiold,phioldold,aux)

    deallocate(tempm1,tau,tautilde,Mm,HR,T,Tnew,Rtilde,Rold,Roldold)

    deallocate(ag,bg,cg,dg,agold,bgold,cgold,dgold,dgoldold,bgoldold)

    deallocate(QM,HM,RT)

    deallocate(qrtau1,qrtau2,work,dwork,ipiv)
deallocate(zsb)

  end subroutine zminresBlock























  



!!$
!!$
!!$
!!$subroutine blockGMRESarnoldi(UPLO,nrow,ncol,m,dsa,isa,jsa,ze,kmax,restarts,Brhs,Xlhs,eps,loops,blockstart)
!!$!n is ncol
!!$implicit none
!!$include 'rundata.f90'
!!$
!!$character :: UPLO
!!$integer :: nrow,ncol,m,kmax,restarts,loops,blockstart
!!$integer, dimension(*) :: isa,jsa
!!$complex (kind=kind(0.0d0)) :: ze
!!$complex (kind=kind(0.0d0)), dimension(*) :: dsa
!!$complex (kind=kind(0.0d0)), dimension(min(nrow,ncol),*) :: Brhs,Xlhs
!!$double precision :: eps !target residual error
!!$
!!$!!!!Arnoldi stuff
!!$complex (kind=kind(0.0d0)), dimension(:,:),allocatable :: V,H,Bsm,ym,Htmp,R,R2,Xtmp
!!$
!!$!!!BLAS and lapack:
!!$character, dimension(6) :: matdescra
!!$integer :: lwork,info
!!$complex (kind=kind(0.0d0)), dimension(:),allocatable :: work
!!$double precision, external :: dznrm2
!!$
!!$double precision :: error,error2
!!$
!!$integer :: i,j,jmax,l
!!$
!!$!double precision, external::elapsed_time
!!$integer :: nmin,nmax
!!$
!!$nmin=min(nrow,ncol)
!!$nmax=max(nrow,ncol)
!!$
!!$
!!$
!!$
!!$
!!$
!!$
!!$
!!$
!!$
!!$enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine dfeast_gcsrevit(N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
    !  Purpose 
    !  =======
    !  FEAST SPARSE INTERFACE
    !  Solve the standard Ax=ex eigenvalue problem
    !  
    !  A GENERAL REAL :: SPARSE FORMAT 
    ! 
    !  DOUBLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !
    !  N          (input)        INTEGER: Size system
    !  sa         (input)        REAL DOUBLE PRECISION (isa(N+1)-1):  Matrix A- CSR format 
    !  isa        (input)        INTEGER(N+1): CSR row array of Matrix A
    !  jsa        (input)        INTEGER(isa(N+1)-1): CSR column array of Matrix A
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  fpm (input/output)        INTEGER(*) : FEAST parameters
    !  epsout     (output)       REAL DOUBLE PRECISION : Error on the trace
    !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
    !  Emid       (input/output) COMPLEX DOUBLE PRECISION: middle of the search interval
    !  r          (input/output) REAL DOUBLE PRECISION: horizontal radius of the search interval
    !                            comments: (input)  Emid,r dummy variables
    !                                      (output) Emid,r calculated/estimated from custom contour
    !                            expert comments: (input)  Emid,r provided by user if default contour i.e. fpm(29)=1
    !                                             (output) Emid,r unchanged if fpm(29)=1 
    !  M0         (input/output) INTEGER: Size subspace
    !  E          (output)       COMPLEX DOUBLE PRECISION(M0)   : Eigenvalues -solution
    !  q          (input/output) COMPLEX DOUBLE PRECISION(N,2*M0) : 
    !                                                       On entry: subspace initial guess if fpm(5)=1 
    !                                                       On exit : Right Eigenvectors-solution (1:N,1:mode)
    !                                                                 Left  Eigenvectors-solution (1:N,M0+1:M0+mode)
    !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
    !  res        (output)       REAL DOUBLE PRECISION(2*M0) : Relative Residual of the solution (1-norm)
    !                                                        if option fpm(6)=1 selected                           
    !  info       (output)       INTEGER: Error handling (0: successful exit)
    !=====================================================================
    ! James Kestyn, Eric Polizzi 2015
    ! ====================================================================
    implicit none
    integer :: N,LDA
    double precision,dimension(*) :: sa
    integer,dimension(*) :: isa,jsa
    integer,dimension(*) :: fpm
    complex(kind=(kind(1.0d0))) :: Emid
    double precision :: r 
    double precision :: epsout
    integer :: loop
    integer :: M0
    complex(kind=(kind(1.0d0))),dimension(*):: E
    complex(kind=(kind(1.0d0))),dimension(N,*):: X
    integer :: mode
    double precision,dimension(*)    :: res
    integer :: info

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!    Wrapper Routine to expert routine: dfeast_gcsrgvx
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    complex(kind=(kind(1.0d0))), dimension(1:fpm(8)) :: Zne,Wne 
    double precision,dimension(1) :: sb ! identity
    integer,dimension(1) :: isb
    integer,dimension(1) :: jsb
    integer :: i

    call zfeast_gcontour(Emid,r,fpm(8),fpm(17),fpm(18),fpm(19),Zne,Wne)

    fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
    ! identity B matrix- option for standard eigenvalue problem
!    do i=1,n
!       sb(i)=1.0d0
!       jsb(i)=i
!       isb(i)=i
!    enddo
!    isb(n+1)=n+1
    call dfeast_gcsrgvxit(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

  end subroutine dfeast_gcsrevit







  subroutine dfeast_gcsrgvxit(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
    !  Purpose 
    !  =======
    !  FEAST SPARSE INTERFACE
    !  Solve the generalized Ax=eBx eigenvalue problem
    !  
    !  A,B GENERAL REAL :: SPARSE FORMAT 
    ! 
    !  DOUBLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !
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
    !  fpm (input/output)        INTEGER(*) : FEAST parameters
    !  epsout     (output)       REAL DOUBLE PRECISION : Error on the trace
    !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
    !  Emid       (input/output) COMPLEX DOUBLE PRECISION: middle of the search interval
    !  r          (input/output) REAL DOUBLE PRECISION: horizontal radius of the search interval
    !                            comments: (input)  Emid,r dummy variables
    !                                      (output) Emid,r calculated/estimated from custom contour
    !                            expert comments: (input)  Emid,r provided by user if default contour i.e. fpm(29)=1
    !                                             (output) Emid,r unchanged if fpm(29)=1 
    !  M0         (input/output) INTEGER: Size subspace
    !  E          (output)       COMPLEX DOUBLE PRECISION(M0)   : Eigenvalues -solution
    !  q          (input/output) COMPLEX DOUBLE PRECISION(N,2*M0) : 
    !                                                       On entry: subspace initial guess if fpm(5)=1 
    !                                                       On exit : Right Eigenvectors-solution (1:N,1:mode)
    !                                                                 Left  Eigenvectors-solution (1:N,M0+1:M0+mode)
    !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
    !  res        (output)       REAL DOUBLE PRECISION(2*M0) : Relative Residual of the solution (1-norm)
    !                                                        if option fpm(6)=1 selected                           
    !  info       (output)       INTEGER: Error handling (0: successful exit)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of ARGUMENTS FOR EXPERT ROUTINE 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  Zne        (input)        COMPLEX DOUBLE PRECISION(fpm(8)): Custom Integration nodes
    !  Wne        (input)        COMPLEX DOUBLE PRECISION(fpm(8)): Custom Integration weights
    !           
    !=====================================================================
    ! James Kestyn, Eric Polizzi 2015
    ! ====================================================================
    implicit none
    include 'f90_noruntime_interface.fi'
    integer :: N
    double precision,dimension(*) :: sa,sb
    integer,dimension(*) :: isa,jsa,isb,jsb
    integer,dimension(*) :: fpm
    complex(kind=(kind(1.0d0))) :: Emid
    double precision :: r 
    double precision :: epsout
    integer :: loop
    integer :: M0
    complex(kind=(kind(1.0d0))),dimension(*):: E
    complex(kind=(kind(1.0d0))),dimension(N,*):: X
    integer :: mode
    double precision,dimension(*)    :: res
    integer :: info
    complex(kind=(kind(1.0d0))),dimension(*):: Zne,Wne
!!!!!!!!!!!!!!!!!!!!!!
    integer :: ijob,infoloc,i,s
    complex(kind=(kind(1.0d0))) :: Ze
    complex(kind=(kind(1.0d0))), dimension(:,:),pointer ::work,workc,zAq,zSq
    complex(kind=(kind(1.0d0))),dimension(:,:),pointer :: saz
    integer,dimension(:),pointer :: isaz,jsaz
    double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
    complex(kind=(kind(1.0d0))),parameter :: ONEC=(DONE,DZERO),ZEROC=(DZERO,DZERO)
    logical :: fact
    integer :: nfact,id
    integer :: rank,code,nb_procs,NEW_COMM_WORLD
!!!! transpose matrix
    !double precision,dimension(:),pointer :: tsa,tsb
    !integer,dimension(:), pointer :: tisa,tjsa,tisb,tjsb
    integer :: opt,nnz
!!!!!!!!!! iterative solver
    complex(kind=(kind(1.0d0))) :: ze2,cze2
    complex(kind=(kind(1.0d0))), dimension(:,:), pointer :: ztempmat
    integer :: linresindex !number of RHS to use in measuring linear system error
    integer :: blockits,blockremainder
    double precision :: lintargeterror,linsyserror !goal error for linear system, return linear system error
    integer :: linloops
    double precision :: linepsout
    integer :: k
character(len=1) :: UPLO,TRANS

    double precision, dimension(M0) :: Er
    
    
!!!MKL mat-vec:
character, dimension(6) :: matdescra
    matdescra(1)='G'
    matdescra(2)='F' ! dummy I guess
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
    IF( N<=0 ) THEN
       INFO = -102
    END IF
    IF( INFO.NE.0 ) THEN
       CALL XERBLA( 'DFEAST_GCSRGVit', -INFO+100 )
       RETURN
    END IF


    infoloc=0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call wallocate_2z(zAq,M0,M0,infoloc)
    call wallocate_2z(zSq,M0,M0,infoloc)
    call wallocate_2z(work,N,2*M0,infoloc)
    call wallocate_2z(workc,N,M0,infoloc)
 call wallocate_2z(ztempmat,N,m0,infoloc)



    
    if (infoloc/=0) then
       info=-1
       return
    end if

!!!!!!!!!!!!!!!! find fpm(22) - # total number of contour points
    fpm(22)=fpm(8) ! default
    if (fpm(29)==1) then! default for generated contour -- half-contour capability case 
       if ((aimag(Emid)==0.0d0).and.((2*(fpm(8)/2))==fpm(8))) fpm(22)=fpm(8)/2  
    end if


!!!!!!!!!!!!!!!!!! Factorizations Set-up  
    fact=.true.
    nfact=0
    !! nfact is local (number of total factorization by contour points)
    do i=rank+1,fpm(22),nb_procs
       nfact=nfact+1
    end do

    if (fpm(10)==1) then
       id=0
    else
       id=1
    end if


!!!!!!!!!!!!!!!!! Set up for Az matrix
    nnz=isa(N+1)-1
    call wallocate_2z(saz,nnz,nfact,infoloc)
    do i=1,nfact
    saz(1:nnz,i)=sa(1:nnz)*(1.0d0,0.0d0)   
 enddo
 
    do i=1,n
       do k=isa(i),isa(i+1)-1
         ! if (jsa(k)==i) saz(k,:)=saz(k,:)+(0.0d0,100.0d0)
          enddo
    enddo
    
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    UPLO='F'
    ijob=-1 ! initialization
    do while (ijob/=0) 
     !  call dfeast_grcix(ijob,N,Ze,work,workc,zAq,zSq,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

       fpm(31)=3 ! for real non-symmetric
       call zfeast_mygrcix(ijob,N,Ze,work,workc,zAq,zSq,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

!       call zfeast_hrcix(ijob,N,Ze,work,workc,zAq,zSq,fpm,epsout,loop,dble(Emid)-r,dble(Emid)+r,M0,Er,X,mode,res,info,Zne,Wne)
!print *,ijob
       
       !call zfeast_hrcix(ijob,N,Ze,work,workc,zAq,zSq,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)
       select case(ijob)

       case(10) !! factorize (zeB-A)

         
          !user iterative solver, don't need factorization
           
          ze2=ze!(0.0d0,100.0d0)*abs(aimag(ze))
!print *,ze          

       case(11) !!solve the linear system (ZeB-A)x=workc(1:N,1:fpm(23)) result in to workc

 
  TRANS='N'
          
    if(loop>0 .and. (epsout .ne. 0.0d0)) then
                  if(fpm(54)==0) then
                      lintargeterror=1.0d-2*epsout
                  else
                      !lintargeterror=epsout*10.0d0**(-1.0d0*fpm(54))
                      lintargeterror=epsout/(1.0d0*fpm(54))

                  end if
               else
                  if (fpm(54)==0) then
                     lintargeterror=1d-2
                     else
                  lintargeterror=1.0d0/fpm(54)!1d-1
               end if
              endif          
                

                !use fpm(53) as block size
                blockits=m0/fpm(53)
                blockremainder=m0-blockits*fpm(53)
                
                do i=1,blockits

!if (loop>0) then
!if (E(i)<0.5d0) then ! outside interval
!!if (i<15.or.i>55) then
!fpm(51)=2 
!if (fpm(62)==1) print *,'leave',loop,i,E(i)
!!cycle
!else 
!fpm(51)=100
!endif
!endif

                    linresindex=fpm(53)
                    if(fpm(11)==1) call blockGMRESarnoldi(UPLO,TRANS,N,N,fpm(53),saz(1,id),isa,jsa,ze2,fpm(51),fpm(50),workc(1,1+(i-1)*fpm(53)),ztempmat(1,1+(i-1)*fpm(53)),lintargeterror,linloops,1+(i-1)*fpm(53))

print *,i,linloops
                 end do
                 
                !print *,'i ',i,blockits
                if(blockremainder .ne. 0) then
                    linresindex=blockremainder
                    if(fpm(11)==1) call blockGMRESarnoldi(UPLO,TRANS,N,N,blockremainder,saz(1,id),isa,jsa,ze2,fpm(51),fpm(50),workc(1,1+(i-1)*fpm(53)),ztempmat(1,1+(i-1)*fpm(53)),lintargeterror,linloops,1+(i-1)*fpm(53))
                 end if

           workc=ztempmat          


        case(20)
cze2=conjg(ze)!+(0.0d0,-100.0d0)*abs(aimag(ze))         
!cze2=ze
!print *,ze  


           
        case(21) !!solve the linear system (ZeB-A)x=workc(1:N,1:fpm(23)) result in to workc

 
  TRANS='C'
  
  if(loop>0 .and. (epsout .ne. 0.0d0)) then
                  if(fpm(54)==0) then
                      lintargeterror=1.0d-2*epsout
                  else
                      !lintargeterror=epsout*10.0d0**(-1.0d0*fpm(54))
                      lintargeterror=epsout/(1.0d0*fpm(54))

                  end if
               else
                   if (fpm(54)==0) then
                     lintargeterror=1d-2
                     else
                  lintargeterror=1.0d0/fpm(54)!1d-1
               end if
                  !lintargeterror=1d-1
              endif          
                

                !use fpm(53) as block size
                blockits=m0/fpm(53)
                blockremainder=m0-blockits*fpm(53)
                
                do i=1,blockits

!if (loop>0) then
!if (E(i)<0.5d0) then ! outside interval
!!if (i<15.or.i>55) then
!fpm(51)=2 
!if (fpm(62)==1) print *,'leave',loop,i,E(i)
!!cycle
!else 
!fpm(51)=100
!endif
!endif

                    linresindex=fpm(53)
                    if(fpm(11)==1) call blockGMRESarnoldi(UPLO,TRANS,N,N,fpm(53),saz(1,id),isa,jsa,cze2,fpm(51),fpm(50),workc(1,1+(i-1)*fpm(53)),ztempmat(1,1+(i-1)*fpm(53)),lintargeterror,linloops,1+(i-1)*fpm(53))

print *,i,linloops
                 end do
                 
                !print *,'i ',i,blockits
                if(blockremainder .ne. 0) then
                    linresindex=blockremainder
                    if(fpm(11)==1) call blockGMRESarnoldi(UPLO,TRANS,N,N,blockremainder,saz(1,id),isa,jsa,cze2,fpm(51),fpm(50),workc(1,1+(i-1)*fpm(53)),ztempmat(1,1+(i-1)*fpm(53)),lintargeterror,linloops,1+(i-1)*fpm(53))
                 end if

           workc=ztempmat          




       case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

          !X(1:N,1)=(1.0d0,1.0d0)
!            call mkl_zcsrmm('N', N, m0, N, (1.0d0,0.0d0), matdescra, saz, jsa, isa, isa(2), X, N, (0.0d0,0.0d0), work, N)
!!$         print *,'mat'
!!$            print *,saz
!!$            print *,'x'
!!$            
!!$            print *,X(1:N,1)
!!$            print *,'sol'
!!$            print *,work(1:N,1)

            print *,'FIX LINE 3381 in iterative_feast.f90'
           stop 
          !call zcsrmm('F','N',N,N,fpm(25),ONEC,saz(1,1),isa,jsa,X(1,fpm(24)),ZEROC,work(1,fpm(24)))

            
         case(31) !! perform multiplication A^H*x(1:N,fpm(34):fpm(34)+fpm(35)-1) result in work(1:N,fpm(34)+fpm(35)-1)

            !X(1:N,1)=(1.0d0,1.0d0)
 !  call mkl_zcsrmm('C', N, m0, N, (1.0d0,0.0d0), matdescra, saz, jsa, isa, isa(2), X(1,M0+1), N, (0.0d0,0.0d0), work(1,M0+1), N)
!!$ print *,'mat'
!!$            print *,saz
!!$            print *,'x'
!!$            
!!$            print *,X(1:N,1)
!!$            print *,'sol'
!!$            print *,work(1:N,1)
!!$
!!$stop

    
            print *,'FIX LINE 3401 in iterative_feast.f90'
           stop           
        !call zcsrmm('F','C',N,N,fpm(35),ONEC,saz(1,1),isa,jsa,X(1,fpm(34)),ZEROC,work(1,fpm(34)))


       case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

          work(1:N,1:m0)=X(1:N,1:m0)

          
!          call dzcsrmm('F','N',N,N,fpm(25),ONEC,sb,isb,jsb,X(1,fpm(24)),ZEROC,work(1,fpm(24)))


       case(41)!! perform multiplication B^H*x(1:N,fpm(34):fpm(34)+fpm(35)-1) result in work(1:N,fpm(34)+fpm(35)-1)

 !         call dzcsrmm('F','C',N,N,fpm(35),ONEC,sb,isb,jsb,X(1,fpm(34)),ZEROC,work(1,fpm(34)))

          work(1:N,m0+1:2*m0)=X(1:N,m0+1:2*m0)


       end select
    end do

!

    call wdeallocate_2z(zAq)
    call wdeallocate_2z(zSq)
    call wdeallocate_2z(work)
    call wdeallocate_2z(workc)
    call wdeallocate_2z(ztempmat)

    call wdeallocate_2z(saz)
   
!E(1:M0)=Er*(1.0d0,0.0d0)   


  end subroutine dfeast_gcsrgvxit



  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine zfeast_mygrcix(ijob,N,Ze,work,workc,zAq,zSq,fpm,epsout,loop,Emid,r,M0,lambda,Q,mode,res,info,Zne,Wne)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Purpose 
  !  =======
  !  FEAST RCI (Reverse Communication Interfaces) -  Includes option for custom integration nodes/weight
  !  Solve generalized eigenvalue problems and obtain eigenvalues lambda, right qr and left ql eigenvectors
  !      Aqr=lambda Bqr and A^Hql=conjg(lambda) B^Hql                            
  ! 
  !  Expert comments: This is also a Kernel routine (contains several cases depending of fpm(31))
  !                   fpm(31)=0 ! complex general (default)  
  !                   fpm(31)=1 ! Real general (contour symmetric via real axis)           
  !                   fpm(31)=2 ! Complex symmetric          
  !                   fpm(31)=3 ! Real general               
  !
  !  A and B can be NON-HERMITIAN (could be real or complex)  
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  ijob       (input/output) INTEGER :: ID of the RCI
  !                            INPUT on first entry: ijob=-1 
  !                            OUTPUT Return values (0,10,20,21,30,31,40,41)-- see FEAST documentation
  !  N          (input)        INTEGER: Size system
  !  Ze         (output)       COMPLEX DOUBLE PRECISION        : Integration node
  !  work       (input/output) COMPLEX DOUBLE PRECISION (N,2*M0):   Workspace
  !                            expert comment: size of (N,M0) if fpm(31)=2 (complex symmetric case)  
  !  workc      (input/output) COMPLEX DOUBLE PRECISION (N,M0):   Workspace 
  !  zAq,zSq    (input/output) COMPLEX DOUBLE PRECISION (M0,M0) : Workspace for Reduced Eigenvalue System
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  fpm     (input/output)    INTEGER(*) : FEAST parameters (see FEAST documentation)
  !  epsout     (output)       REAL DOUBLE PRECISION : Error on the trace
  !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
  !  Emid       (input/output) COMPLEX DOUBLE PRECISION: middle of the search interval
  !  r          (input/output) REAL DOUBLE PRECISION: horizontal radius of the search interval
  !                            comments: (input) Emid,r dummy variables
  !                                      (output) Emid,r calculated/estimated from custom contour
  !                            expert comments: (input)  Emid,r provided by user if default contour i.e. fpm(29)=1
  !                                             (output) Emid,r unchanged if fpm(29)=1 
  !  M0         (input/output) INTEGER: Size subspace
  !  lambda     (output)       COMPLEX DOUBLE PRECISION(M0)   : Eigenvalues -solution
  !  q          (input/output) COMPLEX DOUBLE PRECISION(N,2*M0) : 
  !                                                       On entry: subspace initial guess if fpm(5)=1 
  !                                                       On exit : Right Eigenvectors-solution (1:N,1:mode)
  !                                                                 Left  Eigenvectors-solution (1:N,M0+1:M0+mode)
  !                            expert comment: size of (N,M0) if fpm(31)=2 (complex symmetric case - only right vectors)
  !  
  !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
  !  res        (output)       REAL DOUBLE PRECISION(2*M0) : Relative Residual of the solution (1-norm) 
  !                                                          contains right (1:mode) and left (M0+1:M0+mode)
  !                          
  !                            expert comment: size of (M0) if fpm(31)=2 (complex symmetric case - only rigth vectors)
  !  info       (output)       INTEGER: Error handling (0: successful exit) -- see FEAST documentation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of ARGUMENTS FOR EXPERT ROUTINE 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Zne        (input)        COMPLEX DOUBLE PRECISION(fpm(8)): Custom Integration nodes
  !  Wne        (input)        COMPLEX DOUBLE PRECISION(fpm(8)): Custom Integration weights
  !
  !                            expert comments: if fpm(29)=1, Wne,Zne are generated using default contour
  !                                             important for testing if an eigenvalue is inside an ellipsoid shape vs polygon contour
  !=========================================================================
  !  James Kestyn - Eric Polizzi - 2015
  !=====================================================================

  implicit none
  !-------------------------------------
#ifdef MPI
  include 'mpif.h'
#endif
  !-------------------------------------
  include "f90_noruntime_interface.fi"
  integer :: ijob,N,M0
  complex(kind=(kind(1.0d0))) :: Ze
  complex(kind=(kind(1.0d0))) :: Emid
  double precision :: r
  complex(kind=(kind(1.0d0))),dimension(N,*):: work,workc
  complex(kind=(kind(1.0d0))),dimension(M0,*):: zAq,zSq
  integer,dimension(*) :: fpm
  double precision :: epsout
  integer :: loop
  complex(kind=(kind(1.0d0))),dimension(*):: lambda,Zne,Wne
  complex(kind=(kind(1.0d0))),dimension(N,*):: Q
  integer :: mode
  double precision,dimension(*) :: res
  integer :: info
  !! parameters
  double precision, Parameter :: pi=3.1415926535897932d0
  double precision, Parameter :: DONE=1.0d0, DZERO=0.0d0
  complex(kind=(kind(1.0d0))),parameter :: ONEC=(DONE,DZERO)
  complex(kind=(kind(1.0d0))),parameter :: ZEROC=(DZERO,DZERO)
  double precision, parameter :: ba=-pi/2.0d0,ab=pi/2.0d0
  integer*8,parameter :: fout =6
  !! variable for FEAST
  integer :: i,e,j,k,jj,m_min
  integer,dimension(4) :: iseed
  double precision :: theta
  complex(kind=(kind(1.0d0))) :: jac,aux
  complex(kind=(kind(1.0d0))), dimension(:,:),pointer :: zSqo
  integer, dimension(:),pointer :: fpm_default,color
  logical :: testconv
  complex(kind=(kind(1.0d0))) :: trace, zalpha
  !! Lapack variable 
  double precision, dimension(:),pointer :: work_locr,work_loc2
  complex(kind=(kind(1.0d0))), dimension(:),pointer :: zwork_loc
  complex(kind=(kind(1.0d0))), dimension(:),pointer :: LALPHA,LBETA
  complex(kind=(kind(1.0d0))), dimension(:,:),pointer :: VL,VR
  integer :: lwork_loc,info_lap,infoloc
  character(len=1) :: JOBVL
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

  if ((fpm(29)==1).and.(fpm(31)==3)) then ! default contour
     !Test if Emid is real and fpm(8) even -> factorization on half-contour -- New ID if contour symmetric via real axis
     if ((aimag(Emid)==0.0d0).and.((2*(fpm(8)/2))==fpm(8))) fpm(31)=1 
  end if
  
  if (fpm(14)==2) fpm(31)=2 ! Compute only right contour if estimate is on


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! Initialization!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (ijob==-1) then 
     info=0 ! default value
     if(fpm(29)==0)then ! expert routines Emid,r needs to be calculated
        ! calculate the center of mass and relative maximum length of custom contour
        Emid=ZEROC
        do i=1,fpm(8)
           Emid = Emid + Zne(i)
        enddo
        Emid = Emid / fpm(8)
        r = abs(maxval(abs(Zne(1:fpm(8))))-Emid)
     endif
     if (fpm(1)==1) then
        call wwrite_n(fout)
        call wwrite_s(fout, '***********************************************')  
        call wwrite_n(fout) 
        call wwrite_s(fout, '*********** FEAST- BEGIN **********************')
        call wwrite_n(fout) 
        call wwrite_s(fout, '***********************************************')  
        call wwrite_n(fout)
        if (fpm(31)==3) then
           call wwrite_s(fout, 'Routine DFEAST_G{}{}') !<<<<<<<<
        elseif (fpm(31)==1) then
           call wwrite_s(fout, 'Routine DFEAST_G{}{}') !<<<<<<<<
        elseif (fpm(31)==2) then
           call wwrite_s(fout, 'Routine ZFEAST_S{}{}') !<<<<<<<<
        elseif (fpm(31)==0) then
           call wwrite_s(fout, 'Routine ZFEAST_G{}{}') !<<<<<<<<
        end if
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


        if(fpm(1)==1 .and. rank==0)then
           if (fpm(29)==1) then
           call wwrite_s(fout, 'Search contour [(Emid), r] (')
           else
           call wwrite_s(fout, 'Center of mass for user custom contour (')
           endif
           call wwrite_d(fout,dble(Emid))
           call wwrite_s(fout, ',')
           call wwrite_d(fout,aimag(Emid))
           call wwrite_s(fout, '); ')
           if (fpm(29)==1) call wwrite_d(fout,r)
           call wwrite_n(fout)
        endif
     end if
     call check_feast_fpm_input(fpm,info)
     call dcheck_feast_grci_input(r,M0,N,info)

     if (info/=0) fpm(21)=100 ! The End

!!!!!!!!!!!!!!!!
     IF (info==0) then
        fpm(22)=fpm(8) ! default
        if (fpm(29)==1) then! default generated contour 
          if (fpm(31)==1) fpm(22)=fpm(8)/2  ! fpm(31)=1 -- Emid real--fpm(8) even-   will only have half-contour
        endif  
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
           call wwrite_s(fout, '#Loop | #Eig  |          |Trace|       |     Error-Trace       |     Max-Residual')
           call wwrite_n(fout)
           call wwrite_s(fout, '-----------------------------------------------------------------------------------')
           call wwrite_n(fout)
        endif
        fpm(23)=min(M0,N) ! 'current M0' size (global value)
        fpm(25)=fpm(23) !! 'current M0' size (by default)
        fpm(24)=1  !! origin first column number of vector q for parallel mat-vec (default)
        fpm(35)=fpm(23) !! 'current M0' size (by default)
        fpm(34)=M0+1  !! origin first column number of vector q for parallel mat-vec (default)
        !----------------------------------------------
#ifdef MPI
        if (fpm(23)/nb_procs>=1) then ! criteria for parallelism of mat-vec
           fpm(25)=fpm(23)/nb_procs ! local size of 'current M0'
           if (rank==nb_procs-1) fpm(25)=fpm(25)+fpm(23)-(fpm(23)/nb_procs)*nb_procs 
           fpm(24)=1+rank*(fpm(23)/nb_procs) ! local origin of first column for vector q for parallel mat-vec 
           fpm(35)=fpm(23)/nb_procs ! local size of 'current M0'
           if (rank==nb_procs-1) fpm(35)=fpm(35)+fpm(23)-(fpm(23)/nb_procs)*nb_procs 
           fpm(34)=M0+1+rank*(fpm(23)/nb_procs) ! local origin of first column for vector q for parallel mat-vec 
        end if
#endif
        !-----------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
        fpm(21)=1 ! prepare reentry
        if (fpm(5)==0) then !!! random vectors
           iseed=(/56,890,3456,2333/)
           if (fpm(14)==2) then
              ! copy work into q for multiplication by B matrix, if stochastic estimate  is "on" 
              call ZLARNV(4,iseed,N*fpm(23),q(1,1))
           else
              call ZLARNV(4,iseed,N*fpm(23),work(1,1))
              iseed=(/6,80,456,5421/)
              if (fpm(31)/=2) call ZLARNV(4,iseed,N*fpm(23),work(1,M0+1))
           endif
        end if

        if ((fpm(5)==1).or.(fpm(14)==2)) then !!!!!! q is the initial guess
           !----------------------------------------
#ifdef MPI
           work(1:N,1:fpm(23))=ZEROC 
           if (fpm(31) /= 2) work(1:N,M0+1:M0+fpm(23))=ZEROC 
#endif
           !------------------------------------------
           ijob=40 !! B*Q(1,1:fpm(23))=>work(1:N,1:fpm(23))
           if (fpm(31) /=2 ) fpm(21)=-1
           return
        end if
     end IF ! info=0
!!!!!!!!!!!!!!
  end if   !ijob=-1

  if (fpm(21)==-1) then !! we need to initialize the left eigenvectors as well
     ijob=41!! B^T*Q(1,M0+1:M0+fpm(23))=>work(1:N,M0+1:M0+fpm(23))
     fpm(21)=1 ! prepare reentry
     return
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! CONTOUR INTEGRATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  IF (fpm(21)==1) then !! we initialize a new contour integration
     !------------------------------------------------------------------------
#ifdef MPI
     if ((loop>0).or.(fpm(5)==1).or.(fpm(14)==2)) then        
        if (fpm(23)/nb_procs>=1) call MPI_ALLREDUCE(MPI_IN_PLACE,work(1,1),N*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,NEW_COMM_WORLD,code)
        if ((fpm(31)/=2) .and. (fpm(23)/nb_procs>=1)) call MPI_ALLREDUCE(MPI_IN_PLACE,work(1,M0+1),N*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,NEW_COMM_WORLD,code)
     end if
#endif
     !------------------------------------------------------------------------
     Q(1:N,1:fpm(23))=ZEROC
     if(fpm(31)/=2) Q(1:N,M0+1:M0+fpm(23))=ZEROC
     fpm(20)=1
     fpm(21)=2
     ijob=-2 ! just initialization 
  end IF


!!!!!!!!!!!!
  IF (fpm(21)==2) then !! we start or pursue the contour integration

     IF (info==0) then !! will end up checking info errors returned by FEAST drivers
        do e=fpm(20)+rank,fpm(22),nb_procs !!!! loop over the full contour  

           if (ijob==-2) then !!Factorize the linear system                   
              Ze = Zne(e)   
              fpm(20)=e-rank
              ijob=10 ! for fact
              if ((loop==0).or.(fpm(22)>nb_procs)) return ! no need to factorize again if one linear system per processor 
           end if

           if (ijob==10) then !!Solve the linear system (complex) (zS-A)qr=workc 
              call ZLACPY( 'F', N, fpm(23), work, N, workc, N )
              ijob=11 ! for solve
              return
           end if

           if (ijob==11) then  !! summation
              aux = Wne(e) 
              if (fpm(33)==1) then
                  aux=conjg(aux)
                  workc(1:N,1:fpm(23))= conjg(workc(1:N,1:fpm(23)))
               end if
              call ZAXPY(N*fpm(23),aux,workc,1,Q(1,1),1)

              if(fpm(31)==1) then
                 if (fpm(33)==0) then
                    call ZLACPY( 'F', N, fpm(23),work(1,1) , N, workc, N )
                    workc(1:N,1:fpm(23))= conjg(workc(1:N,1:fpm(23)))
                    fpm(33)=1 
                    ijob=11 ! reuse fact for second solve - must conjugate the solution qr
                    return
                 elseif (fpm(33)==1) then
                    fpm(33)=0
                    ijob=20
                    if ((loop==0).or.(fpm(22)>nb_procs)) return ! no need to factorize again if one linear system per processor
                 end if
              elseif (fpm(31)==2) then
                 ijob=-2 !no solve for ql if complex symmetric
              elseif ((fpm(31)==0).or.(fpm(31)==3)) then 
                 ijob=20 ! for fact
                 if ((loop==0).or.(fpm(22)>nb_procs)) return ! no need to factorize again if one linear system per processor
              end if
           end if !ijob==11

           if (ijob==20) then !!Solve the linear system (complex) (zS-A)^Hql=workc 
              call ZLACPY( 'F', N, fpm(23),work(1,M0+1) , N, workc, N )
              ijob=21 ! for solve
              return 
           end if

           if (ijob==21) then  
              aux = conjg(Wne(e)) 
              if (fpm(33)==1) then
                  aux=conjg(aux)
                  workc(1:N,1:fpm(23))= conjg(workc(1:N,1:fpm(23)))                   
              end if
              call ZAXPY(N*fpm(23),aux,workc,1,Q(1,M0+1),1)

              ijob=-2 ! default identification

              if (fpm(31)==1) then
                 if (fpm(33)==0) then
                    call ZLACPY( 'F', N, fpm(23),work(1,M0+1) , N, workc, N )
                    workc(1:N,1:fpm(23))= conjg(workc(1:N,1:fpm(23)))
                    fpm(33)=1
                    ijob=21 ! reuse fact for second solve - must conjugate the solution qr
                    return
                 elseif (fpm(33)==1) then
                    fpm(33)=0
                 endif
              end if

           end if
        end do
     end IF !! info=0

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
        if(fpm(31)/=2) call MPI_ALLREDUCE(MPI_IN_PLACE,Q(1,M0+1),N*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,NEW_COMM_WORLD,code)
        call MPI_ALLREDUCE(MPI_IN_PLACE,Q(1,1),N*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,NEW_COMM_WORLD,code)
#endif
        !-----------------------------------------------  
     end if


  end IF   ! fpm(21)==1
!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if ((fpm(21)==4).and.(fpm(14)==1)) then !! only qr,ql vectors has been computed and is returned
     info=4
     if (info/=0) fpm(21)=100 ! The End
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!! Stochastic Estimate
  if ((fpm(21)==4).and.(fpm(14)==2)) then !! only qr is returned with stochastic estimation in res
     iseed=(/56,890,3456,2333/)
     call ZLARNV(4,iseed,N*fpm(23),work(1,1)) 
     call ZGEMM('C','N',fpm(23),fpm(23),N,-ONEC,work,N,q,N,ZEROC,zAq,M0) ! projection
     call ZGEMM('C','N',fpm(23),fpm(23),N,ONEC,work,N,work,N,ZEROC,zSq,M0) ! normalization
     theta=DZERO
     do i=1,fpm(23)
        theta=theta+dble(zAq(i,i)/zSq(i,i))
        res(i)=abs(theta*(N/(DONE*i)))
     enddo
     mode=int(real(res(fpm(23))))+1
     info=5
     if (info/=0) fpm(21)=100 ! The End
  end if

!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Form the reduced eigenvalue problem
!!!!!!! Aq xq=eq Sq xq     with Aq=Q^TAQ Sq=Q^TAQ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!form  Bq=> Bq=Q^T B Q
  if (fpm(21)==4) then 
     !------------------------------------------------
#ifdef MPI
     work(1:N,1:fpm(23))=ZEROC
#endif
     !-------------------------------------------------
     fpm(21)=5 ! preparing reenty
     ijob=40 
     return! mat-vec B.qr => workr
  end  if

  if (fpm(21)==5) then

     ! after first loop zSq(1,1) contains old residual-
     if (loop>0) Ze=zSq(1,1) !  dummy variable

     !------------------------------------------------
#ifdef MPI
     zSq(1:M0,1:fpm(23))=ZEROC
#endif
     !-------------------------------------------------
     if(fpm(31)==2)then
        call ZGEMM('T','N',fpm(23),fpm(25),N,ONEC,Q(1,1),N,work(1,fpm(24)),N,ZEROC,zSq(1,fpm(24)),M0)
     else
        call ZGEMM('C','N',fpm(23),fpm(35),N,ONEC,Q(1,M0+1),N,work(1,fpm(24)),N,ZEROC,zSq(1,fpm(24)),M0)
     endif
     !---------------------------------------
#ifdef MPI
     if (fpm(23)/nb_procs>=1) then
        call MPI_ALLREDUCE(MPI_IN_PLACE,work(1,1),N*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,NEW_COMM_WORLD,code)
        call MPI_ALLREDUCE(MPI_IN_PLACE,zSq(1,1),M0*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,NEW_COMM_WORLD,code)
     end if
#endif
     !       !---------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! Spurious test 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
     mode=0
     if ((fpm(13)==1) .and. (fpm(37)==0)) then
        if (loop>0) then
           do i=fpm(24),fpm(24)+fpm(25)-1
              if (res(i)==DONE) then ! indicator for being inside the search interval
                 if( abs( (sqrt(zSq(i,i)))-lambda(i) ) / sqrt(abs(zSq(i,i))) < 0.1d0) mode=mode+1
              end if
           end do
        endif
     endif
!!!!!!!!!!!!!!!!!!!!
#ifdef MPI
     call MPI_ALLREDUCE(MPI_IN_PLACE,mode,1,MPI_INTEGER,MPI_SUM,NEW_COMM_WORLD,code)
#endif
!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! B-biorthogonalization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


     !-------------------------
     ! Diagonalize Bq=Ql^H B Qr
     !-------------------------

     e=0
     call wallocate_2z(VR,fpm(23),fpm(23),infoloc)
     call wallocate_2z(VL,fpm(23),fpm(23),infoloc)
     call wallocate_1z(LALPHA,fpm(23),infoloc)

     if (rank==0) then
        call wallocate_2z(zSqo,fpm(23),fpm(23),infoloc)
        call ZLACPY( 'F', fpm(23), fpm(23),zSq , M0, zSqo, fpm(23) )

        LWORK_LOC=8*fpm(23) !! for lapack eig reduced system
        call wallocate_1z(zWORK_LOC,LWORK_LOC,infoloc)
        call wallocate_1d(WORK_LOCR,2*fpm(23),infoloc)
        if (infoloc/=0) info=-1

        JOBVL='V'
        if (fpm(31)==2) JOBVL='N'        
        call ZGEEV(JOBVL,'V',fpm(23),zSqo,fpm(23),LALPHA,VL,fpm(23),VR,fpm(23),zWORK_LOC,LWORK_LOC,WORK_LOCR,info_lap)

        theta = abs(lalpha(1)) !! Find max lalpha
        do j=1,fpm(23)
           if(theta<abs(lalpha(j))) theta=abs(lalpha(j))
        enddo

        e=0
        do j=1,fpm(23)  !! find dimension of NULL(Bq)=e
           if( abs(LALPHA(j))/theta < 1.0d-16 ) e = e+1
        enddo

        do i=1,e  !! Order eigenpairs by |lalpha_i|
           m_min=1
           theta = abs(lalpha(1))
           do j=1,fpm(23) - i + 1 !find min lalpha and put at end
              if( abs(LALPHA(j)) < theta )then
                 theta = abs(lalpha(j))
                 m_min = j
              endif
           enddo
           call ZSWAP(fpm(23),VR(1,m_min),1,VR(1,fpm(23)-i+1),1) 
    if (fpm(31)/=2) call ZSWAP(fpm(23),VL(1,m_min),1,VL(1,fpm(23)-i+1),1) 
           jac = LALPHA(fpm(23)-i+1)
           LALPHA(fpm(23)-i+1) = LALPHA(m_min)
           LALPHA(m_min)=jac
        enddo

        call wdeallocate_1z(zWORK_LOC)
        call wdeallocate_1d(WORK_LOCR)
        call wdeallocate_2z(zSqo)
     endif ! rank==0

     !-------------
#ifdef MPI
     call MPI_BCAST(e,1,MPI_INTEGER,0,NEW_COMM_WORLD,code)
     call MPI_BCAST(VR,fpm(23)*fpm(23),MPI_DOUBLE_COMPLEX,0,NEW_COMM_WORLD,code)
if (fpm(31)/=2)  call MPI_BCAST(VL,fpm(23)*fpm(23),MPI_DOUBLE_COMPLEX,0,NEW_COMM_WORLD,code)
     call MPI_BCAST(lalpha,fpm(23),MPI_DOUBLE_COMPLEX,0,NEW_COMM_WORLD,code)
#endif

     !--------------
     ! Create Ur and Ul (bi-B-orthogonal)
     !--------------------
     IF (fpm(15)==1) then  

        call ZLACPY('F',N,fpm(23),Q(1,1),N,workc(1,1),N)
        Q(1:N,1:fpm(23))=ZEROC
        call ZGEMM('N','N',N,fpm(25),fpm(23),ONEC,workc(1,1),N,VR(1,fpm(24)),fpm(23),ZEROC,Q(1,fpm(24)),N)
        call ZLACPY('F',N,fpm(23),work(1,1),N,workc(1,1),N)
        work(1:N,1:fpm(23))=ZEROC
        call ZGEMM('N','N',N,fpm(25),fpm(23),ONEC,workc(1,1),N,VR(1,fpm(24)),fpm(23),ZEROC,work(1,fpm(24)),N)
        if(fpm(31)/=2) then
           call ZLACPY('F',N,fpm(23),Q(1,M0+1),N,workc(1,1),N)
           Q(1:N,M0+1:M0+fpm(23))=ZEROC
           call ZGEMM('N','N',N,fpm(35),fpm(23),ONEC,workc(1,1),N,VL(1,fpm(24)),fpm(23),ZEROC,Q(1,fpm(34)),N)
        endif

        !--------------------
#ifdef MPI
        if (fpm(23)/nb_procs>=1) then
           call MPI_ALLREDUCE(MPI_IN_PLACE,work(1,1),N*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,NEW_COMM_WORLD,code)
           call MPI_ALLREDUCE(MPI_IN_PLACE,Q(1,1),N*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,NEW_COMM_WORLD,code)
           if (fpm(31)/=2) call MPI_ALLREDUCE(MPI_IN_PLACE,Q(1,M0+1),N*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,NEW_COMM_WORLD,code)
        endif
#endif

        !----------------------------------------
        !---------- Normalization of Ur and Ul
        !----------------------------------------          
        call ZGEMM('N','N',fpm(23),fpm(23),fpm(23),ONEC,zSq(1,1),M0,VR(1,1),fpm(23),ZEROC,zAq(1,1),M0)
if (fpm(31)/=2) then
        call ZGEMM('C','N',fpm(23),fpm(23),fpm(23),ONEC,VL(1,1),fpm(23),zAq(1,1),M0,ZEROC,zSq(1,1),M0)
else
        call ZGEMM('T','N',fpm(23),fpm(23),fpm(23),ONEC,VR(1,1),fpm(23),zAq(1,1),M0,ZEROC,zSq(1,1),M0)
endif

        do i=1,fpm(23)
           call ZSCAL(N,ONEC/sqrt(zSq(i,i)),Q(1,i),1)
           call ZSCAL(N,ONEC/sqrt(zSq(i,i)),work(1,i),1)
           if(fpm(31)/=2) call ZSCAL(N,ONEC/sqrt(conjg(zSq(i,i))),Q(1,M0+i),1)
        enddo

     end IF ! fpm(15)==1

     call wdeallocate_1z(LALPHA)
     call wdeallocate_2z(VR)
     call wdeallocate_2z(VL)

     !        if ((e/=0) .and. (loop>0))then
     if (e/=0) then
        fpm(23) = fpm(23)-e
        fpm(25)=fpm(23) ! default - single proc
        fpm(35)=fpm(23) ! default - single proc
        if (fpm(1)==1) then 
           call wwrite_s(fout,'Resize Subspace (Bq)')
           call wwrite_t(fout)
           call wwrite_i(fout,fpm(23))
           call wwrite_n(fout)
        endif


        if (fpm(23)/nb_procs>=1) then ! criteria for parallelism of mat-vec
           fpm(25)=fpm(23)/nb_procs ! local size of 'current M0'
           if (rank==nb_procs-1) fpm(25)=fpm(25)+fpm(23)-(fpm(23)/nb_procs)*nb_procs 
           fpm(24)=1+rank*(fpm(23)/nb_procs) ! local origin of first column for vector q for parallel mat-vec 

           fpm(35)=fpm(23)/nb_procs ! local size of 'current M0'
           if (rank==nb_procs-1) fpm(35)=fpm(35)+fpm(23)-(fpm(23)/nb_procs)*nb_procs 
           fpm(34)=M0+1+rank*(fpm(23)/nb_procs) ! local origin of first column for vector q for parallel mat-vec 
        end if

     endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     !--------------------------------
#ifdef MPI
     zSq(1:M0,1:fpm(23))=ZEROC
#endif
     !--------------------------------
     if(fpm(31)==2)then
        call ZGEMM('T','N',fpm(23),fpm(25),N,ONEC,Q(1,1),N,work(1,fpm(24)),N,ZEROC,zSq(1,fpm(24)),M0)
     else
        call ZGEMM('C','N',fpm(23),fpm(35),N,ONEC,Q(1,M0+1),N,work(1,fpm(24)),N,ZEROC,zSq(1,fpm(24)),M0)
     endif


     fpm(21)=6
  endif


!!!!!!!!!form Aq=> Aq=Q^T A Q 
  if (fpm(21)==6) then
     !------------------------------------------------
#ifdef MPI 
     work(1:N,1:fpm(23))=ZEROC
#endif 
     fpm(21)=7 ! preparing reentry
     ijob=30 
     return  ! mat-vec A.qr => workr
  endif


  if (fpm(21)==7) then 
     !------------------------------------------------
#ifdef MPI 
     zAq(1:M0,1:fpm(23))=ZEROC
#endif
     !-------------------------------------------------
     if(fpm(31)==2)then
        call ZGEMM('T','N',fpm(23),fpm(25),N,ONEC,Q(1,1),N,work(1,fpm(24)),N,ZEROC,zAq(1,fpm(24)),M0)
     else
        call ZGEMM('C','N',fpm(23),fpm(35),N,ONEC,Q(1,M0+1),N,work(1,fpm(24)),N,ZEROC,zAq(1,fpm(24)),M0)
     endif

     fpm(21)=8
  endif



  if (fpm(21)==8) then
     !----------------------------------------!(zAq,zSq known to all processors) 
#ifdef MPI 
     if (fpm(23)/nb_procs>=1) then
        call MPI_ALLREDUCE(MPI_IN_PLACE,zAq(1,1),M0*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,NEW_COMM_WORLD,code)
        call MPI_ALLREDUCE(MPI_IN_PLACE,zSq(1,1),M0*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,NEW_COMM_WORLD,code)
     end if
#endif
     !---------------------------------------
     if (fpm(12)==1) then ! customize eigenvalue solver
        fpm(21)=9 ! preparing reentry - could return new value of M0 in fpm(23) if reduced subspace is needed
        ijob=50
        return
     endif
  endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Solve the reduced eigenvalue problem
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (fpm(21)==8) then
     ! if using FEAST-MPI ==> solve on a single proc
     if (rank==0) then
        info_lap=1 ! initialization
        LWORK_LOC=8*fpm(23) !! for lapack eig reduced system
        call wallocate_1z(zWORK_LOC,LWORK_LOC,infoloc)
        call wallocate_1d(WORK_LOCR,8*fpm(23),infoloc)
        if (infoloc/=0) info=-1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        call wallocate_2z(zSqo,fpm(23),fpm(23),infoloc)
        call wallocate_2z(VL,fpm(23),fpm(23),infoloc)
        call wallocate_2z(VR,fpm(23),fpm(23),infoloc)

        if (infoloc/=0) info=-1

        call ZLACPY( 'F', fpm(23), fpm(23),zSq , M0, zSqo, fpm(23) )

        if (fpm(39)==1) then  !!! standard eigenvalue

           do i=1,fpm(23)
              zAq(1:fpm(23),i) = zAq(1:fpm(23),i) / sqrt(zSq(i,i))
              zAq(i,1:fpm(23)) = zAq(i,1:fpm(23)) / sqrt(zSq(i,i))
           enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
           JOBVL='V'
           if (fpm(31)==2) JOBVL='N'  
           call ZGEEV(JOBVL,'V',fpm(23),zAq,M0,LAMBDA,VL,fpm(23),VR,fpm(23), zWORK_LOC, LWORK_LOC, WORK_LOCR,INFO_lap)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

           if(info_lap==0)then
              do j=1,fpm(23)
                 call ZSCAL(fpm(23),ONEC/sqrt(zSq(j,j)),VR(j,1),fpm(23))
       if (fpm(32)/=2)  call ZSCAL(fpm(23),ONEC/conjg(sqrt(zSq(j,j))),VL(j,1),fpm(23))
              end do
           endif

        else                !!! generalized eigenvalue

           call wallocate_1z(LALPHA,fpm(23),infoloc)
           call wallocate_1z(LBETA,fpm(23),infoloc)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           JOBVL='V'
           if (fpm(31)==2) JOBVL='N'  
           call ZGGEV(JOBVL,'V',fpm(23),zAq,M0,zSqo,fpm(23),LALPHA,LBETA,VL,fpm(23),VR,fpm(23), zWORK_LOC, LWORK_LOC, WORK_LOCR,INFO_lap)
           LAMBDA(1:fpm(23))=LALPHA(1:fpm(23))/LBETA(1:fpm(23))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           call wdeallocate_1z(LALPHA)
           call wdeallocate_1z(LBETA)
        endif

        if (info_lap/=0) info=-3 ! problem with lapack

        IF (info_lap==0) then

           call ZGEMM('N','N',fpm(23),fpm(23),fpm(23),ONEC,zSq(1,1),M0,VR(1,1),fpm(23),ZEROC,zAq(1,1),M0)
           if(fpm(31)==2) then
              call ZGEMM('T','N',fpm(23),fpm(23),fpm(23),ONEC,VR(1,1),fpm(23),zAq(1,1),M0,ZEROC,zSqo(1,1),fpm(23))
           else
              call ZGEMM('C','N',fpm(23),fpm(23),fpm(23),ONEC,VL(1,1),fpm(23),zAq(1,1),M0,ZEROC,zSqo(1,1),fpm(23))
           endif


           fpm(37)=0
           do j=1,fpm(23) !!!!!! Check othogonality of ZGGEV eigenvectors
              zalpha = zSqo(j,j)
              do jj=1,fpm(23)
                 zalpha = zalpha - zSqo(jj,j)
              enddo
              if (log10(abs(zalpha)/(fpm(23)))>-10) then
                 fpm(37) = 1 
              endif
           enddo

           if (fpm(37)==0) then  !!!!!! normalize
              do j=1,fpm(23)
                 call ZSCAL(fpm(23),ONEC/sqrt(zSqo(j,j)),VR(1,j),1)
                 if(fpm(31)/=2) call ZSCAL(fpm(23),ONEC/conjg(sqrt(zSqo(j,j))),VL(1,j),1)
              enddo
           endif

           call ZLACPY( 'F', fpm(23), fpm(23),VR , fpm(23), zAq, M0 )
    if (fpm(32)/=2) call ZLACPY( 'F', fpm(23), fpm(23),VL , fpm(23), zSq, M0 )

        end IF


        call wdeallocate_2z(zSqo)
        call wdeallocate_2z(VL)
        call wdeallocate_2z(VR)

        call wdeallocate_1z(zwork_loc)
        call wdeallocate_1d(work_locr)
     end if !(rank 0)
     !-------------------------------- !(info common to all processors) 
#ifdef MPI
     call MPI_BCAST(info,1,MPI_INTEGER,0,NEW_COMM_WORLD,code)
#endif
     !--------------------------------
     if (info/=0) fpm(21)=100 ! the end


     if (info==0) then

        !---------------------------------------- !(zAq==> vectors, lambda and fpm(23), known to all processors) 
#ifdef MPI
        call MPI_BCAST(zAq,M0*( fpm(23) ),MPI_DOUBLE_COMPLEX,0,NEW_COMM_WORLD,code)
        call MPI_BCAST(zSq,M0*( fpm(23) ),MPI_DOUBLE_COMPLEX,0,NEW_COMM_WORLD,code)
        call MPI_BCAST(lambda,fpm(23),MPI_DOUBLE_COMPLEX,0,NEW_COMM_WORLD,code)
        call MPI_BCAST(fpm(37),1,MPI_INTEGER,0,NEW_COMM_WORLD,code)
#endif
        !-----------------------------

        fpm(21)=9
     end if
  end if !! fpm(21)=8



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Compute Ritz vectors Xr=Qr*xr  Xl=Ql*xl 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(21)==9 ) then
     call ZCOPY( N*fpm(23), Q(1,1), 1, workc, 1 ) 
     !! option - non shifted
#ifdef MPI 
     Q(1:N,1:fpm(23))=ZEROC
#endif
     call ZGEMM('N','N',N,fpm(25),fpm(23),ONEC,workc(1,1),N,zAq(1,fpm(24)),M0,ZEROC,Q(1,fpm(24)),N)

     if(fpm(31)/=2)  then
        call ZCOPY( N*fpm(23), Q(1,M0+1), 1, workc, 1 ) 
#ifdef MPI 
        Q(1:N,M0+1:M0+fpm(23))=ZEROC
#endif
        call ZGEMM('N','N',N,fpm(35),fpm(23),ONEC,workc(1,1),N,zSq(1,fpm(24)),M0,ZEROC,Q(1,fpm(34)),N)
     endif
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Compute Residuals ||AXr-lambda*BXr|| ||A^H Xl-conjg(lambda)*B^H Xl|| 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (fpm(21)==9) then
     fpm(21)=10 ! preparing reentry
     ijob=30 
     return  ! mat-vec A*qr => workr
  endif

  if (fpm(21)==10) then
     call ZLACPY( 'F', N, fpm(25),work(1,fpm(24)), N, workc(1,fpm(24)), N ) ! workc = Ax
     fpm(21)=11 ! preparing reentry
     ijob=40 
     !----------------------------------------
#ifdef MPI
     work(1:N,1:fpm(23))=ZEROC  !! work is also needed for outer-loop if any 
#endif
     !------------------------------------------
     return  ! mat-vec S*qr => workr
  endif

  if (fpm(21)==11) then
     !----------------------------------------
#ifdef MPI
     res(1:fpm(23))=DZERO
#endif
     !-------- Absolute residual
     do i=fpm(24),fpm(24)+fpm(25)-1
        res(i)=sum(abs(workc(1:N,i)-lambda(i)*work(1:N,i)))/sum(abs((abs(Emid)+r)*work(1:N,i)))
     end  do
     !----------------------------------------
#ifdef MPI 
     if (fpm(23)/nb_procs>=1) call MPI_ALLREDUCE(MPI_IN_PLACE,res,fpm(23),MPI_DOUBLE_PRECISION,MPI_SUM,NEW_COMM_WORLD,code)
#endif
     !-----------------------------------------
     fpm(21)=-9 
     if(fpm(31)==2) fpm(21)=-12
  end if


  if (fpm(21)==-9) then
     fpm(21)=-10 ! preparing reentry
     ijob=31 
     return  ! mat-vec A^C.ql => work(1,fpm(34))
  endif

  if (fpm(21)==-10) then
     call ZLACPY( 'F', N, fpm(35),work(1,fpm(34)), N, workc(1,fpm(24)), N ) ! workc = A^C.XL
     fpm(21)=-11 ! preparing reentry
     ijob=41 
     !----------------------------------------
#ifdef MPI
     work(1:N,M0+1:M0+fpm(23))=ZEROC  !! work is also needed for outer-loop if any 
#endif
     !------------------------------------------
     return  ! mat-vec S^C.ql => work
  endif

  if (fpm(21)==-11) then
     !----------------------------------------
#ifdef MPI
     res(M0+1:M0+fpm(23))=DZERO
#endif
     !-------- Absolute residual
     do i=fpm(34),fpm(34)+fpm(35)-1
        res(i)=sum(abs(workc(1:N,i-M0)-conjg(lambda(i-M0))*work(1:N,i)))/sum(abs((abs(Emid)+r)*work(1:N,i)))
     end  do
     !----------------------------------------
#ifdef MPI 
     if (fpm(23)/nb_procs>=1) call MPI_ALLREDUCE(MPI_IN_PLACE,res(M0+1),fpm(23),MPI_DOUBLE_PRECISION,MPI_SUM,NEW_COMM_WORLD,code)
#endif
     !-----------------------------------------
     fpm(21)=-12
  end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Compute Trace / count eigenvalues (remove spurious if any)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  IF (fpm(21)==-12) then

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !---color mapping --> check if eigenvalue inside the contour
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     call wallocate_1i(color,fpm(23),infoloc)

     if (fpm(29)==0) then ! expert custom contour
        call zfeast_inside_contourx(Zne,fpm(8),lambda,fpm(23),color)
     else
        call zfeast_inside_contour(Emid,r,fpm(18),fpm(19),lambda,fpm(23),color)
     endif
!!!!!!!!!!!!

     k=0
     trace=ZEROC
     theta=DZERO
!!! count how many eigenvalues have converged + trace + max residual
     do i=1,fpm(23)
        if (color(i)==1) then ! inside the search interval
           k=k+1
           trace=trace+lambda(i)
           if (res(i)>theta) theta=res(i) ! max residual (right residual)
           if (fpm(31)/=2) then
              if (res(M0+i)>theta) theta=res(M0+i) ! max residual (left residual)
           end if
        endif
     enddo
!!!!!!!!! remove spurious if any  ---- (work only with the right subspace) ----
     if ((mode==0).and.(k>0)) mode=k ! wait before looking into mode 
     ! Rq: if all eigenvalues k are spurious...FEAST will not converge
     if (mode<k) then
        do j=1,k-mode
           theta=DZERO
           e=1
           do i=1,fpm(23)
              if (color(i)==1) then ! inside the search interval
                 if (res(i)>theta) then
                    e=i
                    theta=res(i) !max
                 endif
              end if
           enddo
           trace=trace-lambda(e) ! update trace
           res(e)=-DONE !spurious
           if (fpm(31)/=2) res(e+M0)=-DONE !spurious
        end do
!!! max residual
        theta=DZERO
        do i=1,fpm(23)
           if (color(i)==1) then ! inside the search interval
              if (res(i)>theta) theta=res(i) ! max residual (right residual)
              if (fpm(31)/=2) then
                 if (res(M0+i)>theta) theta=res(M0+i) ! max residual (left residual)
              end if
           end if
        enddo
     end if

!     if (mode==0) info=1  ! no eigenvalue detected in the interval   !!!<<<<< here
     if (loop>1) then! wait second iteration (spurious related)
        if ((mode==M0).and.(mode/=N)) info=3 !size subspace too small
     endif
     if (info/=0) then
        fpm(21)=100 ! The End
        call wdeallocate_1i(color)
     end if

  end IF


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Check! FEAST Convergence
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       

  if (fpm(21)==-12) then
     testconv=.false. ! initialization

     ! here Ze contains the old trace         
     if (loop==0) epsout=DONE
     if (loop>0) then !! trace convergence
        epsout=abs(trace-Ze)/(abs(Emid) + r) 
        if ((fpm(6)==0).and.(log10((epsout))<(-fpm(3)))) testconv=.true.
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
        ZSq(1,1)=trace ! dummy
        if (loop==fpm(4)) then
           info=2 ! FEAST did not converge (#loop reaches maximum)
           testconv=.true. ! return final eigenvector anyway
        endif
     endif
  end if

  
    if (mode==0) testconv=.false. !assume we'll find something, don't stop until we do !!!!<<<<<<<<<<<<<<



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! FEAST exit IF Convergence - FEAST iteration IF NOT 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!            

  if (fpm(21)==-12) then

     if (.not.testconv) then !!! need FEAST refinement loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! rational function for the inner eigenvalues
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        res(1:fpm(23))=DZERO ! temp indicator for being inside the contour
        call wallocate_1z(zwork_loc,fpm(23),infoloc)
        call zfeast_grationalx(Zne,Wne,fpm(8),lambda,fpm(23),zwork_loc)
        do i=1,fpm(23)
           if (color(i)==1)  res(i)=DONE
           lambda(i)=zwork_loc(i)
           !if (fpm(31) == 1) lambda(i) = (2.0d0,0.0d0)*dble(zwork_loc(i))
        enddo
        call wdeallocate_1z(zwork_loc)
        call wdeallocate_1i(color)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! Remark: q and work=S*q is known already,  using FEAST-MPI work and q are already distributed
        fpm(21)=1   ! prepare reentry- reloop (with contour integration)
        !fpm(21)=4 ! reloop (without contour integration) -in this case work=q (actually does not need "work")
        loop=loop+1
        ijob=-2 ! do nothing
        return  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     else    ! final eigenvectors/eigenvalues
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
        
#ifdef MPI       
        if (fpm(23)/nb_procs>=1) call MPI_ALLREDUCE(MPI_IN_PLACE,q(1,1),N*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,NEW_COMM_WORLD,code)
        if (fpm(31)/=2) then
         if (fpm(23)/nb_procs>=1) call MPI_ALLREDUCE(MPI_IN_PLACE,q(1,M0+1),N*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,NEW_COMM_WORLD,code)
        endif
#endif

!! reorder (shift) lambda, eigenvector and residual
          call ZCOPY( N*fpm(23),Q(1,1), 1, work(1,1), 1 )
          if(fpm(31)/=2) call ZCOPY( N*fpm(23),Q(1,M0+1), 1, work(1,fpm(23)+1), 1 )
          call wallocate_1z(zwork_loc,fpm(23),infoloc)
          if (fpm(31)/=2) then
          call wallocate_1d(work_loc2,2*fpm(23),infoloc) 
          else
          call wallocate_1d(work_loc2,fpm(23),infoloc) 
          end if
          call ZCOPY(fpm(23),lambda , 1, zwork_loc, 1 )
          call DCOPY(fpm(23),res , 1, work_loc2, 1 )
          if(fpm(31)/=2) call DCOPY(fpm(23), res(M0+1), 1, work_loc2(fpm(23)+1), 1 )

          M0=fpm(23)  ! update value of M0 (new subspace)

          Q(1:N,1:fpm(23))=ZEROC
          if(fpm(31)/=2) Q(1:N,M0+1:M0+fpm(23))=ZEROC
          e=0
          j=0
          k=0
          do i=1,fpm(23)  
          if (color(i)==1) then ! inside the search interval
             if (work_loc2(i)/=-DONE) then ! not spurious 
                    e=e+1
                    q(1:N,e)=work(1:N,i)
                    lambda(e)=zwork_loc(i)
                    res(e)=work_loc2(i)
                    if(fpm(31)/=2) then
                    q(1:N,M0+e)=work(1:N,M0+i)
                    res(M0+e)=work_loc2(M0+i)
                 endif
              endif
           else
              j=j+1
              q(1:N,mode+j)=work(1:N,i)
              lambda(mode+j)=zwork_loc(i)
              res(mode+j)=work_loc2(i)
              if(fpm(31)/=2) then
                 q(1:N,M0+mode+j)=work(1:N,M0+i)
                 res(M0+mode+j)=work_loc2(M0+i)
              endif
           end if
           if (work_loc2(i)==-DONE) then ! spurious at the end
              k=k+1
              q(1:N,fpm(23)-k+1)=work(1:N,i)
              lambda(fpm(23)-k+1)=zwork_loc(i)
              res(fpm(23)-k+1)=work_loc2(i)
              if(fpm(31)/=2) then
                 q(1:N,M0+fpm(23)-k+1)=work(1:N,M0+i)
                 res(M0+fpm(23)-k+1)=work_loc2(M0+i)
              end if
           endif
        end do

        call wdeallocate_1z(zwork_loc)
        call wdeallocate_1d(work_loc2)
        call wdeallocate_1i(color)

!!!!!!!!!!!!!!!!!!!!!!!!!!
        fpm(21)=100 ! The End

        if ((info==0).and.(fpm(37)==1)) info=6

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

          if (info==6) then
             call wwrite_s(fout, '==>WARNING: FEAST converge but subspace is not biorthonormal')
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

  end subroutine zfeast_mygrcix







  
