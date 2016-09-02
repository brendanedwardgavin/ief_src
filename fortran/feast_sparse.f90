  !=========================================================================================
  !Copyright (c) 2009-2013, The Regents of the University of Massachusetts, Amherst.
  !Developed by E. Polizzi
  !All rights reserved.
  !
  !Redistribution and use in source and binary forms, with or without modification, 
  !are permitted provided that the following conditions are met:
  !
  !1. Redistributions of source code must retain the above copyright notice, this list of conditions 
  !   and the following disclaimer.
  !2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions 
  !   and the following disclaimer in the documentation and/or other materials provided with the distribution.
  !3. Neither the name of the University nor the names of its contributors may be used to endorse or promote
  !    products derived from this software without specific prior written permission.
  !
  !THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, 
  !BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
  !ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
  !EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
  !SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
  !LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING 
  !IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
  !==========================================================================================
  
  
  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!! FEAST PREDEFINED SPARSE INTERFACES !!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! List of routines:
!-------------------

!{S,D,C,Z}FEAST_{SCSR,HCSR}{EV,GV}

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Symmetric eigenvalue
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! single/double precision
! real symmetric 
!{s,d}feast_scsrgv ! generalized
!{s,d}feast_scsrev ! standard

! complex Hermitian
!{c,z}feast_hcsrgv ! generalized
!{c,z}feast_hcsrev ! standard
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  include 'lsprim.f90' !! Sparse primitives


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  
  subroutine dfeast_scsrgv(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
    !  Purpose 
    !  =======
    !  FEAST SPARSE INTERFACE (SPARSE CSR FORMAT)
    !  Solve the generalized Ax=eBx eigenvalue problem
    !  
    !  A REAL SYMMETRIC, B SYMMETRIC POSITIVE DEFINITE 
    ! 
    !  DOUBLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  UPLO       (input)        CHARACTER: specifies whether the full part 'F', or the upper 'U' or lower 'L'
    !                                       triangular part of the matrix(ces) is(are) being supplied.
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
    !  fpm (input/output)     INTEGER(*) : FEAST parameters (see FEAST documentation)
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
    !                                                        if option fpm(6)=1 selected                           
    !  info       (output)       INTEGER: Error handling (0: successful exit) -- see FEAST documentation
    !=====================================================================
    ! Eric Polizzi 2009-2013
    ! ====================================================================

    implicit none
    include 'f90_noruntime_interface.fi'
    character(len=1) :: UPLO
    integer :: N
    double precision,dimension(*),target:: sa,sb
    integer,dimension(*),target:: isa,jsa,isb,jsb
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer :: ijob,infoloc
    complex(kind=(kind(1.0d0))) :: Ze
    complex(kind=(kind(1.0d0))),dimension(:),pointer :: saz
    integer,dimension(:),pointer :: isaz,jsaz
    complex(kind=(kind(1.0d0))),dimension(:,:),pointer ::workc,caux,vv,r
    double precision, dimension(:,:),pointer ::work,Aq,Sq
    double precision,parameter :: DONE=1.0d0, DZERO=0.0d0
    complex(kind=(kind(1.0d0))),parameter :: ONEC=(DONE,DZERO)
!!!!!for pardiso
    integer(8),dimension(64) :: pt
    integer,dimension(64) :: iparm
    integer :: mtype
    integer :: MAXFCT,MNUM,PHASE,MSGLVL
    integer :: idum,nrhs
    double precision :: ddum1,err
!!!!! csr-upper format
    double precision,dimension(:),pointer :: ssa,ssb,nres
    integer,dimension(:), pointer :: sisa,sjsa,sisb,sjsb
    integer :: i,k 
    integer :: opt,nnza,nnzb,nnz
    integer :: itmax,im

!!!!!!!!!!!!!! Check INPUT PARAMETERS
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

    info=-1 ! by default

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! FORMAT CONVERSION TO CSR-UPPER for PARDISO !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if ((UPLO=='F').or.(UPLO=='f')) then !!! full csr to upper-csr
       call wallocate_1i(sisa,n+1,infoloc)
if (infoloc/=0) return
       call wallocate_1i(sisb,n+1,infoloc)
if (infoloc/=0) return
 
!!<<<
       call wallocate_1d(ssa,1,infoloc) ! dummy
if (infoloc/=0) return
       call wallocate_1i(sjsa,1,infoloc) !dummy
if (infoloc/=0) return
       call wallocate_1d(ssb,1,infoloc) !dummy
if (infoloc/=0) return
       call wallocate_1i(sjsb,1,infoloc) !dummy
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
if (infoloc/=0) return
       call wallocate_1i(sjsa,nnza,infoloc)
if (infoloc/=0) return
       call wallocate_1d(ssb,nnzb,infoloc)
if (infoloc/=0) return
       call wallocate_1i(sjsb,nnzb,infoloc)
if (infoloc/=0) return

       opt=2
       call dcsr2csr_up(opt,N,sa,isa,jsa,ssa,sisa,sjsa)
       call dcsr2csr_up(opt,N,sb,isb,jsb,ssb,sisb,sjsb)

!!!!!!!!!!!!!!!!!!!!!!!!!
    elseif ((UPLO=='U').or.(UPLO=='u')) then !! upper-csr already 

       nnza=isa(n+1)-1
       ssa => sa(1:nnza)
       sisa => isa(1:n+1)
       sjsa => jsa(1:nnza)

       nnzb=isb(n+1)-1
       ssb =>  sb(1:nnzb)
       sisb => isb(1:n+1)
       sjsb => jsb(1:nnzb)



    elseif ((UPLO=='L').or.(UPLO=='l')) then !!! lower-csr to upper-csr

       nnza=isa(n+1)-1
       nnzb=isb(n+1)-1
       call wallocate_1d(ssa,nnza,infoloc)
if (infoloc/=0) return
       call wallocate_1i(sjsa,nnza,infoloc)
if (infoloc/=0) return
       call wallocate_1d(ssb,nnzb,infoloc)
if (infoloc/=0) return
       call wallocate_1i(sjsb,nnzb,infoloc)
if (infoloc/=0) return
       call wallocate_1i(sisa,n+1,infoloc)
if (infoloc/=0) return
       call wallocate_1i(sisb,n+1,infoloc)
if (infoloc/=0) return

       call dcsr_transpose(N,sa,isa,jsa,ssa,sisa,sjsa)
       call dcsr_transpose(N,sb,isb,jsb,ssb,sisb,sjsb)

    endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! FEAST INITIALIZATION 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call wallocate_2d(Aq,M0,M0,infoloc)
    if (infoloc/=0) return
    call wallocate_2d(Sq,M0,M0,infoloc)
    if (infoloc/=0) return
    call wallocate_2d(work,N,M0,infoloc)
    if (infoloc/=0) return
    call wallocate_2z(workc,N,M0,infoloc)
    if (infoloc/=0) return
    call wallocate_1i(isaz,n+1,infoloc)
    if (infoloc/=0) return
    !!<<<
    call wallocate_1z(saz,1,infoloc) ! dummy
    if (infoloc/=0) return
    call wallocate_1i(jsaz,1,infoloc)! dummy
    if (infoloc/=0) return
    !!>>>
    opt=1
    call zdaddcsr(N,opt,ONEC,ssa,sisa,sjsa,ONEC,ssb,sisb,sjsb,saz,isaz,jsaz) !! get isaz
    nnz=isaz(n+1)-1
    !!<<<
    call wdeallocate_1z(saz)
    call wdeallocate_1i(jsaz)
    !!>>>
    call wallocate_1z(saz,nnz,infoloc)
    if (infoloc/=0) return
    call wallocate_1i(jsaz,nnz,infoloc)
    if (infoloc/=0) return

    opt=2
    call zdaddcsr(N,opt,ONEC,ssa,sisa,sjsa,ONEC,ssb,sisb,sjsb,saz,isaz,jsaz) !! get jsaz
!!!!!!!!!!!!!!!
    call wallocate_2z(caux,N,M0,infoloc)
    if (infoloc/=0) return
    call wallocate_2z(r,N,M0,infoloc)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  analysis step (symbolic factorizarion for pardiso)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    MAXFCT=1 ! one factorization to consider normal+transpose
    MTYPE=6  ! complex and symmetric 
    call pardisoinit(PT,MTYPE,IPARM)

!!!!!!!!!!!!
if (fpm(64)==1) then
do i=1,64
if (fpm(64+i)/=-9999) iparm(i)=fpm(64+i)
enddo
endif
!!!!!!!!!!!!

    !IPARM(3)=fpm(10) !omp_num_threads !! openmp -number of threads
    IPARM(6)=1 ! solution and rhs are input/output, attention caux is always used
    MSGLVL=0 !0- no output, 1- output
    PHASE=11
    MNUM=1
    call PARDISO(PT,MAXFCT,MNUM,MTYPE,PHASE,n,saz,isaz,jsaz,idum,M0,IPARM,MSGLVL,workc(:,1:M0),caux(:,1:M0),infoloc) 
    if (infoloc/=0) then
       info=-2
       return
    end if


im=32!50!6
allocate(vv(N,im+1)) 
itmax=5
err=1d-14 

       call wallocate_1d(nres,M0,infoloc) ! dummy


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!! FEAST !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

    ijob=-1 ! initialization

    do while (ijob/=0)

       call dfeast_srci(ijob,N,Ze,work,workc,Aq,Sq,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
       select case(ijob)

       case(10) !! Factorize (ZeB-A)
          opt=3
          call zdaddcsr(N,opt,-ONEC,ssa,sisa,sjsa,Ze,ssb,sisb,sjsb,saz,isaz,jsaz) !! get saz

if (fpm(43)==0) then
          PHASE=22
          call PARDISO(PT,MAXFCT,MNUM,MTYPE,PHASE,n,saz,isaz,jsaz,idum,M0,IPARM,MSGLVL,workc(:,1:M0),caux(:,1:M0),infoloc)
          if (infoloc/=0) info=-2
endif

       case(11) !! Solve (ZeB-A)x=work2(1:N,1:M0) result in to work2


select case(fpm(43))
case(0)
          PHASE=33
          call PARDISO(PT,MAXFCT,MNUM,MTYPE,PHASE,n,saz,isaz,jsaz,idum,M0,IPARM,MSGLVL,workc(:,1:M0),caux(:,1:M0),infoloc)
          if (infoloc/=0) info=-2


 case(10) !!!!! gmrez
 caux(1:N,1:M0)=(0.0d0,0.0d0) !!! zero initial guess
 do i=1,M0
!  call Bgmrez(n,im,workc(1,i),caux(1,i),vv,err,itmax,-6,saz,isaz,jsaz,ssb*(1.0d0,0.0d0),'U')
call gmrez(n,im,workc(1,i),caux(1,i),vv,err,itmax,-6,saz,isaz,jsaz,'U')

  enddo

 ! new residual (b-A*x(k))
         r(1:N,1:M0)=workc(1:N,1:M0)
         call zcsrmm('U',N,N,M0,(-1.0d0,0.0d0),saz,isaz,jsaz,caux,(1.0d0,0.0d0),r)
do i=1,M0
      nres(i)=sum(abs(r(:,i)))/sum(abs(workc(:,i))) ! norm relative residual 
 enddo
     
 print *,maxval(nres),minval(nres) 
workc(1:N,1:M0)=caux(1:N,1:M0)
     




end select






       case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

          call dcsrmm('U',N,N,fpm(25),DONE,ssa,sisa,sjsa,X(1,fpm(24)),DZERO,work(1,fpm(24)))

       case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

          call dcsrmm('U',N,N,fpm(25),DONE,ssb,sisb,sjsb,X(1,fpm(24)),DZERO,work(1,fpm(24)))

       end select
    end do
!!!!!!!!!!!!!!!!!!!!!!!!
!!!! release memory
!!!!!!!!!!!!!!!!!!!!!!!!
    PHASE=-1 

    call PARDISO(PT,MAXFCT,MNUM,MTYPE,PHASE,n,saz,isaz,jsaz,idum,M0,IPARM,MSGLVL,workc,caux,infoloc)
    if (infoloc/=0) then
       info=-2
       return
    end if

 

    call wdeallocate_2d(Aq)
    call wdeallocate_2d(Sq)
    call wdeallocate_2d(work)
    call wdeallocate_2z(workc)
    call wdeallocate_2z(caux)

    call wdeallocate_1z(saz)
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



  end subroutine dfeast_scsrgv


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  



subroutine dfeast_scsrev(UPLO,N,sa,isa,jsa,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
    !  Purpose 
    !  =======
    !  FEAST SPARSE INTERFACE (SPARSE CSR FORMAT)
    !  Solve the standard Ax=ex eigenvalue problem
    !  
    !  A REAL SYMMETRIC
    ! 
    !  DOUBLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  UPLO       (input)        CHARACTER: specifies whether the full part 'F', or the upper 'U' or lower 'L'
    !                                       triangular part of the matrix(ces) is(are) being supplied.
    !  N          (input)        INTEGER: Size system 
    !  sa         (input)        REAL DOUBLE PRECISION (isa(N+1)-1):  Matrix A- CSR format 
    !  isa        (input)        INTEGER(N+1): CSR row array of Matrix A   
    !  jsa        (input)        INTEGER(isa(N+1)-1): CSR column array of Matrix A
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  fpm (input/output)     INTEGER(*) : FEAST parameters (see FEAST documentation)
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
    !                                                        if option fpm(6)=1 selected                           
    !  info       (output)       INTEGER: Error handling (0: successful exit) -- see FEAST documentation
    !=====================================================================
    ! Eric Polizzi 2009-2013
    ! ====================================================================

    implicit none
    include 'f90_noruntime_interface.fi'
    character(len=1) :: UPLO
    integer :: N
    double precision,dimension(*),target:: sa
    integer,dimension(*),target:: isa,jsa
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer :: ijob,infoloc
    complex(kind=(kind(1.0d0))) :: Ze
    complex(kind=(kind(1.0d0))),dimension(:),pointer :: saz
    integer,dimension(:),pointer :: isaz,jsaz
    complex(kind=(kind(1.0d0))),dimension(:,:),pointer ::workc,caux
    double precision, dimension(:,:),pointer ::work,Aq,Sq
    double precision,parameter :: DONE=1.0d0, DZERO=0.0d0
    complex(kind=(kind(1.0d0))),parameter :: ONEC=(DONE,DZERO)
!!!!!for pardiso
    integer(8),dimension(64) :: pt
    integer,dimension(64) :: iparm
    integer :: mtype
    integer :: MAXFCT,MNUM,PHASE,MSGLVL
    integer :: idum,nrhs
    double precision :: ddum1
!!!!! csr-upper format
    double precision,dimension(:),pointer :: ssa,ssb
    integer,dimension(:), pointer :: sisa,sjsa,sisb,sjsb
    integer :: i,k 
    integer :: opt,nnza,nnzb,nnz


!!!!!!!!!!!!!! Check INPUT PARAMETERS
    INFO = 0
    IF ((UPLO/='F').and.(UPLO/='f').and.(UPLO/='L').and.(UPLO/='l').and.(UPLO/='U').and.(UPLO/='u')) THEN
       INFO=-101
    ELSE IF ( N<=0 ) THEN
       INFO = -102
    END IF
    IF( INFO.NE.0 ) THEN
       CALL XERBLA( 'DFEAST_SCSREV', -INFO+100 )
       RETURN
    END IF

    info=-1 ! by default

infoloc=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! identity B matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!
    nnzb=n
    call wallocate_1i(sisb,n+1,infoloc)
    call wallocate_1d(ssb,nnzb,infoloc)
    call wallocate_1i(sjsb,nnzb,infoloc)
    if (infoloc/=0) then
       info=-1
       return
    end if

    do i=1,n
       sisb(i)=i
       sjsb(i)=i
       ssb(i)=DONE
    enddo
    sisb(n+1)=nnzb+1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! FORMAT CONVERSION TO CSR-UPPER for PARDISO !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if ((UPLO=='F').or.(UPLO=='f')) then !!! full csr to upper-csr
       call wallocate_1i(sisa,n+1,infoloc)
if (infoloc/=0) return
      
!!<<<
       call wallocate_1d(ssa,1,infoloc) ! dummy
if (infoloc/=0) return
       call wallocate_1i(sjsa,1,infoloc) !dummy
if (infoloc/=0) return
      
!!>>>
       opt=1
       call dcsr2csr_up(opt,N,sa,isa,jsa,ssa,sisa,sjsa)
       nnza=sisa(n+1)-1 
!!<<<
       call wdeallocate_1d(ssa)
       call wdeallocate_1i(sjsa)
!!>>>

       call wallocate_1d(ssa,nnza,infoloc)
if (infoloc/=0) return
       call wallocate_1i(sjsa,nnza,infoloc)
if (infoloc/=0) return

       opt=2
       call dcsr2csr_up(opt,N,sa,isa,jsa,ssa,sisa,sjsa)
      
!!!!!!!!!!!!!!!!!!!!!!!!!
    elseif ((UPLO=='U').or.(UPLO=='u')) then !! upper-csr already 

       nnza=isa(n+1)-1
       ssa => sa(1:nnza)
       sisa => isa(1:n+1)
       sjsa => jsa(1:nnza)

    elseif ((UPLO=='L').or.(UPLO=='l')) then !!! lower-csr to upper-csr

       nnza=isa(n+1)-1
       call wallocate_1d(ssa,nnza,infoloc)
if (infoloc/=0) return
       call wallocate_1i(sjsa,nnza,infoloc)
if (infoloc/=0) return
       call wallocate_1i(sisa,n+1,infoloc)
if (infoloc/=0) return

       call dcsr_transpose(N,sa,isa,jsa,ssa,sisa,sjsa)
    endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! FEAST INITIALIZATION 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call wallocate_2d(Aq,M0,M0,infoloc)
    if (infoloc/=0) return
    call wallocate_2d(Sq,M0,M0,infoloc)
    if (infoloc/=0) return
    call wallocate_2d(work,N,M0,infoloc)
    if (infoloc/=0) return
    call wallocate_2z(workc,N,M0,infoloc)
    if (infoloc/=0) return
    call wallocate_1i(isaz,n+1,infoloc)
    if (infoloc/=0) return
    !!<<<
    call wallocate_1z(saz,1,infoloc) ! dummy
    if (infoloc/=0) return
    call wallocate_1i(jsaz,1,infoloc)! dummy
    if (infoloc/=0) return
    !!>>>
    opt=1
    call zdaddcsr(N,opt,ONEC,ssa,sisa,sjsa,ONEC,ssb,sisb,sjsb,saz,isaz,jsaz) !! get isaz
    nnz=isaz(n+1)-1
    !!<<<
    call wdeallocate_1z(saz)
    call wdeallocate_1i(jsaz)
    !!>>>
    call wallocate_1z(saz,nnz,infoloc)
    if (infoloc/=0) return
    call wallocate_1i(jsaz,nnz,infoloc)
    if (infoloc/=0) return

    opt=2
    call zdaddcsr(N,opt,ONEC,ssa,sisa,sjsa,ONEC,ssb,sisb,sjsb,saz,isaz,jsaz) !! get jsaz
!!!!!!!!!!!!!!!
    call wallocate_2z(caux,N,M0,infoloc)
    if (infoloc/=0) return

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  analysis step (symbolic factorizarion for pardiso)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    MAXFCT=1 ! one factorization to consider normal+transpose
    MTYPE=6  ! complex and symmetric 
    call pardisoinit(PT,MTYPE,IPARM)

!!!!!!!!!!!!
if (fpm(64)==1) then
do i=1,64
if (fpm(64+i)/=-9999) iparm(i)=fpm(64+i)
enddo
endif
!!!!!!!!!!!!

    !IPARM(3)=fpm(10) !omp_num_threads !! openmp -number of threads
    IPARM(6)=1 ! solution and rhs are input/output, attention caux is always used
    MSGLVL=0 !0- no output, 1- output
    PHASE=11
    MNUM=1
    call PARDISO(PT,MAXFCT,MNUM,MTYPE,PHASE,n,saz,isaz,jsaz,idum,M0,IPARM,MSGLVL,workc(:,1:M0),caux(:,1:M0),infoloc) 
    if (infoloc/=0) then
       info=-2
       return
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!! FEAST !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

    ijob=-1 ! initialization

    do while (ijob/=0)

       call dfeast_srci(ijob,N,Ze,work,workc,Aq,Sq,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
       select case(ijob)

       case(10) !! Factorize (ZeB-A)
          opt=3
          call zdaddcsr(N,opt,-ONEC,ssa,sisa,sjsa,Ze,ssb,sisb,sjsb,saz,isaz,jsaz) !! get saz

          PHASE=22
          call PARDISO(PT,MAXFCT,MNUM,MTYPE,PHASE,n,saz,isaz,jsaz,idum,M0,IPARM,MSGLVL,workc(:,1:M0),caux(:,1:M0),infoloc)
          if (infoloc/=0) info=-2

       case(11) !! Solve (ZeB-A)x=work2(1:N,1:M0) result in to work2
          PHASE=33
          call PARDISO(PT,MAXFCT,MNUM,MTYPE,PHASE,n,saz,isaz,jsaz,idum,M0,IPARM,MSGLVL,workc(:,1:M0),caux(:,1:M0),infoloc)
          if (infoloc/=0) info=-2

       case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

          call dcsrmm('U',N,N,fpm(25),DONE,ssa,sisa,sjsa,X(1,fpm(24)),DZERO,work(1,fpm(24)))

       case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

        call DLACPY( 'F', N, fpm(25),X(1,fpm(24)), N, work(1,fpm(24)), N )
       end select
    end do
!!!!!!!!!!!!!!!!!!!!!!!!
!!!! release memory
!!!!!!!!!!!!!!!!!!!!!!!!
    PHASE=-1 

    call PARDISO(PT,MAXFCT,MNUM,MTYPE,PHASE,n,saz,isaz,jsaz,idum,M0,IPARM,MSGLVL,workc,caux,infoloc)
    if (infoloc/=0) then
       info=-2
       return
    end if

 

    call wdeallocate_2d(Aq)
    call wdeallocate_2d(Sq)
    call wdeallocate_2d(work)
    call wdeallocate_2z(workc)
    call wdeallocate_2z(caux)

    call wdeallocate_1z(saz)
    call wdeallocate_1i(isaz)
    call wdeallocate_1i(jsaz)


    if ((UPLO/='U').and.(UPLO/='u')) then
       call wdeallocate_1d(ssa)
       call wdeallocate_1i(sisa)
       call wdeallocate_1i(sjsa)
end if

       call wdeallocate_1d(ssb)
       call wdeallocate_1i(sisb)
       call wdeallocate_1i(sjsb)
   

  end subroutine dfeast_scsrev




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




subroutine dfeast_iscsrev(UPLO,N,sa,isa,jsa,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
    !  Purpose 
    !  =======
    !  FEAST SPARSE INTERFACE (SPARSE CSR FORMAT)
    !  Solve the standard Ax=ex eigenvalue problem
    !  
    !  A REAL SYMMETRIC
    ! 
    !  DOUBLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  UPLO       (input)        CHARACTER: specifies whether the full part 'F', or the upper 'U' or lower 'L'
    !                                       triangular part of the matrix(ces) is(are) being supplied.
    !  N          (input)        INTEGER: Size system 
    !  sa         (input)        REAL DOUBLE PRECISION (isa(N+1)-1):  Matrix A- CSR format 
    !  isa        (input)        INTEGER(N+1): CSR row array of Matrix A   
    !  jsa        (input)        INTEGER(isa(N+1)-1): CSR column array of Matrix A
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  fpm (input/output)     INTEGER(*) : FEAST parameters (see FEAST documentation)
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
    !                                                        if option fpm(6)=1 selected                           
    !  info       (output)       INTEGER: Error handling (0: successful exit) -- see FEAST documentation
    !=====================================================================
    ! Eric Polizzi 2009-2013
    ! ====================================================================

    implicit none
    include 'f90_noruntime_interface.fi'
    character(len=1) :: UPLO
    integer :: N
    double precision,dimension(*),target:: sa
    integer,dimension(*),target:: isa,jsa
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer :: ijob,infoloc
    complex(kind=(kind(1.0d0))) :: Ze
    complex(kind=(kind(1.0d0))),dimension(:),pointer :: saz,zdummy,ztemp,usaz,alpha,beta
    integer,dimension(:),pointer :: isaz,jsaz
    complex(kind=(kind(1.0d0))),dimension(:,:),pointer ::workc,caux,vv,r,oworkc,d,ztemp2,ztemp3,cX,ctemp
    double precision, dimension(:,:),pointer ::work,Aq,Sq,vX
    double precision,parameter :: DONE=1.0d0, DZERO=0.0d0
    complex(kind=(kind(1.0d0))),parameter :: ONEC=(DONE,DZERO)
!!!!! csr-upper format
    double precision,dimension(:),pointer :: ssa,ssb,nres,vE
    integer,dimension(:), pointer :: sisa,sjsa,sisb,sjsb,uisa,ujsa
    integer :: i,k,im,it,itmax,j 
    integer :: opt,nnza,nnzb,nnz
double precision :: err,omega,ares
logical :: comb
integer :: infob
integer :: mode1,rank

rank=0



!!!!!!!!!!!!!! Check INPUT PARAMETERS
    INFO = 0
    IF ((UPLO/='F').and.(UPLO/='f').and.(UPLO/='L').and.(UPLO/='l').and.(UPLO/='U').and.(UPLO/='u')) THEN
       INFO=-101
    ELSE IF ( N<=0 ) THEN
       INFO = -102
    END IF
    IF( INFO.NE.0 ) THEN
       CALL XERBLA( 'DFEAST_SCSREV', -INFO+100 )
       RETURN
    END IF

    info=-1 ! by default

infoloc=0

nnz=isa(n+1)-1



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! FEAST INITIALIZATION 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call wallocate_2d(vX,N,M0,infoloc) 

   call wallocate_1d(nres,M0,infoloc)
    call wallocate_2d(Aq,M0,M0,infoloc)
    call wallocate_2d(Sq,M0,M0,infoloc)
    call wallocate_2d(work,N,M0,infoloc)
    call wallocate_2z(workc,N,M0,infoloc)
call wallocate_2z(oworkc,N,M0,infoloc)

    call wallocate_2z(r,N,M0,infoloc)
    call wallocate_1z(saz,nnz,infoloc) 
    call wallocate_1z(zdummy,n,infoloc) 
 call wallocate_1z(ztemp,n,infoloc) 
 
!!!!!!!!!!!!!!!
    call wallocate_2z(caux,N,M0,infoloc)
call wallocate_2z(ztemp2,N,M0,infoloc)
call wallocate_2z(ztemp3,N,M0,infoloc)
call wallocate_2z(d,N,M0,infoloc)
call wallocate_1z(alpha,M0,infoloc)
call wallocate_1z(beta,M0,infoloc)


call wallocate_1z(usaz,nnz,infoloc) 
call wallocate_1i(ujsa,nnz,infoloc) 
call wallocate_1i(uisa,n+1,infoloc) 

    allocate(vE(M0))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
im=32!50!6
    allocate(vv(N,im+1))


itmax=4000

!im=itmax

print *,'itmax',itmax
err=1d-10!16


omega=0.1d0

mode1=0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!! FEAST !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

    ijob=-1 ! initialization

    do while (ijob/=0)

       call dfeast_srci(ijob,N,Ze,work,workc,Aq,Sq,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
       select case(ijob)

       case(10) !! Factorize (ZeB-A)
       !   opt=3
       !   call zdaddcsr(N,opt,-ONEC,ssa,sisa,sjsa,Ze,ssb,sisb,sjsb,saz,isaz,jsaz) !! get saz

!saz(1:nnz)=-(1.0d0,0.0d0)*sa(1:nnz)
!do  i=1,n
!do k=isa(i),isa(i+1)-1
!if (jsa(k)==i) saz(k)=saz(k)+Ze
!enddo
!enddo





         
       case(11) !! Solve (ZeB-A)x=work2(1:N,1:M0) result in to work2
saz(1:nnz)=-(1.0d0,0.0d0)*sa(1:nnz)
do  i=1,n
do k=isa(i),isa(i+1)-1
if (jsa(k)==i) saz(k)=saz(k)+Ze
enddo
enddo


!!!! extract diagonal matrix
     do i=1,N
        do k=isa(i),isa(i+1)-1
           if (jsa(k)==i) zdummy(i)=saz(k)
        enddo
     enddo
!!!! tranpose of the matrix
 call zcsr_htranspose(N,saz,isa,jsa,usaz,uisa,ujsa)
usaz=conjg(usaz)


 

caux(1:N,1:M0)=(0.0d0,0.0d0)
if (loop>0) then
do i=1,M0
caux(1:N,i)=workc(1:N,i)/(Ze-E(i))
enddo
endif
         



oworkc(1:N,1:M0)=workc(1:N,1:M0)






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
select case(fpm(43))


case(10) !!!!! gmrez
caux(1:N,1:M0)=(0.0d0,0.0d0) !!! zero initial guess

do i=1,M0
call gmrez(n,im,workc(1,i),caux(1,i),vv,err,itmax,-6,saz,isa,jsa,UPLO)
enddo
!itmax=itmax-2
!im=im-2
!stop

case(1) !! jacobi (with relaxation)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  Iterative  Jacobi !!!!!!!!!!!!!!!!!!!!!!
!!at iteration k: x(k)=x(k-1)+D^-1*(b-A*x(k-1)) !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! residual
 ! new residual (b-A*x(k))
        r(1:N,1:M0)=workc(1:N,1:M0)
        call zcsrmm('L',N,N,M0,(-1.0d0,0.0d0),saz,isa,jsa,caux,(1.0d0,0.0d0),r)
do i=1,M0
     nres(i)=sum(abs(r(:,i)))/sum(abs(workc(:,i))) ! norm relative residual 
enddo

it=0
     do while ((maxval(nres)>err).and.(it<itmax))
        it=it+1 ! # of iterations

        ! x(k)=x(k-1)+D^-1*r
        do i=1,N
           caux(i,1:M0)=caux(i,1:M0)+0.1d0*r(i,1:M0)/zdummy(i)
        enddo

        ! new residual (b-A*x(k))
        r(1:N,1:M0)=workc(1:N,1:M0)
        call zcsrmm('L',N,N,M0,(-1.0d0,0.0d0),saz,isa,jsa,caux,(1.0d0,0.0d0),r)
do i=1,M0
     nres(i)=sum(abs(r(:,i)))/sum(abs(workc(:,i))) ! norm relative residual 
enddo    
!print *,it,maxval(nres),minval(nres)
     end do
!stop

case(11) !! Simple iteration

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  Iterative  Jacobi !!!!!!!!!!!!!!!!!!!!!!
!!at iteration k: x(k)=(A*x(k-1)+b)/z !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! residual
 ! new residual (b-A*x(k))
 !       r(1:N,1:M0)=workc(1:N,1:M0)
 !       call zcsrmm('L',N,N,M0,(-1.0d0,0.0d0),saz,isa,jsa,caux,(1.0d0,0.0d0),r)
!do i=1,M0
!     nres(i)=sum(abs(r(:,i)))/sum(abs(workc(:,i))) ! norm relative residual 
!enddo

it=0
     do while ((it<itmax))
        it=it+1 ! # of iterations
r(1:N,1:M0)=workc(1:N,1:M0)
call zcsrmm('L',N,N,M0,(1.0d0,0.0d0),saz,isa,jsa,caux,(1.0d0,0.0d0),r)

        ! x(k)=x(k-1)+D^-1*r
        do i=1,N
           caux(i,1:M0)=r(i,1:M0)/Ze
        enddo

        ! new residual (b-A*x(k))
       ! r(1:N,1:M0)=workc(1:N,1:M0)
       ! call zcsrmm('L',N,N,M0,(-1.0d0,0.0d0),saz,isa,jsa,caux,(1.0d0,0.0d0),r)
!do i=1,M0
!     nres(i)=sum(abs(r(:,i)))/sum(abs(workc(:,i))) ! norm relative residual 
!enddo    
!print *,it,maxval(nres),minval(nres)
     end do

case(12) !! jacobi (with preconditioenr)
      
!caux(1:N,1:M0)=(0.0d0,0.0d0)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!  Iterative  Jacobi !!!!!!!!!!!!!!!!!!!!!!
  !!at iteration k: x(k)=x(k-1)+D^-1*(b-A*x(k-1)) !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (rank==0) print *,'jacobi'
if (mode1/=0) then      
if (rank==0) print *,'jacobi with preconditioner',mode1
allocate(cX(N,mode1))
cX(1:N,1:mode1)=(1.0d0,0.0d0)*vX(1:N,1:mode1)
allocate(ctemp(mode1,M0))
endif
  !!! residual
   ! new residual (b-A*x(k))
          r(1:N,1:M0)=workc(1:N,1:M0)
          call zcsrmm('L',N,N,M0,(-1.0d0,0.0d0),saz,isa,jsa,caux,(1.0d0,0.0d0),r)
  do i=1,M0
       nres(i)=sum(abs(r(:,i)))/sum(abs(workc(:,i))) ! norm relative residual 
  enddo
      
  it=0
       do while ((maxval(nres)>err).and.(it<itmax))
          it=it+1 ! # of iterations
      
          if (mode1==0) then
!             x(k)=x(k-1)+D^-1*r

          do i=1,N
             caux(i,1:M0)=caux(i,1:M0)+0.1d0*r(i,1:M0)/zdummy(i)
          enddo

else



 call ZGEMM('C','N',mode1,M0,N,(1.0d0,0.0d0),cX,N,r,N,(0.0d0,0.0d0),ctemp,mode)
do i=1,mode1
ctemp(i,1:M0)=ctemp(i,1:M0)/(Ze-vE(i)*(1.0d0,0.0d0))
enddo
call ZGEMM('N','N',N,M0,mode1,(1.0d0,0.0d0),cX,N,ctemp,mode,(1.0d0,0.0d0),caux,N)

endif


          ! new residual (b-A*x(k))
          r(1:N,1:M0)=workc(1:N,1:M0)
          call zcsrmm('L',N,N,M0,(-1.0d0,0.0d0),saz,isa,jsa,caux,(1.0d0,0.0d0),r)
  do i=1,M0
       nres(i)=sum(abs(r(:,i)))/sum(abs(workc(:,i))) ! norm relative residual 
  enddo
  !print *,it,maxval(nres),minval(nres)
       end do
  !stop

if (mode1/=0) then
 deallocate(ctemp)
deallocate(cX)
endif
!stop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
case(2)  !!! Jacobi on (A^H.A)x=A^Hf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! new diagonal
do i=1,N
r(1:N,1)=(0.0d0,0.0d0)
r(i,1)=(1.0d0,0.0d0)
call zcsrmm('L',N,N,1,(1.0d0,0.0d0),saz,isa,jsa,r,(0.0d0,0.0d0),ztemp)
zdummy(i)=sum(abs(ztemp)**2)
enddo



r(1:N,1:M0)=workc(1:N,1:M0)
call zcsrmm('L',N,N,M0,(1.0d0,0.0d0),conjg(saz),isa,jsa,r,(0.0d0,0.0d0),workc)



!!! residual
 ! new residual (b-A^H*A*x(k))
        r(1:N,1:M0)=caux(1:N,1:M0)
        call zcsrmm('L',N,N,M0,(1.0d0,0.0d0),saz,isa,jsa,r,(0.0d0,0.0d0),caux)

        r(1:N,1:M0)=workc(1:N,1:M0)
        call zcsrmm('L',N,N,M0,(-1.0d0,0.0d0),conjg(saz),isa,jsa,caux,(1.0d0,0.0d0),r)

do i=1,M0
     nres(i)=sum(abs(r(:,i)))/sum(abs(workc(:,i))) ! norm relative residual 
enddo
    

it=0
     do while ((maxval(nres)>err).and.(it<itmax))
        it=it+1 ! # of iterations

        ! x(k)=x(k-1)+D^-1*r
        do i=1,N
           caux(i,1:M0)=caux(i,1:M0)+0.1d0*r(i,1:M0)/zdummy(i)
        enddo

        ! new residual (b-A*x(k))
        r(1:N,1:M0)=caux(1:N,1:M0)
        call zcsrmm('L',N,N,M0,(1.0d0,0.0d0),saz,isa,jsa,r,(0.0d0,0.0d0),caux)

        r(1:N,1:M0)=workc(1:N,1:M0)
        call zcsrmm('L',N,N,M0,(-1.0d0,0.0d0),conjg(saz),isa,jsa,caux,(1.0d0,0.0d0),r)
do i=1,M0
     nres(i)=sum(abs(r(:,i)))/sum(abs(workc(:,i))) ! norm relative residual 
enddo    

print *,it,maxval(nres),minval(nres)

     end do



case(3) !!!!!!!!!!!!! Gauss-Seidel
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  Iterative  Gauss-Seidel !!!!!!!!!!!!!!!!
     !!at iteration k: x(k)=x(k-1)+(D+L)^-1*(b-A*x(k-1)) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


 r(1:N,1:M0)=workc(1:N,1:M0)
        call zcsrmm('L',N,N,M0,(-1.0d0,0.0d0),saz,isa,jsa,caux,(1.0d0,0.0d0),r)
do i=1,M0
     nres(i)=sum(abs(r(:,i)))/sum(abs(workc(:,i))) ! norm relative residual 
enddo
  

     it=0
     do while ((maxval(nres)>err).and.(it<itmax))
        it=it+1 ! # of iterations

        ! x(k)=x(k-1)+(D+L)^-1*r
        call zcsrsv('L',N,M0,saz,isa,jsa,r)
        caux=caux+r


        ! new residual (b-A*x(k))
 r(1:N,1:M0)=workc(1:N,1:M0)
        call zcsrmm('L',N,N,M0,(-1.0d0,0.0d0),saz,isa,jsa,caux,(1.0d0,0.0d0),r)
do i=1,M0
     nres(i)=sum(abs(r(:,i)))/sum(abs(workc(:,i))) ! norm relative residual 
enddo
    

      !  print *,it,nres
     end do


case(4)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  Iterative  SSOR (omega=1)!!!!!!!!!!!!!!!!!!!!!!!!
     !!at iteration k: x(k)=x(k-1)+[(D+U)^-1]*D*[(D+L)^-1]*(b-A*x(k-1)) !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


do i=1,N
do k=isa(i),isa(i+1)-1
if (i==jsa(k)) saz(k)=saz(k)/omega
enddo
do k=uisa(i),uisa(i+1)-1
if (i==ujsa(k)) usaz(k)=usaz(k)/omega
enddo
enddo

zdummy=zdummy*(2.0d0-omega)/omega


 r(1:N,1:M0)=workc(1:N,1:M0)
        call zcsrmm('L',N,N,M0,(-1.0d0,0.0d0),saz,isa,jsa,caux,(1.0d0,0.0d0),r)
do i=1,M0
     nres(i)=sum(abs(r(:,i)))/sum(abs(workc(:,i))) ! norm relative residual 
enddo
  

     it=0
     do while ((maxval(nres)>err).and.(it<itmax))
        it=it+1 ! # of iterations

        ! x(k)=x(k-1)+(D+L)^-1*r
        call zcsrsv('L',N,M0,saz,isa,jsa,r)
        do i=1,N
           r(i,1:M0)=r(i,1:M0)*zdummy(i) !D*(D+L)^-1*r
        enddo
        call zcsrsv('U',N,M0,usaz,uisa,ujsa,r)
        caux=caux+r


        ! new residual (b-A*x(k))
 r(1:N,1:M0)=workc(1:N,1:M0)
        call zcsrmm('L',N,N,M0,(-1.0d0,0.0d0),saz,isa,jsa,caux,(1.0d0,0.0d0),r)
do i=1,M0
     nres(i)=sum(abs(r(:,i)))/sum(abs(workc(:,i))) ! norm relative residual 
enddo
    




     
!print *,maxval(nres),minval(nres)

     end do

do i=1,N
do k=isa(i),isa(i+1)-1
if (i==jsa(k)) saz(k)=saz(k)*omega
enddo
enddo



!!!!!!!!!!!!!!!!!!!!!!!!!!1
case(5) !! bicgstab


caux(1:N,1:M0)=(0.0d0,0.0d0)

comb=.false.!.true.   ! comments/timing yes or no         
it=itmax
call zbicgstab('L',N,saz,isa,jsa,M0,workc,caux,nres,ares,it,err,comb,infob) 

!print *,'bicgstab',it,ares


!!!!!!!!!!!!!!!!!!!!!!!!!!1
case(51) !! bicgstab- spectral preconditioner


caux(1:N,1:M0)=(0.0d0,0.0d0)

comb=.false.!.true.   ! comments/timing yes or no         
it=itmax
!if (loop<=1) then
if (mode1==0) then
if (rank==0) print *,'bicgstab'
call zbicgstab('L',N,saz,isa,jsa,M0,workc,caux,nres,ares,it,err,comb,infob) 
else
if (rank==0) print *,'bicgstabi- with preconditioner',mode1
!comb=.true.
!call zbicgstab('L',N,saz,isa,jsa,M0,workc,caux,nres,ares,it,err,comb,infob) 
call zbicgstabp('L',N,saz,isa,jsa,vX,vE,Ze,mode1,M0,workc,caux,nres,ares,it,err,comb,infob) 
endif
!print *,'bicgstab',it,ares






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
case(6) !! CG on normal equation





! update rhs
r(1:N,1:M0)=workc(1:N,1:M0)
call zcsrmm('L',N,N,M0,(1.0d0,0.0d0),conjg(saz),isa,jsa,r,(0.0d0,0.0d0),workc)

!! residual
 ! new residual (b-A^H*A*x(k))
        r(1:N,1:M0)=caux(1:N,1:M0)
        call zcsrmm('L',N,N,M0,(1.0d0,0.0d0),saz,isa,jsa,r,(0.0d0,0.0d0),caux)

        r(1:N,1:M0)=workc(1:N,1:M0)
        call zcsrmm('L',N,N,M0,(-1.0d0,0.0d0),conjg(saz),isa,jsa,caux,(1.0d0,0.0d0),r)

do i=1,M0
     nres(i)=sum(abs(r(:,i)))/sum(abs(workc(:,i))) ! norm relative residual 
enddo



     ! search direction initial
     allocate(d(1:N,1:M0))

     d(1:N,1:M0)=r(1:N,1:M0)

     it=0
     do while ((maxval(nres)>err).and.(it<itmax))
        it=it+1 ! # of iterations

       ! A*d
!        call dcsrmm('F',N,N,1,1.0d0,sa,isa,jsa,d,0.0d0,dummy)

        call zcsrmm('L',N,N,M0,(1.0d0,0.0d0),saz,isa,jsa,d,(0.0d0,0.0d0),ztemp3)
        call zcsrmm('L',N,N,M0,(1.0d0,0.0d0),conjg(saz),isa,jsa,ztemp3,(0.0d0,0.0d0),ztemp2)


       ! alpha=d^t.d/(d^t.A.d)
do i=1,M0 
        alpha(i)=sum(abs(r(:,i))**2)/sum(conjg(d(:,i))*ztemp2(:,i))
enddo

        ! x(k)=x(k-1)+alpha*d(k-1) (new iterate)
do i=1,M0
           caux(:,i)=caux(:,i)+alpha(i)*d(:,i)
enddo

        ! beta
do i=1,M0
        beta(i)=(1.0d0,0.0d0)/sum(abs(r(:,i))**2)
enddo

        ! new residual r=r-alpha*A*p
        r(:,i)=r(:,i)-alpha(i)*ztemp2(:,i)


do i=1,M0
     nres(i)=sum(abs(r(:,i)))/sum(abs(workc(:,i))) ! norm relative residual 
enddo

        !Update search direction
do i=1,M0
        beta(i)=beta(i)*sum(abs(r(:,i))**2)
        d(:,i)=r(:,i)+beta(i)*d(:,i)
     end do


enddo





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
case(7) !! CG




!! residual
r(1:N,1:M0)=workc(1:N,1:M0)
        call zcsrmm('L',N,N,M0,(-1.0d0,0.0d0),saz,isa,jsa,caux,(1.0d0,0.0d0),r)
do i=1,M0
     nres(i)=sum(abs(r(:,i)))/sum(abs(workc(:,i))) ! norm relative residual 
enddo



     ! search direction initial
     allocate(d(1:N,1:M0))

     d(1:N,1:M0)=r(1:N,1:M0)

     it=0
     do while ((maxval(nres)>err).and.(it<itmax))
        it=it+1 ! # of iterations

       ! A*d
!        call dcsrmm('F',N,N,1,1.0d0,sa,isa,jsa,d,0.0d0,dummy)

        call zcsrmm('L',N,N,M0,(1.0d0,0.0d0),saz,isa,jsa,d,(0.0d0,0.0d0),ztemp2)
!        call zcsrmm('L',N,N,M0,(1.0d0,0.0d0),conjg(saz),isa,jsa,ztemp3,(0.0d0,0.0d0),ztemp2)


       ! alpha=r^t.r/(d^t.A.d)
do i=1,M0 
        alpha(i)=sum(abs(r(:,i))**2)/sum(conjg(d(:,i))*ztemp2(:,i))
enddo

        ! x(k)=x(k-1)+alpha*d(k-1) (new iterate)
do i=1,M0
           caux(:,i)=caux(:,i)+alpha(i)*d(:,i)
enddo

        ! beta
do i=1,M0
        beta(i)=(1.0d0,0.0d0)/sum(abs(r(:,i))**2)
enddo

        ! new residual r=r-alpha*A*p
        r(:,i)=r(:,i)-alpha(i)*ztemp2(:,i)


do i=1,M0
     nres(i)=sum(abs(r(:,i)))/sum(abs(workc(:,i))) ! norm relative residual 
enddo

        !Update search direction
do i=1,M0
        beta(i)=beta(i)*sum(abs(r(:,i))**2)
        d(:,i)=r(:,i)+beta(i)*d(:,i)
     end do


enddo



end select



        ! new residual (b-A*x(k))
        r(1:N,1:M0)=oworkc(1:N,1:M0)
        call zcsrmm('L',N,N,M0,(-1.0d0,0.0d0),saz,isa,jsa,caux,(1.0d0,0.0d0),r)
do i=1,M0
     nres(i)=sum(abs(r(:,i)))/sum(abs(oworkc(:,i))) ! norm relative residual 
enddo
    



print *,maxval(nres),minval(nres)

!stop

workc(1:N,1:M0)=caux(1:N,1:M0)


       case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

          call dcsrmm(UPLO,N,N,fpm(25),DONE,sa,isa,jsa,X(1,fpm(24)),DZERO,work(1,fpm(24)))

       case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

        call DLACPY( 'F', N, fpm(25),X(1,fpm(24)), N, work(1,fpm(24)), N )




       case(60) 
!             vX(1:N,1:M0)=X(1:N,1:M0)

!Vx(1:N,1:M0)=0.0d0
mode1=0
do i=1,M0
if (res(i)<=1d-3) then
mode1=mode1+1
vX(1:N,mode1)=X(1:N,i)
vE(mode1)=E(i)
endif
enddo

!print *,'here',j,mode,M0,E(j+i-1)

       end select
    end do
!



    call wdeallocate_2d(Aq)
    call wdeallocate_2d(Sq)
    call wdeallocate_2d(work)
    call wdeallocate_2z(workc)
    call wdeallocate_2z(caux)

    call wdeallocate_1z(saz)
    !call wdeallocate_1i(isaz)
    !call wdeallocate_1i(jsaz)


       !call wdeallocate_1d(ssb)
       !call wdeallocate_1i(sisb)
       !call wdeallocate_1i(sjsb)
   

  end subroutine dfeast_iscsrev




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!










!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  
  subroutine zfeast_hcsrgv(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
    !  Purpose 
    !  =======
    !  FEAST SPARSE INTERFACE (SPARSE CSR FORMAT)
    !  Solve the generalized Ax=eBx eigenvalue problem
    !  
    !  A COMPLEX HERMITIAN, B HERMITIAN POSITIVE DEFINITE 
    ! 
    !  DOUBLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  UPLO       (input)        CHARACTER: specifies whether the full part 'F', or the upper 'U' or lower 'L'
    !                                       triangular part of the matrix(ces) is(are) being supplied.
    !  N          (input)        INTEGER: Size system 
    !  sa         (input)        COMPLEX DOUBLE PRECISION (isa(N+1)-1):  Matrix A- CSR format 
    !  isa        (input)        INTEGER(N+1): CSR row array of Matrix A   
    !  jsa        (input)        INTEGER(isa(N+1)-1): CSR column array of Matrix A
    !  sb         (input)        COMPLEX DOUBLE PRECISION (isb(N+1)-1):  Matrix B- CSR format 
    !  isb        (input)        INTEGER(N+1): CSR row array of Matrix B   
    !  jsb        (input)        INTEGER(isb(N+1)-1): CSR column array of Matrix B
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  fpm (input/output) INTEGER(*) : FEAST parameters (see FEAST documentation)
    !  epsout     (output)       REAL DOUBLE PRECISION : Error on the trace
    !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
    !  Emin,Emax  (input)        REAL DOUBLE PRECISION: search interval
    !  M0         (input/output) INTEGER: Size subspace
    !  lambda     (output)       REAL DOUBLE PRECISION(M0)   : Eigenvalues -solution
    !  q          (input/output) COMPLEX DOUBLE PRECISION(N,M0) : 
    !                                                       On entry: subspace initial guess if fpm(5)=1 
    !                                                       On exit : Eigenvectors-solution
    !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
    !  res        (output)       REAL DOUBLE PRECISION(M0) : Relative Residual of the solution (1-norm)
    !                                                        if option fpm(6)=1 selected                           
    !  info       (output)       INTEGER: Error handling (0: successful exit) -- see FEAST documentation
    !=====================================================================
    ! Eric Polizzi 2009-2013
    ! ====================================================================

    implicit none
    include 'f90_noruntime_interface.fi'
    character(len=1) :: UPLO
    integer :: N
    complex(kind=kind(1.0d0)),dimension(*),target:: sa,sb
    integer,dimension(*),target:: isa,jsa,isb,jsb
    integer,dimension(*) :: fpm
    double precision :: epsout 
    integer :: loop
    double precision :: Emin,Emax
    integer :: M0
    double precision,dimension(*)  :: E
    complex(kind=kind(1.0d0)),dimension(N,*):: X
    integer :: mode
    double precision,dimension(*)    :: res
    integer :: info
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer :: ijob,infoloc
    complex(kind=(kind(1.0d0))) :: Ze
    complex(kind=(kind(1.0d0))),dimension(:),pointer :: saz,tsaz
    integer,dimension(:),pointer :: isaz,jsaz
    complex(kind=(kind(1.0d0))),dimension(:,:),pointer ::work1,work2,caux,zAq,zSq
    double precision,parameter :: DONE=1.0d0, DZERO=0.0d0
    complex(kind=(kind(1.0d0))),parameter :: ONEC=(DONE,DZERO), ZEROC=(DZERO,DZERO)
!!!!!for pardiso
    integer(8),dimension(64) :: pt1,pt2
    integer,dimension(64) :: iparm
    integer :: mtype
    integer :: MAXFCT,MNUM,PHASE,MSGLVL
    integer :: idum,nrhs
    double precision :: ddum1
!!!!! full csr format
    complex(kind=kind(1.0d0)),dimension(:),pointer :: ssa,ssb
    integer,dimension(:), pointer :: sisa,sjsa,sisb,sjsb
    integer :: i,k 
    integer :: opt,nnza,nnzb,nnz


!!!!!!!!!!!!!! Check INPUT PARAMETERS
    INFO = 0
    IF ((UPLO/='F').and.(UPLO/='f').and.(UPLO/='L').and.(UPLO/='l').and.(UPLO/='U').and.(UPLO/='u')) THEN
       INFO=-101
    ELSE IF ( N<=0 ) THEN
       INFO = -102
    END IF
    IF( INFO.NE.0 ) THEN
       CALL XERBLA( 'ZFEAST_HCSRGV', -INFO+100 )
       RETURN
    END IF

    info=-1 ! by default
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! FORMAT CONVERSION TO full CSR for PARDISO !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Remark: everything should be csr since we work with zS-H which is actually unsymmetric 


    if ((UPLO=='F').or.(UPLO=='f')) then !!! full csr already

       nnza=isa(n+1)-1
       ssa => sa(1:nnza)
       sisa => isa(1:n+1)
       sjsa => jsa(1:nnza)

       nnzb=isb(n+1)-1
       ssb => sb(1:nnza)
       sisb => isb(1:n+1)
       sjsb => jsb(1:nnza)

    else !! upper-csr or lower-csr to full csr

       nnza=2*(isa(n+1)-1) ! may be overestimated
       nnzb=2*(isb(n+1)-1)
       call wallocate_1z(ssa,nnza,infoloc)
       if (infoloc/=0) return
       call wallocate_1i(sjsa,nnza,infoloc)
       if (infoloc/=0) return
       call wallocate_1z(ssb,nnzb,infoloc)
       if (infoloc/=0) return
       call wallocate_1i(sjsb,nnzb,infoloc)
       if (infoloc/=0) return
       call wallocate_1i(sisa,n+1,infoloc)
       if (infoloc/=0) return
       call wallocate_1i(sisb,n+1,infoloc)
       if (infoloc/=0) return

       call zhcsr_uplo_to_csr(N,sa,isa,jsa,ssa,sisa,sjsa)
       call zhcsr_uplo_to_csr(N,sb,isb,jsb,ssb,sisb,sjsb)

    endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! FEAST INITIALIZATION 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call wallocate_2z(zAq,M0,M0,infoloc)
    if (infoloc/=0) return
    call wallocate_2z(zSq,M0,M0,infoloc)
    if (infoloc/=0) return
    call wallocate_2z(work1,N,M0,infoloc)
    if (infoloc/=0) return
    call wallocate_2z(work2,N,M0,infoloc)
    if (infoloc/=0) return
    call wallocate_1i(isaz,n+1,infoloc)
    if (infoloc/=0) return
    !!<<<
    call wallocate_1z(saz,1,infoloc) ! dummy
    if (infoloc/=0) return
    call wallocate_1i(jsaz,1,infoloc)! dummy
    if (infoloc/=0) return
    !!>>>
    opt=1
    call zaddcsr(N,opt,ONEC,ssa,sisa,sjsa,ONEC,ssb,sisb,sjsb,saz,isaz,jsaz) !! get isaz
    nnz=isaz(n+1)-1
    !!<<<
    call wdeallocate_1z(saz)
    call wdeallocate_1i(jsaz)
    !!>>>
    call wallocate_1z(saz,nnz,infoloc)
    if (infoloc/=0) return
    call wallocate_1i(jsaz,nnz,infoloc)
    if (infoloc/=0) return


    opt=2
    call zaddcsr(N,opt,ONEC,ssa,sisa,sjsa,ONEC,ssb,sisb,sjsb,saz,isaz,jsaz) !! get jsaz
!!!!!!!!!!!!!!!
    call wallocate_2z(caux,N,M0,infoloc)
    if (infoloc/=0) return


 MAXFCT=1 ! One factorization to consider (default)

if (fpm(11)==0) then!!! needs the transpose option 
     MAXFCT=2 ! two factorization to consider normal+transpose
    call wallocate_1z(tsaz,nnz,infoloc)
    if (infoloc/=0) return
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  analysis step (symbolic factorizarion for pardiso)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
    MTYPE=3 ! complex and structurally symmetric 
    call pardisoinit(PT1,MTYPE,IPARM)
if (fpm(11)==0) call pardisoinit(PT2,MTYPE,IPARM)

!!!!!!!!!!!!
if (fpm(64)==1) then
do i=1,64
if (fpm(64+i)/=-9999) iparm(i)=fpm(64+i)
enddo
endif
!!!!!!!!!!!!


    !IPARM(3)=fpm(10) !omp_num_threads !! openmp -number of threads
    IPARM(6)=1 ! solution and rhs are input/output, attention caux is always used
    MSGLVL=0 !0- no output, 1- output
    PHASE=11
    MNUM=1 
    call PARDISO(PT1,MAXFCT,MNUM,MTYPE,PHASE,n,saz,isaz,jsaz,idum,M0,IPARM,MSGLVL,work2(:,1:M0),caux(:,1:M0),infoloc)
if (fpm(11)==0) then 
    MNUM=2 
    call PARDISO(PT2,MAXFCT,MNUM,MTYPE,PHASE,n,tsaz,isaz,jsaz,idum,M0,IPARM,MSGLVL,work2(:,1:M0),caux(:,1:M0),infoloc)
endif
    if (infoloc/=0) then
       info=-2
       return
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!! FEAST !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

    ijob=-1 ! initialization

    do while (ijob/=0)

       call zfeast_hrci(ijob,N,Ze,work1,work2,zAq,zSq,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
       select case(ijob)

       case(10) !! Factorize (ZeB-A)
          opt=3
          call zaddcsr(N,opt,-ONEC,ssa,sisa,sjsa,Ze,ssb,sisb,sjsb,saz,isaz,jsaz) !! get saz

          PHASE=22
          MNUM=1
          call PARDISO(PT1,MAXFCT,MNUM,MTYPE,PHASE,n,saz,isaz,jsaz,idum,M0,IPARM,MSGLVL,work2(:,1:M0),caux(:,1:M0),infoloc)
          if (infoloc/=0) info=-2

       case(11) !! Solve (ZeB-A)x=work2(1:N,1:M0) result in to work2
          IPARM(12)=0 ! normal solve
          PHASE=33
          MNUM=1
          call PARDISO(PT1,MAXFCT,MNUM,MTYPE,PHASE,n,saz,isaz,jsaz,idum,M0,IPARM,MSGLVL,work2(:,1:M0),caux(:,1:M0),infoloc)
          if (infoloc/=0) info=-2

       case(20) !! Factorize (ZeB-A)^T   
if (fpm(11)==0) then
          opt=3
          call zaddcsr(N,opt,-ONEC,ssa,sisa,sjsa,conjg(Ze),ssb,sisb,sjsb,tsaz,isaz,jsaz) !! get saz

          PHASE=22
          MNUM=2
          call PARDISO(PT2,MAXFCT,MNUM,MTYPE,PHASE,n,tsaz,isaz,jsaz,idum,M0,IPARM,MSGLVL,work2(:,1:M0),caux(:,1:M0),infoloc)
          if (infoloc/=0) info=-2
endif
       case(21) !!  Solve  (ZeB-A)^Tx=work2(1:N,1:M0) result in to work2
if (fpm(11)==0) then
          IPARM(12)=0 ! normal solve
          PHASE=33
          MNUM=2
          call PARDISO(PT2,MAXFCT,MNUM,MTYPE,PHASE,n,tsaz,isaz,jsaz,idum,M0,IPARM,MSGLVL,work2(:,1:M0),caux(:,1:M0),infoloc)
          if (infoloc/=0) info=-2
else
IPARM(12)=1 ! transpose conjugate solve
  PHASE=33
          MNUM=1
          call PARDISO(PT1,MAXFCT,MNUM,MTYPE,PHASE,n,saz,isaz,jsaz,idum,M0,IPARM,MSGLVL,work2(:,1:M0),caux(:,1:M0),infoloc)
          if (infoloc/=0) info=-2
endif

       case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

          call zhcsrmm('F',N,N,fpm(25),ONEC,ssa,sisa,sjsa,X(1,fpm(24)),ZEROC,work1(1,fpm(24)))

       case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

          call zhcsrmm('F',N,N,fpm(25),ONEC,ssb,sisb,sjsb,X(1,fpm(24)),ZEROC,work1(1,fpm(24)))

       end select
    end do
!!!!!!!!!!!!!!!!!!!!!!!!
!!!! release memory
!!!!!!!!!!!!!!!!!!!!!!!!
    PHASE=-1 

    MNUM=1
    call PARDISO(PT1,MAXFCT,MNUM,MTYPE,PHASE,n,saz,isaz,jsaz,idum,M0,IPARM,MSGLVL,work2,caux,infoloc)
    if (infoloc/=0) then
       info=-2
       return
    end if

if (fpm(11)==0) then
    MNUM=2
    call PARDISO(PT2,MAXFCT,MNUM,MTYPE,PHASE,n,tsaz,isaz,jsaz,idum,M0,IPARM,MSGLVL,work2,caux,infoloc)
    if (infoloc/=0) then
       info=-2
       return
    end if
 end if

    call wdeallocate_2z(zAq)
    call wdeallocate_2z(zSq)
    call wdeallocate_2z(work1)
    call wdeallocate_2z(work2)
    call wdeallocate_2z(caux)

    call wdeallocate_1z(saz)
    call wdeallocate_1i(isaz)
    call wdeallocate_1i(jsaz)


if (fpm(11)==0)  call wdeallocate_1z(tsaz) !! transpose option


    if ((UPLO/='F').and.(UPLO/='f')) then
       call wdeallocate_1z(ssa)
       call wdeallocate_1i(sisa)
       call wdeallocate_1i(sjsa)

       call wdeallocate_1z(ssb)
       call wdeallocate_1i(sisb)
       call wdeallocate_1i(sjsb)
    endif



  end subroutine zfeast_hcsrgv



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine zfeast_hcsrev(UPLO,N,sa,isa,jsa,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
    !  Purpose 
    !  =======
    !  FEAST SPARSE INTERFACE (SPARSE CSR FORMAT)
    !  Solve the standard Ax=ex eigenvalue problem
    !  
    !  A COMPLEX HERMITIAN 
    ! 
    !  DOUBLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  UPLO       (input)        CHARACTER: specifies whether the full part, or the upper or lower
    !                                       triangular part of the matrix(ces) is(are) being supplied.
    !  N          (input)        INTEGER: Size system
    !  sa         (input)        COMPLEX DOUBLE PRECISION (isa(N+1)-1):  Matrix A- CSR format 
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
    !  lambda     (output)       REAL DOUBLE PRECISION(M0)   : Eigenvalues -solution
    !  q          (input/output) COMPLEX DOUBLE PRECISION(N,M0) : 
    !                                                       On entry: subspace initial guess if fpm(5)=1 
    !                                                       On exit : Eigenvectors-solution
    !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
    !  res        (output)       REAL DOUBLE PRECISION(M0) : Relative Residual of the solution (1-norm)
    !                                                        if option fpm(6)=1 selected                           
    !  info       (output)       INTEGER: Error handling (0: successful exit)
    !=====================================================================
    ! Eric Polizzi 2009-2013
    ! ====================================================================


    implicit none
    include 'f90_noruntime_interface.fi'
    character(len=1) :: UPLO
    integer :: N
    complex(kind=kind(1.0d0)),dimension(*),target:: sa
    integer,dimension(*),target:: isa,jsa
    integer,dimension(*) :: fpm
    double precision :: epsout 
    integer :: loop
    double precision :: Emin,Emax
    integer :: M0
    double precision,dimension(*)  :: E
    complex(kind=kind(1.0d0)),dimension(N,*):: X
    integer :: mode
    double precision,dimension(*)    :: res
    integer :: info
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer :: ijob,infoloc
    complex(kind=(kind(1.0d0))) :: Ze
    complex(kind=(kind(1.0d0))),dimension(:),pointer :: saz,tsaz
    integer,dimension(:),pointer :: isaz,jsaz
    complex(kind=(kind(1.0d0))),dimension(:,:),pointer ::work1,work2,caux,zAq,zSq
    double precision,parameter :: DONE=1.0d0, DZERO=0.0d0
    complex(kind=(kind(1.0d0))),parameter :: ONEC=(DONE,DZERO), ZEROC=(DZERO,DZERO)
!!!!!for pardiso
    integer(8),dimension(64) :: pt1,pt2
    integer,dimension(64) :: iparm
    integer :: mtype
    integer :: MAXFCT,MNUM,PHASE,MSGLVL
    integer :: idum,nrhs
    double precision :: ddum1
!!!!! full csr format
    complex(kind=kind(1.0d0)),dimension(:),pointer :: ssa,ssb
    integer,dimension(:), pointer :: sisa,sjsa,sisb,sjsb
    integer :: i,k 
    integer :: opt,nnza,nnzb,nnz


!!!!!!!!!!!!!! Check INPUT PARAMETERS
    INFO = 0
    IF ((UPLO/='F').and.(UPLO/='f').and.(UPLO/='L').and.(UPLO/='l').and.(UPLO/='U').and.(UPLO/='u')) THEN
       INFO=-101
    ELSE IF ( N<=0 ) THEN
       INFO = -102
    END IF
    IF( INFO.NE.0 ) THEN
       CALL XERBLA( 'ZFEAST_HCSREV', -INFO+100 )
       RETURN
    END IF

    infoloc=0
    info=-1 ! by default
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! identity B matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!
    nnzb=n
    call wallocate_1i(sisb,n+1,infoloc)
    call wallocate_1z(ssb,nnzb,infoloc)
    call wallocate_1i(sjsb,nnzb,infoloc)
    if (infoloc/=0) then
       info=-1
       return
    end if

    do i=1,n
       sisb(i)=i
       sjsb(i)=i
       ssb(i)=ONEC
    enddo
    sisb(n+1)=nnzb+1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! FORMAT CONVERSION TO full CSR for PARDISO !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Remark: everything should be csr since we work with zS-H which is actually unsymmetric 


    if ((UPLO=='F').or.(UPLO=='f')) then !!! full csr already

       nnza=isa(n+1)-1
       ssa => sa(1:nnza)
       sisa => isa(1:n+1)
       sjsa => jsa(1:nnza)


    else !! upper-csr or lower-csr to full csr

       nnza=2*(isa(n+1)-1) ! may be overestimated
       call wallocate_1z(ssa,nnza,infoloc)
       if (infoloc/=0) return
       call wallocate_1i(sjsa,nnza,infoloc)
       if (infoloc/=0) return
       call wallocate_1i(sisa,n+1,infoloc)
       if (infoloc/=0) return
      
       call zhcsr_uplo_to_csr(N,sa,isa,jsa,ssa,sisa,sjsa)
      
    endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! FEAST INITIALIZATION 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call wallocate_2z(zAq,M0,M0,infoloc)
    if (infoloc/=0) return
    call wallocate_2z(zSq,M0,M0,infoloc)
    if (infoloc/=0) return
    call wallocate_2z(work1,N,M0,infoloc)
    if (infoloc/=0) return
    call wallocate_2z(work2,N,M0,infoloc)
    if (infoloc/=0) return
    call wallocate_1i(isaz,n+1,infoloc)
    if (infoloc/=0) return
    !!<<<
    call wallocate_1z(saz,1,infoloc) ! dummy
    if (infoloc/=0) return
    call wallocate_1i(jsaz,1,infoloc)! dummy
    if (infoloc/=0) return
    !!>>>
    opt=1
    call zaddcsr(N,opt,ONEC,ssa,sisa,sjsa,ONEC,ssb,sisb,sjsb,saz,isaz,jsaz) !! get isaz
    nnz=isaz(n+1)-1
    !!<<<
    call wdeallocate_1z(saz)
    call wdeallocate_1i(jsaz)
    !!>>>
    call wallocate_1z(saz,nnz,infoloc)
    if (infoloc/=0) return
    call wallocate_1i(jsaz,nnz,infoloc)
    if (infoloc/=0) return


    opt=2
    call zaddcsr(N,opt,ONEC,ssa,sisa,sjsa,ONEC,ssb,sisb,sjsb,saz,isaz,jsaz) !! get jsaz
!!!!!!!!!!!!!!!
    call wallocate_2z(caux,N,M0,infoloc)
    if (infoloc/=0) return


 MAXFCT=1 ! One factorization to consider (default)

if (fpm(11)==0) then!!!  transpose option not possible 
     MAXFCT=2 ! two factorization to consider normal+transpose
    call wallocate_1z(tsaz,nnz,infoloc)
    if (infoloc/=0) return
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  analysis step (symbolic factorizarion for pardiso)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
    MTYPE=3 ! complex and structurally symmetric 
    call pardisoinit(PT1,MTYPE,IPARM)
if (fpm(11)==0) call pardisoinit(PT2,MTYPE,IPARM)

!!!!!!!!!!!!
if (fpm(64)==1) then
do i=1,64
if (fpm(64+i)/=-9999) iparm(i)=fpm(64+i)
enddo
endif
!!!!!!!!!!!!


    !IPARM(3)=fpm(10) !omp_num_threads !! openmp -number of threads
    IPARM(6)=1 ! solution and rhs are input/output, attention caux is always used
    MSGLVL=0 !0 !0- no output, 1- output
    PHASE=11
    MNUM=1 
    call PARDISO(PT1,MAXFCT,MNUM,MTYPE,PHASE,n,saz,isaz,jsaz,idum,M0,IPARM,MSGLVL,work2(:,1:M0),caux(:,1:M0),infoloc)
if (fpm(11)==0) then 
    MNUM=2 
    call PARDISO(PT2,MAXFCT,MNUM,MTYPE,PHASE,n,tsaz,isaz,jsaz,idum,M0,IPARM,MSGLVL,work2(:,1:M0),caux(:,1:M0),infoloc)
endif
    if (infoloc/=0) then
       info=-2
       return
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!! FEAST !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

    ijob=-1 ! initialization

    do while (ijob/=0)

       call zfeast_hrci(ijob,N,Ze,work1,work2,zAq,zSq,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
       select case(ijob)

       case(10) !! Factorize (ZeB-A)
          opt=3
          call zaddcsr(N,opt,-ONEC,ssa,sisa,sjsa,Ze,ssb,sisb,sjsb,saz,isaz,jsaz) !! get saz

          PHASE=22
          MNUM=1
          call PARDISO(PT1,MAXFCT,MNUM,MTYPE,PHASE,n,saz,isaz,jsaz,idum,M0,IPARM,MSGLVL,work2(:,1:M0),caux(:,1:M0),infoloc)
          if (infoloc/=0) info=-2

       case(11) !! Solve (ZeB-A)x=work2(1:N,1:M0) result in to work2
          IPARM(12)=0 ! normal solve
          PHASE=33
          MNUM=1
          call PARDISO(PT1,MAXFCT,MNUM,MTYPE,PHASE,n,saz,isaz,jsaz,idum,M0,IPARM,MSGLVL,work2(:,1:M0),caux(:,1:M0),infoloc)
          if (infoloc/=0) info=-2

       case(20) !! Factorize (ZeB-A)^T   
if (fpm(11)==0) then
          opt=3
          call zaddcsr(N,opt,-ONEC,ssa,sisa,sjsa,conjg(Ze),ssb,sisb,sjsb,tsaz,isaz,jsaz) !! get saz

          PHASE=22
          MNUM=2
          call PARDISO(PT2,MAXFCT,MNUM,MTYPE,PHASE,n,tsaz,isaz,jsaz,idum,M0,IPARM,MSGLVL,work2(:,1:M0),caux(:,1:M0),infoloc)
          if (infoloc/=0) info=-2
endif
       case(21) !!  Solve  (ZeB-A)^Tx=work2(1:N,1:M0) result in to work2
if (fpm(11)==0) then
          IPARM(12)=0 ! normal solve
          PHASE=33
          MNUM=2
          call PARDISO(PT2,MAXFCT,MNUM,MTYPE,PHASE,n,tsaz,isaz,jsaz,idum,M0,IPARM,MSGLVL,work2(:,1:M0),caux(:,1:M0),infoloc)
          if (infoloc/=0) info=-2
else
IPARM(12)=1 ! transpose conjugate solve
  PHASE=33
          MNUM=1
          call PARDISO(PT1,MAXFCT,MNUM,MTYPE,PHASE,n,saz,isaz,jsaz,idum,M0,IPARM,MSGLVL,work2(:,1:M0),caux(:,1:M0),infoloc)
          if (infoloc/=0) info=-2
endif

       case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

          call zhcsrmm('F',N,N,fpm(25),ONEC,ssa,sisa,sjsa,X(1,fpm(24)),ZEROC,work1(1,fpm(24)))

       case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

          call ZLACPY( 'F', N, fpm(25),X(1,fpm(24)) , N, work1(1,fpm(24)), N )

       end select
    end do
!!!!!!!!!!!!!!!!!!!!!!!!
!!!! release memory
!!!!!!!!!!!!!!!!!!!!!!!!
    PHASE=-1 

    MNUM=1
    call PARDISO(PT1,MAXFCT,MNUM,MTYPE,PHASE,n,saz,isaz,jsaz,idum,M0,IPARM,MSGLVL,work2,caux,infoloc)
    if (infoloc/=0) then
       info=-2
       return
    end if

if (fpm(11)==0) then
    MNUM=2
    call PARDISO(PT2,MAXFCT,MNUM,MTYPE,PHASE,n,tsaz,isaz,jsaz,idum,M0,IPARM,MSGLVL,work2,caux,infoloc)
    if (infoloc/=0) then
       info=-2
       return
    end if
 end if

    call wdeallocate_2z(zAq)
    call wdeallocate_2z(zSq)
    call wdeallocate_2z(work1)
    call wdeallocate_2z(work2)
    call wdeallocate_2z(caux)

    call wdeallocate_1z(saz)
    call wdeallocate_1i(isaz)
    call wdeallocate_1i(jsaz)


if (fpm(11)==0)  call wdeallocate_1z(tsaz) !! transpose option


    if ((UPLO/='F').and.(UPLO/='f')) then
       call wdeallocate_1z(ssa)
       call wdeallocate_1i(sisa)
       call wdeallocate_1i(sjsa)
    endif


       call wdeallocate_1z(ssb)
       call wdeallocate_1i(sisb)
       call wdeallocate_1i(sjsb)

  end subroutine zfeast_hcsrev


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  
  subroutine sfeast_scsrgv(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
    !  Purpose 
    !  =======
    !  FEAST SPARSE INTERFACE (SPARSE CSR FORMAT)
    !  Solve the generalized Ax=eBx eigenvalue problem
    !  
    !  A REAL SYMMETRIC, B SYMMETRIC POSITIVE DEFINITE 
    ! 
    !  SINGLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  UPLO       (input)        CHARACTER: specifies whether the full part 'F', or the upper 'U' or lower 'L'
    !                                       triangular part of the matrix(ces) is(are) being supplied.
    !  N          (input)        INTEGER: Size system
    !  sa         (input)        REAL SINGLE PRECISION (isa(N+1)-1):  Matrix A- CSR format 
    !  isa        (input)        INTEGER(N+1): CSR row array of Matrix A   
    !  jsa        (input)        INTEGER(isa(N+1)-1): CSR column array of Matrix A
    !  sb         (input)        REAL SINGLE PRECISION (isb(N+1)-1):  Matrix B- CSR format 
    !  isb        (input)        INTEGER(N+1): CSR row array of Matrix B   
    !  jsb        (input)        INTEGER(isb(N+1)-1): CSR column array of Matrix B
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  fpm (input/output)     INTEGER(*) : FEAST parameters (see FEAST documentation)
    !  epsout     (output)       REAL SINGLE PRECISION : Error on the trace
    !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
    !  Emin,Emax  (input)        REAL SINGLE PRECISION: search interval
    !  M0         (input/output) INTEGER: Size subspace
    !  lambda     (output)       REAL SINGLE PRECISION(M0)   : Eigenvalues -solution
    !  q          (input/output) REAL SINGLE PRECISION(N,M0) : 
    !                                                       On entry: subspace initial guess if fpm(5)=1 
    !                                                       On exit : Eigenvectors-solution
    !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
    !  res        (output)       REAL SINGLE PRECISION(M0) : Relative Residual of the solution (1-norm)
    !                                                        if option fpm(6)=1 selected                           
    !  info       (output)       INTEGER: Error handling (0: successful exit) -- see FEAST documentation
    !=====================================================================
    ! Eric Polizzi 2009-2013
    ! ====================================================================

    implicit none
    include 'f90_noruntime_interface.fi'
    character(len=1) :: UPLO
    integer :: N
    real,dimension(*),target:: sa,sb
    integer,dimension(*),target:: isa,jsa,isb,jsb
    integer,dimension(*) :: fpm
    real :: epsout 
    integer :: loop
    real :: Emin,Emax
    integer :: M0
    real,dimension(*)  :: E
    real,dimension(N,*):: X
    integer :: mode
    real,dimension(*)    :: res
    integer :: info
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer :: ijob,infoloc
    complex :: Ze
    complex,dimension(:),pointer :: saz
    integer,dimension(:),pointer :: isaz,jsaz
    complex,dimension(:,:),pointer ::workc,caux
    real, dimension(:,:),pointer ::work,Aq,Sq
    real,parameter :: SONE=1.0e0, SZERO=0.0e0
    complex,parameter :: ONEC=(SONE,SZERO)
!!!!!for pardiso
    integer(8),dimension(64) :: pt
    integer,dimension(64) :: iparm
    integer :: mtype
    integer :: MAXFCT,MNUM,PHASE,MSGLVL
    integer :: idum,nrhs
    real :: ddum1
!!!!! csr-upper format
    real,dimension(:),pointer :: ssa,ssb
    integer,dimension(:), pointer :: sisa,sjsa,sisb,sjsb
    integer :: i,k 
    integer :: opt,nnza,nnzb,nnz


!!!!!!!!!!!!!! Check INPUT PARAMETERS
    INFO = 0
    IF ((UPLO/='F').and.(UPLO/='f').and.(UPLO/='L').and.(UPLO/='l').and.(UPLO/='U').and.(UPLO/='u')) THEN
       INFO=-101
    ELSE IF ( N<=0 ) THEN
       INFO = -102
    END IF
    IF( INFO.NE.0 ) THEN
       CALL XERBLA( 'SFEAST_SCSRGV', -INFO+100 )
       RETURN
    END IF

    info=-1 ! by default

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! FORMAT CONVERSION TO CSR-UPPER for PARDISO !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if ((UPLO=='F').or.(UPLO=='f')) then !!! full csr to upper-csr
       call wallocate_1i(sisa,n+1,infoloc)
if (infoloc/=0) return
       call wallocate_1i(sisb,n+1,infoloc)
if (infoloc/=0) return
 
!!<<<
       call wallocate_1s(ssa,1,infoloc) ! dummy
if (infoloc/=0) return
       call wallocate_1i(sjsa,1,infoloc) !dummy
if (infoloc/=0) return
       call wallocate_1s(ssb,1,infoloc) !dummy
if (infoloc/=0) return
       call wallocate_1i(sjsb,1,infoloc) !dummy
!!>>>
       opt=1
       call scsr2csr_up(opt,N,sa,isa,jsa,ssa,sisa,sjsa)
       nnza=sisa(n+1)-1 
       call scsr2csr_up(opt,N,sb,isb,jsb,ssb,sisb,sjsb)
       nnzb=sisb(n+1)-1 
!!<<<
       call wdeallocate_1s(ssa)
       call wdeallocate_1i(sjsa)
       call wdeallocate_1s(ssb)
       call wdeallocate_1i(sjsb)
!!>>>

       call wallocate_1s(ssa,nnza,infoloc)
if (infoloc/=0) return
       call wallocate_1i(sjsa,nnza,infoloc)
if (infoloc/=0) return
       call wallocate_1s(ssb,nnzb,infoloc)
if (infoloc/=0) return
       call wallocate_1i(sjsb,nnzb,infoloc)
if (infoloc/=0) return

       opt=2
       call scsr2csr_up(opt,N,sa,isa,jsa,ssa,sisa,sjsa)
       call scsr2csr_up(opt,N,sb,isb,jsb,ssb,sisb,sjsb)

!!!!!!!!!!!!!!!!!!!!!!!!!
    elseif ((UPLO=='U').or.(UPLO=='u')) then !! upper-csr already 

       nnza=isa(n+1)-1
       ssa => sa(1:nnza)
       sisa => isa(1:n+1)
       sjsa => jsa(1:nnza)

       nnzb=isb(n+1)-1
       ssb =>  sb(1:nnzb)
       sisb => isb(1:n+1)
       sjsb => jsb(1:nnzb)



    elseif ((UPLO=='L').or.(UPLO=='l')) then !!! lower-csr to upper-csr

       nnza=isa(n+1)-1
       nnzb=isb(n+1)-1
       call wallocate_1s(ssa,nnza,infoloc)
if (infoloc/=0) return
       call wallocate_1i(sjsa,nnza,infoloc)
if (infoloc/=0) return
       call wallocate_1s(ssb,nnzb,infoloc)
if (infoloc/=0) return
       call wallocate_1i(sjsb,nnzb,infoloc)
if (infoloc/=0) return
       call wallocate_1i(sisa,n+1,infoloc)
if (infoloc/=0) return
       call wallocate_1i(sisb,n+1,infoloc)
if (infoloc/=0) return

       call scsr_transpose(N,sa,isa,jsa,ssa,sisa,sjsa)
       call scsr_transpose(N,sb,isb,jsb,ssb,sisb,sjsb)

    endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! FEAST INITIALIZATION 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call wallocate_2s(Aq,M0,M0,infoloc)
    if (infoloc/=0) return
    call wallocate_2s(Sq,M0,M0,infoloc)
    if (infoloc/=0) return
    call wallocate_2s(work,N,M0,infoloc)
    if (infoloc/=0) return
    call wallocate_2c(workc,N,M0,infoloc)
    if (infoloc/=0) return
    call wallocate_1i(isaz,n+1,infoloc)
    if (infoloc/=0) return
    !!<<<
    call wallocate_1c(saz,1,infoloc) ! dummy
    if (infoloc/=0) return
    call wallocate_1i(jsaz,1,infoloc)! dummy
    if (infoloc/=0) return
    !!>>>
    opt=1
    call csaddcsr(N,opt,ONEC,ssa,sisa,sjsa,ONEC,ssb,sisb,sjsb,saz,isaz,jsaz) !! get isaz
    nnz=isaz(n+1)-1
    !!<<<
    call wdeallocate_1c(saz)
    call wdeallocate_1i(jsaz)
    !!>>>
    call wallocate_1c(saz,nnz,infoloc)
    if (infoloc/=0) return
    call wallocate_1i(jsaz,nnz,infoloc)
    if (infoloc/=0) return

    opt=2
    call csaddcsr(N,opt,ONEC,ssa,sisa,sjsa,ONEC,ssb,sisb,sjsb,saz,isaz,jsaz) !! get jsaz
!!!!!!!!!!!!!!!
    call wallocate_2c(caux,N,M0,infoloc)
    if (infoloc/=0) return

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  analysis step (symbolic factorizarion for pardiso)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    MAXFCT=1 ! one factorization to consider normal+transpose
    MTYPE=6  ! complex and symmetric 
    call pardisoinit(PT,MTYPE,IPARM)
!!!!!!!!!!!!
if (fpm(64)==1) then
do i=1,64
if (fpm(64+i)/=-9999) iparm(i)=fpm(64+i)
enddo
endif
!!!!!!!!!!!!

    !IPARM(3)=fpm(10) !omp_num_threads !! openmp -number of threads
    IPARM(6)=1 ! solution and rhs are input/output, attention caux is always used
    MSGLVL=0 !0- no output, 1- output
    PHASE=11
!!!!!!!! single precision pardiso (MKL only)
             IPARM(28)=1 ! pardiso single precision
!!!!!!!!!!!!!!!
    MNUM=1
    call PARDISO(PT,MAXFCT,MNUM,MTYPE,PHASE,n,saz,isaz,jsaz,idum,M0,IPARM,MSGLVL,workc(:,1:M0),caux(:,1:M0),infoloc) 
    if (infoloc/=0) then
       info=-2
       return
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!! FEAST !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

    ijob=-1 ! initialization

    do while (ijob/=0)

       call sfeast_srci(ijob,N,Ze,work,workc,Aq,Sq,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
       select case(ijob)

       case(10) !! Factorize (ZeB-A)
          opt=3
          call csaddcsr(N,opt,-ONEC,ssa,sisa,sjsa,Ze,ssb,sisb,sjsb,saz,isaz,jsaz) !! get saz

          PHASE=22
          call PARDISO(PT,MAXFCT,MNUM,MTYPE,PHASE,n,saz,isaz,jsaz,idum,M0,IPARM,MSGLVL,workc(:,1:M0),caux(:,1:M0),infoloc)
          if (infoloc/=0) info=-2

       case(11) !! Solve (ZeB-A)x=work2(1:N,1:M0) result in to work2
          PHASE=33
          call PARDISO(PT,MAXFCT,MNUM,MTYPE,PHASE,n,saz,isaz,jsaz,idum,M0,IPARM,MSGLVL,workc(:,1:M0),caux(:,1:M0),infoloc)
          if (infoloc/=0) info=-2

       case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

          call scsrmm('U',N,N,fpm(25),SONE,ssa,sisa,sjsa,X(1,fpm(24)),SZERO,work(1,fpm(24)))

       case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

          call scsrmm('U',N,N,fpm(25),SONE,ssb,sisb,sjsb,X(1,fpm(24)),SZERO,work(1,fpm(24)))

       end select
    end do
!!!!!!!!!!!!!!!!!!!!!!!!
!!!! release memory
!!!!!!!!!!!!!!!!!!!!!!!!
    PHASE=-1 

    call PARDISO(PT,MAXFCT,MNUM,MTYPE,PHASE,n,saz,isaz,jsaz,idum,M0,IPARM,MSGLVL,workc,caux,infoloc)
    if (infoloc/=0) then
       info=-2
       return
    end if

 

    call wdeallocate_2s(Aq)
    call wdeallocate_2s(Sq)
    call wdeallocate_2s(work)
    call wdeallocate_2c(workc)
    call wdeallocate_2c(caux)

    call wdeallocate_1c(saz)
    call wdeallocate_1i(isaz)
    call wdeallocate_1i(jsaz)


    if ((UPLO/='U').and.(UPLO/='u')) then
       call wdeallocate_1s(ssa)
       call wdeallocate_1i(sisa)
       call wdeallocate_1i(sjsa)

       call wdeallocate_1s(ssb)
       call wdeallocate_1i(sisb)
       call wdeallocate_1i(sjsb)
    endif

  end subroutine sfeast_scsrgv



subroutine sfeast_scsrev(UPLO,N,sa,isa,jsa,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
    !  Purpose 
    !  =======
    !  FEAST SPARSE INTERFACE (SPARSE CSR FORMAT)
    !  Solve the standard Ax=ex eigenvalue problem
    !  
    !  A REAL SYMMETRIC
    ! 
    !  SINGLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  UPLO       (input)        CHARACTER: specifies whether the full part 'F', or the upper 'U' or lower 'L'
    !                                       triangular part of the matrix(ces) is(are) being supplied.
    !  N          (input)        INTEGER: Size system
    !  sa         (input)        REAL SINGLE PRECISION (isa(N+1)-1):  Matrix A- CSR format 
    !  isa        (input)        INTEGER(N+1): CSR row array of Matrix A   
    !  jsa        (input)        INTEGER(isa(N+1)-1): CSR column array of Matrix A
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  fpm (input/output)     INTEGER(*) : FEAST parameters (see FEAST documentation)
    !  epsout     (output)       REAL SINGLE PRECISION : Error on the trace
    !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
    !  Emin,Emax  (input)        REAL SINGLE PRECISION: search interval
    !  M0         (input/output) INTEGER: Size subspace
    !  lambda     (output)       REAL SINGLE PRECISION(M0)   : Eigenvalues -solution
    !  q          (input/output) REAL SINGLE PRECISION(N,M0) : 
    !                                                       On entry: subspace initial guess if fpm(5)=1 
    !                                                       On exit : Eigenvectors-solution
    !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
    !  res        (output)       REAL SINGLE PRECISION(M0) : Relative Residual of the solution (1-norm)
    !                                                        if option fpm(6)=1 selected                           
    !  info       (output)       INTEGER: Error handling (0: successful exit) -- see FEAST documentation
    !=====================================================================
    ! Eric Polizzi 2009-2013
    ! ====================================================================


    implicit none
    include 'f90_noruntime_interface.fi'
    character(len=1) :: UPLO
    integer :: N
    real,dimension(*),target:: sa
    integer,dimension(*),target:: isa,jsa
    integer,dimension(*) :: fpm
    real :: epsout 
    integer :: loop
    real :: Emin,Emax
    integer :: M0
    real,dimension(*)  :: E
    real,dimension(N,*):: X
    integer :: mode
    real,dimension(*)    :: res
    integer :: info
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer :: ijob,infoloc
    complex :: Ze
    complex,dimension(:),pointer :: saz
    integer,dimension(:),pointer :: isaz,jsaz
    complex,dimension(:,:),pointer ::workc,caux
    real, dimension(:,:),pointer ::work,Aq,Sq
    real,parameter :: SONE=1.0e0, SZERO=0.0e0
    complex,parameter :: ONEC=(SONE,SZERO)
!!!!!for pardiso
    integer(8),dimension(64) :: pt
    integer,dimension(64) :: iparm
    integer :: mtype
    integer :: MAXFCT,MNUM,PHASE,MSGLVL
    integer :: idum,nrhs
    real :: ddum1
!!!!! csr-upper format
    real,dimension(:),pointer :: ssa,ssb
    integer,dimension(:), pointer :: sisa,sjsa,sisb,sjsb
    integer :: i,k 
    integer :: opt,nnza,nnzb,nnz


!!!!!!!!!!!!!! Check INPUT PARAMETERS
    INFO = 0
    IF ((UPLO/='F').and.(UPLO/='f').and.(UPLO/='L').and.(UPLO/='l').and.(UPLO/='U').and.(UPLO/='u')) THEN
       INFO=-101
    ELSE IF ( N<=0 ) THEN
       INFO = -102
    END IF
    IF( INFO.NE.0 ) THEN
       CALL XERBLA( 'SFEAST_SCSREV', -INFO+100 )
       RETURN
    END IF

    info=-1 ! by default
 infoloc=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! identity B matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!
    nnzb=n
    call wallocate_1i(sisb,n+1,infoloc)
    call wallocate_1s(ssb,nnzb,infoloc)
    call wallocate_1i(sjsb,nnzb,infoloc)
    if (infoloc/=0) then
       info=-1
       return
    end if

    do i=1,n
       sisb(i)=i
       sjsb(i)=i
       ssb(i)=SONE
    enddo
    sisb(n+1)=nnzb+1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! FORMAT CONVERSION TO CSR-UPPER for PARDISO !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if ((UPLO=='F').or.(UPLO=='f')) then !!! full csr to upper-csr
       call wallocate_1i(sisa,n+1,infoloc)
if (infoloc/=0) return
     
!!<<<
       call wallocate_1s(ssa,1,infoloc) ! dummy
if (infoloc/=0) return
       call wallocate_1i(sjsa,1,infoloc) !dummy
if (infoloc/=0) return
!!>>>
       opt=1
       call scsr2csr_up(opt,N,sa,isa,jsa,ssa,sisa,sjsa)
       nnza=sisa(n+1)-1 
!!<<<
       call wdeallocate_1s(ssa)
       call wdeallocate_1i(sjsa)
!!>>>

       call wallocate_1s(ssa,nnza,infoloc)
if (infoloc/=0) return
       call wallocate_1i(sjsa,nnza,infoloc)
if (infoloc/=0) return

       opt=2
       call scsr2csr_up(opt,N,sa,isa,jsa,ssa,sisa,sjsa)
      
!!!!!!!!!!!!!!!!!!!!!!!!!
    elseif ((UPLO=='U').or.(UPLO=='u')) then !! upper-csr already 

       nnza=isa(n+1)-1
       ssa => sa(1:nnza)
       sisa => isa(1:n+1)
       sjsa => jsa(1:nnza)

    elseif ((UPLO=='L').or.(UPLO=='l')) then !!! lower-csr to upper-csr

       nnza=isa(n+1)-1
       call wallocate_1s(ssa,nnza,infoloc)
if (infoloc/=0) return
       call wallocate_1i(sjsa,nnza,infoloc)
if (infoloc/=0) return
       call wallocate_1i(sisa,n+1,infoloc)
if (infoloc/=0) return
    
       call scsr_transpose(N,sa,isa,jsa,ssa,sisa,sjsa)

    endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! FEAST INITIALIZATION 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call wallocate_2s(Aq,M0,M0,infoloc)
    if (infoloc/=0) return
    call wallocate_2s(Sq,M0,M0,infoloc)
    if (infoloc/=0) return
    call wallocate_2s(work,N,M0,infoloc)
    if (infoloc/=0) return
    call wallocate_2c(workc,N,M0,infoloc)
    if (infoloc/=0) return
    call wallocate_1i(isaz,n+1,infoloc)
    if (infoloc/=0) return
    !!<<<
    call wallocate_1c(saz,1,infoloc) ! dummy
    if (infoloc/=0) return
    call wallocate_1i(jsaz,1,infoloc)! dummy
    if (infoloc/=0) return
    !!>>>
    opt=1
    call csaddcsr(N,opt,ONEC,ssa,sisa,sjsa,ONEC,ssb,sisb,sjsb,saz,isaz,jsaz) !! get isaz
    nnz=isaz(n+1)-1
    !!<<<
    call wdeallocate_1c(saz)
    call wdeallocate_1i(jsaz)
    !!>>>
    call wallocate_1c(saz,nnz,infoloc)
    if (infoloc/=0) return
    call wallocate_1i(jsaz,nnz,infoloc)
    if (infoloc/=0) return

    opt=2
    call csaddcsr(N,opt,ONEC,ssa,sisa,sjsa,ONEC,ssb,sisb,sjsb,saz,isaz,jsaz) !! get jsaz
!!!!!!!!!!!!!!!
    call wallocate_2c(caux,N,M0,infoloc)
    if (infoloc/=0) return

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  analysis step (symbolic factorizarion for pardiso)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    MAXFCT=1 ! one factorization to consider normal+transpose
    MTYPE=6  ! complex and symmetric 
    call pardisoinit(PT,MTYPE,IPARM)
!!!!!!!!!!!!
if (fpm(64)==1) then
do i=1,64
if (fpm(64+i)/=-9999) iparm(i)=fpm(64+i)
enddo
endif
!!!!!!!!!!!!

    !IPARM(3)=fpm(10) !omp_num_threads !! openmp -number of threads
    IPARM(6)=1 ! solution and rhs are input/output, attention caux is always used
    MSGLVL=0 !0- no output, 1- output
    PHASE=11
!!!!!!!! single precision pardiso (MKL only)
             IPARM(28)=1 ! pardiso single precision
!!!!!!!!!!!!!!!
    MNUM=1
    call PARDISO(PT,MAXFCT,MNUM,MTYPE,PHASE,n,saz,isaz,jsaz,idum,M0,IPARM,MSGLVL,workc(:,1:M0),caux(:,1:M0),infoloc) 
    if (infoloc/=0) then
       info=-2
       return
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!! FEAST !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

    ijob=-1 ! initialization

    do while (ijob/=0)

       call sfeast_srci(ijob,N,Ze,work,workc,Aq,Sq,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
       select case(ijob)

       case(10) !! Factorize (ZeB-A)
          opt=3
          call csaddcsr(N,opt,-ONEC,ssa,sisa,sjsa,Ze,ssb,sisb,sjsb,saz,isaz,jsaz) !! get saz

          PHASE=22
          call PARDISO(PT,MAXFCT,MNUM,MTYPE,PHASE,n,saz,isaz,jsaz,idum,M0,IPARM,MSGLVL,workc(:,1:M0),caux(:,1:M0),infoloc)
          if (infoloc/=0) info=-2

       case(11) !! Solve (ZeB-A)x=work2(1:N,1:M0) result in to work2
          PHASE=33
          call PARDISO(PT,MAXFCT,MNUM,MTYPE,PHASE,n,saz,isaz,jsaz,idum,M0,IPARM,MSGLVL,workc(:,1:M0),caux(:,1:M0),infoloc)
          if (infoloc/=0) info=-2

       case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

          call scsrmm('U',N,N,fpm(25),SONE,ssa,sisa,sjsa,X(1,fpm(24)),SZERO,work(1,fpm(24)))

       case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

          call SLACPY( 'F', N, fpm(25),X(1,fpm(24)) , N, work(1,fpm(24)), N )

       end select
    end do
!!!!!!!!!!!!!!!!!!!!!!!!
!!!! release memory
!!!!!!!!!!!!!!!!!!!!!!!!
    PHASE=-1 

    call PARDISO(PT,MAXFCT,MNUM,MTYPE,PHASE,n,saz,isaz,jsaz,idum,M0,IPARM,MSGLVL,workc,caux,infoloc)
    if (infoloc/=0) then
       info=-2
       return
    end if

 

    call wdeallocate_2s(Aq)
    call wdeallocate_2s(Sq)
    call wdeallocate_2s(work)
    call wdeallocate_2c(workc)
    call wdeallocate_2c(caux)

    call wdeallocate_1c(saz)
    call wdeallocate_1i(isaz)
    call wdeallocate_1i(jsaz)


    if ((UPLO/='U').and.(UPLO/='u')) then
       call wdeallocate_1s(ssa)
       call wdeallocate_1i(sisa)
       call wdeallocate_1i(sjsa)
endif
       call wdeallocate_1s(ssb)
       call wdeallocate_1i(sisb)
       call wdeallocate_1i(sjsb)
    

  end subroutine sfeast_scsrev

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  
  subroutine cfeast_hcsrgv(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
    !  Purpose 
    !  =======
    !  FEAST SPARSE INTERFACE (SPARSE CSR FORMAT)
    !  Solve the generalized Ax=eBx eigenvalue problem
    !  
    !  A COMPLEX HERMITIAN, B HERMITIAN POSITIVE DEFINITE 
    ! 
    !  SINGLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  UPLO       (input)        CHARACTER: specifies whether the full part 'F', or the upper 'U' or lower 'L'
    !                                       triangular part of the matrix(ces) is(are) being supplied.
    !  N          (input)        INTEGER: Size system
    !  sa         (input)        COMPLEX SINGLE PRECISION (isa(N+1)-1):  Matrix A- CSR format 
    !  isa        (input)        INTEGER(N+1): CSR row array of Matrix A   
    !  jsa        (input)        INTEGER(isa(N+1)-1): CSR column array of Matrix A
    !  sb         (input)        COMPLEX SINGLE PRECISION (isb(N+1)-1):  Matrix B- CSR format 
    !  isb        (input)        INTEGER(N+1): CSR row array of Matrix B   
    !  jsb        (input)        INTEGER(isb(N+1)-1): CSR column array of Matrix B
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  fpm (input/output) INTEGER(*) : FEAST parameters (see FEAST documentation)
    !  epsout     (output)       REAL SINGLE PRECISION : Error on the trace
    !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
    !  Emin,Emax  (input)        REAL SINGLE PRECISION: search interval
    !  M0         (input/output) INTEGER: Size subspace
    !  lambda     (output)       REAL SINGLE PRECISION(M0)   : Eigenvalues -solution
    !  q          (input/output) COMPLEX SINGLE PRECISION(N,M0) : 
    !                                                       On entry: subspace initial guess if fpm(5)=1 
    !                                                       On exit : Eigenvectors-solution
    !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
    !  res        (output)       REAL SINGLE PRECISION(M0) : Relative Residual of the solution (1-norm)
    !                                                        if option fpm(6)=1 selected                           
    !  info       (output)       INTEGER: Error handling (0: successful exit) -- see FEAST documentation
    !=====================================================================
    ! Eric Polizzi 2009-2013
    ! ====================================================================

    implicit none
    include 'f90_noruntime_interface.fi'
    character(len=1) :: UPLO
    integer :: N
    complex,dimension(*),target:: sa,sb
    integer,dimension(*),target:: isa,jsa,isb,jsb
    integer,dimension(*) :: fpm
    real :: epsout 
    integer :: loop
    real :: Emin,Emax
    integer :: M0
    real,dimension(*)  :: E
    complex,dimension(N,*):: X
    integer :: mode
    real,dimension(*)    :: res
    integer :: info
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer :: ijob,infoloc
    complex :: Ze
    complex,dimension(:),pointer :: saz,tsaz
    integer,dimension(:),pointer :: isaz,jsaz
    complex,dimension(:,:),pointer ::work1,work2,caux,zAq,zSq
    real,parameter :: SONE=1.0e0, SZERO=0.0e0
    complex,parameter :: ONEC=(SONE,SZERO), ZEROC=(SZERO,SZERO)
!!!!!for pardiso
    integer(8),dimension(64) :: pt1,pt2
    integer,dimension(64) :: iparm
    integer :: mtype
    integer :: MAXFCT,MNUM,PHASE,MSGLVL
    integer :: idum,nrhs
    real :: ddum1
!!!!! full csr format
    complex,dimension(:),pointer :: ssa,ssb
    integer,dimension(:), pointer :: sisa,sjsa,sisb,sjsb
    integer :: i,k 
    integer :: opt,nnza,nnzb,nnz


!!!!!!!!!!!!!! Check INPUT PARAMETERS
    INFO = 0
    IF ((UPLO/='F').and.(UPLO/='f').and.(UPLO/='L').and.(UPLO/='l').and.(UPLO/='U').and.(UPLO/='u')) THEN
       INFO=-101
    ELSE IF ( N<=0 ) THEN
       INFO = -102
    END IF
    IF( INFO.NE.0 ) THEN
       CALL XERBLA( 'CFEAST_HCSRGV', -INFO+100 )
       RETURN
    END IF

    info=-1 ! by default
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! FORMAT CONVERSION TO full CSR for PARDISO !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Remark: everything should be csr since we work with zS-H which is actually unsymmetric 


    if ((UPLO=='F').or.(UPLO=='f')) then !!! full csr already

       nnza=isa(n+1)-1
       ssa => sa(1:nnza)
       sisa => isa(1:n+1)
       sjsa => jsa(1:nnza)

       nnzb=isb(n+1)-1
       ssb => sb(1:nnza)
       sisb => isb(1:n+1)
       sjsb => jsb(1:nnza)

    else !! upper-csr or lower-csr to full csr

       nnza=2*(isa(n+1)-1) ! may be overestimated
       nnzb=2*(isb(n+1)-1)
       call wallocate_1c(ssa,nnza,infoloc)
       if (infoloc/=0) return
       call wallocate_1i(sjsa,nnza,infoloc)
       if (infoloc/=0) return
       call wallocate_1c(ssb,nnzb,infoloc)
       if (infoloc/=0) return
       call wallocate_1i(sjsb,nnzb,infoloc)
       if (infoloc/=0) return
       call wallocate_1i(sisa,n+1,infoloc)
       if (infoloc/=0) return
       call wallocate_1i(sisb,n+1,infoloc)
       if (infoloc/=0) return

       call chcsr_uplo_to_csr(N,sa,isa,jsa,ssa,sisa,sjsa)
       call chcsr_uplo_to_csr(N,sb,isb,jsb,ssb,sisb,sjsb)

    endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! FEAST INITIALIZATION 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call wallocate_2c(zAq,M0,M0,infoloc)
    if (infoloc/=0) return
    call wallocate_2c(zSq,M0,M0,infoloc)
    if (infoloc/=0) return
    call wallocate_2c(work1,N,M0,infoloc)
    if (infoloc/=0) return
    call wallocate_2c(work2,N,M0,infoloc)
    if (infoloc/=0) return
    call wallocate_1i(isaz,n+1,infoloc)
    if (infoloc/=0) return
    !!<<<
    call wallocate_1c(saz,1,infoloc) ! dummy
    if (infoloc/=0) return
    call wallocate_1i(jsaz,1,infoloc)! dummy
    if (infoloc/=0) return
    !!>>>
    opt=1
    call caddcsr(N,opt,ONEC,ssa,sisa,sjsa,ONEC,ssb,sisb,sjsb,saz,isaz,jsaz) !! get isaz
    nnz=isaz(n+1)-1
    !!<<<
    call wdeallocate_1c(saz)
    call wdeallocate_1i(jsaz)
    !!>>>
    call wallocate_1c(saz,nnz,infoloc)
    if (infoloc/=0) return
    call wallocate_1i(jsaz,nnz,infoloc)
    if (infoloc/=0) return


    opt=2
    call caddcsr(N,opt,ONEC,ssa,sisa,sjsa,ONEC,ssb,sisb,sjsb,saz,isaz,jsaz) !! get jsaz
!!!!!!!!!!!!!!!
    call wallocate_2c(caux,N,M0,infoloc)
    if (infoloc/=0) return




 MAXFCT=1 ! One factorization to consider (default)

if (fpm(11)==0) then!!  needs the transpose option
    MAXFCT=2 ! two factorization to consider normal+transpose
    call wallocate_1c(tsaz,nnz,infoloc)
    if (infoloc/=0) return
 end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  analysis step (symbolic factorizarion for pardiso)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    MTYPE=3 ! complex and structurally symmetric 
    call pardisoinit(PT1,MTYPE,IPARM)
if (fpm(11)==0)  call pardisoinit(PT2,MTYPE,IPARM)

!!!!!!!!!!!!
if (fpm(64)==1) then
do i=1,64
if (fpm(64+i)/=-9999) iparm(i)=fpm(64+i)
enddo
endif
!!!!!!!!!!!!

    !IPARM(3)=fpm(10) !omp_num_threads !! openmp -number of threads
    IPARM(6)=1 ! solution and rhs are input/output, attention caux is always used
    MSGLVL=0 !0- no output, 1- output
    PHASE=11
    IPARM(28)=1 ! pardiso single precision (MKL)
    MNUM=1 
    call PARDISO(PT1,MAXFCT,MNUM,MTYPE,PHASE,n,saz,isaz,jsaz,idum,M0,IPARM,MSGLVL,work2(:,1:M0),caux(:,1:M0),infoloc) 
if (fpm(11)==0) then
    MNUM=2 
    call PARDISO(PT2,MAXFCT,MNUM,MTYPE,PHASE,n,tsaz,isaz,jsaz,idum,M0,IPARM,MSGLVL,work2(:,1:M0),caux(:,1:M0),infoloc)
endif
    if (infoloc/=0) then
       info=-2
       return
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!! FEAST !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

    ijob=-1 ! initialization

    do while (ijob/=0)

       call cfeast_hrci(ijob,N,Ze,work1,work2,zAq,zSq,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
       select case(ijob)

       case(10) !! Factorize (ZeB-A)
          opt=3
          call caddcsr(N,opt,-ONEC,ssa,sisa,sjsa,Ze,ssb,sisb,sjsb,saz,isaz,jsaz) !! get saz

          PHASE=22
          MNUM=1
          call PARDISO(PT1,MAXFCT,MNUM,MTYPE,PHASE,n,saz,isaz,jsaz,idum,M0,IPARM,MSGLVL,work2(:,1:M0),caux(:,1:M0),infoloc)
          if (infoloc/=0) info=-2

       case(11) !! Solve (ZeB-A)x=work2(1:N,1:M0) result in to work2
          IPARM(12)=0 ! normal solve
          PHASE=33
          MNUM=1
          call PARDISO(PT1,MAXFCT,MNUM,MTYPE,PHASE,n,saz,isaz,jsaz,idum,M0,IPARM,MSGLVL,work2(:,1:M0),caux(:,1:M0),infoloc)
          if (infoloc/=0) info=-2

       case(20) !! Factorize (ZeB-A)^T   
if (fpm(11)==0) then
          opt=3
          call caddcsr(N,opt,-ONEC,ssa,sisa,sjsa,conjg(Ze),ssb,sisb,sjsb,tsaz,isaz,jsaz) !! get saz

          PHASE=22
          MNUM=2
          call PARDISO(PT2,MAXFCT,MNUM,MTYPE,PHASE,n,tsaz,isaz,jsaz,idum,M0,IPARM,MSGLVL,work2(:,1:M0),caux(:,1:M0),infoloc)
          if (infoloc/=0) info=-2
endif

       case(21) !!  Solve  (ZeB-A)^Tx=work2(1:N,1:M0) result in to work2
if (fpm(11)==0) then
          IPARM(12)=0 ! normal solve
          PHASE=33
          MNUM=2
          call PARDISO(PT2,MAXFCT,MNUM,MTYPE,PHASE,n,tsaz,isaz,jsaz,idum,M0,IPARM,MSGLVL,work2(:,1:M0),caux(:,1:M0),infoloc)
          if (infoloc/=0) info=-2
else
 IPARM(12)=1 ! transpose conjugate option solve
          PHASE=33
          MNUM=1
          call PARDISO(PT1,MAXFCT,MNUM,MTYPE,PHASE,n,saz,isaz,jsaz,idum,M0,IPARM,MSGLVL,work2(:,1:M0),caux(:,1:M0),infoloc)
          if (infoloc/=0) info=-2
endif

       case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

          call chcsrmm('F',N,N,fpm(25),ONEC,ssa,sisa,sjsa,X(1,fpm(24)),ZEROC,work1(1,fpm(24)))

       case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

          call chcsrmm('F',N,N,fpm(25),ONEC,ssb,sisb,sjsb,X(1,fpm(24)),ZEROC,work1(1,fpm(24)))

       end select
    end do
!!!!!!!!!!!!!!!!!!!!!!!!
!!!! release memory
!!!!!!!!!!!!!!!!!!!!!!!!
    PHASE=-1 

    MNUM=1
    call PARDISO(PT1,MAXFCT,MNUM,MTYPE,PHASE,n,saz,isaz,jsaz,idum,M0,IPARM,MSGLVL,work2,caux,infoloc)
    if (infoloc/=0) then
       info=-2
       return
    end if

if (fpm(11)==0) then
    MNUM=2
    call PARDISO(PT2,MAXFCT,MNUM,MTYPE,PHASE,n,tsaz,isaz,jsaz,idum,M0,IPARM,MSGLVL,work2,caux,infoloc)
    if (infoloc/=0) then
       info=-2
       return
    end if
endif

    call wdeallocate_2c(zAq)
    call wdeallocate_2c(zSq)
    call wdeallocate_2c(work1)
    call wdeallocate_2c(work2)
    call wdeallocate_2c(caux)

    call wdeallocate_1c(saz)
    call wdeallocate_1i(isaz)
    call wdeallocate_1i(jsaz)

if (fpm(11)==0)  call wdeallocate_1c(tsaz) !! transpose option


    if ((UPLO/='F').and.(UPLO/='f')) then
       call wdeallocate_1c(ssa)
       call wdeallocate_1i(sisa)
       call wdeallocate_1i(sjsa)

       call wdeallocate_1c(ssb)
       call wdeallocate_1i(sisb)
       call wdeallocate_1i(sjsb)
    endif



  end subroutine cfeast_hcsrgv


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine cfeast_hcsrev(UPLO,N,sa,isa,jsa,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
    !  Purpose 
    !  =======
    !  FEAST SPARSE INTERFACE (SPARSE CSR FORMAT)
    !  Solve the standard Ax=ex eigenvalue problem
    !  
    !  A COMPLEX HERMITIAN 
    ! 
    !  SINGLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  UPLO       (input)        CHARACTER: specifies whether the full part, or the upper or lower
    !                                       triangular part of the matrix(ces) is(are) being supplied.
    !  N          (input)        INTEGER: Size system
    !  sa         (input)        COMPLEX SINGLE PRECISION (isa(N+1)-1):  Matrix A- CSR format 
    !  isa        (input)        INTEGER(N+1): CSR row array of Matrix A   
    !  jsa        (input)        INTEGER(isa(N+1)-1): CSR column array of Matrix A
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  fpm (input/output) INTEGER(*) : FEAST parameters
    !  epsout     (output)       REAL SINGLE PRECISION : Error on the trace
    !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
    !  Emin,Emax  (input)        REAL SINGLE PRECISION: search interval
    !  M0         (input/output) INTEGER: Size subspace
    !  lambda     (output)       REAL SINGLE PRECISION(M0)   : Eigenvalues -solution
    !  q          (input/output) COMPLEX SINGLE PRECISION(N,M0) : 
    !                                                       On entry: subspace initial guess if fpm(5)=1 
    !                                                       On exit : Eigenvectors-solution
    !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
    !  res        (output)       REAL SINGLE PRECISION(M0) : Relative Residual of the solution (1-norm)
    !                                                        if option fpm(6)=1 selected                           
    !  info       (output)       INTEGER: Error handling (0: successful exit)
    !=====================================================================
    ! Eric Polizzi 2009-2013
    ! ====================================================================


    implicit none
    include 'f90_noruntime_interface.fi'
    character(len=1) :: UPLO
    integer :: N
    complex,dimension(*),target:: sa
    integer,dimension(*),target:: isa,jsa
    integer,dimension(*) :: fpm
    real :: epsout 
    integer :: loop
    real :: Emin,Emax
    integer :: M0
    real,dimension(*)  :: E
    complex,dimension(N,*):: X
    integer :: mode
    real,dimension(*)    :: res
    integer :: info
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer :: ijob,infoloc
    complex :: Ze
    complex,dimension(:),pointer :: saz,tsaz
    integer,dimension(:),pointer :: isaz,jsaz
    complex,dimension(:,:),pointer ::work1,work2,caux,zAq,zSq
    real,parameter :: SONE=1.0e0, SZERO=0.0e0
    complex,parameter :: ONEC=(SONE,SZERO), ZEROC=(SZERO,SZERO)
!!!!!for pardiso
    integer(8),dimension(64) :: pt1,pt2
    integer,dimension(64) :: iparm
    integer :: mtype
    integer :: MAXFCT,MNUM,PHASE,MSGLVL
    integer :: idum,nrhs
    real :: ddum1
!!!!! full csr format
    complex,dimension(:),pointer :: ssa,ssb
    integer,dimension(:), pointer :: sisa,sjsa,sisb,sjsb
    integer :: i,k 
    integer :: opt,nnza,nnzb,nnz


!!!!!!!!!!!!!! Check INPUT PARAMETERS
    INFO = 0
    IF ((UPLO/='F').and.(UPLO/='f').and.(UPLO/='L').and.(UPLO/='l').and.(UPLO/='U').and.(UPLO/='u')) THEN
       INFO=-101
    ELSE IF ( N<=0 ) THEN
       INFO = -102
    END IF
    IF( INFO.NE.0 ) THEN
       CALL XERBLA( 'CFEAST_HCSREV', -INFO+100 )
       RETURN
    END IF

    info=-1 ! by default
  infoloc=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! identity B matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!
    nnzb=n
    call wallocate_1i(sisb,n+1,infoloc)
    call wallocate_1c(ssb,nnzb,infoloc)
    call wallocate_1i(sjsb,nnzb,infoloc)
    if (infoloc/=0) then
       info=-1
       return
    end if

    do i=1,n
       sisb(i)=i
       sjsb(i)=i
       ssb(i)=ONEC
    enddo
    sisb(n+1)=nnzb+1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! FORMAT CONVERSION TO full CSR for PARDISO !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Remark: everything should be csr since we work with zS-H which is actually unsymmetric 


    if ((UPLO=='F').or.(UPLO=='f')) then !!! full csr already

       nnza=isa(n+1)-1
       ssa => sa(1:nnza)
       sisa => isa(1:n+1)
       sjsa => jsa(1:nnza)

    else !! upper-csr or lower-csr to full csr

       nnza=2*(isa(n+1)-1) ! may be overestimated
       call wallocate_1c(ssa,nnza,infoloc)
       if (infoloc/=0) return
       call wallocate_1i(sjsa,nnza,infoloc)
       if (infoloc/=0) return
       call wallocate_1i(sisa,n+1,infoloc)
       if (infoloc/=0) return
     
       call chcsr_uplo_to_csr(N,sa,isa,jsa,ssa,sisa,sjsa)

    endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! FEAST INITIALIZATION 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call wallocate_2c(zAq,M0,M0,infoloc)
    if (infoloc/=0) return
    call wallocate_2c(zSq,M0,M0,infoloc)
    if (infoloc/=0) return
    call wallocate_2c(work1,N,M0,infoloc)
    if (infoloc/=0) return
    call wallocate_2c(work2,N,M0,infoloc)
    if (infoloc/=0) return
    call wallocate_1i(isaz,n+1,infoloc)
    if (infoloc/=0) return
    !!<<<
    call wallocate_1c(saz,1,infoloc) ! dummy
    if (infoloc/=0) return
    call wallocate_1i(jsaz,1,infoloc)! dummy
    if (infoloc/=0) return
    !!>>>
    opt=1
    call caddcsr(N,opt,ONEC,ssa,sisa,sjsa,ONEC,ssb,sisb,sjsb,saz,isaz,jsaz) !! get isaz
    nnz=isaz(n+1)-1
    !!<<<
    call wdeallocate_1c(saz)
    call wdeallocate_1i(jsaz)
    !!>>>
    call wallocate_1c(saz,nnz,infoloc)
    if (infoloc/=0) return
    call wallocate_1i(jsaz,nnz,infoloc)
    if (infoloc/=0) return


    opt=2
    call caddcsr(N,opt,ONEC,ssa,sisa,sjsa,ONEC,ssb,sisb,sjsb,saz,isaz,jsaz) !! get jsaz
!!!!!!!!!!!!!!!
    call wallocate_2c(caux,N,M0,infoloc)
    if (infoloc/=0) return




 MAXFCT=1 ! One factorization to consider (default)

if (fpm(11)==0) then!!  needs the transpose option
    MAXFCT=2 ! two factorization to consider normal+transpose
    call wallocate_1c(tsaz,nnz,infoloc)
    if (infoloc/=0) return
 end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  analysis step (symbolic factorizarion for pardiso)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    MTYPE=3 ! complex and structurally symmetric 
    call pardisoinit(PT1,MTYPE,IPARM)
if (fpm(11)==0)  call pardisoinit(PT2,MTYPE,IPARM)

!!!!!!!!!!!!
if (fpm(64)==1) then
do i=1,64
if (fpm(64+i)/=-9999) iparm(i)=fpm(64+i)
enddo
endif
!!!!!!!!!!!!

    !IPARM(3)=fpm(10) !omp_num_threads !! openmp -number of threads
    IPARM(6)=1 ! solution and rhs are input/output, attention caux is always used
    MSGLVL=0 !0- no output, 1- output
    PHASE=11
    IPARM(28)=1 ! pardiso single precision (MKL)
    MNUM=1 
    call PARDISO(PT1,MAXFCT,MNUM,MTYPE,PHASE,n,saz,isaz,jsaz,idum,M0,IPARM,MSGLVL,work2(:,1:M0),caux(:,1:M0),infoloc) 
if (fpm(11)==0) then
    MNUM=2 
    call PARDISO(PT2,MAXFCT,MNUM,MTYPE,PHASE,n,tsaz,isaz,jsaz,idum,M0,IPARM,MSGLVL,work2(:,1:M0),caux(:,1:M0),infoloc)
endif
    if (infoloc/=0) then
       info=-2
       return
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!! FEAST !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

    ijob=-1 ! initialization

    do while (ijob/=0)

       call cfeast_hrci(ijob,N,Ze,work1,work2,zAq,zSq,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
       select case(ijob)

       case(10) !! Factorize (ZeB-A)
          opt=3
          call caddcsr(N,opt,-ONEC,ssa,sisa,sjsa,Ze,ssb,sisb,sjsb,saz,isaz,jsaz) !! get saz

          PHASE=22
          MNUM=1
          call PARDISO(PT1,MAXFCT,MNUM,MTYPE,PHASE,n,saz,isaz,jsaz,idum,M0,IPARM,MSGLVL,work2(:,1:M0),caux(:,1:M0),infoloc)
          if (infoloc/=0) info=-2

       case(11) !! Solve (ZeB-A)x=work2(1:N,1:M0) result in to work2
          IPARM(12)=0 ! normal solve
          PHASE=33
          MNUM=1
          call PARDISO(PT1,MAXFCT,MNUM,MTYPE,PHASE,n,saz,isaz,jsaz,idum,M0,IPARM,MSGLVL,work2(:,1:M0),caux(:,1:M0),infoloc)
          if (infoloc/=0) info=-2

       case(20) !! Factorize (ZeB-A)^T   
if (fpm(11)==0) then
          opt=3
          call caddcsr(N,opt,-ONEC,ssa,sisa,sjsa,conjg(Ze),ssb,sisb,sjsb,tsaz,isaz,jsaz) !! get saz

          PHASE=22
          MNUM=2
          call PARDISO(PT2,MAXFCT,MNUM,MTYPE,PHASE,n,tsaz,isaz,jsaz,idum,M0,IPARM,MSGLVL,work2(:,1:M0),caux(:,1:M0),infoloc)
          if (infoloc/=0) info=-2
endif

       case(21) !!  Solve  (ZeB-A)^Tx=work2(1:N,1:M0) result in to work2
if (fpm(11)==0) then
          IPARM(12)=0 ! normal solve
          PHASE=33
          MNUM=2
          call PARDISO(PT2,MAXFCT,MNUM,MTYPE,PHASE,n,tsaz,isaz,jsaz,idum,M0,IPARM,MSGLVL,work2(:,1:M0),caux(:,1:M0),infoloc)
          if (infoloc/=0) info=-2
else
 IPARM(12)=1 ! transpose conjugate option solve
          PHASE=33
          MNUM=1
          call PARDISO(PT1,MAXFCT,MNUM,MTYPE,PHASE,n,saz,isaz,jsaz,idum,M0,IPARM,MSGLVL,work2(:,1:M0),caux(:,1:M0),infoloc)
          if (infoloc/=0) info=-2
endif

       case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

          call chcsrmm('F',N,N,fpm(25),ONEC,ssa,sisa,sjsa,X(1,fpm(24)),ZEROC,work1(1,fpm(24)))

       case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

          call CLACPY( 'F', N, fpm(25),X(1,fpm(24)) , N, work1(1,fpm(24)), N )

       end select
    end do
!!!!!!!!!!!!!!!!!!!!!!!!
!!!! release memory
!!!!!!!!!!!!!!!!!!!!!!!!
    PHASE=-1 

    MNUM=1
    call PARDISO(PT1,MAXFCT,MNUM,MTYPE,PHASE,n,saz,isaz,jsaz,idum,M0,IPARM,MSGLVL,work2,caux,infoloc)
    if (infoloc/=0) then
       info=-2
       return
    end if

if (fpm(11)==0) then
    MNUM=2
    call PARDISO(PT2,MAXFCT,MNUM,MTYPE,PHASE,n,tsaz,isaz,jsaz,idum,M0,IPARM,MSGLVL,work2,caux,infoloc)
    if (infoloc/=0) then
       info=-2
       return
    end if
endif

    call wdeallocate_2c(zAq)
    call wdeallocate_2c(zSq)
    call wdeallocate_2c(work1)
    call wdeallocate_2c(work2)
    call wdeallocate_2c(caux)

    call wdeallocate_1c(saz)
    call wdeallocate_1i(isaz)
    call wdeallocate_1i(jsaz)

if (fpm(11)==0)  call wdeallocate_1c(tsaz) !! transpose option


    if ((UPLO/='F').and.(UPLO/='f')) then
       call wdeallocate_1c(ssa)
       call wdeallocate_1i(sisa)
       call wdeallocate_1i(sjsa)
    endif
       call wdeallocate_1c(ssb)
       call wdeallocate_1i(sisb)
       call wdeallocate_1i(sjsb)
    

  end subroutine cfeast_hcsrev



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine zfeast_gcsrev(N,ssa,sisa,sjsa,fpm,epsout,loop,Emid,r,M0,E,XR,XL,mode,resR,resL,info)
    !  Purpose 
    !  =======
    !  FEAST SPARSE INTERFACE (SPARSE CSR FORMAT)
    !  Solve the standard Ax=ex eigenvalue problem
    !  
    !  A COMPLEX HERMITIAN 
    ! 
    !  DOUBLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  UPLO       (input)        CHARACTER: specifies whether the full part, or the upper or lower
    !                                       triangular part of the matrix(ces) is(are) being supplied.
    !  N          (input)        INTEGER: Size system
    !  sa         (input)        COMPLEX DOUBLE PRECISION (isa(N+1)-1):  Matrix A- CSR format 
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
    !  lambda     (output)       REAL DOUBLE PRECISION(M0)   : Eigenvalues -solution
    !  q          (input/output) COMPLEX DOUBLE PRECISION(N,M0) : 
    !                                                       On entry: subspace initial guess if fpm(5)=1 
    !                                                       On exit : Eigenvectors-solution
    !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
    !  res        (output)       REAL DOUBLE PRECISION(M0) : Relative Residual of the solution (1-norm)
    !                                                        if option fpm(6)=1 selected                           
    !  info       (output)       INTEGER: Error handling (0: successful exit)
    !=====================================================================
    ! Eric Polizzi 2009-2013
    ! ====================================================================


    implicit none
    include 'f90_noruntime_interface.fi'
    !character(len=1) :: UPLO
    integer :: N
    complex(kind=kind(1.0d0)),dimension(*) :: ssa
    integer,dimension(*) :: sisa,sjsa
    integer,dimension(*) :: fpm
    complex(kind=kind(1.0d0)) :: epsout,Emid 
    integer :: loop
    double precision :: r
    integer :: M0
    complex(kind=kind(1.0d0)),dimension(*)  :: E
    complex(kind=kind(1.0d0)),dimension(N,*):: XR,XL
    integer :: mode
    double precision,dimension(*)    :: resR,resL
    integer :: info
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer :: ijob,infoloc
    complex(kind=(kind(1.0d0))) :: Ze
    complex(kind=(kind(1.0d0))),dimension(:),pointer :: saz,tsaz
    integer,dimension(:),pointer :: isaz,jsaz
    complex(kind=(kind(1.0d0))),dimension(:,:),pointer ::workr,workl,workc,caux,zAq,zSq
    double precision,parameter :: DONE=1.0d0, DZERO=0.0d0
    complex(kind=(kind(1.0d0))),parameter :: ONEC=(DONE,DZERO), ZEROC=(DZERO,DZERO)
!!!!!for pardiso
    integer(8),dimension(64) :: pt1,pt2
    integer,dimension(64) :: iparm
    integer :: mtype
    integer :: MAXFCT,MNUM,PHASE,MSGLVL
    integer :: idum,nrhs
    double precision :: ddum1
!!!!! full csr format
    complex(kind=kind(1.0d0)),dimension(:),pointer :: ssb,tssa
    integer,dimension(:),pointer :: sisb,sjsb,tsisa,tsjsa
    integer :: i,k 
    integer :: opt,nnza,nnzb,nnz


!!!!!!!!!!!!!! Check INPUT PARAMETERS
    INFO = 0
    !IF ((UPLO/='F').and.(UPLO/='f').and.(UPLO/='L').and.(UPLO/='l').and.(UPLO/='U').and.(UPLO/='u')) THEN
     !  INFO=-101
    !ELSE 
    IF ( N<=0 ) THEN
       INFO = -102
    END IF
    IF( INFO.NE.0 ) THEN
       CALL XERBLA( 'ZFEAST_GCSREV', -INFO+100 )
       RETURN
    END IF

    infoloc=0
    info=-1 ! by default
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! identity B matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!
    nnzb=n
    call wallocate_1i(sisb,n+1,infoloc)
    call wallocate_1z(ssb,nnzb,infoloc)
    call wallocate_1i(sjsb,nnzb,infoloc)
    if (infoloc/=0) then
       info=-1
       return
    end if

    do i=1,n
       sisb(i)=i
       sjsb(i)=i
       ssb(i)=ONEC
    enddo
    sisb(n+1)=nnzb+1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! FORMAT CONVERSION TO full CSR for PARDISO !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Remark: everything should be csr since we work with zS-H which is actually unsymmetric 


!!! create the transpose for case 31
 call wallocate_1z(tssa,sisa(n+1)-1,infoloc) 
 call wallocate_1i(tsjsa,sisa(n+1)-1,infoloc)
 call wallocate_1i(tsisa,n+1,infoloc)
call zcsr_htranspose(n,ssa,sisa,sjsa,tssa,tsisa,tsjsa)
!do i=1,3
!print *,i,sjsa(sisa(i):sisa(i+1)-1),ssa(sisa(i):sisa(i+1)-1)
!enddo

!do i=1,3
!print *,i,tsjsa(tsisa(i):tsisa(i+1)-1),tssa(tsisa(i):tsisa(i+1)-1)
!enddo
!stop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! FEAST INITIALIZATION 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call wallocate_2z(zAq,M0,M0,infoloc)
    if (infoloc/=0) return
    call wallocate_2z(zSq,M0,M0,infoloc)
    if (infoloc/=0) return
    call wallocate_2z(workr,N,M0,infoloc)
    if (infoloc/=0) return
    call wallocate_2z(workl,N,M0,infoloc)
    if (infoloc/=0) return
    call wallocate_2z(workc,N,M0,infoloc)
    if (infoloc/=0) return
    call wallocate_1i(isaz,n+1,infoloc)
    if (infoloc/=0) return
    !!<<<
    call wallocate_1z(saz,1,infoloc) ! dummy
    if (infoloc/=0) return
    call wallocate_1i(jsaz,1,infoloc)! dummy
    if (infoloc/=0) return
    !!>>>
    opt=1
    call zaddcsr(N,opt,ONEC,ssa,sisa,sjsa,ONEC,ssb,sisb,sjsb,saz,isaz,jsaz) !! get isaz
    nnz=isaz(n+1)-1
    !!<<<
    call wdeallocate_1z(saz)
    call wdeallocate_1i(jsaz)
    !!>>>
    call wallocate_1z(saz,nnz,infoloc)
    if (infoloc/=0) return
    call wallocate_1i(jsaz,nnz,infoloc)
    if (infoloc/=0) return


    opt=2
    call zaddcsr(N,opt,ONEC,ssa,sisa,sjsa,ONEC,ssb,sisb,sjsb,saz,isaz,jsaz) !! get jsaz
!!!!!!!!!!!!!!!
    call wallocate_2z(caux,N,M0,infoloc)
    if (infoloc/=0) return


 MAXFCT=1 ! One factorization to consider (default)

if (fpm(11)==0) then!!!  transpose option not possible 
     MAXFCT=2 ! two factorization to consider normal+transpose
    call wallocate_1z(tsaz,nnz,infoloc)
    if (infoloc/=0) return
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  analysis step (symbolic factorizarion for pardiso)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
    MTYPE=13  ! complex and unsymmetric 
    call pardisoinit(PT1,MTYPE,IPARM)
if (fpm(11)==0) call pardisoinit(PT2,MTYPE,IPARM)

!!!!!!!!!!!!
if (fpm(64)==1) then
do i=1,64
if (fpm(64+i)/=-9999) iparm(i)=fpm(64+i)
enddo
endif
!!!!!!!!!!!!
    !IPARM(3)=fpm(10) !omp_num_threads !! openmp -number of threads
    IPARM(6)=1 ! solution and rhs are input/output, attention caux is always used
    MSGLVL=0 !0- no output, 1- output
    PHASE=11
    MNUM=1 
    call PARDISO(PT1,MAXFCT,MNUM,MTYPE,PHASE,n,saz,isaz,jsaz,idum,M0,IPARM,MSGLVL,workc(:,1:M0),caux(:,1:M0),infoloc)
if (fpm(11)==0) then 
    MNUM=2 
    call PARDISO(PT2,MAXFCT,MNUM,MTYPE,PHASE,n,tsaz,isaz,jsaz,idum,M0,IPARM,MSGLVL,workc(:,1:M0),caux(:,1:M0),infoloc)
endif
    if (infoloc/=0) then
       info=-2
       return
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!! FEAST !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

    ijob=-1 ! initialization

    do while (ijob/=0)

       call zfeast_grci(ijob,N,Ze,workr,workl,workc,zAq,zSq,fpm,epsout,loop,Emid,r,M0,E,XR,XL,mode,resR,resL,info)

       !call zfeast_srci(ijob,N,Ze,workr,workl,workc,zAq,zSq,fpm,epsout,loop,Emid,r,M0,E,XR,XL,mode,resR,resL,info)
       select case(ijob)

       case(10) !! Factorize (ZeB-A)
          opt=3
!print *,'before'
          call zaddcsr(N,opt,-ONEC,ssa,sisa,sjsa,Ze,ssb,sisb,sjsb,saz,isaz,jsaz) !! get saz
!print *,'after'
          PHASE=22
          MNUM=1
          call PARDISO(PT1,MAXFCT,MNUM,MTYPE,PHASE,n,saz,isaz,jsaz,idum,M0,IPARM,MSGLVL,workc(:,1:M0),caux(:,1:M0),infoloc)
          if (infoloc/=0) info=-2

       case(11) !! Solve (ZeB-A)x=workc(1:N,1:M0) result in to workc
          IPARM(12)=0 ! normal solve
          PHASE=33
          MNUM=1
          call PARDISO(PT1,MAXFCT,MNUM,MTYPE,PHASE,n,saz,isaz,jsaz,idum,M0,IPARM,MSGLVL,workc(:,1:M0),caux(:,1:M0),infoloc)
          if (infoloc/=0) info=-2

       case(20) !! Factorize (ZeB-A)^T   
if (fpm(11)==0) then
          opt=3
          call zaddcsr(N,opt,-ONEC,ssa,sisa,sjsa,conjg(Ze),ssb,sisb,sjsb,tsaz,isaz,jsaz) !! get saz

          PHASE=22
          MNUM=2
          call PARDISO(PT2,MAXFCT,MNUM,MTYPE,PHASE,n,tsaz,isaz,jsaz,idum,M0,IPARM,MSGLVL,workc(:,1:M0),caux(:,1:M0),infoloc)
          if (infoloc/=0) info=-2
endif
       case(21) !!  Solve  (ZeB-A)^Tx=workc(1:N,1:M0) result in to workc
if (fpm(11)==0) then
          IPARM(12)=0 ! normal solve
          PHASE=33
          MNUM=2
          call PARDISO(PT2,MAXFCT,MNUM,MTYPE,PHASE,n,tsaz,isaz,jsaz,idum,M0,IPARM,MSGLVL,workc(:,1:M0),caux(:,1:M0),infoloc)
          if (infoloc/=0) info=-2
else
IPARM(12)=1 ! transpose conjugate solve
  PHASE=33
          MNUM=1
          call PARDISO(PT1,MAXFCT,MNUM,MTYPE,PHASE,n,saz,isaz,jsaz,idum,M0,IPARM,MSGLVL,workc(:,1:M0),caux(:,1:M0),infoloc)
          if (infoloc/=0) info=-2
endif

       case(30) !! perform multiplication A*xr(1:N,fpm(24):fpm(24)+fpm(25)-1) result in workr(1:N,fpm(24)+fpm(25)-1)

          call zhcsrmm('F',N,N,fpm(25),ONEC,ssa,sisa,sjsa,XR(1,fpm(24)),ZEROC,workr(1,fpm(24)))

      case(31) !! perform multiplication A^c*xl(1:N,fpm(24):fpm(24)+fpm(25)-1) result in workl(1:N,fpm(24)+fpm(25)-1)

          call zhcsrmm('F',N,N,fpm(25),ONEC,tssa,tsisa,tsjsa,XL(1,fpm(24)),ZEROC,workl(1,fpm(24)))



       case(40) !! perform multiplication B*xr(1:N,fpm(24):fpm(24)+fpm(25)-1) result in workr(1:N,fpm(24)+fpm(25)-1)

          call ZLACPY( 'F', N, fpm(25),XR(1,fpm(24)) , N, workr(1,fpm(24)), N )

       case(41) !! perform multiplication B^c*xl(1:N,fpm(24):fpm(24)+fpm(25)-1) result in workl(1:N,fpm(24)+fpm(25)-1)

          call ZLACPY( 'F', N, fpm(25),XL(1,fpm(24)) , N, workl(1,fpm(24)), N )



       end select
    end do
!!!!!!!!!!!!!!!!!!!!!!!!
!!!! release memory
!!!!!!!!!!!!!!!!!!!!!!!!
    PHASE=-1 

    MNUM=1
    call PARDISO(PT1,MAXFCT,MNUM,MTYPE,PHASE,n,saz,isaz,jsaz,idum,M0,IPARM,MSGLVL,workc,caux,infoloc)
    if (infoloc/=0) then
       info=-2
       return
    end if

if (fpm(11)==0) then
    MNUM=2
    call PARDISO(PT2,MAXFCT,MNUM,MTYPE,PHASE,n,tsaz,isaz,jsaz,idum,M0,IPARM,MSGLVL,workc,caux,infoloc)
    if (infoloc/=0) then
       info=-2
       return
    end if
 end if

    call wdeallocate_2z(zAq)
    call wdeallocate_2z(zSq)
    call wdeallocate_2z(workr)
    call wdeallocate_2z(workl)
    call wdeallocate_2z(workc)
    call wdeallocate_2z(caux)

    call wdeallocate_1z(saz)
    call wdeallocate_1i(isaz)
    call wdeallocate_1i(jsaz)

    call wdeallocate_1z(tssa)
    call wdeallocate_1i(tsisa)
    call wdeallocate_1i(tsjsa)


if (fpm(11)==0)  call wdeallocate_1z(tsaz) !! transpose option

  


       call wdeallocate_1z(ssb)
       call wdeallocate_1i(sisb)
       call wdeallocate_1i(sjsb)

  end subroutine zfeast_gcsrev


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  





      subroutine  gmrez (n,im,rhs,sol,vv,eps,maxits,iout,saz,isa,jsa,UPLO)
       implicit none
       integer n,im,maxits,iout
       real*8 eps 
       complex*16 vv(n,*),sol(n),rhs(n)
       complex*16 saz(*)
       integer isa(*),jsa(*)
       character(len=1) UPLO
!-----------------------------------------------------------------------
! complex gmres algorithm . simple version .  (may 23,1985)
!----------------------------------------------------------------------- 
! parameter list:
! n     == size of problem
! im    == size of krylov subspa!e:  should not ex!eed 50 in this
!          version (can be reset in code. looking at comment below)
! rhs   == right hand side
! sol   == initial guess on input, approximate solution on output
! vv    == work space of size n x (im+1)
! eps   == tolerance for stopping criterion. process is stopped
!          as soon as ( ||.|| is the euclidean norm):
!          || current residual||/||initial residual|| <= eps
! maxits== maximum number of iterations allowed
! iout  == output unit number number for printing intermediate results
!          if (iout .le. 0) no statistics are printed.
!-----------------------------------------------------------------------
! subroutines used =
! ope(n,x,y)  ==  matrix by vector multiplication delivers y=ax, given x.
! ddot        == dot product function.  replace by blas routine ddot
! daxpy       == y <-- y+ax  routine. replace by blas routine daxpy
!----------------------------------------------------------------------- 
      integer ih
      parameter (ih = 50)
      complex*16  hh(ih+1,ih), c(ih), s(ih), rs(ih+1), dconjg,tz 
      real*8 ro, t, gam, epsmac, eps1, znrm2 
      integer n1, j, its, k, k1, i, i1,ii
      complex*16 zdotc 
!-------------------------------------------------------------
!     arnoldi size should not exceed immax in this version..
!     to reset modify parameter instruction accordingly. 
!-------------------------------------------------------------
      data epsmac/1.d-16/
      n1 = n + 1
      its = 0
!
!     outer loop starts here. -compute initial residual  
!
 10   continue
!      call ope (n, sol, vv)

!call zhcsrmm(UPLO,n,n,1,(1.0d0,0.0d0),saz,isa,jsa,sol,(0.0d0,0.0d0),vv)

call zcsrmm(UPLO,n,n,1,(1.0d0,0.0d0),saz,isa,jsa,sol,(0.0d0,0.0d0),vv)



      do 21 j=1,n
         vv(j,1) = rhs(j) - vv(j,1)
 21   continue
!     compute norm of r0 
      ro = znrm2(n, vv) 
      if (ro .eq. 0.0d0) return
      tz = dcmplx(1.0d0/ro,0.0d0) 
      do 210 j=1, n
         vv(j,1) = vv(j,1)*tz
 210  continue
      if (its .eq. 0) eps1=eps*ro
!      if (iout .gt. 0) write(iout, 199) its, ro
!
!     initialize 1-st term  of rhs of hessenberg system..
!     
      rs(1) = dcmplx(ro,0.0d0) 
      i = 0
 4    i=i+1
      its = its + 1
      i1 = i + 1
!      call  ope (n, vv(1,i), vv(1,i1))

!call zhcsrmm(UPLO,n,n,1,(1.0d0,0.0d0),saz,isa,jsa,vv(1,i),(0.0d0,0.0d0),vv(1,i1))

call zcsrmm(UPLO,n,n,1,(1.0d0,0.0d0),saz,isa,jsa,vv(1,i),(0.0d0,0.0d0),vv(1,i1))
!-----------------------------------------
!     modified gram - schmidt...
!-----------------------------------------
      do 55 j=1, i
!         tz = zdot(n, vv(1,i1),1, vv(1,j),1)
          tz = zdotc(n, vv(1,j),1, vv(1,i1),1)

         hh(j,i) = tz
         call zaxpy(n,-tz,vv(1,j),1,vv(1,i1),1)


 55   continue
      t = znrm2(n, vv(1,i1)) 
      hh(i1,i) = dcmplx(t)
      tz = dcmplx(1.0d0/t,0.0d0)
      do 57  k=1,n
         vv(k,i1) = vv(k,i1)*tz
 57   continue 
!
!     done with modified gram schimd and arnoldi step..
!     now  update factorization of hh
!     
      if (i .eq. 1) goto 121
!-------- apply previous transformations to i-th column of h
      do 66 k=2,i
         k1 = k-1
         tz = hh(k1,i)
         hh(k1,i) = dconjg(c(k1))*tz + s(k1)*hh(k,i)
         hh(k,i) = -s(k1)*tz + c(k1)*hh(k,i)
 66   continue
 121  gam=dsqrt(cdabs(hh(i,i))**2 + cdabs(hh(i1,i))**2)
      if (gam .eq. 0.0d0) gam = epsmac
!--------  determine next plane rotation  
      tz = dcmplx(1.0d0/gam , 0.0d0)
      c(i) = hh(i,i)*tz
      s(i) = hh(i1,i)*tz
      rs(i1) = -s(i)*rs(i)
      rs(i) = dconjg(c(i))*rs(i)
!-------- determine residual norm and test for convergence-
      hh(i,i) = dconjg(c(i)) *hh(i,i) + s(i)*hh(i1,i)
      ro = cdabs(rs(i1))
      if (iout .gt. 0) write(iout, 199) its, ro
      if (i .lt. im .and. (ro .gt. eps1))  goto 4

!      if (iout .gt. 0) write(iout, 199) its, ro

!     compute solution. first solve upper triangular system.
      
      rs(i) = rs(i)/hh(i,i)
      do 30 ii=2,i
         k=i-ii+1
         tz=rs(k)
         do 40 j=k+1,i
            tz = tz-hh(k,j)*rs(j)
 40      continue 
         rs(k) = tz/hh(k,k)
 30   continue
!
!     done with back substitution..
!     now form linear combination to get solution
!     
      do 16 j=1, i
         tz = rs(j)
         call zaxpy(n,tz,vv(1,j),1,sol,1)
 16   continue
!---------  restart outer loop  when necessary
      if (ro .gt. eps1 .and. its .lt. maxits) goto 10
 199  format('   its =', i4, ' res. norm =', d20.6)
      return
!-----------------------------------------------------------------------
!-------------------------------end-of-gmrez --------------------------- 
      end
!---------------replace the following 2 routines by the corresponding--- 
! blas routines -------------------------------------------------------- 
!----------------------------------------------------------------------- 
!      subroutine zaxpy(n,t,x,indx,y,indy)
!      complex*16 x(*), y(*), t
!-------------------------------------------------------------------
!     does the following operation --- replce by blas equivalent...
!     y <--- y + t * x ,   (replace by the blas routine daxpy )
!     indx and indy are supposed to be one here
!-------------------------------------------------------------------
!      do 1 k=1,n
!         y(k) = y(k) + x(k)*t
! 1    continue
!      return
!--------------end of daxpy ----------------------------------------
!      end
!--------------------------------------------------------------------
      function zdot(n,x,ix,y,iy)
      complex*16 zdot, x(*), y(*), t
!--------------------------------------------------------------------
!     computes the inner product t=(x,y) -- replace by blas routine..
!--------------------------------------------------------------------
      t = (0.0d0,0.0d0)
      do 10 j=1,n
         t = t + x(j)*dconjg(y(j))

 10   continue
      zdot=t
      return
!*******end of zdot   *******************************************
      end
!
      function znrm2(n,x)
      complex*16 x(*) 
      real*8 znrm2, t
!-------------------------------------------------------------------
!     computes the inner product t=(x,y) -- replace by blas routine..
!-------------------------------------------------------------------
      t = 0.0d0
      do 10 j=1,n
         t = t + cdabs(x(j))**2
 10   continue
      znrm2 = dsqrt(t)
      return
!*******end of znrm2   *******************************************
      end










      subroutine  Bgmrez (n,im,rhs,sol,vv,eps,maxits,iout,saz,isa,jsa,sbz,UPLO)
       implicit none
       integer n,im,maxits,iout
       real*8 eps 
       complex*16 vv(n,*),sol(n),rhs(n)
       complex*16 saz(*),sbz(*)
       integer isa(*),jsa(*)
       character(len=1) UPLO
!-----------------------------------------------------------------------
! complex gmres algorithm . simple version .  (may 23,1985)
!----------------------------------------------------------------------- 
! parameter list:
! n     == size of problem
! im    == size of krylov subspa!e:  should not ex!eed 50 in this
!          version (can be reset in code. looking at comment below)
! rhs   == right hand side
! sol   == initial guess on input, approximate solution on output
! vv    == work space of size n x (im+1)
! eps   == tolerance for stopping criterion. process is stopped
!          as soon as ( ||.|| is the euclidean norm):
!          || current residual||/||initial residual|| <= eps
! maxits== maximum number of iterations allowed
! iout  == output unit number number for printing intermediate results
!          if (iout .le. 0) no statistics are printed.
!-----------------------------------------------------------------------
! subroutines used =
! ope(n,x,y)  ==  matrix by vector multiplication delivers y=ax, given x.
! ddot        == dot product function.  replace by blas routine ddot
! daxpy       == y <-- y+ax  routine. replace by blas routine daxpy
!----------------------------------------------------------------------- 
      integer ih,p,pk
      parameter (ih = 50)
      complex*16  hh(ih+1,ih), c(ih), s(ih), rs(ih+1), dconjg,tz 
      complex*16 ztemp(n)
      real*8 ro, t, gam, epsmac, eps1, znrm2 
      integer n1, j, its, k, k1, i, i1,ii
      complex*16 zdotc 
!-------------------------------------------------------------
!     arnoldi size should not exceed immax in this version..
!     to reset modify parameter instruction accordingly. 
!-------------------------------------------------------------
      data epsmac/1.d-16/
      n1 = n + 1
      its = 0
!
!     outer loop starts here. -compute initial residual  
!
 10   continue
!      call ope (n, sol, vv)

!call zhcsrmm(UPLO,n,n,1,(1.0d0,0.0d0),saz,isa,jsa,sol,(0.0d0,0.0d0),vv)

call zcsrmm(UPLO,n,n,1,(1.0d0,0.0d0),saz,isa,jsa,sol,(0.0d0,0.0d0),vv)



      do 21 j=1,n
         vv(j,1) = rhs(j) - vv(j,1)
 21   continue
!     compute norm of r0 
!      ro = znrm2(n, vv)
call zcsrmm(UPLO,n,n,1,(1.0d0,0.0d0),sbz,isa,jsa,vv,(0.0d0,0.0d0),ztemp)
ro = zdotc(n, vv,1, ztemp,1)
ro=dsqrt(abs(ro))

 
      if (ro .eq. 0.0d0) return
      tz = dcmplx(1.0d0/ro,0.0d0) 
      do 210 j=1, n
         vv(j,1) = vv(j,1)*tz
 210  continue
      if (its .eq. 0) eps1=eps*ro
!      if (iout .gt. 0) write(iout, 199) its, ro
!
!     initialize 1-st term  of rhs of hessenberg system..
!     
      rs(1) = dcmplx(ro,0.0d0) 
      i = 0
 4    i=i+1
      its = its + 1
      i1 = i + 1
!      call  ope (n, vv(1,i), vv(1,i1))

!call zhcsrmm(UPLO,n,n,1,(1.0d0,0.0d0),saz,isa,jsa,vv(1,i),(0.0d0,0.0d0),vv(1,i1))

call zcsrmm(UPLO,n,n,1,(1.0d0,0.0d0),saz,isa,jsa,vv(1,i),(0.0d0,0.0d0),vv(1,i1))
!-----------------------------------------
!     modified gram - schmidt...
!-----------------------------------------
      do 55 j=1, i
!!!!         tz = zdot(n, vv(1,i1),1, vv(1,j),1)

!          tz = zdotc(n, vv(1,j),1, vv(1,i1),1)
!sbz(1:isa(n)-1)=(0.0d0,0.0d0)
!do p=1,n
!sbz(isa(p))=(1.0d0,0.0d0)
!enddo
call zcsrmm(UPLO,n,n,1,(1.0d0,0.0d0),sbz,isa,jsa,vv(1,i1),(0.0d0,0.0d0),ztemp)
tz = zdotc(n, vv(1,j),1, ztemp,1)


         hh(j,i) = tz
         call zaxpy(n,-tz,vv(1,j),1,vv(1,i1),1)


 55   continue
!      t = znrm2(n, vv(1,i1)) 
call zcsrmm(UPLO,n,n,1,(1.0d0,0.0d0),sbz,isa,jsa,vv(1,i1),(0.0d0,0.0d0),ztemp)
t = zdotc(n, vv(1,i1),1, ztemp,1)
t=dsqrt(abs(t))



      hh(i1,i) = dcmplx(t)
      tz = dcmplx(1.0d0/t,0.0d0)
      do 57  k=1,n
         vv(k,i1) = vv(k,i1)*tz
 57   continue 
!
!     done with modified gram schimd and arnoldi step..
!     now  update factorization of hh
!     
      if (i .eq. 1) goto 121
!-------- apply previous transformations to i-th column of h
      do 66 k=2,i
         k1 = k-1
         tz = hh(k1,i)
         hh(k1,i) = dconjg(c(k1))*tz + s(k1)*hh(k,i)
         hh(k,i) = -s(k1)*tz + c(k1)*hh(k,i)
 66   continue
 121  gam=dsqrt(cdabs(hh(i,i))**2 + cdabs(hh(i1,i))**2)
      if (gam .eq. 0.0d0) gam = epsmac
!--------  determine next plane rotation  
      tz = dcmplx(1.0d0/gam , 0.0d0)
      c(i) = hh(i,i)*tz
      s(i) = hh(i1,i)*tz
      rs(i1) = -s(i)*rs(i)
      rs(i) = dconjg(c(i))*rs(i)
!-------- determine residual norm and test for convergence-
      hh(i,i) = dconjg(c(i)) *hh(i,i) + s(i)*hh(i1,i)
      ro = cdabs(rs(i1))
      if (iout .gt. 0) write(iout, 199) its, ro
      if (i .lt. im .and. (ro .gt. eps1))  goto 4

!      if (iout .gt. 0) write(iout, 199) its, ro

!     compute solution. first solve upper triangular system.
      
      rs(i) = rs(i)/hh(i,i)
      do 30 ii=2,i
         k=i-ii+1
         tz=rs(k)
         do 40 j=k+1,i
            tz = tz-hh(k,j)*rs(j)
 40      continue 
         rs(k) = tz/hh(k,k)
 30   continue
!
!     done with back substitution..
!     now form linear combination to get solution
!     
      do 16 j=1, i
         tz = rs(j)
         call zaxpy(n,tz,vv(1,j),1,sol,1)
 16   continue
!---------  restart outer loop  when necessary
      if (ro .gt. eps1 .and. its .lt. maxits) goto 10
 199  format('   its =', i4, ' res. norm =', d20.6)
      return
!-----------------------------------------------------------------------
!-------------------------------end-of-gmrez --------------------------- 
      end















SUBROUTINE zcsrsv(UPLO,N,rhs,a,ia,ja,B)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      Purpose: This subroutine solve
  !
  !      b=A^-1*x      (Ax=b A=L or A=U)
  !  
  !      where A is a NxN triangular matrix stored in CSR format (a,ia,ja), 
  !      and X is Mxrhs and B is Nxrhs 
  !           
  !
  !      Arguments: 
  !
  !  UPLO (input) CHARACTER*1
  !       if UPLO='L'  A is lower triangular
  !       if UPLO='U'  A is upper triangular 
  !
  !  N    (input) INTEGER
  !        The number of row/column of the matrix A and row of matrix B.  N >= 0.
  !  rhs  (input) INTEGER
  !        The number of column of matrices B and C
  ! 
  !
  !  a,ia,ja (input) CSR format for the matrix A (a DOUBLE PRECISION)
  !
  !  B     (input/output) DOUBLE PRECISION
  !        matrix of size (Nxrhs). On exit contains the solution.
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  implicit none
  character(len=1) :: UPLO
  integer :: N,rhs
  complex(kind=kind(1.0d0)),dimension(*) ::a
  integer,dimension(*) :: ia
  integer,dimension(*) :: ja
  complex(kind=kind(1.0d0)),dimension(N,*) ::b
!!!!!!!!        

integer ::i,j,k
complex(kind=kind(1.0d0)),dimension(rhs) :: summ



  if (UPLO=='L') then


     DO  i=1,N
        summ=b(i,1:rhs)
        DO  k=ia(i), ia(i+1)-1
           if (ja(k)<i)  summ=summ-a(k)*b(ja(k),1:rhs)
           if (ja(k)==i) j=k
        end DO
        b(i,1:rhs)=summ/a(j)
     end DO

  elseif (UPLO=='U') then

     DO  i=N,1,-1
        summ=b(i,1:rhs)
        DO  k=ia(i), ia(i+1)-1
           if (ja(k)>i)  summ=summ-a(k)*b(ja(k),1:rhs)
           if (ja(k)==i) j=k
        end DO
        b(i,1:rhs)=summ/a(j)
     end DO



  end if


end SUBROUTINE zcsrsv





SUBROUTINE zcsrmm(UPLO,N,M,rhs,alpha,a,ia,ja,X,beta,B)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      Purpose: This subroutine peforms
  !
  !      b=alpha*A*x+beta*b      
  !  
  !      where the NxM matrix A is stored in CSR format (a,ia,ja), 
  !      and X is Mxrhs and B is Nxrhs 
  !           
  !
  !      Arguments: 
  !
  !  UPLO (input) CHARACTER*1
  !       if UPLO='F' all the matrix A is provided in CSR format ('F'ull format)
  !       if UPLO='L' only the lower part of the matrix A is provided- symmetric CSR (N must be equal to M)
  !       if UPLO='U' only the upper part of the matrix A is provided- symmetric CSR (N must be equal to M) 
  !
  !  N    (input) INTEGER
  !        The number of row of the matrix A and row of matrix B.  N >= 0.
  !  M    (input) INTEGER 
  !        The number of column of matrix A and row of matrix X. M>=0; M=N (square matrix) for UPLO=L,U
  !  rhs  (input) INTEGER
  !        The number of column of matrices B and C
  ! 
  !  alpha (input) DOUBLE PRECISION
  !
  !  a,ia,ja (input) CSR format for the matrix A (a DOUBLE PRECISION)
  !
  !  X     (input) DOUBLE PRECISION
  !        matrix of size (Mxrhs)
  !
  !  beta  (input) DOUBLE PRECISION
  !
  !  B     (input/output) DOUBLE PRECISION
  !        matrix of size (Nxrhs). On exit contains the solution.
 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  implicit none
  character(len=1) :: UPLO
  integer :: N,M,rhs
  complex(kind=kind(1.0d0)) :: alpha,beta
  complex(kind=kind(1.0d0)),dimension(*) ::a
  integer,dimension(*) :: ia
  integer,dimension(*) :: ja
  complex(kind=kind(1.0d0)),dimension(M,*):: x
  complex(kind=kind(1.0d0)),dimension(N,*) ::b
!!!!!!!!        

  integer ::i,j,k,s

  complex(kind=kind(1.0d0)),parameter :: ZERO=(0.0d0,0.0d0)
  complex(kind=kind(1.0d0)),dimension(rhs) :: summ

  if (UPLO=='F') then

     do i=1,N
        summ=ZERO
        do k=ia(i),ia(i+1)-1
           j=ja(k)
           summ(1:rhs)=summ(1:rhs)+alpha*a(k)*x(j,1:rhs)
        end do
        b(i,1:rhs)=beta*b(i,1:rhs)+summ(1:rhs)
     end do

  elseif ((UPLO=='L').or.(UPLO=='U')) then

     do s=1,rhs

        if (beta/=ZERO) then
           b(1:N,s)=beta*b(1:N,s)
        else
           b(1:N,s)=ZERO
        endif
        do i=1,N
           do k=ia(i),ia(i+1)-1
              j=ja(k)
              b(i,s)=b(i,s)+alpha*a(k)*x(j,s)
              if (j/=i) b(j,s)=b(j,s)+alpha*a(k)*x(i,s)
           end do
        end do

     enddo

  end if

end SUBROUTINE zcsrmm








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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
    ONEC=(1.0d0,0.0d0)
    ZEROC=(0.0d0,0.0d0)


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
 call zcsrmm(UPLO,N,N,M0,-ONEC,saz,isa,jsa,xj,(1.0d0,0.0d0),rj)

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
       call zcsrmm(UPLO,N,N,M0,ONEC,saz,isa,jsa,pb,(0.0d0,0.0d0),v)


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
          call zcsrmm(UPLO,N,N,M0,ONEC,saz,isa,jsa,sb,(0.0d0,0.0d0),t)


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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine  zbicgstabp(UPLO,N,saz,isa,jsa,X,lambda,Ze,M,M0,fj,xj,res,ares,nbit_out,epso,comd,info) 
    implicit none
    character(len=1) :: UPLO
    integer :: N,M0,M
    complex(kind=kind(1.0d0)),dimension(*) :: saz
    complex(kind=kind(1.0d0)) :: Ze
    double precision, dimension(N,*) :: X
    double precision,dimension(M) :: lambda
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
    complex(kind=kind(1.0d0)), dimension(:,:),allocatable :: rb,pp,pb,v,ss,sb,t,rj,cX,ztemp
    complex(kind=kind(1.0d0)),dimension(:),allocatable ::rho_1,prho_1,pbeta,palpha,pomega,prho_2,aux0,paux0,aux1,paux1
    double precision,dimension(:),allocatable :: resf,nss

    complex(kind=kind(1.0d0)) :: ONEC,ZEROC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
    ONEC=(1.0d0,0.0d0)
    ZEROC=(0.0d0,0.0d0)


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


    allocate(ztemp(1:M,1:M0))
    allocate(cX(1:N,1:M))

    cX(1:N,1:M)=(1.0d0,0.0d0)*X(1:N,1:M)
!!!!!!!!!!!!!


    s=M0
    sizeAj=N

    info=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Initial rj,xj,res and resf !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!! SOLVE M*xj=fj ==> xj
 !!   call ZCOPY(sizeAj*s,fj(1,1),1,xj(1,1),1)   
 !!   call zSOLVER_for_BICGSTAB(pre,xj,info)

 call ZGEMM('C','N',M,M0,N,(1.0d0,0.0d0),cX,N,fj,N,(0.0d0,0.0d0),ztemp,M)
do i=1,M
ztemp(i,1:M0)=ztemp(i,1:M0)/(Ze-lambda(i)*(1.0d0,0.0d0))
enddo
call ZGEMM('N','N',N,M0,M,(1.0d0,0.0d0),cX,N,ztemp,M,(0.0d0,0.0d0),xj,N)


      
!!!!!!!!!!!!!!! MATVEC for the initial residual r=f-Ax
    call ZCOPY(sizeAj*s,fj(1,1),1,rj(1,1),1)

!    call zMATVEC_for_BICGSTAB(-ONEC,mat,xj,ONEC,rj)
 call zcsrmm(UPLO,N,N,M0,-ONEC,saz,isa,jsa,xj,(1.0d0,0.0d0),rj)

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
!       call ZCOPY(sizeAj*s,pp(1,1),1,pb(1,1),1)
       if (comd) call system_clock(tc1,tim) 
       !call zSOLVER_for_BICGSTAB(pre,pb,info)


call ZGEMM('C','N',M,M0,N,(1.0d0,0.0d0),cX,N,pp,N,(0.0d0,0.0d0),ztemp,M)
do i=1,M
ztemp(i,1:M0)=ztemp(i,1:M0)/(Ze-lambda(i)*(1.0d0,0.0d0))
enddo
call ZGEMM('N','N',N,M0,M,(1.0d0,0.0d0),cX,N,ztemp,M,(0.0d0,0.0d0),pb,N)





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
       call zcsrmm(UPLO,N,N,M0,ONEC,saz,isa,jsa,pb,(0.0d0,0.0d0),v)


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
          !call ZCOPY(sizeAj*s,ss(1,1),1,sb(1,1),1)
          if (comd) call system_clock(tc1,tim)
          !call zSOLVER_for_BICGSTAB(pre,sb,info)

call ZGEMM('C','N',M,M0,N,(1.0d0,0.0d0),cX,N,ss,N,(0.0d0,0.0d0),ztemp,M)
do i=1,M
ztemp(i,1:M0)=ztemp(i,1:M0)/(Ze-lambda(i)*(1.0d0,0.0d0))
enddo
call ZGEMM('N','N',N,M0,M,(1.0d0,0.0d0),cX,N,ztemp,M,(0.0d0,0.0d0),sb,N)




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
          call zcsrmm(UPLO,N,N,M0,ONEC,saz,isa,jsa,sb,(0.0d0,0.0d0),t)


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

  end subroutine zbicgstabp


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!







