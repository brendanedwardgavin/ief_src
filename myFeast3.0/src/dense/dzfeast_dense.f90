!=========================================================================================
!Copyright (c) 2009-2015, The Regents of the University of Massachusetts, Amherst.
!E. Polizzi research lab
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
!!!!!!!!!!!!!!!!! FEAST PREDEFINED DENSE INTERFACES !!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!! DOUBLE PRECISION VERSION

! List of routines:
!-------------------

!{D,Z}FEAST_{SY,HE,GE}{EV,GV}


!111111111111!!!!!! EXPERT ROUTINE for GENERALIZED PROBLEM (CORE)
! dfeast_sygvx
! zfeast_hegvx
! dfeast_gegvx
! zfeast_gegvx
! zfeast_sygvx

!222222222222!!!!  DEFAULT ROUTINES for GENERALIZED PROBLEM (Wrappers to expert generalized)
! dfeast_sygv
! zfeast_hegv
! dfeast_gegv
! zfeast_gegv
! zfeast_sygv

!333333333333!!!! EXPERT ROUTINE for STANDARD PROBLEM (Wrappers to expert generalized)
! dfeast_syevx
! zfeast_heevx
! dfeast_geevx
! zfeast_geevx
! zfeast_syevx

!44444444444!!!!! DEFAULT ROUTINES for STANDARD PROBLEM (Wrappers to expert generalized)
! dfeast_syev
! zfeast_heev
! dfeast_geev
! zfeast_geev
! zfeast_syev




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!11111111111111111111111111111111111111111111111111111111111
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine dfeast_sygvx(UPLO,N,A,LDA,B,LDB,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the generalized Ax=eBx eigenvalue problem
  !  
  !  A REAL SYMMETRIC, B SYMMETRIC POSITIVE DEFINITE:: DENSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
  !                           triangular part of the matrix(ces) is(are) being supplied.
  !  N          (input)        INTEGER: Size system
  !  A          (input)        REAL DOUBLE PRECISION (LDA,N):  Matrix A 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=N)
  !  B          (input)        REAL DOUBLE PRECISION (LDB,N):  Matrix B 
  !  LDB        (input)        INTEGER: Leading dimension of matrix B (LDB>=N)
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
  integer :: N,LDA,LDB
  double precision,dimension(LDA,*):: A
  double precision,dimension(LDB,*):: B
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
  complex(kind=(kind(1.0d0))), dimension(:,:),pointer ::workc
  double precision, dimension(:,:),pointer ::work,Aq,Sq
  complex(kind=(kind(1.0d0))), dimension(:,:,:),pointer :: Az
  integer, dimension(:,:),pointer ::ipivloc
  double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
  complex(kind=(kind(1.0d0))),parameter :: ONEC=(DONE,DZERO)
  logical :: fact
  integer :: nfact,id
  character(len=1) :: UPLO2
  integer :: rank,code,nb_procs,NEW_COMM_WORLD

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
  ELSE IF(LDA<N ) THEN
     INFO = -104
  ELSE IF((LDB<N ).and.(LDB/=-1)) THEN
     INFO = -106
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'DFEAST_SYGV', -INFO+100 )
     RETURN
  END IF
  infoloc=0

  call wallocate_2d(Aq,M0,M0,infoloc)
  call wallocate_2d(Sq,M0,M0,infoloc)
  call wallocate_2d(work,N,M0,infoloc)
  call wallocate_2z(workc,N,M0,infoloc)

  if (infoloc/=0) then
     info=-1
     return
  end if

!!! Format conversion
  UPLO2=UPLO
  IF ((UPLO=='F').or.(UPLO=='f')) UPLO2='L'
  IF ((UPLO=='L').or.(UPLO=='l')) UPLO2='L'
  IF ((UPLO=='U').or.(UPLO=='u')) UPLO2='U'

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
!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call wallocate_3z(Az,N,N,nfact,infoloc)
  call wallocate_2i(ipivloc,N,nfact,infoloc)
  if (infoloc/=0) then
     info=-1
     return
  end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ijob=-1 ! initialization 
  do while (ijob/=0)
     call dfeast_srcix(ijob,N,Ze,work,workc,Aq,Sq,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)    

     select case(ijob)
     case(10) !! factorize (zeB-A)

        if (fpm(10)==1) then
           id=id+1
           if (id==nfact+1) then
              id=1
              fact=.false.
           endif
        endif

        if (fact) then
           if (LDB==-1) then !! standard
              Az(1:N,1:N,id)=-A(1:N,1:N)*ONEC
              do i=1,N
                 Az(i,i,id)=Az(i,i,id)+Ze
              enddo
           else        
              Az(1:N,1:N,id)=Ze*B(1:N,1:N)-A(1:N,1:N)*ONEc 
           end if

           if(fpm(11)==0) then
           call ZSYTRF(UPLO2,N,Az(1,1,id),N,IPIVloc(1,id),workc,N*M0,INFOloc)     
           if (infoloc/=0) then
              info=-2
              return
           end if
           end if
        endif ! fact true

     case(11) !!solve the linear system (ZeB-A)x=workc(1:N,1:fpm(23)) result in to workc

        !if(fpm(11)==0) then
        call ZSYTRS( UPLO2, N, fpm(23), Az(1,1,id), N, IPIVloc(1,id), workc, N, INFOloc )
        if (infoloc/=0) then
           info=-2
           return
        end if
        !else
        !    !print *,'starting rfd'
        !    if(loop>3) then
        !        call refineSolveDense('N',Az(1,1,id),workc,n,fpm(23),fpm(50),fpm(51),0)
        !    else
        !        call refineSolveDense('N',Az(1,1,id),workc,n,fpm(23),fpm(50),fpm(51),0)
        !    end if
        !    !print *,'ending rfd'
        !end if

     case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

        call DSYMM ('L', UPLO2, N, fpm(25), DONE, A, LDA, X(1,fpm(24)), N, DZERO,work(1,fpm(24)), N)

     case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

        if (LDB==-1) then ! standard
           call DLACPY( 'F', N, fpm(25),X(1,fpm(24)) , N, work(1,fpm(24)), N )
        else
           call DSYMM ('L', UPLO2, N, fpm(25), DONE, B, LDB, X(1,fpm(24)), N, DZERO,work(1,fpm(24)), N)
        end if

     end select
  end do



  call wdeallocate_2d(Aq)
  call wdeallocate_2d(Sq)
  call wdeallocate_2d(work)
  call wdeallocate_2z(workc)
  call wdeallocate_3z(Az)
  call wdeallocate_2i(ipivloc)

end subroutine dfeast_sygvx





subroutine zfeast_hegvx(UPLO,N,A,LDA,B,LDB,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the generalized Ax=eBx eigenvalue problem
  !  
  !  A COMPLEX HERMITIAN, B HERMITIAN POSITIVE DEFINITE:: DENSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
  !                           triangular part of the matrix(ces) is(are) being supplied.
  !  N          (input)        INTEGER: Size system
  !  A          (input)        COMPLEX DOUBLE PRECISION (LDA,N):  Matrix A 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=N)
  !  B          (input)        COMPLEX DOUBLE PRECISION (LDB,N):  Matrix B 
  !  LDB        (input)        INTEGER: Leading dimension of matrix B (LDB>=N)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  fpm (input/output) INTEGER(*) : FEAST parameters
  !  epsout     (output)       REAL DOUBLE PRECISION : Error on the trace
  !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
  !  Emin,Emax  (input)        REAL DOUBLE PRECISION: search interval
  !  M0         (input/output) INTEGER: Size subspace
  !  E          (output)       REAL DOUBLE PRECISION(M0)   : Eigenvalues -solution
  !  q          (input/output) COMPLEX DOUBLE PRECISION(N,M0) : 
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
  integer :: N,LDA,LDB
  complex(kind=(kind(1.0d0))),dimension(LDA,*):: A
  complex(kind=(kind(1.0d0))),dimension(LDB,*):: B
  integer,dimension(*) :: fpm
  double precision :: epsout 
  integer :: loop
  double precision :: Emin,Emax
  integer :: M0
  double precision,dimension(*)  :: E
  complex(kind=(kind(1.0d0))),dimension(N,*):: X
  integer :: mode
  double precision,dimension(*)    :: res
  integer :: info
  complex(kind=(kind(1.0d0))),dimension(*) :: Zne,Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer :: ijob,infoloc,i,s
  complex(kind=(kind(1.0d0))) :: Ze
  complex(kind=(kind(1.0d0))), dimension(:,:),pointer ::work,workc,zAq,zSq
  complex(kind=(kind(1.0d0))), dimension(:,:,:),pointer :: Az,cAz
  complex(kind=(kind(1.0d0))), dimension(:),pointer ::ztmp
  integer, dimension(:,:),pointer ::ipivloc
  double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
  complex(kind=(kind(1.0d0))),parameter :: ONEC=(DONE,DZERO),ZEROC=(DZERO,DZERO)
  logical :: fact
  integer :: nfact,id
  integer :: rank,code,nb_procs,NEW_COMM_WORLD

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
  ELSE IF(LDA<N ) THEN
     INFO = -104
  ELSE IF((LDB<N ).and.(LDB/=-1)) THEN
     INFO = -106
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'ZFEAST_HEGV', -INFO+100 )
     RETURN
  END IF
  infoloc=0

  call wallocate_2z(zAq,M0,M0,infoloc)
  call wallocate_2z(zSq,M0,M0,infoloc)
  call wallocate_2z(work,N,M0,infoloc)
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call wallocate_3z(Az,N,N,nfact,infoloc)
  call wallocate_3z(cAz,N,N,nfact,infoloc)

  call wallocate_2i(ipivloc,N,nfact,infoloc)
  call wallocate_1z(ztmp,N,infoloc)
  if (infoloc/=0) then
     info=-1
     return
  end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ijob=-1 ! initialization
  do while (ijob/=0) 
     call zfeast_hrcix(ijob,N,Ze,work,workc,zAq,zSq,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)

     select case(ijob)
     case(10) !! factorize (zeB-A)

        if (fpm(10)==1) then
           id=id+1
           if (id==nfact+1) then
              id=1
              fact=.false.
           endif
        endif

        if (fact) then
           if (LDB==-1) then !! standard
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Format CONVERSION

              if ((UPLO=='L').or.(UPLO=='l')) then
                 do i=1,N-1
                    s=N-i
                    Az(i:N,i,id)=-A(i:N,i)
                    Az(i,i,id)=Az(i,i,id)+Ze
                    call ZCOPY(s,Az(i+1,i,id),1,Az(i,i+1,id),N)
                    call ZLACGV(s, Az(i,i+1,id), N )
                 enddo
                 Az(N,N,id)=Ze-A(N,N)*ONEC

              elseif ((UPLO=='U').or.(UPLO=='u')) then
                 do i=1,N-1
                    s=N-i
                    Az(i,i:N,id)=-A(i,i:N)
                    Az(i,i,id)=Az(i,i,id)+Ze
                    call ZCOPY(s,Az(i,i+1,id),N,Az(i+1,i,id),1)
                    call ZLACGV(s, Az(i+1,i,id), 1 )
                 enddo
                 Az(N,N,id)=Ze-A(N,N)*ONEC

              else

                 Az(1:N,1:N,id)=-A(1:N,1:N)
                 do i=1,N
                    Az(i,i,id)=Az(i,i,id)+Ze
                 enddo

              end if


           else 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Format CONVERSION
              if ((UPLO=='L').or.(UPLO=='l')) then
                 do i=1,N-1
                    s=N-i
                    Az(i:N,i,id)=Ze*B(i:N,i)
                    call ZCOPY(s,Az(i+1,i,id),1,Az(i,i+1,id),N)
                    call ZLACGV(s, Az(i,i+1,id), N )
                    call ZSCAL(s,Ze/conjg(Ze),Az(i,i+1,id),N)
                    Az(i:N,i,id)=Az(i:N,i,id)-A(i:N,i)
                 enddo
                 Az(N,N,id)=Ze*B(N,N)-A(N,N)*ONEC
                 do i=1,N-1
                    s=N-i
                    call ZCOPY(s,A(i+1,i),1,ztmp(1),1)
                    call ZLACGV(s,ztmp(1),1)
                    call ZAXPY(s,-ONEC,ztmp(1),1,Az(i,i+1,id),N)
                 enddo


              elseif ((UPLO=='U').or.(UPLO=='u')) then

                 do i=1,N-1
                    s=N-i
                    Az(i,i:N,id)=Ze*B(i,i:N)
                    call ZCOPY(s,Az(i,i+1,id),N,Az(i+1,i,id),1)
                    call ZLACGV(s, Az(i+1,i,id), 1 )
                    call ZSCAL(s,Ze/conjg(Ze),Az(i+1,i,id),1)
                    Az(i,i:N,id)=Az(i,i:N,id)-A(i,i:N)
                 enddo
                 Az(N,N,id)=Ze*B(N,N)-A(N,N)*ONEC
                 do i=1,N-1
                    s=N-i
                    call ZCOPY(s,A(i,i+1),N,ztmp(1),1)
                    call ZLACGV(s,ztmp(1),1)
                    call ZAXPY(s,-ONEC,ztmp(1),1,Az(i+1,i,id),1)
                 enddo
              else ! full 
                 Az(1:N,1:N,id)=Ze*B(1:N,1:N)-A(1:N,1:N)
              end if
           end if

           if(fpm(11)==0) then
           call ZGETRF(N,N,Az(1,1,id),N,IPIVloc(1,id),INFOloc)     
           if (infoloc/=0) then
              info=-2
              return
           end if
           end if
        end if ! fact true

     case(11) !!solve the linear system (ZeB-A)x=workc(1:N,1:fpm(23)) result in to workc
        if(fpm(11)==0) then
        call ZGETRS( 'N', N, fpm(23), Az(1,1,id), N, IPIVloc(1,id), workc, N, INFOloc )
        if (infoloc/=0) then
           info=-2
           return
        end if
        else
            !call refineSolveDense('N',Az(1,1,id),workc,n,fpm(23),fpm(50),fpm(51))
            if(loop>3) then
                call refineSolveDense('N',Az(1,1,id),workc,n,fpm(23),fpm(50),fpm(51),0)
            else
                call refineSolveDense('N',Az(1,1,id),workc,n,fpm(23),fpm(50),fpm(51),0)
            end if
        end if

     case(20) 
        !cAz(1:N,1:N,id)=conjg(Ze)*B(1:N,1:N)-A(1:N,1:N)
        cAz(1:N,1:N,id)=transpose(conjg(Az(1:n,1:n,id)))

     case(21) !!solve the linear system (ZeB-A)^H x=workc(1:N,1:fpm(23)) result in to workc
        !cAz(1:N,1:N,id)=transpose(conjg(Az(1:n,1:n,id)))

        if(fpm(11)==0) then
        call ZGETRS( 'C', N, fpm(23), Az(1,1,id), N, IPIVloc(1,id), workc, N, INFOloc )
        if (infoloc/=0) then
           info=-2
           return
        end if
        else

            !call refineSolveDense('C',Az(1,1,id),workc,n,fpm(23),fpm(50),fpm(51))
            if(loop>3) then
                call refineSolveDense('C',cAz(1,1,id),workc,n,fpm(23),fpm(50),fpm(51),0)
            else
                call refineSolveDense('C',cAz(1,1,id),workc,n,fpm(23),fpm(50),fpm(51),0)
            end if
        end if

     case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

        if ((UPLO=='F').or.(UPLO=='f')) then
           call ZGEMM('N','N',N,fpm(25),N,ONEC,A,LDA,X(1,fpm(24)),N,ZEROC,work(1,fpm(24)),N)
        else
           call ZHEMM ('L', UPLO, N, fpm(25), ONEC, A, LDA, X(1,fpm(24)), N, ZEROC,work(1,fpm(24)), N)
        endif

     case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

        if (LDB==-1) then ! standard
           call ZLACPY( 'F', N, fpm(25),X(1,fpm(24)) , N, work(1,fpm(24)), N )
        else
           if ((UPLO=='F').or.(UPLO=='f')) then
              call ZGEMM('N','N',N,fpm(25),N,ONEC,B,LDB,X(1,fpm(24)),N,ZEROC,work(1,fpm(24)),N)
           else
              call ZHEMM ('L', UPLO, N, fpm(25), ONEC, B, LDB, X(1,fpm(24)), N, ZEROC,work(1,fpm(24)), N)
           endif
        end if


     end select
  end do

  call wdeallocate_2z(zAq)
  call wdeallocate_2z(zSq)
  call wdeallocate_2z(work)
  call wdeallocate_2z(workc)
  call wdeallocate_3z(Az)
  call wdeallocate_2i(ipivloc)
  call wdeallocate_1z(ztmp)


end subroutine zfeast_hegvx







subroutine dfeast_gegvx(N,A,LDA,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the generalized Ax=eBx eigenvalue problem
  !  
  !  A,B GENERAL REAL :: DENSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !
  !  N          (input)        INTEGER: Size system
  !  A          (input)        REAL DOUBLE PRECISION (LDA,N):  Matrix A 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=N)
  !  B          (input)        REAL DOUBLE PRECISION (LDB,N):  Matrix B 
  !  LDB        (input)        INTEGER: Leading dimension of matrix B (LDB>=N)
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
  integer :: N,LDA,LDB
  double precision,dimension(LDA,*):: A
  double precision,dimension(LDB,*):: B
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
  complex(kind=(kind(1.0d0))), dimension(:,:,:),pointer :: Az
  integer, dimension(:,:),pointer ::ipivloc
  double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
  complex(kind=(kind(1.0d0))),parameter :: ONEC=(DONE,DZERO),ZEROC=(DZERO,DZERO)
  logical :: fact
  integer :: nfact,id
  integer :: rank,code,nb_procs,NEW_COMM_WORLD

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
  ELSE IF(LDA<N ) THEN
     INFO = -104
  ELSE IF((LDB<N ).and.(LDB/=-1)) THEN
     INFO = -106
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'DFEAST_GEGV', -INFO+100 )
     RETURN
  END IF
  infoloc=0

  call wallocate_2z(zAq,M0,M0,infoloc)
  call wallocate_2z(zSq,M0,M0,infoloc)
  call wallocate_2z(work,N,2*M0,infoloc)
  call wallocate_2z(workc,N,M0,infoloc)

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
!!!!!!!!!!!!!!!!!!!!!!!!!!
  call wallocate_3z(Az,N,N,nfact,infoloc)
  call wallocate_2i(ipivloc,N,nfact,infoloc)
  if (infoloc/=0) then
     info=-1
     return
  end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ijob=-1 ! initialization
  do while (ijob/=0) 
     call dfeast_grcix(ijob,N,Ze,work,workc,zAq,zSq,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

     select case(ijob)
     case(10) !! factorize (zeB-A)

        if (fpm(10)==1) then
           id=id+1
           if (id==nfact+1) then
              id=1
              fact=.false.
           endif
        endif

        if (fact) then
           if (LDB==-1) then !! standard
              Az(1:N,1:N,id)=-A(1:N,1:N)*ONEC
              do i=1,N
                 Az(i,i,id)=Az(i,i,id)+Ze
              enddo
           else     
              Az(1:N,1:N,id)=Ze*B(1:N,1:N)-A(1:N,1:N)*ONEC
           end if


           call ZGETRF(N,N,Az(1,1,id),N,IPIVloc(1,id),INFOloc)     
           if (infoloc/=0) then
              info=-2
              return
           end if
        end if ! fact true


     case(11) !!solve the linear system (ZeB-A)x=workc(1:N,1:fpm(23)) result in to workc
        call ZGETRS( 'N', N, fpm(23), Az(1,1,id), N, IPIVloc(1,id), workc, N, INFOloc )
        if (infoloc/=0) then
           info=-2
           return
        end if

     case(21) !!solve the linear system (ZeB-A)x=workc(1:N,1:fpm(23)) result in to workc

        call ZGETRS( 'C', N, fpm(23), Az(1,1,id), N, IPIVloc(1,id), workc, N, INFOloc )
        if (infoloc/=0) then
           info=-2
           return
        end if

     case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

        call ZGEMM('N','N',N,fpm(25),N,ONEC,A(1:N,1:N)*ONEC,LDA,X(1,fpm(24)),N,ZEROC,work(1,fpm(24)),N)
        ! With MKL-blas- one should use 
        !call DZGEMM('N','N',N,fpm(25),N,ONEC,A(1:N,1:N),LDA,X(1,fpm(24)),N,ZEROC,work(1,fpm(24)),N)


     case(31) !! perform multiplication A^T*x(1:N,fpm(34):fpm(34)+fpm(35)-1) result in work(1:N,fpm(34)+fpm(35)-1)

        call ZGEMM('C','N',N,fpm(35),N,ONEC,A(1:N,1:N)*ONEC,LDA,X(1,fpm(34)),N,ZEROC,work(1,fpm(34)),N)
        ! With MKL-blas- one should use
        ! call DZGEMM('T','N',N,fpm(35),N,ONEC,A(1:N,1:N),LDA,X(1,fpm(34)),N,ZEROC,work(1,fpm(34)),N)

     case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

        if (LDB==-1) then ! standard
           call ZLACPY( 'F', N, fpm(25),X(1,fpm(24)),N,work(1,fpm(24)),N)
        else
           call ZGEMM('N','N',N,fpm(25),N,ONEC,B(1:N,1:N)*ONEC,LDB,X(1,fpm(24)),N,ZEROC,work(1,fpm(24)),N)
           ! With MKL-blas- one should use
           ! call DZGEMM('N','N',N,fpm(25),N,ONEC,B(1:N,1:N),LDB,X(1,fpm(24)),N,ZEROC,work(1,fpm(24)),N)
        end if

     case(41)!! perform multiplication B^T*x(1:N,fpm(34):fpm(34)+fpm(35)-1) result in work(1:N,fpm(34)+fpm(35)-1)

        if (LDB==-1) then ! standard
           call ZLACPY( 'F', N, fpm(35),X(1,fpm(34)),N,work(1,fpm(34)),N)
        else
           call ZGEMM('C','N',N,fpm(35),N,ONEC,B(1:N,1:N)*ONEC,LDB,X(1,fpm(34)),N,ZEROC,work(1,fpm(34)),N)
           ! With MKL-blas- one should use
           !call DZGEMM('T','N',N,fpm(35),N,ONEC,B(1:N,1:N),LDB,X(1,fpm(34)),N,ZEROC,work(1,fpm(34)),N)
        end if

     end select
  end do

  call wdeallocate_2z(zAq)
  call wdeallocate_2z(zSq)
  call wdeallocate_2z(work)
  call wdeallocate_2z(workc)
  call wdeallocate_3z(Az)
  call wdeallocate_2i(ipivloc)

end subroutine dfeast_gegvx




subroutine zfeast_gegvx(N,A,LDA,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the generalized Ax=eBx eigenvalue problem
  !  
  !  A,B GENERAL COMPLEX :: DENSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !
  !  N          (input)        INTEGER: Size system
  !  A          (input)        COMPLEX DOUBLE PRECISION (LDA,N):  Matrix A 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=N)
  !  B          (input)        COMPLEX DOUBLE PRECISION (LDB,N):  Matrix B 
  !  LDB        (input)        INTEGER: Leading dimension of matrix B (LDB>=N)
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
  integer :: N,LDA,LDB
  complex(kind=(kind(1.0d0))),dimension(LDA,*):: A
  complex(kind=(kind(1.0d0))),dimension(LDB,*):: B
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
  complex(kind=(kind(1.0d0))), dimension(:,:,:),pointer :: Az
  integer, dimension(:,:),pointer ::ipivloc
  double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
  complex(kind=(kind(1.0d0))),parameter :: ONEC=(DONE,DZERO),ZEROC=(DZERO,DZERO)
  logical :: fact
  integer :: nfact,id
  integer :: rank,code,nb_procs,NEW_COMM_WORLD

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
  ELSE IF(LDA<N ) THEN
     INFO = -104
  ELSE IF((LDB<N ).and.(LDB/=-1)) THEN
     INFO = -106
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'ZFEAST_GEGV', -INFO+100 )
     RETURN
  END IF
  infoloc=0

  call wallocate_2z(zAq,M0,M0,infoloc)
  call wallocate_2z(zSq,M0,M0,infoloc)
  call wallocate_2z(work,N,2*M0,infoloc)
  call wallocate_2z(workc,N,M0,infoloc)

  if (infoloc/=0) then
     info=-1
     return
  end if


!!!!!!!!!!!!!!!! find fpm(22) - # total number of contour points
  fpm(22)=fpm(8) ! default


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
!!!!!!!!!!!!!!!!!!!!!!!!!!
  call wallocate_3z(Az,N,N,nfact,infoloc)
  call wallocate_2i(ipivloc,N,nfact,infoloc)
  if (infoloc/=0) then
     info=-1
     return
  end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ijob=-1 ! initialization
  do while (ijob/=0) 
     call zfeast_grcix(ijob,N,Ze,work,workc,zAq,zSq,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

     select case(ijob)
     case(10) !! factorize (zeB-A)

        if (fpm(10)==1) then
           id=id+1
           if (id==nfact+1) then
              id=1
              fact=.false.
           endif
        endif

        if (fact) then
           if (LDB==-1) then !! standard
              Az(1:N,1:N,id)=-A(1:N,1:N)
              do i=1,N
                 Az(i,i,id)=Az(i,i,id)+Ze
              enddo
           else     
              Az(1:N,1:N,id)=Ze*B(1:N,1:N)-A(1:N,1:N)
           end if


           call ZGETRF(N,N,Az(1,1,id),N,IPIVloc(1,id),INFOloc)     
           if (infoloc/=0) then
              info=-2
              return
           end if
        end if ! fact true


     case(11) !!solve the linear system (ZeB-A)x=workc(1:N,1:fpm(23)) result in to workc
        call ZGETRS( 'N', N, fpm(23), Az(1,1,id), N, IPIVloc(1,id), workc, N, INFOloc )
        if (infoloc/=0) then
           info=-2
           return
        end if

     case(21) !!solve the linear system (ZeB-A)x=workc(1:N,1:fpm(23)) result in to workc

        call ZGETRS( 'C', N, fpm(23), Az(1,1,id), N, IPIVloc(1,id), workc, N, INFOloc )
        if (infoloc/=0) then
           info=-2
           return
        end if

     case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

        call ZGEMM('N','N',N,fpm(25),N,ONEC,A(1:N,1:N),LDA,X(1,fpm(24)),N,ZEROC,work(1,fpm(24)),N)

     case(31) !! perform multiplication A^H*x(1:N,fpm(34):fpm(34)+fpm(35)-1) result in work(1:N,fpm(34)+fpm(35)-1)

        call ZGEMM('C','N',N,fpm(35),N,ONEC,A(1:N,1:N),LDA,X(1,fpm(34)),N,ZEROC,work(1,fpm(34)),N)

     case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

        if (LDB==-1) then ! standard
           call ZLACPY( 'F', N, fpm(25),X(1,fpm(24)),N,work(1,fpm(24)),N)
        else
           call ZGEMM('N','N',N,fpm(25),N,ONEC,B(1:N,1:N),LDB,X(1,fpm(24)),N,ZEROC,work(1,fpm(24)),N)
        end if

     case(41)!! perform multiplication B^H*x(1:N,fpm(34):fpm(34)+fpm(35)-1) result in work(1:N,fpm(34)+fpm(35)-1)

        if (LDB==-1) then ! standard
           call ZLACPY( 'F', N, fpm(35),X(1,fpm(34)),N,work(1,fpm(34)),N)
        else
           call ZGEMM('C','N',N,fpm(35),N,ONEC,B(1:N,1:N),LDB,X(1,fpm(34)),N,ZEROC,work(1,fpm(34)),N)
        end if

     end select
  end do

  call wdeallocate_2z(zAq)
  call wdeallocate_2z(zSq)
  call wdeallocate_2z(work)
  call wdeallocate_2z(workc)
  call wdeallocate_3z(Az)
  call wdeallocate_2i(ipivloc)


end subroutine zfeast_gegvx







subroutine zfeast_sygvx(UPLO,N,A,LDA,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the generalized Ax=eBx eigenvalue problem
  !  
  !  A,B COMPLEX SYMMETRIC :: DENSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
  !                           triangular part of the matrix(ces) is(are) being supplied.
  !  N          (input)        INTEGER: Size system
  !  A          (input)        COMPLEX DOUBLE PRECISION (LDA,N):  Matrix A 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=N)
  !  B          (input)        COMPLEX DOUBLE PRECISION (LDB,N):  Matrix B 
  !  LDB        (input)        INTEGER: Leading dimension of matrix B (LDB>=N)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  fpm (input/output) INTEGER(*) : FEAST parameters
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
  !  q          (input/output) COMPLEX DOUBLE PRECISION(N,M0) : 
  !                                                       On entry: subspace initial guess if fpm(5)=1 
  !                                                       On exit : Right Eigenvectors-solution
  !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
  !  res        (output)       REAL DOUBLE PRECISION(M0) : Relative Residual of the solution (1-norm)
  !                                                        if option fpm(6)=1 selected                           
  !  info       (output)       INTEGER: Error handling (0: successful exit)
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of ARGUMENTS FOR EXPERT ROUTINE 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Zne        (input)        COMPLEX DOUBLE PRECISION(fpm(8)): Custom Integration nodes
  !  Wne        (input)        COMPLEX DOUBLE PRECISION(fpm(8)): Custom Integration weights
  !           
  !=====================================================================
  ! James Kestyn, Eric Polizzi 2015
  !=====================================================================
  implicit none
  include 'f90_noruntime_interface.fi'
  character(len=1) :: UPLO
  integer :: N,LDA,LDB
  complex(kind=(kind(1.0d0))),dimension(LDA,*):: A
  complex(kind=(kind(1.0d0))),dimension(LDB,*):: B
  integer,dimension(*) :: fpm
  double precision :: epsout 
  integer :: loop
  complex(kind=(kind(1.0d0))) :: Emid
  double precision :: r
  integer :: M0
  complex(kind=(kind(1.0d0))),dimension(*)  :: E
  complex(kind=(kind(1.0d0))),dimension(N,*):: X
  integer :: mode
  double precision,dimension(*)    :: res
  integer :: info
  complex(kind=(kind(1.0d0))),dimension(*) :: Zne,Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer :: ijob,infoloc,i,s
  complex(kind=(kind(1.0d0))) :: Ze
  complex(kind=(kind(1.0d0))), dimension(:,:),pointer ::workc,work,zAq,zSq
  complex(kind=(kind(1.0d0))), dimension(:,:,:),pointer :: Az
  integer, dimension(:,:),pointer ::ipivloc
  double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
  complex(kind=(kind(1.0d0))),parameter :: ONEC=(DONE,DZERO), ZEROC=(DZERO,DZERO)
  logical :: fact
  integer :: nfact,id
  character(len=1) :: UPLO2
  integer :: rank,code,nb_procs,NEW_COMM_WORLD

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
  ELSE IF(LDA<N ) THEN
     INFO = -104
  ELSE IF((LDB<N ).and.(LDB/=-1)) THEN
     INFO = -106
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'ZFEAST_SYGV', -INFO+100 )
     RETURN
  END IF
  infoloc=0

  call wallocate_2z(zAq,M0,M0,infoloc)
  call wallocate_2z(zSq,M0,M0,infoloc)
  call wallocate_2z(work,N,M0,infoloc)
  call wallocate_2z(workc,N,M0,infoloc)

  if (infoloc/=0) then
     info=-1
     return
  end if

!!! Format conversion
  UPLO2=UPLO
  IF ((UPLO=='F').or.(UPLO=='f')) UPLO2='L'
  IF ((UPLO=='L').or.(UPLO=='l')) UPLO2='L'
  IF ((UPLO=='U').or.(UPLO=='u')) UPLO2='U'



!!!!!!!!!!!!!!!! find fpm(22) - # total number of contour points
  fpm(22)=fpm(8) ! default
 

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
!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call wallocate_3z(Az,N,N,nfact,infoloc)
  call wallocate_2i(ipivloc,N,nfact,infoloc)
  if (infoloc/=0) then
     info=-1
     return
  end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ijob=-1 ! initialization 
  do while (ijob/=0)
     call zfeast_srcix(ijob,N,Ze,work,workc,zAq,zSq,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)    

     select case(ijob)
     case(10) !! factorize (zeB-A)

        if (fpm(10)==1) then
           id=id+1
           if (id==nfact+1) then
              id=1
              fact=.false.
           endif
        endif

        if (fact) then
           if (LDB==-1) then !! standard
              Az(1:N,1:N,id)=-A(1:N,1:N)
              do i=1,N
                 Az(i,i,id)=Az(i,i,id)+Ze
              enddo
           else        
              Az(1:N,1:N,id)=Ze*B(1:N,1:N)-A(1:N,1:N) 
           end if

            if (fpm(11)==0) then
           call ZSYTRF(UPLO2,N,Az(1,1,id),N,IPIVloc(1,id),workc,N*M0,INFOloc)     
           if (infoloc/=0) then
              info=-2
              return
           end if
           end if
        endif ! fact true

     case(11) !!solve the linear system (ZeB-A)x=workc(1:N,1:fpm(23)) result in to workc

        if( fpm(11)==0) then
        call ZSYTRS( UPLO2, N, fpm(23), Az(1,1,id), N, IPIVloc(1,id), workc, N, INFOloc )
        if (infoloc/=0) then
           info=-2
           return
        end if
        else
                call refineSolveDense('N',Az(1,1,id),workc,n,fpm(23),fpm(50),fpm(51),0)
        end if

     case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

        call ZSYMM ('L', UPLO2, N, fpm(25), ONEC, A, LDA, X(1,fpm(24)), N, ZEROC,work(1,fpm(24)), N)

     case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

        if (LDB==-1) then ! standard
           call ZLACPY( 'F', N, fpm(25),X(1,fpm(24)) , N, work(1,fpm(24)), N )
        else
           call ZSYMM ('L', UPLO2, N, fpm(25), ONEC, B, LDB, X(1,fpm(24)), N, ZEROC,work(1,fpm(24)), N)
        end if

     end select
  end do



  call wdeallocate_2z(zAq)
  call wdeallocate_2z(zSq)
  call wdeallocate_2z(work)
  call wdeallocate_2z(workc)
  call wdeallocate_3z(Az)
  call wdeallocate_2i(ipivloc)

end subroutine zfeast_sygvx







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!22222222222222222222222222222222222222222222222222222222222
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



subroutine dfeast_sygv(UPLO,N,A,LDA,B,LDB,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the generalized Ax=eBx eigenvalue problem
  !  
  !  A REAL SYMMETRIC, B SYMMETRIC POSITIVE DEFINITE:: DENSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
  !                           triangular part of the matrix(ces) is(are) being supplied.
  !  N          (input)        INTEGER: Size system
  !  A          (input)        REAL DOUBLE PRECISION (LDA,N):  Matrix A 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=N)
  !  B          (input)        REAL DOUBLE PRECISION (LDB,N):  Matrix B 
  !  LDB        (input)        INTEGER: Leading dimension of matrix B (LDB>=N)
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
  integer :: N,LDA,LDB
  double precision,dimension(LDA,*):: A
  double precision,dimension(LDB,*):: B
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
!!!    Wrapper Routine to expert routine: dfeast_sygvx
!!!    Create Zne, Wne arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  complex(kind=(kind(1.0d0))),dimension(fpm(2)) :: Zne,Wne 

  call zfeast_contour(Emin,Emax,fpm(2),fpm(16),fpm(18),Zne,Wne)

  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  call dfeast_sygvx(UPLO,N,A,LDA,B,LDB,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)

end subroutine dfeast_sygv




subroutine zfeast_hegv(UPLO,N,A,LDA,B,LDB,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the generalized Ax=eBx eigenvalue problem
  !  
  !  A COMPLEX HERMITIAN, B HERMITIAN POSITIVE DEFINITE:: DENSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
  !                           triangular part of the matrix(ces) is(are) being supplied.
  !  N          (input)        INTEGER: Size system
  !  A          (input)        COMPLEX DOUBLE PRECISION (LDA,N):  Matrix A 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=N)
  !  B          (input)        COMPLEX DOUBLE PRECISION (LDB,N):  Matrix B 
  !  LDB        (input)        INTEGER: Leading dimension of matrix B (LDB>=N)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  fpm (input/output) INTEGER(*) : FEAST parameters
  !  epsout     (output)       REAL DOUBLE PRECISION : Error on the trace
  !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
  !  Emin,Emax  (input)        REAL DOUBLE PRECISION: search interval
  !  M0         (input/output) INTEGER: Size subspace
  !  E          (output)       REAL DOUBLE PRECISION(M0)   : Eigenvalues -solution
  !  q          (input/output) COMPLEX DOUBLE PRECISION(N,M0) : 
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
  integer :: N,LDA,LDB
  complex(kind=(kind(1.0d0))),dimension(LDA,*):: A
  complex(kind=(kind(1.0d0))),dimension(LDB,*):: B
  integer,dimension(*) :: fpm
  double precision :: epsout 
  integer :: loop
  double precision :: Emin,Emax
  integer :: M0
  double precision,dimension(*)  :: E
  complex(kind=(kind(1.0d0))),dimension(N,*):: X
  integer :: mode
  double precision,dimension(*)    :: res
  integer :: info

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!    Wrapper Routine to expert routine: zfeast_hegvx
!!!    Create Zne, Wne arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  complex(kind=(kind(1.0d0))),dimension(fpm(2)) :: Zne,Wne 

  call zfeast_contour(Emin,Emax,fpm(2),fpm(16),fpm(18),Zne,Wne)

  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  call zfeast_hegvx(UPLO,N,A,LDA,B,LDB,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)


end subroutine zfeast_hegv



subroutine dfeast_gegv(N,A,LDA,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the generalized Ax=eBx eigenvalue problem
  !  
  !  A,B GENERAL REAL:: DENSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !
  !  N          (input)        INTEGER: Size system
  !  A          (input)        REAL DOUBLE PRECISION (LDA,N):  Matrix A 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=N)
  !  B          (input)        REAL DOUBLE PRECISION (LDB,N):  Matrix B 
  !  LDB        (input)        INTEGER: Leading dimension of matrix B (LDB>=N)
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
  integer :: N,LDA,LDB
  double precision,dimension(LDA,*):: A
  double precision,dimension(LDB,*):: B
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
!!!    Wrapper Routine to expert routine: dfeast_gegvx
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  complex(kind=(kind(1.0d0))), dimension(1:fpm(8)) :: Zne,Wne 

  call zfeast_gcontour(Emid,r,fpm(8),fpm(17),fpm(18),fpm(19),Zne,Wne)

  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  call dfeast_gegvx(N,A,LDA,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

end subroutine dfeast_gegv




subroutine zfeast_gegv(N,A,LDA,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the generalized Ax=eBx eigenvalue problem
  !  
  !  A,B GENERAL COMPLEX :: DENSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !
  !  N          (input)        INTEGER: Size system
  !  A          (input)        COMPLEX DOUBLE PRECISION (LDA,N):  Matrix A 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=N)
  !  B          (input)        COMPLEX DOUBLE PRECISION (LDB,N):  Matrix B 
  !  LDB        (input)        INTEGER: Leading dimension of matrix B (LDB>=N)
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
  integer :: N,LDA,LDB
  complex(kind=(kind(1.0d0))),dimension(LDA,*):: A
  complex(kind=(kind(1.0d0))),dimension(LDB,*):: B
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
!!!    Wrapper Routine to expert routine: zfeast_gegvx
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  complex(kind=(kind(1.0d0))), dimension(1:fpm(8)) :: Zne,Wne 

  call zfeast_gcontour(Emid,r,fpm(8),fpm(17),fpm(18),fpm(19),Zne,Wne)

  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  call zfeast_gegvx(N,A,LDA,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

end subroutine zfeast_gegv




subroutine zfeast_sygv(UPLO,N,A,LDA,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the generalized Ax=eBx eigenvalue problem
  !  
  !  A,B GENERAL COMPLEX SYMMETRIC:: DENSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
  !                           triangular part of the matrix(ces) is(are) being supplied.
  !                           Rq: It is a dummy parameter here since matrix is general
  !
  !  N          (input)        INTEGER: Size system
  !  A          (input)        COMPLEX DOUBLE PRECISION (LDA,N):  Matrix A 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=N)
  !  B          (input)        COMPLEX DOUBLE PRECISION (LDB,N):  Matrix B 
  !  LDB        (input)        INTEGER: Leading dimension of matrix B (LDB>=N)
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
  !  q          (input/output) COMPLEX DOUBLE PRECISION(N,M0) : 
  !                                                       On entry: subspace initial guess if fpm(5)=1 
  !                                                       On exit : Right Eigenvectors-solution 
  !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
  !  res        (output)       REAL DOUBLE PRECISION(M0) : Relative Residual of the solution (1-norm)
  !                                                        if option fpm(6)=1 selected                           
  !  info       (output)       INTEGER: Error handling (0: successful exit)
  !=====================================================================
  ! James Kestyn, Eric Polizzi 2015
  ! ====================================================================
  implicit none
  character(len=1) :: UPLO
  integer :: N,LDA,LDB
  complex(kind=(kind(1.0d0))),dimension(LDA,*):: A
  complex(kind=(kind(1.0d0))),dimension(LDB,*):: B
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
!!!    Wrapper Routine to expert routine: zfeast_sygvx
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  complex(kind=(kind(1.0d0))), dimension(1:fpm(8)) :: Zne,Wne 

  call zfeast_gcontour(Emid,r,fpm(8),fpm(17),fpm(18),fpm(19),Zne,Wne)

  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  call zfeast_sygvx(UPLO,N,A,LDA,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

end subroutine zfeast_sygv





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!33333333333333333333333333333333333333333333333333333333333
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





subroutine dfeast_syevx(UPLO,N,A,LDA,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the standard Ax=ex eigenvalue problem
  !  
  !  A REAL SYMMETRIC :: DENSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
  !                           triangular part of the matrix(ces) is(are) being supplied.
  !  N          (input)        INTEGER: Size system
  !  A          (input)        REAL DOUBLE PRECISION (LDA,N):  Matrix A 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=N)
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
  integer :: N,LDA
  double precision,dimension(LDA,*):: A
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
  double precision,dimension(1) :: B ! dummy
  integer :: LDB

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!    Wrapper Routine to expert routine: dfeast_sygvx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  LDB=-1 ! B is a dummy- option for standard eigenvalue problem
  call dfeast_sygvx(UPLO,N,A,LDA,B,LDB,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)

end subroutine dfeast_syevx





subroutine zfeast_heevx(UPLO,N,A,LDA,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the standard Ax=ex eigenvalue problem
  !  
  !  A COMPLEX HERMITIAN :: DENSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
  !                           triangular part of the matrix(ces) is(are) being supplied.
  !  N          (input)        INTEGER: Size system
  !  A          (input)        COMPLEX DOUBLE PRECISION (LDA,N):  Matrix A 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=N)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  fpm (input/output) INTEGER(*) : FEAST parameters
  !  epsout     (output)       REAL DOUBLE PRECISION : Error on the trace
  !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
  !  Emin,Emax  (input)        REAL DOUBLE PRECISION: search interval
  !  M0         (input/output) INTEGER: Size subspace
  !  E          (output)       REAL DOUBLE PRECISION(M0)   : Eigenvalues -solution
  !  q          (input/output) COMPLEX DOUBLE PRECISION(N,M0) : 
  !                                                       On entry: subspace initial guess if fpm(5)=1 
  !                                                       On exit : Eigenvectors-solution
  !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
  !  res        (output)       REAL DOUBLE PRECISION(M0) : Relative Residual of the solution (1-norm)
  !                                                        if option fpm(6)=1 selected                           
  !  info       (output)       INTEGER: Error handling (0: successful exit)
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
  integer :: N,LDA
  complex(kind=(kind(1.0d0))),dimension(LDA,*):: A
  integer,dimension(*) :: fpm
  double precision :: epsout 
  integer :: loop
  double precision :: Emin,Emax
  integer :: M0
  double precision,dimension(*)  :: E
  complex(kind=(kind(1.0d0))),dimension(N,*):: X
  integer :: mode
  double precision,dimension(*)    :: res
  integer :: info
  complex(kind=(kind(1.0d0))),dimension(*) :: Zne,Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!    Wrapper Routine to expert routine: zfeast_hegvx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  complex(kind=(kind(1.0d0))),dimension(1) :: B ! dummy
  integer :: LDB
  LDB=-1 ! B is a dummy- option for standard eigenvalue problem
  call zfeast_hegvx(UPLO,N,A,LDA,B,LDB,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)

end subroutine zfeast_heevx



subroutine dfeast_geevx(N,A,LDA,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the standard Ax=ex eigenvalue problem
  !  
  !  A GENERAL REAL :: DENSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  N          (input)        INTEGER: Size system
  !  A          (input)        REAL DOUBLE PRECISION (LDA,N):  Matrix A 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=N)
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
  integer :: N,LDA
  double precision,dimension(LDA,*):: A
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!    Wrapper Routine to expert routine: dfeast_gegvx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  double precision,dimension(1) :: B ! dummy
  integer :: LDB
  LDB=-1 ! B is a dummy- option for standard eigenvalue problem
  call dfeast_gegvx(N,A,LDA,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)


end subroutine dfeast_geevx




subroutine zfeast_geevx(N,A,LDA,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the standard Ax=ex eigenvalue problem
  !  
  !  A GENERAL COMPLEX :: DENSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !
  !  N          (input)        INTEGER: Size system
  !  A          (input)        COMPLEX DOUBLE PRECISION (LDA,N):  Matrix A 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=N)
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
  integer :: N,LDA
  complex(kind=(kind(1.0d0))),dimension(LDA,*):: A
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!    Wrapper Routine to expert routine: zfeast_gegvx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  complex(kind=(kind(1.0d0))),dimension(1) :: B ! dummy
  integer :: LDB
  LDB=-1 ! B is a dummy- option for standard eigenvalue problem
  call zfeast_gegvx(N,A,LDA,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)


end subroutine zfeast_geevx



subroutine zfeast_syevx(UPLO,N,A,LDA,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the standard Ax=ex eigenvalue problem
  !  
  !  A GENERAL COMPLEX SYMMETRIC:: DENSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
  !                           triangular part of the matrix(ces) is(are) being supplied.
  !                           Rq: It is a dummy parameter here since matrix is general
  !
  !  N          (input)        INTEGER: Size system
  !  A          (input)        COMPLEX DOUBLE PRECISION (LDA,N):  Matrix A 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=N)
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
  !  q          (input/output) COMPLEX DOUBLE PRECISION(N,M0) : 
  !                                                       On entry: subspace initial guess if fpm(5)=1 
  !                                                       On exit : Right Eigenvectors-solution 
  !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
  !  res        (output)       REAL DOUBLE PRECISION(M0) : Relative Residual of the solution (1-norm)
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
  character(len=1) :: UPLO
  integer :: N,LDA
  complex(kind=(kind(1.0d0))),dimension(LDA,*):: A
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!    Wrapper Routine to expert routine: zfeast_sygvx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  complex(kind=(kind(1.0d0))),dimension(1) :: B ! dummy
  integer :: LDB
  LDB=-1 ! B is a dummy- option for standard eigenvalue problem
  call zfeast_sygvx(UPLO,N,A,LDA,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)


end subroutine zfeast_syevx




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!44444444444444444444444444444444444444444444444444444444444
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





subroutine dfeast_syev(UPLO,N,A,LDA,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the standard Ax=ex eigenvalue problem
  !  
  !  A REAL SYMMETRIC:: DENSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
  !                           triangular part of the matrix(ces) is(are) being supplied.
  !  N          (input)        INTEGER: Size system
  !  A          (input)        REAL DOUBLE PRECISION (LDA,N):  Matrix A 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=N)
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
  integer :: N,LDA
  double precision,dimension(LDA,*):: A
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
!!!    Wrapper Routine to expert routine: dfeast_sygvx
!!!    Create Zne, Wne arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  complex(kind=(kind(1.0d0))),dimension(fpm(2)) :: Zne,Wne 
  double precision,dimension(1) :: B ! dummy
  integer :: LDB

  call zfeast_contour(Emin,Emax,fpm(2),fpm(16),fpm(18),Zne,Wne)

  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  LDB=-1 ! B is a dummy- option for standard eigenvalue problem
  call dfeast_sygvx(UPLO,N,A,LDA,B,LDB,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)

end subroutine dfeast_syev








subroutine zfeast_heev(UPLO,N,A,LDA,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the standard Ax=ex eigenvalue problem
  !  
  !  A COMPLEX HERMITIAN :: DENSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
  !                           triangular part of the matrix(ces) is(are) being supplied.
  !  N          (input)        INTEGER: Size system
  !  A          (input)        COMPLEX DOUBLE PRECISION (LDA,N):  Matrix A 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=N)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  fpm (input/output) INTEGER(*) : FEAST parameters
  !  epsout     (output)       REAL DOUBLE PRECISION : Error on the trace
  !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
  !  Emin,Emax  (input)        REAL DOUBLE PRECISION: search interval
  !  M0         (input/output) INTEGER: Size subspace
  !  E          (output)       REAL DOUBLE PRECISION(M0)   : Eigenvalues -solution
  !  q          (input/output) COMPLEX DOUBLE PRECISION(N,M0) : 
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
  include 'f90_noruntime_interface.fi'
  character(len=1) :: UPLO
  integer :: N,LDA
  complex(kind=(kind(1.0d0))),dimension(LDA,*):: A
  integer,dimension(*) :: fpm
  double precision :: epsout 
  integer :: loop
  double precision :: Emin,Emax
  integer :: M0
  double precision,dimension(*)  :: E
  complex(kind=(kind(1.0d0))),dimension(N,*):: X
  integer :: mode
  double precision,dimension(*)    :: res
  integer :: info
!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!    Wrapper Routine to expert routine: zfeast_hegvx
!!!    Create Zne, Wne arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  complex(kind=(kind(1.0d0))),dimension(fpm(2)) :: Zne,Wne 
  complex(kind=(kind(1.0d0))),dimension(1) :: B ! dummy
  integer :: LDB

  call zfeast_contour(Emin,Emax,fpm(2),fpm(16),fpm(18),Zne,Wne)

  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  LDB=-1 ! B is a dummy- option for standard eigenvalue problem
  call zfeast_hegvx(UPLO,N,A,LDA,B,LDB,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)

end subroutine zfeast_heev




subroutine dfeast_geev(N,A,LDA,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the standard Ax=ex eigenvalue problem
  !  
  !  A GENERAL REAL :: DENSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !
  !  N          (input)        INTEGER: Size system
  !  A          (input)        REAL DOUBLE PRECISION (LDA,N):  Matrix A 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=N)
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
  double precision,dimension(LDA,*):: A
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
!!!    Wrapper Routine to expert routine: dfeast_gegvx
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  complex(kind=(kind(1.0d0))), dimension(1:fpm(8)) :: Zne,Wne 
  double precision,dimension(1) :: B ! dummy
  integer :: LDB

  call zfeast_gcontour(Emid,r,fpm(8),fpm(17),fpm(18),fpm(19),Zne,Wne)

  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  LDB=-1 ! B is a dummy- option for standard eigenvalue problem
  call dfeast_gegvx(N,A,LDA,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

end subroutine dfeast_geev




subroutine zfeast_geev(N,A,LDA,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the standard Ax=ex eigenvalue problem
  !  
  !  A GENERAL COMPLEX :: DENSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  N          (input)        INTEGER: Size system
  !  A          (input)        COMPLEX DOUBLE PRECISION (LDA,N):  Matrix A 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=N)
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
  complex(kind=(kind(1.0d0))),dimension(LDA,*):: A
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
!!!    Wrapper Routine to expert routine: zfeast_gegvx
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  complex(kind=(kind(1.0d0))), dimension(1:fpm(8)) :: Zne,Wne 
  complex(kind=(kind(1.0d0))),dimension(1) :: B ! dummy
  integer :: LDB

  call zfeast_gcontour(Emid,r,fpm(8),fpm(17),fpm(18),fpm(19),Zne,Wne)

  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  LDB=-1 ! B is a dummy- option for standard eigenvalue problem
  call zfeast_gegvx(N,A,LDA,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

end subroutine zfeast_geev




subroutine zfeast_syev(UPLO,N,A,LDA,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the standard Ax=ex eigenvalue problem
  !  
  !  A GENERAL COMPLEX SYMMETRIC:: DENSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
  !                           triangular part of the matrix(ces) is(are) being supplied.
  !                           Rq: It is a dummy parameter here since matrix is general
  !
  !  N          (input)        INTEGER: Size system
  !  A          (input)        COMPLEX DOUBLE PRECISION (LDA,N):  Matrix A 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=N)
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
  !  q          (input/output) COMPLEX DOUBLE PRECISION(N,M0) : 
  !                                                       On entry: subspace initial guess if fpm(5)=1 
  !                                                       On exit : Right Eigenvectors-solution 
  !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
  !  res        (output)       REAL DOUBLE PRECISION(M0) : Relative Residual of the solution (1-norm)
  !                                                        if option fpm(6)=1 selected                           
  !  info       (output)       INTEGER: Error handling (0: successful exit)
  !=====================================================================
  ! James Kestyn, Eric Polizzi 2015
  ! ====================================================================
  implicit none
  character(len=1) :: UPLO
  integer :: N,LDA
  complex(kind=(kind(1.0d0))),dimension(LDA,*):: A
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
!!!    Wrapper Routine to expert routine: zfeast_sygvx
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  complex(kind=(kind(1.0d0))), dimension(1:fpm(8)) :: Zne,Wne 
  complex(kind=(kind(1.0d0))),dimension(1) :: B ! dummy
  integer :: LDB

  call zfeast_gcontour(Emid,r,fpm(8),fpm(17),fpm(18),fpm(19),Zne,Wne)

  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  LDB=-1 ! B is a dummy- option for standard eigenvalue problem
  call zfeast_sygvx(UPLO,N,A,LDA,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

end subroutine zfeast_syev


subroutine refineSolveDense(transA,A,B,n,m,r,m0,usex)
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

