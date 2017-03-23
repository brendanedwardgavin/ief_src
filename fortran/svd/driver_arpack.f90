      program driver_arpack
        implicit none
!     dseupd   ARPACK routine that returns Ritz values and (optionally)
!             Ritz vectors. 
!     daxpy    Level 1 BLAS that computes y <- alpha*x+y.
!     dnrm2   Level 1 BLAS that computes the norm of a complex vector.

!     %-----------------------------%
!     | Define leading dimensions   |
!     | for all arrays.             |
!     | MAXN:   Maximum dimension   |
!     |         of the A allowed.   |
!     | MAXNEV: Maximum NEV allowed |
!     | MAXNCV: Maximum NCV allowed |
!     %-----------------------------%
!     %--------------%
!     | Local Arrays |
!     %--------------%
      integer,dimension(11)  ::    iparam
      integer,dimension(11)  ::    ipntr(11)
      logical,dimension(:),allocatable ::          select
      double precision,allocatable,dimension(:) :: workd,resid,workl,ax,bx
      double precision,allocatable,dimension(:,:) :: v !! eigenvectors
      double precision,allocatable,dimension(:,:) :: d   !! eigenvalues
!     %---------------%
!     | Local Scalars |
!     %---------------%
      character ::        bmat*1, which*2
      integer ::          ido, nev, ncv, lworkl, info, ierr,nconv, maxitr, ishfts, mode
      Double precision :: tol
      logical       ::    rvec 
!     %-----------------------------%
!     | BLAS & LAPACK routines used |
!     %-----------------------------%
      Double precision ::  dnrm2
      external         dnrm2 , daxpy 
      integer :: i,j,k
      character(len=1) :: cc,uplo
      integer :: n,nnza,ncol,nrow,nmin,nmax
      integer,dimension(:), allocatable :: ica,jca,isa,jsa
      double precision,dimension(:),allocatable :: dca,dsa
      character(len=1),dimension(4) ::   matdescra
      double precision :: sigma

        integer :: mmcounter
      integer :: ti,tf,tim
        
        double precision :: resnorm
        
      character(len=100) :: name,nevstr,ncvstr

!     %-----------------------%
!     | Executable statements |
!     %-----------------------%
!     %----------------------------------------------------%
!     | The number N is the dimension of the matrix.  A   |
      !     | a standard eigenvalue problem is solved.
!      NEV is the number of eigenvalues (closest   |
!     | to SIGMAR) to be approximated.  Since the          |
!     | shift-invert mode is used,  WHICH is set to 'LM'.  |
!     | The user can modify NEV, NCV, SIGMA to solve       |
!     | problems of different sizes, and to get different  |
!     | parts of the spectrum.  However, The following     |
!     | conditions must be satisfied:                      |
!     |                     N <= MAXN,                     |
!     |                   NEV <= MAXNEV,                   |
!     |               NEV + 2 <= NCV <= MAXNCV             |
!     %----------------------------------------------------%
        mmcounter=0
      
      if (iargc()<5) then
            print *,'Too few arguments'
            print *,'usage: svd_arpack Matrix_name UPLO Which(LM,LA,SA) N_EV Subspace'
      end if

      call getarg(1,name)
      call getarg(2,uplo)
      call getarg(3,which)
      call getarg(4,nevstr)
      call getarg(5,ncvstr)
      read(nevstr,*) nev
      read(ncvstr,*) ncv      

      !read *,nev ! # eigen
      !read *,ncv ! M0

    print *,'Reading matrix...' 
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
  read(10,*) nrow,ncol,nnza
  allocate(ica(nnza))
  allocate(jca(nnza))
     allocate(dca(nnza))
     do i=1,nnza
        read(10,*) ica(i),jca(i),dca(i)
     end do
  close(10)
    print *,'Done reading. Converting to CSR'
 
  !! create csr format
  allocate(isa(1:nrow+1))
  allocate(jsa(1:nnza))
  allocate(dsa(1:nnza))
  call dcoo2csr(nrow,nnza,ica,jca,dca,isa,jsa,dsa)
 
    print *,'Done converting.'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!! INFORMATION ABOUT MATRIX !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  print *,'matrix name ',trim(name)
  print *,'matrix -coordinate format- size',nrow,ncol
  print *,'sparse matrix A- nnz',nnza

 
  print *, 'ARPACK: nev, nsv',nev,ncv

     
      
      bmat  = 'I'
      !which = 'LA'


      matdescra(1) = 'S'  !SYMMETRIC or GENERAL
      matdescra(2) = UPLO  !UPLO INDICATOR

if (UPLO=="F") matdescra(1) = 'G'  !SYMMETRIC or GENERAL


      matdescra(3) = 'N'  !unit or non-unit diagonal
      matdescra(4) = 'F'  !1 based indexing or 0 based


nmin=min(nrow,ncol)
nmax=max(nrow,ncol)



      
call system_clock(ti,tim)

     

!     %-----------------------------------------------------%
!     | The work array WORKL is used in ZNAUPD  as           |
!     | workspace.  Its dimension LWORKL is set as          |
!     | illustrated below.  The parameter TOL determines    |
!     | the stopping criterion. If TOL<=0, machine          |
!     | precision is used.  The variable IDO is used for    |
!     | reverse communication, and is initially set to 0.   |
!     | Setting INFO=0 indicates that a random vector is    |
!     | generated in ZNAUPD  to start the Arnoldi iteration. |
!     %-----------------------------------------------------%
      lworkl = ncv*(ncv+8) 
      tol    = 1d-10
      ido    = 0
      info   = 0

!!!! allocate
 allocate(workd(1:3*nmin))
      allocate(resid(nmin))
      allocate(workl(lworkl))
      allocate(d(1:ncv,1:2))
allocate(v(1:nmin,1:ncv))
      allocate(select(1:ncv))
allocate(ax(nmax))       
allocate(bx(nmin))       
      
      
      
!     %---------------------------------------------------%
!     | This program uses exact shifts with respect to    |
!     | the current Hessenberg matrix (IPARAM(1) = 1).    |
!     | IPARAM(3) specifies the maximum number of Arnoldi |
!     | iterations allowed. Mode 3 of DSAUPD  is used      |
!     | (IPARAM(7) = 1).  All these options can be        |
!     | changed by the user. For details see the          |
!     | documentation in DSAUPD .                          |
!     %---------------------------------------------------%

      ishfts = 1
      maxitr = 100000
      mode   = 1
      
      iparam(1) = ishfts
      iparam(3) = maxitr 
      iparam(7) = mode


      
!     %------------------------------------------%
!     | M A I N   L O O P(Reverse communication) | 
!     %------------------------------------------%
 20   continue
!        %---------------------------------------------%
!        | Repeatedly call the routine ZNAUPD  and take |
!        | actions indicated by parameter IDO until    |
!        | either convergence is indicated or maxitr   |
!        | has been exceeded.                          |
!        %---------------------------------------------%
           call dsaupd ( ido, bmat, nmin, which, nev, tol, resid,ncv, v, nmin, iparam, ipntr, workd, workl,lworkl, info )

      if (ido .eq. -1 .or. ido .eq. 1) then
            mmcounter=mmcounter+1
            print *,'Matmul!',mod(mmcounter,ncv),mmcounter/ncv
            !print *,'       ',dnrm2(nmin,resid,1)
!         c
!c           %--------------------------------------%
!c           | Perform matrix vector multiplication |
!c           |              y <--- OP*x             |
!c           | The user should supply his/her own   |
!c           | matrix vector multiplication routine |
!c           | here that takes workd(ipntr(1)) as   |
!c           | the input, and return the result to  |
!c           | workd(ipntr(2)).                     |
!c           %--------------------------------------%
!c
!c            call av (nx, workd(ipntr(1)), workd(ipntr(2)))
         !c         print*,'MAIN LOOP MULTIPLICATION'
if (ncol==nrow) then
   call mkl_dcsrmm('N', nrow, 1, nrow, 1.0d0, matdescra, dsa,jsa,isa,isa(2), workd(ipntr(1)), nrow, 0.0d0, workd(ipntr(2)), nrow)
   elseif (nrow>ncol) then
            call mkl_dcsrmm('N', nrow, 1, ncol, 1.0d0, matdescra, dsa,jsa,isa,isa(2), workd(ipntr(1)), ncol, 0.0d0, ax, nrow) 
            call mkl_dcsrmm('T', nrow, 1, ncol, 1.0d0, matdescra, dsa,jsa,isa,isa(2), ax, nrow, 0.0d0, workd(ipntr(2)), ncol)
    elseif (nrow<ncol) then
            call mkl_dcsrmm('T', nrow, 1, ncol, 1.0d0, matdescra, dsa,jsa,isa,isa(2), workd(ipntr(1)), nrow, 0.0d0, ax, ncol) 
            call mkl_dcsrmm('N', nrow, 1, ncol, 1.0d0, matdescra, dsa,jsa,isa,isa(2), ax, ncol, 0.0d0, workd(ipntr(2)), nrow)
            
endif

!c         print*,'OUT OF MAIN LOOP MULTIPLICATION'    
!c
!c           %-----------------------------------------%
!c           | L O O P   B A C K to call DSAUPD again. |
!c           %-----------------------------------------%
!c

            go to 20
         
         end if
!     %-----------------------------------------%
!     | Either we have convergence, or there is |
!     | an error.                               |
!     %-----------------------------------------%
      if ( info .lt. 0 ) then
!        %----------------------------%
!        |  Error message, check the  |
!        |  documentation in DSAUPD    |
!        %----------------------------%
         print *, ' '
         print *, ' Error with _saupd, info = ', info
         print *, ' Check the documentation of _naupd.'
         print *, ' ' 
      else 
!        %-------------------------------------------%
!        | No fatal errors occurred.                 |
!        | Post-Process using DSEUPD .                |
!        |                                           |
!        | Computed eigenvalues may be extracted.    |
!        |                                           |
!        | Eigenvectors may also be computed now if  |
!        | desired.  (indicated by rvec = .true.)    |
!        %-------------------------------------------%
         rvec = .true.
         !call zneupd  (rvec, 'A', select, d, v, n, sigma, &
         !    workev, bmat, n, which, nev, tol, resid, ncv, v, & 
!             n, iparam, ipntr, workd, workl, lworkl, rwork, &
 !            ierr)
   call dseupd ( rvec, 'A', select, d, v, nmin, sigma, bmat, nmin, which, nev, tol, resid, ncv, v, nmin, iparam, ipntr, workd, workl, lworkl, ierr )

!c        %----------------------------------------------%
!c        | Eigenvalues are returned in the first column |
!c        | of the two dimensional array D and the       |
!c        | corresponding eigenvectors are returned in   |
!c        | the first NEV columns of the two dimensional |
!c        | array V if requested.  Otherwise, an         |
!c        | orthogonal basis for the invariant subspace  |
!c        | corresponding to the eigenvalues in D is     |
!c        | returned in V.                               |
!c        %----------------------------------------------%
!c
 
         if ( ierr .ne. 0) then
!           %------------------------------------%
!           | Error condition:                   |
!           | Check the documentation of ZNEUPD . |
!           %------------------------------------%
             print *, ' '
             print *, ' Error with _seupd, info = ', ierr
             print *, ' Check the documentation of _neupd. '
             print *, ' '
         else
             nconv = iparam(5)
             do 80 j=1, nconv

!                c
!c               %---------------------------%
!c               | Compute the residual norm |
!c               |                           |
!c               |   ||  A*x - lambda*x ||   |
!c               |                           |
!c               | for the NCONV accurately  |
!c               | computed eigenvalues and  |
!c               | eigenvectors.  (iparam(5) |
!c               | indicates how many are    |
!c               | accurate to the requested |
!c               | tolerance)                |
!c               %---------------------------%
!c
!c                call av(nx, v(1,j), ax)
!c               print*,"LOOP=",j
!c               print*,"MULTIPLICATION TEST AFTER=",sum(dsa),sum(v(1:N,j)
              
!c     &         ), sum(ax)                     
if (ncol==nrow) then
   call mkl_dcsrmm('N', nrow, 1, nrow, 1.0d0, matdescra, dsa,jsa,isa,isa(2), v(1,j), nrow, 0.0d0, bx, nrow)
   elseif (nrow>ncol) then
            call mkl_dcsrmm('N', nrow, 1, ncol, 1.0d0, matdescra, dsa,jsa,isa,isa(2), v(1,j), ncol, 0.0d0, ax, nrow) 
            call mkl_dcsrmm('T', nrow, 1, ncol, 1.0d0, matdescra, dsa,jsa,isa,isa(2), ax, nrow, 0.0d0, bx, ncol)
    elseif (nrow<ncol) then
            call mkl_dcsrmm('T', nrow, 1, ncol, 1.0d0, matdescra, dsa,jsa,isa,isa(2), v(1,j), nrow, 0.0d0, ax, ncol) 
            call mkl_dcsrmm('N', nrow, 1, ncol, 1.0d0, matdescra, dsa,jsa,isa,isa(2), ax, ncol, 0.0d0, bx, nrow)            
endif



!c               print*,"MULTIPLICATION TEST AFTER=",sum(dsa),sum(v(1:N,j)
!c     &         ), sum(ax)                

                call daxpy(nmin, -d(j,1), v(1,j), 1, bx, 1)
                d(j,2) = dnrm2(nmin, bx, 1)
                d(j,2) = d(j,2) / abs(d(j,1))
                
80           end do
!            %-----------------------------%
!            | Display computed residuals. |
!            %-----------------------------%
            ! call dmout (6, nconv, 2, d, ncv, -6,  &
            !     'Ritz values (Real, Imag) and direct residuals')
if (nrow/=ncol) then
print *,'singular values'
do i=1,nconv
print *,i,sqrt(d(i,1)),d(i,2)
enddo
else
print *,'eigenvalues'
do i=1,nconv
print *,i,d(i,1),d(i,2)
enddo
endif




!             call dmout (6, nconv, 2, d, ncv, -6,  &
!                 'Ritz values (Real, Imag) and direct residuals')

          end if

call system_clock(tf,tim)

!        %-------------------------------------------%
!        | Print additional convergence information. |
!        %-------------------------------------------%
         print *, ' Simulation Time', 1.0d0*(tf-ti)/(1.0d0*tim)
         if ( info .eq. 1) then
             print *, ' '
             print *, ' Maximum number of iterations reached.'
             print *, ' '
         else if ( info .eq. 3) then
             print *, ' '
             print *, ' No shifts could be applied during implicit',  &
     &                ' Arnoldi update, try increasing NCV.'
             print *, ' '
         end if
         print *, ' '
         print *, 'ARPACK_DRIVER '
         print *, '============='
         print *, ' '
         print *, ' Size of the matrix is ', nrow,ncol
         print *, ' The number of Ritz values requested is ', nev 
         print *, ' The number of Arnoldi vectors generated',&
                 ' (NCV) is ', ncv
         print *, ' What portion of the spectrum: ', which
         print *, ' The number of converged Ritz values is ', &
                   nconv
         print *, ' The number of Implicit Arnoldi update', &
                 ' iterations taken is ', iparam(3)
         print *, ' The number of OP*x is ', iparam(9)
         print *, ' The convergence criterion is ', tol
         print *, ' '
!
      end if
!
 9000 continue
!
    end program driver_arpack

   

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
