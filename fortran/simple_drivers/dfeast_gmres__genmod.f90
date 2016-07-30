        !COMPILER-GENERATED INTERFACE MODULE: Thu Jul 28 16:03:25 2016
        MODULE DFEAST_GMRES__genmod
          INTERFACE 
            SUBROUTINE DFEAST_GMRES(IJOB,STATEVARS,BRHS,X,V,AV,AX,ZE,N,M&
     &,RESTARTS,M0,XWORK,WORKIN,AV2)
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: IJOB
              INTEGER(KIND=4) :: STATEVARS(3)
              COMPLEX(KIND=8) :: BRHS(N,*)
              COMPLEX(KIND=8) :: X(N,*)
              COMPLEX(KIND=8) :: V(N,*)
              COMPLEX(KIND=8) :: AV(N,*)
              COMPLEX(KIND=8) :: AX(N,*)
              COMPLEX(KIND=8) :: ZE
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: RESTARTS
              INTEGER(KIND=4) :: M0
              REAL(KIND=8) :: XWORK(N,*)
              REAL(KIND=8) :: WORKIN(N,*)
              COMPLEX(KIND=8) :: AV2(N,*)
            END SUBROUTINE DFEAST_GMRES
          END INTERFACE 
        END MODULE DFEAST_GMRES__genmod
