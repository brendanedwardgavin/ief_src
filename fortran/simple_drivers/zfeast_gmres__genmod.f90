        !COMPILER-GENERATED INTERFACE MODULE: Thu Jul 28 16:03:26 2016
        MODULE ZFEAST_GMRES__genmod
          INTERFACE 
            SUBROUTINE ZFEAST_GMRES(IJOB,STATEVARS,BRHS,X,V,AV,AX,ZE,N,M&
     &,MAXM,EPS,RESTARTS,M0,XWORK,WORKIN,AV2,TIMES)
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: IJOB
              INTEGER(KIND=4) :: STATEVARS(3)
              COMPLEX(KIND=8) :: BRHS(N,*)
              COMPLEX(KIND=8) :: X(N,*)
              COMPLEX(KIND=8) :: V(N,*)
              COMPLEX(KIND=8) :: AV(N,*)
              COMPLEX(KIND=8) :: AX(N,*)
              COMPLEX(KIND=8) :: ZE
              INTEGER(KIND=4) :: MAXM
              REAL(KIND=8) :: EPS
              INTEGER(KIND=4) :: RESTARTS
              INTEGER(KIND=4) :: M0
              COMPLEX(KIND=8) :: XWORK(N,*)
              COMPLEX(KIND=8) :: WORKIN(N,*)
              COMPLEX(KIND=8) :: AV2(N,*)
              REAL(KIND=8) :: TIMES
            END SUBROUTINE ZFEAST_GMRES
          END INTERFACE 
        END MODULE ZFEAST_GMRES__genmod
