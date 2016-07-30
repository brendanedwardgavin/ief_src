        !COMPILER-GENERATED INTERFACE MODULE: Thu Jul 28 16:03:26 2016
        MODULE PRE_DMULT__genmod
          INTERFACE 
            SUBROUTINE PRE_DMULT(N,M,A,X,Y,MOUT,ZE)
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              COMPLEX(KIND=8), INTENT(IN) :: A(N,N)
              COMPLEX(KIND=8), INTENT(IN) :: X(N,M)
              COMPLEX(KIND=8), INTENT(IN) :: Y(N,M)
              COMPLEX(KIND=8), INTENT(OUT) :: MOUT(N,M)
              COMPLEX(KIND=8) :: ZE
            END SUBROUTINE PRE_DMULT
          END INTERFACE 
        END MODULE PRE_DMULT__genmod
