      SUBROUTINE FORCE(N, R, F, U, BX, BY, BZ, DE, BELTA, R0, R_CUTOFF)
      INTEGER N
      DOUBLE PRECISION R(3*N)
      DOUBLE PRECISION F(3*N)
      DOUBLE PRECISION U(1)
      DOUBLE PRECISION BX, BY, BZ
      DOUBLE PRECISION DE
      DOUBLE PRECISION BELTA
      DOUBLE PRECISION R0
      DOUBLE PRECISION R_CUTOFF, RCsq
      DOUBLE PRECISION R_IJ_X, R_IJ_Y, R_IJ_Z, R_IJ
      DOUBLE PRECISION UU, DUU, UU_SHIFT, DUMMY
      INTEGER I,J
      DUMMY=DBLE(2.0)*BELTA*DE
      RCsq=R_CUTOFF*R_CUTOFF
      UU_SHIFT=DEXP(BELTA*(R0-R_CUTOFF))
      UU_SHIFT=UU_SHIFT*UU_SHIFT-DBLE(2.0)*UU_SHIFT
      U(1)=DBLE(0.0)
      DO I=1, N-1
        DO J=I+1, N
          R_IJ_X=R(3*I-2)-R(3*J-2)
          R_IJ_Y=R(3*I-1)-R(3*J-1)
          R_IJ_Z=R(3*I)-R(3*J)
          R_IJ_X=R_IJ_X-BX*DNINT(R_IJ_X/BX)
          R_IJ_Y=R_IJ_Y-BY*DNINT(R_IJ_Y/BY)
          R_IJ_Z=R_IJ_Z-BZ*DNINT(R_IJ_Z/BZ)
          R_IJ=R_IJ_X*R_IJ_X+R_IJ_Y*R_IJ_Y+R_IJ_Z*R_IJ_Z
          IF(R_IJ<RCsq) THEN
            R_IJ=DSQRT(R_IJ)
            UU=DEXP(BELTA*(R0 - R_IJ))
            DUU=UU*UU-UU
            F(3*I-2)=F(3*I-2)+DUU*R_IJ_X/R_IJ
            F(3*I-1)=F(3*I-1)+DUU*R_IJ_Y/R_IJ
            F(3*I)=F(3*I)+DUU*R_IJ_Z/R_IJ
            F(3*J-2)=F(3*J-2)-DUU*R_IJ_X/R_IJ
            F(3*J-1)=F(3*J-1)-DUU*R_IJ_Y/R_IJ
            F(3*J)=F(3*J)-DUU*R_IJ_Z/R_IJ
            U(1)=U(1)+ DUU-UU-UU_SHIFT
          END IF
        END DO
      END DO
      DO I=1, 3*N
        F(I)=F(I)*DUMMY
      END DO
      U(1)=U(1)*DE
      END SUBROUTINE FORCE

