      SUBROUTINE FORCE(N, R, F, U, BX, BY, BZ, U0, EPS, SIGMA, R_CUTOFF)
      INTEGER N
      DOUBLE PRECISION R(3*N)
      DOUBLE PRECISION F(3*N)
      DOUBLE PRECISION U(1)
      DOUBLE PRECISION BX, BY, BZ
      DOUBLE PRECISION U0
      DOUBLE PRECISION EPS
      DOUBLE PRECISION SIGMA
      DOUBLE PRECISION R_CUTOFF, RCsq
      DOUBLE PRECISION R_IJ_X, R_IJ_Y, R_IJ_Z, R_IJ
      DOUBLE PRECISION iRsq, iRs2, iRs6, iRs12
      DOUBLE PRECISION UU, DUU,UU_SHIFT
      INTEGER I,J
      RCsq=R_CUTOFF*R_CUTOFF
      iRs2=SIGMA*SIGMA/(RCsq)
      iRs6=iRs2*iRs2*iRs2
      iRs12=iRs6*iRs6
      UU_SHIFT=(iRs12-iRs6)
!      UU_SHIFT=UU_SHIFT*UU_SHIFT-DBLE(2.0)*UU_SHIFT
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
          IF(R_IJ<=RCsq) THEN
            iRsq=DBLE(1.0)/R_IJ
            iRs2=SIGMA*SIGMA*iRsq
            iRs6=iRs2*iRs2*iRs2
            iRs12=iRs6*iRs6
            UU=(iRs12-iRs6)
!            R_IJ=DSQRT(R_IJ)
            DUU=(DBLE(12.0)*iRs12-DBLE(6.0)*iRs6)*iRsq
            F(3*I-2)=F(3*I-2)+DUU*R_IJ_X
            F(3*I-1)=F(3*I-1)+DUU*R_IJ_Y
            F(3*I)=F(3*I)+DUU*R_IJ_Z
            F(3*J-2)=F(3*J-2)-DUU*R_IJ_X
            F(3*J-1)=F(3*J-1)-DUU*R_IJ_Y
            F(3*J)=F(3*J)-DUU*R_IJ_Z
!            U(1)=U(1)+ DUU-UU-UU_SHIFT
            U(1)=U(1)+UU-UU_SHIFT
          END IF
        END DO
      END DO
      DO I=1, 3*N
        F(I)=F(I)*U0
      END DO
      U(1)=U(1)*U0
      END SUBROUTINE FORCE

