      FUNCTION SLENGTH(RADA,NSPS,NANG,CK)
C=======================================================================
C  Scattering length in the one-channel R problem
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A,B,D-H,O-Z)
      IMPLICIT COMPLEX*16 (C)
      DIMENSION CK(*)
C
      CTMP=0.D0
      DO n=1,NSPS
        CTMP=CTMP+(0.D0,1.D0)/CK(n)
      ENDDO
      SLENGTH=DBLE(NANG)*RADA+DREAL(CTMP)
C
      RETURN
      END
C=======================================================================