      SUBROUTINE SPRO(RADA,NSPS,CK,CAK,CS)
C=======================================================================
C  S-matrix in the one-channel R problem for L.GE.0
C  the 'product' formula
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A,B,D-H,O-Z)
      IMPLICIT COMPLEX*16 (C)
      DIMENSION CK(*)
C
      CS=CDEXP(-(0.D0,2.D0)*CAK*RADA)
      DO n=1,NSPS
        CS=CS*(CK(n)+CAK)/(CK(n)-CAK)
      ENDDO
C
      RETURN
      END
C=======================================================================