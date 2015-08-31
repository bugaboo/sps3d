      FUNCTION CBESPOLN(N,CZ)
C=======================================================================
C  The Bessel polynomial y_N(z)
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A,B,D-H,O-Z)
      IMPLICIT COMPLEX*16 (C)
C
      CP=1.D0
      CQ=1.D0
      DO i=1,N
        CTMP=CP
        CP=DBLE(2*i-1)*CZ*CP+CQ
        CQ=CTMP
      ENDDO
      CBESPOLN=CP
C
      RETURN 
      END
C=======================================================================