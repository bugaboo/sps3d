      SUBROUTINE SSUML(LAN,RADA,NSPS,CK,CPHI,CAK,CS)
C=======================================================================
C  S-matrix in the one-channel R problem for L.GE.0
C  the 'sum' formula 
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A,B,D-H,O-Z)
      IMPLICIT COMPLEX*16 (C)
      DIMENSION CK(*),CPHI(*)
C
      CS=0.D0
      DO n=1,NSPS
        CS=CS+CPHI(n)*CPHI(n)/CK(n)/(CK(n)-CAK)
      ENDDO
      CTMP=CBESPOLN(LAN,(0.D0,1.D0)/CAK/RADA)
      CTMP1=CBESPOLN(LAN,-(0.D0,1.D0)/CAK/RADA)
      CS=CTMP1+(0.D0,1.D0)*CAK*CS/CTMP
      CS=CDEXP(-(0.D0,2.D0)*CAK*RADA)*(-1)**LAN*CS/CTMP
C
      RETURN
      END
C=======================================================================