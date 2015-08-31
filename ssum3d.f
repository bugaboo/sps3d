      SUBROUTINE SSUML(LAN,RADA,NSPS,CK,CPHI,AK,CS)
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
        CS=CS+CPHI(n)*CPHI(n)/CK(n)/(CK(n)-AK)
      ENDDO
      CTMP=CBESPOLN(LAN,(0.D0,1.D0)/AK/RADA)
      CS=DCONJG(CTMP)+(0.D0,1.D0)*AK*CS/CTMP
      CS=CDEXP(-(0.D0,2.D0)*AK*RADA)*(-1)**LAN*CS/CTMP
C
      RETURN
      END
C=======================================================================