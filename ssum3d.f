      SUBROUTINE SSUM3D(L,RADA,NSPS,NANG,CK,CPHI,CAK,CS)
C=======================================================================
C  S-matrix in the one-channel R problem for L.GE.0
C  the 'sum' formula 
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A,B,D-H,O-Z)
      IMPLICIT COMPLEX*16 (C)
      DIMENSION CK(*),CPHI(NANG,*),CS(NANG,*), L(NANG)
C
      DO nu=1,NANG
	DO mu=1,NANG
	  CTS=0.D0
	  DO n=1,NSPS
	    CTS=CTS+CPHI(nu,n)*CPHI(mu,n)/CK(n)/(CK(n)-CAK)
	  ENDDO
	  CTMP=CBESPOLN(L(nu),(0.D0,1.D0)/CAK/RADA)
	  CTMP1=CBESPOLN(L(mu),(0.D0,1.D0)/CAK/RADA)
	  CTS=(0.D0,1.D0)**(L(mu)+L(nu)+1)*CAK*CTS/CTMP
	  IF (nu.EQ.mu) CTS=CTS+CBESPOLN(L(nu),-(0.D0,1.D0)/CAK/RADA)
     &			*(-1)**L(nu)
	  CTS=CDEXP(-(0.D0,2.D0)*CAK*RADA)*CTS/CTMP1
	  CS(nu,mu)=CTS
	ENDDO
      ENDDO
      
C
      RETURN
      END
C=======================================================================