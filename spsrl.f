      SUBROUTINE SPSRL(KOUT,KPOL,LAN,RADA,NDVR,X,W,T,DVR,PIR,RAD,
     &                 NSPS,CK,CE,CVEC,CPHI,NB,NA,NOI)
C=======================================================================
C  SPS basis for the R problem with L.GE.0 [PRA 75, 062704 (2007)]
C  KPOL.EQ.0 - Legendre polynomials
C  KPOL.NE.0 - Jacobi(0,2*LAN) polynomials
C  'product' formula for surface amplitudes
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A,B,D-H,O-Z)
      IMPLICIT COMPLEX*16 (C)
      DIMENSION CZER(28)
      DIMENSION X(*),W(*),T(NDVR+1,*),DVR(*),PIR(*),RAD(*)
      DIMENSION CK(*),CE(*),CVEC(NDVR,*),CPHI(*)
      ALLOCATABLE HAM(:,:),VKR(:),VKI(:),VK(:),WK(:,:),CWK(:,:)
C
C  Setting up DVR arrays
C
      SR=0.5D0*RADA
      NDVR2=2*NDVR
      CALL DVRR(KPOL,LAN,NDVR,X,W,T,DVR,PIR)
      TMP=0.5D0*DBLE(LAN*(LAN+1))
      ij=0
      DO j=1,NDVR
        RAD(j)=SR*(1.D0+X(j))
        PIR(j)=PIR(j)/RAD(j)
        DO i=1,j
          ij=ij+1
          DVR(ij)=DVR(ij)/(RAD(i)*RAD(j))+PIR(i)*PIR(j)
        ENDDO
        DVR(ij)=DVR(ij)+TMP/RAD(j)/RAD(j) +POTR(RAD(j))
      ENDDO

C
C  Zeros of the reverse Bessel polynomial
C
      CALL BESZER(LAN,CZER)
C
C  Setting up SPS Hamiltonian
C
      ALLOCATE(HAM(NSPS,NSPS))
      DO j=1,NSPS
        DO i=1,NSPS
          HAM(i,j)=0.D0
        ENDDO
      ENDDO
      DO i=1,NDVR
        HAM(i,NDVR+i)=1.D0
      ENDDO
      ij=0
      DO j=1,NDVR
        DO i=1,j
          ij=ij+1
          TMP=-2.D0*DVR(ij)
          HAM(NDVR+i,j)=TMP
          HAM(NDVR+j,i)=TMP
          TMP=2.D0*RADA*PIR(i)*PIR(j)
          HAM(NDVR+i,NDVR+j)=TMP
          HAM(NDVR+j,NDVR+i)=TMP
        ENDDO
      ENDDO
      TMP=DSQRT(2.D0/RADA)
      DO ip=1,LAN/2
        ip1=NDVR2+2*ip-1
        ip2=ip1+1
        ZER1=DREAL(CZER(ip))
        ZER2=DIMAG(CZER(ip))
        DO i=1,NDVR
          TMPi=TMP*PIR(i)
          HAM(NDVR+i,ip1)=2.D0*TMPi
          HAM(ip1,i)=-ZER1*TMPi
          HAM(ip2,i)=-ZER2*TMPi
        ENDDO
        HAM(ip1,ip1)=-ZER1/RADA
        HAM(ip1,ip2)=+ZER2/RADA
        HAM(ip2,ip1)=-ZER2/RADA
        HAM(ip2,ip2)=-ZER1/RADA
      ENDDO
      IF(MOD(LAN,2).EQ.0) GOTO 10
      ZER=CZER(LAN/2+1)
      DO i=1,NDVR
        TMPi=TMP*PIR(i)
        HAM(NDVR+i,NSPS)=TMPi
        HAM(NSPS,i)=-ZER*TMPi
        IF (LAN.EQ.1) THEN
C	  PRINT *, PIR(i)
	ENDIF
      ENDDO
      HAM(NSPS,NSPS)=-ZER/RADA
 10   CONTINUE
      SELECT CASE(LAN)
	CASE(0)
	OPEN(2,FILE='haml0')
	CASE(1)
	OPEN(2,FILE='haml1')
	CASE(2)
	OPEN(2,FILE='haml2')
	CASE(10)
	OPEN(2,FILE='haml10')
      END SELECT
      DO i=1,NSPS
	DO j=1,NSPS
	  write(2,77) HAM(i,j)
	ENDDO
      ENDDO
      CLOSE(2)

C  Solution of the SPS EVP
C
C --- IMSL
c     ALLOCATE(CWK(NSPS,NSPS))
c     CALL DEVCRG(NSPS,HAM,NSPS,CK,CWK,NSPS)
c      DO n=1,NSPS
c        CK(n)=-(0.D0,1.D0)*CK(n)
c        DO i=1,NDVR
c          CVEC(i,n)=CWK(i,n)
c        ENDDO
c      ENDDO
c      DEALLOCATE(CWK)
C --- LAPACK
      NVK=10*NSPS
      ALLOCATE(VKR(NSPS),VKI(NSPS),VK(NVK),WK(NSPS,NSPS))
      CALL DGEEV('N','V',NSPS,HAM,NSPS,VKR,VKI,ERR,1,WK,NSPS,
     &           VK,NVK,info)
      DO n=1,NSPS
        CK(n)=-(0.D0,1.D0)*DCMPLX(VKR(n),VKI(n))
        IF(VKI(n).EQ.0.D0) THEN
          DO i=1,NDVR
            CVEC(i,n)=DCMPLX(WK(i,n),0.D0)
          ENDDO
        ELSE IF(VKI(n).GT.0.D0) THEN
          DO i=1,NDVR
            CVEC(i,n)  =DCMPLX(WK(i,n),+WK(i,n+1))
            CVEC(i,n+1)=DCMPLX(WK(i,n),-WK(i,n+1))
          ENDDO
        ENDIF
      ENDDO

      DEALLOCATE(VKR,VKI,VK,WK)
      DEALLOCATE(HAM)
C
C
C  Ordering the SPS
C
      CALL SPSODR(NDVR,NSPS,CK,CE,CVEC,NB,NA,NOI)
      
	
C
C  Normalization and surface amplitudes
C
      DO n=1,NSPS
        CSUM=0.D0
        DO i=1,NDVR
          CSUM=CSUM+CVEC(i,n)*CVEC(i,n)
        ENDDO
        CSUM=SR*CSUM
        CTMP=CBESPOLN(LAN,(0.D0,1.D0)/CK(n)/RADA)
        CPHI2=(0.D0,2.D0)*CK(n)*CTMP*CTMP
        DO m=1,NSPS
          IF(m.NE.n) CPHI2=CPHI2*(CK(n)+CK(m))/(CK(n)-CK(m))
        ENDDO
        CPHI(n)=CDSQRT(CPHI2)
        CTMP=1.D0
        DO ip=1,LAN
          CTMP=CTMP+CZER(ip)/((0.D0,1.D0)*CK(n)*RADA+CZER(ip))**2
        ENDDO
        CTMP=CDSQRT((1.D0-(0.D0,0.5D0)*CPHI2*CTMP/CK(n))/CSUM)
        DO i=1,NDVR
          CVEC(i,n)=CVEC(i,n)*CTMP
        ENDDO
      ENDDO
C
C  Orthogonality
C
      DO n=1,NSPS
	DO k=1,NSPS
	  CSUMC=0.D0
	  DO i=1,NDVR
	    CSUMC = CSUMC + CVEC(i,n)*CVEC(i,k)
	  ENDDO
	  CSUMXI = 0.D0
	  
	  DO j=1,LAN
	    CSUMXI=CSUMXI+CZER(j)**2/((0.D0,1.D0)*CK(n)*RADA+CZER(j))
     &			/((0.D0,1.D0)*CK(k)*RADA+CZER(j))
	  ENDDO
	  CSUMC=CSUMC+(0.D0,1.D0)*CPHI(n)*CPHI(m)*(1.D0+CSUMXI)
     & 			/(CK(n)+CK(m))
c	 IF (n.NE.k.AND.CDABS(CSUMC).GT.1.D-6) PRINT *,CDABS(CSUMC),n,k
	ENDDO
      ENDDO
C
C  The sign of surface amplitudes from the 'sum' formula
C
      DO n=1,NSPS
        CSUM=0.D0
        DO i=1,NDVR
          CSUM=CSUM+CVEC(i,n)*PIR(i)
        ENDDO
        CSUM=RADA*CSUM
        CTMP=CSUM/CPHI(n)
        IF(DREAL(CTMP).LT.0.D0) CPHI(n)=-CPHI(n)
      ENDDO
      PRINT *, "CPHI L"
      DO n=1,NSPS
	CSUM=0.D0
	DO i=1,NDVR
	  CSUM = CSUM + PIR(i)*CVEC(i,n)
	ENDDO
	CSUM=CSUM*RADA
	PRINT *, CPHI(n), CSUM-CPHI(n)
      ENDDO
      IF(KOUT.EQ.0) RETURN
C
C  Writing out the eigenvalues
C
      OPEN(1,FILE='spseig')
      DO n=1,NSPS
        write(1,77) CK(n),CE(n)
      ENDDO
      CLOSE(1)
C
C  Writing out hamiltonian
C
      
C
 77   FORMAT(4(E19.12,1X))
      RETURN
      END
C=======================================================================