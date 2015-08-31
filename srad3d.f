C***********************************************************************
      PROGRAM SRADL
C
C  SPS expansion for the R problem with L.GE.0 [PRA 75, 062704 (2007)]
C  Partial-wave scattering characteristics
C
C***********************************************************************
      IMPLICIT REAL*8 (A,B,D-H,O-Z)
      IMPLICIT COMPLEX*16 (C)
      PARAMETER(PI=3.141592653589793238462643D0)
      PARAMETER(LANMAX=28)
      ALLOCATABLE X(:),W(:),T(:,:),DVR(:),PIR(:),RAD(:)
      ALLOCATABLE CK(:),CE(:),CVEC(:,:),CPHI(:)
      NAMELIST /INF/MODEL,RADA,NDVR,KPOL,LAN,KEYA,AMIN,AMAX,NUMA
      COMMON /POT_C/MODEL
C
C  Input parameters
C
      OPEN(1,file='inf')
      READ(1,INF)
      CLOSE(1)
      IF(LAN.GT.LANMAX) STOP ' *** LANMAX ERROR'
      NDVRS=(NDVR*(NDVR+1))/2
      ALLOCATE(X(NDVR),W(NDVR),T(NDVR+1,NDVR+1),DVR(NDVRS),
     &         PIR(NDVR),RAD(NDVR))
      NSPS=2*NDVR+LAN
      ALLOCATE(CK(NSPS),CE(NSPS),CVEC(NDVR,NSPS),CPHI(NSPS))
C
C  SPS basis
C
      CALL SPSR3D(1,KPOL,LAN,RADA,NDVR,X,W,T,DVR,PIR,RAD,
     &           NSPS,CK,CE,CVEC,CPHI,NB,NA,NOI)
      IF(NB.GT.0) WRITE(*,70) DREAL(CE(NB))
      CALPHA = DCMPLX(RADA, 0)
      DO i = 1,NSPS
	CALPHA = CALPHA + DCMPLX(0, 1)/CK(i)
      ENDDO
      print *, "alpha = ", DREAL(CALPHA)
	
         stop
C
C  Scattering calculations
C
      KEY=1
      SELECT CASE(KEY)
C === S-matrix - real momenta
      CASE(1)
      FAC=4.D0*PI*DBLE(2*LAN+1)
      CALL IKESTEP(KEYA,AMIN,AMAX,NUMA)
      DO 10 iA=0,NUMA
        CALL KESTEP(iA,AK,EE)
        IF(AK.EQ.0.D0) GOTO 10
C --- Product formula
        CALL SPRO(RADA,NSPS,CK,AK,CSP)
        DELP=(0.D0,-0.5D0)*CDLOG(CSP)
        SIGP=FAC*(DSIN(DELP)/AK)**2
        write(11,77) AK,EE,CSP,DELP,SIGP
C --- Sum formula
        CALL SSUML(LAN,RADA,NSPS,CK,CPHI,AK,CSS)
        DELS=(0.D0,-0.5D0)*CDLOG(CSS)
        SIGS=FAC*(DSIN(DELS)/AK)**2
        write(12,77) AK,EE,CSS,DELS,SIGS
C --- test
        IF(CDABS(CSS/CSP-1.D0).GT.1.D-6) write(*,*) AK,SIGP,SIGS
 10   CONTINUE
C === S-matrix - complex momenta
      CASE(2)
      DKR=0.01D0
      NKR=500
      DKI=0.001D0
      NKI=500
      DO ikr=0,NKR
        DO iki=0,NKI
          CAK=DCMPLX(DKR*DBLE(ikr),DKI*DBLE(iki))
          CS=CDEXP(-(0.D0,2.D0)*CAK*RADA)
c          CS=1.D0
          DO n=1,NSPS
            CS=CS*(CK(n)+CAK)/(CK(n)-CAK)
          ENDDO
c          write(21,77) CAK,CS
          write(22,77) CAK,CS
        ENDDO
        write(22,77)
      ENDDO
      write(*,77) CS
      END SELECT
C
 70   FORMAT(' ground state energy = ',E19.12)
 77   FORMAT(99(E19.12,1X))
      STOP
      END
C***********************************************************************