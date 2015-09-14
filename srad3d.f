C***********************************************************************
      PROGRAM MAIN
C
C  SPS expansion for the R problem with L.GE.0 [PRA 75, 062704 (2007)]
C  Partial-wave scattering characteristics
C
C***********************************************************************
      USE ANGSYM
      IMPLICIT REAL*8 (A,B,D-H,O-Z)
      IMPLICIT COMPLEX*16 (C)
      PARAMETER(PI=3.141592653589793238462643D0)
      PARAMETER(LANMAX=28)
      ALLOCATABLE CK(:),CE(:),CPHI(:,:),CVEC(:,:),CS(:,:), CEL(:)
      ALLOCATABLE L(:),M(:),PIR(:)
      NAMELIST /INF3D/MODEL,RADA,NDVR,KPOL,LMAX,LAN,MAN,KSYM,NTET,NPHI,
     & 		KEYA,AMIN,AMAX,NUMA
      COMMON /POT_C/MODEL
      CHARACTER(LEN=20) :: SPSNAME
C
C  Input parameters
C
      OPEN(1,file='inf')
      READ(1,INF3D)
      CLOSE(1)
      IF(LAN.GT.LANMAX) STOP ' *** LANMAX ERROR'
      
      CALL ANGBAS(KSYM,L,M,LMAX,NANG,LAN,MAN)
      NSPS = 0
      DO i= 1, NANG
	NSPS = NSPS + 2 * NDVR + L(i)
      ENDDO
      ALLOCATE(CK(NSPS),CE(NSPS),CVEC(NDVR*NANG,NSPS))
      ALLOCATE(PIR(NDVR),CPHI(NANG,NSPS),CS(NANG,NANG),CEL(NANG))
      CALL SPS3D(LMAX,NANG,L,M,RADA,NDVR,NTET,NPHI,PIR,CK,
     &	    CE,CVEC,CPHI,NB,NA,NOI) 
      IF(NB.GT.0) WRITE(*,70) DREAL(CE(NB))
      CAK=5.D0
      CALL SSUM3D(L,RADA,NSPS,NANG,CK,CPHI,CAK,CS)
      CALL CDETS(CS,CEL,CDET,NANG)
      PRINT *, "Sum det", CDET
      CALL SPRO(RADA,NSPS,CK,CAK,CDETP)
      CDETP = CDETP*CDEXP(-(0.D0,2.D0)*CAK*RADA*(NANG-1))
      PRINT *, "Prod det", CDETP
      
C  Printing eigenvalues

      WRITE(SPSNAME,"(A7,I0,A1,I0,A1,I2.2,A1,I2.2)")
     & "spseig_",KSYM,'_',NDVR,'_',LMAX,'_',INT(RADA)
      OPEN(1,FILE=SPSNAME)
      DO n=1,NSPS
        write(1,77) CK(n),CE(n)
      ENDDO
      CLOSE(1)
C  END
      DEALLOCATE(PIR,CK,CE,CVEC,CPHI,CS,CEL)
 70   FORMAT(' ground state energy = ',E19.12)
 77   FORMAT(4(E19.12,1X))
      END PROGRAM