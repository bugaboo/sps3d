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
      ALLOCATABLE CK(:,:),CE(:),CPHI(:,:,:),CVEC(:,:),CS(:,:), CEL(:)
      ALLOCATABLE CKT(:), CET(:), CPHIT(:,:), CVECT(:,:), CST(:,:)
      ALLOCATABLE L(:,:),M(:,:),PIR(:), NANG(:), NSPS(:)
      ALLOCATABLE KSYMA(:),LANA(:),MANA(:),LT(:),MT(:)
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
      
      SELECT CASE(KSYM)
      CASE (0)
	NSYM = 1
	ALLOCATE(KSYMA(NSYM))
	KSYMA(1) = KSYM
      CASE (10,20)
	NSYM = 2
	ALLOCATE (KSYMA(NSYM))
	KSYMA(1) = 10
	KSYMA(2) = 20
      CASE (1, 2)
	NSYM = 2*LAN + 1
	ALLOCATE (KSYMA(NSYM))
	DO i = 1, NSYM
	  KSYMA = KSYM
	ENDDO
      END SELECT
      ALLOCATE(LANA(NSYM), MANA(NSYM), NANG(NSYM), NSPS(NSYM))
      DO i = 1, NSYM
	LANA(i) = LAN
	MANA(i) = MAN
      ENDDO
      
      NSPSM = 0
      NANGMAX = 0
      DO i = 1, NSYM
        CALL ANGBAS(KSYMA(i),LT,MT,LMAX,NANG(i),LANA(i),MANA(i))
        IF (NANG(i).GT.NANGMAX) NANGMAX = NANG(i)
        NSPS(i) = 0
	DO j= 1, NANG(i)
	  NSPS(i) = NSPS(i) + 2 * NDVR + LT(j)
	ENDDO
	PRINT *, NSPS(i), NSPSM
	IF (NSPS(i).GT.NSPSM) NSPSM = NSPS(i)
	PRINT *, NSPS(i), NSPSM, 1
	DEALLOCATE(LT, MT)
      ENDDO
      
      ALLOCATE (L(NANGMAX,NSYM), M(NANGMAX,NSYM))
      L = 0
      M = 0
      ALLOCATE(CK(NSPSM,NSYM))
      ALLOCATE(PIR(NDVR),CPHI(NANGMAX,NSPSM,NSYM))
      DO i = 1, NSYM
	CALL ANGBAS(KSYMA(i),LT,MT,LMAX,NANG(i),LANA(i),MANA(i))
	L(1:NANG(i),i) = LT
	M(1:NANG(i),i) = MT
	ALLOCATE(CKT(NSPS(i)),CE(NSPS(i)),CVEC(NDVR*NANG(i),NSPS(i)))
	ALLOCATE(CPHIT(NANG(i),NSPS(i)))
	NVEC = NDVR * NANG(i)
	CALL SPS3D(LMAX,NANG(i),LT,MT,RADA,NDVR,NTET,NPHI,PIR,CKT,CE
     &	    	,CVEC,CPHIT,NB,NA,NOI) 
	CK(:NSPS(i),i) = CKT
	CPHI(:NANG(i),:NSPS(i),i) = CPHIT
	DEALLOCATE(LT, MT,CKT,CE,CVEC,CPHIT)
      ENDDO
      CAK=5.D0
      ALLOCATE (CS(NANG(1), NANG(1)),CEL(NANG(1)))
      CALL SSUM3D(L(:NANG(1),1),RADA,NSPS(1),NANG(1),CK(:NSPS(1),1),
     &		 CPHI(:NANG(1),:NSPS(1),1),CAK,CS)
      CALL CDETS(CS,CEL,CDET,NANG(1))
      PRINT *, "DET = ", CDET
c      IF (.NOT. UCHECK(CS,NANG)) PRINT *, "UNITARITY PROBLEM"

C  END
C      DEALLOCATE(PIR,CPHI,CS,CEL)
 70   FORMAT(' ground state energy = ',E19.12)
 77   FORMAT(4(E19.12,1X))
 
      CONTAINS 
        FUNCTION UCHECK(CS, NANG)
          LOGICAL UCHECK
          DIMENSION CS(NANG, NANG)


          UCHECK = .FALSE.
          DO i = 1, NANG
            DO j = 1, NANG
              CTMP = 0.D0
              DO k = 1, NANG
                CTMP = CTMP + CS(i,k) * DCONJG(CS(j,k))
              ENDDO
              IF (i.EQ.j) CTMP = CTMP - 1.D0
              IF (CDABS(CTMP) .GT. 1.D-10) THEN
                PRINT *, "NONUNITARY S MATRIX", i, j, CTMP
                RETURN 
              ENDIF
            ENDDO
          ENDDO
          UCHECK = .TRUE.
          RETURN 
        END FUNCTION
        END PROGRAM