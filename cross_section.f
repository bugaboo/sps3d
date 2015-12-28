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
      ALLOCATABLE KSYMA(:),LANA(:),MANA(:),LT(:),MT(:),MULTIP(:)
      ALLOCATABLE DELTA(:),CDETP(:), CDELTA(:)
      NAMELIST /INF3D/MODEL,RADA,NDVR,KPOL,LMAX,LAN,MAN,KSYM,NTET,NPHI,
     & 		KEYA,AMIN,AMAX,NUMA
      COMMON /POT_C/MODEL
      CHARACTER(LEN=30) :: DELTANAME, FILENAME
      LOGICAL CONS_INPUT, EIG_FROM_FILE
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
      CASE (90,-90)
	NSYM = 2
	ALLOCATE (KSYMA(NSYM))
	KSYMA(1) = 90
	KSYMA(2) = -90
      CASE (1)
	NSYM = 2 * LMAX + 1
	ALLOCATE (KSYMA(NSYM))
	DO i = 1, NSYM
	  KSYMA(i) = KSYM
	ENDDO
      CASE (2)
	NSYM = LMAX + 1
	ALLOCATE (KSYMA(NSYM))
	KSYMA = 2
      CASE (10,20,30,40,50,60,70,80)
	NSYM = 8
	ALLOCATE (KSYMA(NSYM))
	DO i = 1, 8
	  KSYMA(i) = i * 10
	ENDDO
      CASE (100)
	NSYM = 1
	ALLOCATE (KSYMA(NSYM))
	KSYMA = 10
      END SELECT
      ALLOCATE(LANA(NSYM),MANA(NSYM),NANG(NSYM),NSPS(NSYM),MULTIP(NSYM))
      SELECT CASE(KSYM)
      CASE (2)
	DO i = 1, NSYM
	  LANA(i) = i - 1
	  MULTIP(i) = 2 * i - 1
	  MANA(i) = 0
	ENDDO
      CASE (-90, 0, 10,20,30,40,50,60,70,80, 90, 100)
	DO i = 1, NSYM
	  MULTIP(i) = 1
	  LANA(i) = LAN
	  MANA(i) = MAN
	ENDDO
      CASE (1)
	DO i = 1, NSYM
	  MULTIP(i) = 1
	  LANA(i) = LAN
	  MANA(i) = -LMAX - 1 + i
	ENDDO
      END SELECT
      
      NSPSM = 0
      NANGMAX = 0
      DO i = 1, NSYM
        CALL ANGBAS(KSYMA(i),LT,MT,LMAX,NANG(i),LANA(i),MANA(i))
        IF (NANG(i).GT.NANGMAX) NANGMAX = NANG(i)
        NSPS(i) = 0
	DO j= 1, NANG(i)
	  NSPS(i) = NSPS(i) + 2 * NDVR + LT(j)
	ENDDO
	IF (NSPS(i).GT.NSPSM) NSPSM = NSPS(i)
	DEALLOCATE(LT, MT)
      ENDDO
      
      ALLOCATE (L(NANGMAX,NSYM), M(NANGMAX,NSYM))
      L = 0
      M = 0
      DO i = 1, NSYM
	CALL ANGBAS(KSYMA(i),LT,MT,LMAX,NANG(i),LANA(i),MANA(i))
	L(1:NANG(i),i) = LT
	M(1:NANG(i),i) = MT
	DEALLOCATE (LT, MT)
      ENDDO
      ALLOCATE(CK(NSPSM,NSYM))
      ALLOCATE(PIR(NDVR),CPHI(NANGMAX,NSPSM,NSYM))
      
      CONS_INPUT = CONSOLE_INPUT(EMAX, NS, EIG_FROM_FILE)
      IF (EIG_FROM_FILE) THEN
	DO i = 1, NSYM
	  ALLOCATE(CKT(NSPS(i)), CPHIT(NANG(i),NSPS(i)))
	  WRITE(FILENAME,"(A8,I0,A1,I0,A1,I0,A1,I2.2,A1,I0)")
     & "eigenval",MODEL,'_',KSYMA(i),'_',NDVR,'_',LMAX,'_',INT(RADA)
	  OPEN(UNIT = 1, FILE=FILENAME, FORM = 'unformatted', 
     &		ACCESS = 'direct', RECL = NSPS(i) * 16)
	  READ(1, REC = 1) CKT
c	  do j = 1, nsps(i)
c	    write (*, 77) real(ckt(j)), aimag(ckt(j))
c	  enddo
	  WRITE(FILENAME,"(A8,I0,A1,I0,A1,I0,A1,I2.2,A1,I0)")
     & "eigenvec",MODEL,'_',KSYMA(i),'_',NDVR,'_',LMAX,'_',INT(RADA)
	  OPEN(UNIT = 1, FILE=FILENAME, FORM = 'unformatted', 
     &		ACCESS = 'direct', RECL = NSPS(i) * NANG(i) * 16)
	  READ(1, REC = 1) CPHIT
	  CK(:NSPS(i),i) = CKT
	  CPHI(:NANG(i),:NSPS(i),i) = CPHIT
	  DEALLOCATE(CKT,CPHIT)
	ENDDO
      ELSE
	DO i = 1, NSYM
	  ALLOCATE(CKT(NSPS(i)),CE(NSPS(i)),CVEC(NDVR*NANG(i),NSPS(i)))
	  ALLOCATE(CPHIT(NANG(i),NSPS(i)))
	  CALL SPS3D(LMAX,NANG(i),L(1:NANG(i),i),M(1:NANG(i),i),RADA,NDVR,
     &	    	NTET,NPHI,PIR,CKT,CE,CVEC,CPHIT,NB,NA,NOI) 
	  WRITE(FILENAME,"(A8,I0,A1,I0,A1,I0,A1,I2.2,A1,I0)")
     & "eigenval",MODEL,'_',KSYMA(i),'_',NDVR,'_',LMAX,'_',INT(RADA)
	  OPEN(UNIT = 1, FILE=FILENAME, FORM = 'unformatted', 
     &		ACCESS = 'direct', RECL = NSPS(i) * 16)
	  WRITE(1, REC = 1) CKT
	  WRITE(FILENAME,"(A8,I0,A1,I0,A1,I0,A1,I2.2,A1,I0)")
     & "eigenvec",MODEL,'_',KSYMA(i),'_',NDVR,'_',LMAX,'_',INT(RADA)
	  OPEN(UNIT = 1, FILE=FILENAME, FORM = 'unformatted', 
     &		ACCESS = 'direct', RECL = NSPS(i) * NANG(i) * 16)
	  WRITE(1, REC = 1) CPHIT
	  CK(:NSPS(i),i) = CKT
	  CPHI(:NANG(i),:NSPS(i),i) = CPHIT
	  DEALLOCATE(CKT,CE,CVEC,CPHIT)
	ENDDO
      ENDIF

      IF (CONS_INPUT) THEN
	ALLOCATE(CDELTA(NS), DELTA(NS), CDETP(NS))
	CDELTA = 1.D0
	DELTA = 0.D0
	CDETP = 1.D0
	DO i = 1, NSYM
	  ALLOCATE (CS(NANG(i), NANG(i)),CEL(NANG(i)))
	  WRITE(DELTANAME, "(A5,I0)"), "delta", KSYMA(i)
	  OPEN(1, FILE=DELTANAME)
	  DO k = 1, NS
	    CAK =  DBLE(k) * EMAX / NS
C	    CALL SSUM3D(L(:NANG(i),i),RADA,NSPS(i),NANG(i),
C     &		 CK(:NSPS(i),i),CPHI(:NANG(i),:NSPS(i),i),CAK,CS)
C	    CALL CDETS(CS,CEL,CDET,NANG(i))
C	    WRITE(1, "(E19.12,1X)", ADVANCE="NO") DBLE(k)*EMAX/DBLE(NS)
C	    DO n = 1, NANG(i)
C	      WRITE(1, "(E19.12,1X)", ADVANCE="NO") CEL(n)
C	    ENDDO
C	    WRITE(1,*)
C	    CDELTA(k) = CDELTA(k) * CDET
	    DELTA(k) = DELTA(k) + MULTIP(i)*AIMAG(CDLOG(CDET)) / 2.D0
	    CALL SPRO(RADA,NSPS(i),CK(:NSPS(i),i),CAK,CDETTP)
	    CTMP = CDEXP(-(0.D0,2.D0)*CAK*RADA*(NANG(i)-1))	    
	    CDETP(k) = CDETP(k)*(CDETTP*CTMP)**MULTIP(i)
	  ENDDO
	  CLOSE(1)
	  DEALLOCATE(CS, CEL)
	ENDDO
	OPEN(2, FILE='delta')
	PREVDEL = 0.0
	OFFSET = 0.0
	EPS = 1
	DO k = 1, NS
	  DO WHILE (DELTA(k).GT.PI/2.D0)
	    DELTA(k) = DELTA(k) - PI
	  ENDDO
	  DO WHILE (DELTA(k).LT.-PI/2.D0)
	    DELTA(k) = DELTA(k) + PI
	  ENDDO
	  DEL = AIMAG(CDLOG(CDETP(k)))/2.D0 + OFFSET
	  IF (DEL - PREVDEL .GT. (PI - EPS)) THEN
	    OFFSET = OFFSET - PI
	    DEL = DEL - PI
	  ENDIF
	  IF (PREVDEL - DEL .GT. (PI - EPS)) THEN
	    OFFSET = OFFSET + PI
	    DEL = DEL + PI
	  ENDIF
	  PREVDEL = DEL
	  WRITE(2, "(E19.12,1X,E19.12,1X,E19.12)")
     &	DBLE(k)*EMAX/DBLE(NS), DELTA(k), DEL
	ENDDO
	DEALLOCATE(DELTA, CDETP)
	CLOSE(2)
      ENDIF
      
      
      CAK=5.D-2
      CDETF = 1.D0
      DO i = 1, NSYM
	ALLOCATE (CS(NANG(i), NANG(i)),CEL(NANG(i)))
	CALL SSUM3D(L(:NANG(i),i),RADA,NSPS(i),NANG(i),CK(:NSPS(i),i),
     &		 CPHI(:NANG(i),:NSPS(i),i),CAK,CS)
	CALL CDETS(CS,CEL,CDET,NANG(i))
	PRINT *, "DET = ", CDET
	CDETF = CDETF * CDET
	DEALLOCATE(CS, CEL)
      ENDDO
      PRINT *, "CDETFULL=", AIMAG(CDLOG(CDETF))/2.D0
           
      
c      IF (.NOT. UCHECK(CS,NANG)) PRINT *, "UNITARITY PROBLEM"

C  END
C      DEALLOCATE(PIR,CPHI,CS,CEL)
 70   FORMAT(' ground state energy = ',E19.12)
 76   FORMAT(3(E19.12,1X))
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
        
        FUNCTION CONSOLE_INPUT(EMAX, NS, EIG_FROM_FILE)
	  LOGICAL, INTENT(OUT) :: EIG_FROM_FILE
	  REAL*8, INTENT(OUT) :: EMAX
	  INTEGER, INTENT(OUT) :: NS
	  LOGICAL CONSOLE_INPUT
	  CHARACTER(LEN=20) :: STR
	  
	  CONSOLE_INPUT = .FALSE.
	  IF (IARGC() .LT. 2) RETURN
	  CALL GETARG(1, STR)
	  READ(STR,*) EMAX
	  CALL GETARG(2, STR)
	  READ(STR,*) NS
	  IF (IARGC() .EQ. 2) THEN
	    EIG_FROM_FILE = .FALSE.
	  ELSE 
	    CALL GETARG(3, STR)
	    EIG_FROM_FILE = (STR .EQ. '-file')
	  ENDIF
	  CONSOLE_INPUT = .TRUE.
	  RETURN
        END FUNCTION
        END PROGRAM