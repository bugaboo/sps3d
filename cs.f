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
      ALLOCATABLE DELTA(:),CDETP(:), SIGMAT(:), X(:, :)
      NAMELIST /INF3D/MODEL,RADA,NDVR,LMAX,LAN,MAN,KSYM,NTET,NPHI,
     &             KEYA,AMIN,AMAX,NUMA
      COMMON /POT_C/MODEL
      CHARACTER(LEN=30) :: SIGMANAME, FILENAME, TMPNAME
      LOGICAL CONS_INPUT, EIG_FROM_FILE, OX, OY
      REAL*8 :: KMAX
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
      
      CONS_INPUT = CONSOLE_INPUT(KMAX, NS, EIG_FROM_FILE, OX, OY)
      IF (EIG_FROM_FILE) THEN
        DO i = 1, NSYM
          ALLOCATE(CKT(NSPS(i)), CPHIT(NANG(i),NSPS(i)))
          TMPNAME = SYMNAME(MODEL,KSYMA(i),NDVR,LMAX,INT(RADA),MANA(i))
          FILENAME = 'eigenval' // TMPNAME
          PRINT *, FILENAME
          OPEN(UNIT = 1, FILE=FILENAME, FORM = 'unformatted', 
     &            ACCESS = 'direct', RECL = NSPS(i) * 16)
          READ(1, REC = 1) CKT
          CLOSE(1)

          FILENAME = 'eigenvec' // TMPNAME
          OPEN(UNIT = 1, FILE=FILENAME, FORM = 'unformatted', 
     &      ACCESS = 'direct', RECL = NSPS(i) * NANG(i) * 16)
          READ(1, REC = 1) CPHIT
          CLOSE(1)
          CK(:NSPS(i),i) = CKT
          CPHI(:NANG(i),:NSPS(i),i) = CPHIT
          DEALLOCATE(CKT,CPHIT)
        ENDDO
      ELSE
      DO i = 1, NSYM
        ALLOCATE(CKT(NSPS(i)),CE(NSPS(i)),CVEC(NDVR*NANG(i),NSPS(i)))
        ALLOCATE(CPHIT(NANG(i),NSPS(i)))
        CALL SPS3D(LMAX,NANG(i),L(1:NANG(i),i),M(1:NANG(i),i),RADA,NDVR,  
     &                NTET,NPHI,PIR,CKT,CE,CVEC,CPHIT,NB,NA,NOI) 
        IF (KSYM .NE. 2) THEN
          TMPNAME = SYMNAME(MODEL,KSYMA(i),NDVR,LMAX,INT(RADA),MANA(i))
          FILENAME = 'eigenval' // TMPNAME
          PRINT *, FILENAME
          OPEN(UNIT = 1, FILE=FILENAME, FORM = 'unformatted', 
     &            ACCESS = 'direct', RECL = NSPS(i) * 16)
          WRITE(1, REC = 1) CKT
          CLOSE(1)
          FILENAME = 'eigenvec' // TMPNAME
          OPEN(UNIT = 1, FILE=FILENAME, FORM = 'unformatted', 
     &            ACCESS = 'direct', RECL = NSPS(i) * NANG(i) * 16)
          WRITE(1, REC = 1) CPHIT
          CLOSE(1)
        ENDIF
        CK(:NSPS(i),i) = CKT
        CPHI(:NANG(i),:NSPS(i),i) = CPHIT
        DEALLOCATE(CKT,CE,CVEC,CPHIT)
      ENDDO
      ENDIF

      IF (CONS_INPUT) THEN
        
        ALLOCATE(X(NSYM, NANGMAX))
        IF (OX .OR. OY) THEN
          IF (OX) PHI = 0
          IF (OY) PHI = PI / 2
          DO i = 1, NSYM
            DO nu = 1, NANG(i)
              X(i, nu) = REHARM(L(nu, i), M(nu, i), PI / 4.D0, PHI)
            ENDDO
          ENDDO
        ENDIF
        
        ALLOCATE(SIGMAT(NS))
        SIGMAT = 0.D0
        DELTA = 0.D0
        CDETP = 1.D0
          OPEN(3, FILE = 'scat')
        DO i = 1, NSYM
          print *, 'isym = ', i
          ALLOCATE (CS(NANG(i), NANG(i)),CEL(NANG(i)))
          DO k = 1, NS
            AK =  DBLE(k) * KMAX / NS
            CALL SSUM3D(L(:NANG(i),i),RADA,NSPS(i),NANG(i),
     &          CK(:NSPS(i),i),CPHI(:NANG(i),:NSPS(i),i),DCMPLX(AK),CS)
c            IF (KSYM.EQ.2) THEN
c              DO m2 = -L(1, i), L(1, i)
c                SIGMAT(k) = SIGMAT(k) + SIGMA(L(:NANG(i),i),
c     &                            m2,NANG(i),CS,AK,0.D0,0.D0)
c              ENDDO
c            ELSE
c              SIGMAT(k) = SIGMAT(k) + SIGMA(L(:NANG(i),i),M(:NANG(i),i),
c     &                                    NANG(i),CS,AK,0.D0,0.D0)
            DO nu = 1, NANG(i)
              DO mu = 1, NANG(i)
                WRITE(1, *) CS(nu, mu)
              ENDDO
            ENDDO
            DO nu = 1, NANG(i)
              CF = 0
              DO mu = 1, NANG(i) 
                CT = -CS(nu, mu)
                IF (nu.EQ.mu) CT = CT + 1.D0
                CT = CT * (0.D0, 2.D0) * PI / AK * 
     &                CDEXP(DCMPLX(0.D0,PI*DBLE(L(mu, i))/2.D0))
                IF (OX .OR. OY) THEN
                  CF = CF + CT * X(i, mu)
                ELSE
                  IF (M(mu, i).EQ.0) THEN
                    CF = CF + CT*DSQRT(DBLE(2 * L(mu, i) + 1)/4.D0/PI)
                  ENDIF
                ENDIF
              ENDDO
              SIGMAT(k) = SIGMAT(k) + CDABS(CF)**2
            ENDDO
          ENDDO

          DEALLOCATE(CS, CEL)
        ENDDO
        CLOSE(3)        
        IF (OX) THEN
          TMPNAME = 'sigmaox'
        ELSE IF (OY) THEN
          TMPNAME = 'sigmaoy'
        ELSE
          TMPNAME = 'sigmaoz'
        ENDIF
      ENDIF
      
      WRITE(SIGMANAME,"(A7,I0,'_',I0,'_',I0,'_',I0)") 
     &          TMPNAME, MODEL, INT(RADA), NDVR, LMAX
c      SIGMANAME = TMPNAME // SIGMANAME
      PRINT *, SIGMANAME
      OPEN(1, FILE=SIGMANAME)
      DO k = 1, NS
        WRITE (1, 76) DBLE(k) * KMAX / NS, SIGMAT(k)
      ENDDO
      CLOSE(1)           
      
      DEALLOCATE(X)
c      IF (.NOT. UCHECK(CS,NANG)) PRINT *, "UNITARITY PROBLEM"

C  END
C      DEALLOCATE(PIR,CPHI,CS,CEL)
 70   FORMAT(' ground state energy = ',E19.12)
 76   FORMAT(3(E19.12,1X))
 77   FORMAT(4(E19.12,1X))
 78   FORMAT((E19.12,1X,E19.12,1X,E19.12))
 
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
        
        FUNCTION CONSOLE_INPUT(KMAX, NS, EIG_FROM_FILE, OX, OY)
        LOGICAL, INTENT(OUT) :: EIG_FROM_FILE, OX, OY
        REAL*8, INTENT(OUT) :: KMAX
        INTEGER, INTENT(OUT) :: NS
        LOGICAL CONSOLE_INPUT
        CHARACTER(LEN=20) :: STR
        
        CONSOLE_INPUT = .FALSE.
        OX = .FALSE.
        OY = .FALSE.
        EIG_FROM_FILE = .FALSE.
        IF (IARGC() .LT. 2) RETURN
        CALL GETARG(1, STR)
        READ(STR,*) KMAX
        CALL GETARG(2, STR)
        READ(STR,*) NS
        IF (IARGC() .EQ. 2) THEN
          EIG_FROM_FILE = .FALSE.
        ELSE 
          CALL GETARG(IARGC(), STR)
          EIG_FROM_FILE = (STR .EQ. '-file')
          CALL GETARG(3, STR)
          IF (STR .EQ. 'ox') OX = .TRUE.
          IF (STR .EQ. 'oy') OY = .TRUE.
        ENDIF
        CONSOLE_INPUT = .TRUE.
        RETURN
        END FUNCTION
        FUNCTION SYMNAME(MODEL, KSYM, N, L, R, MAN) RESULT(RES)
        CHARACTER(20) :: RES
        INTEGER, INTENT(IN) :: MODEL, KSYM, N, L, R, MAN
        SELECT CASE(KSYM)
        CASE (1)
        WRITE(RES, 78) MODEL, N, L, R, MAN
        CASE (-90, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
        WRITE(RES, 78) MODEL, KSYM, N, L, R
        END SELECT
        RETURN
 78   FORMAT(I0,'_',I0,'_',I0,'_',I0,'_',I0)
        END FUNCTION
        END PROGRAM
