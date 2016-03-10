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
      ALLOCATABLE CKT(:), CET(:), CPHIT(:,:), CVECT(:,:), CSA(:,:,:)
      ALLOCATABLE L(:,:),M(:,:),PIR(:), NANG(:), NSPS(:), XK(:,:)
      ALLOCATABLE KSYMA(:),LANA(:),MANA(:),LT(:),MT(:),MULTIP(:)
      ALLOCATABLE DELTA(:),CDETP(:), SIGMAT(:), X(:, :), SIGMA2T(:,:)
      NAMELIST /INF3D/MODEL,RADA,NDVR,LMAX,LAN,MAN,KSYM,NTET,NPHI,
     &             KEYA,AMIN,AMAX,NUMA
      COMMON /POT_C/MODEL
      CHARACTER(LEN=30) :: SIGMANAME, FILENAME, TMPNAME
      LOGICAL CONS_INPUT, EIG_FROM_FILE, OX, OY
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
      KSYMA = 1
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
      
      NJ = 1
      NS = 1
      CONS_INPUT = CONSOLE_INPUT(NJ, NS, EIG_FROM_FILE, OX, OY)
      print *, 'eig from file', EIG_FROM_FILE
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
      
      SELECT CASE (NJ)
      CASE (1)
        ALLOCATE(SIGMAT(NS))
        SIGMAT = 0.D0
        STEP = PI / DBLE(NS)
        PHI = 0.D0
        AK = 0.2D0
        DO i = 1, NSYM
          print *, 'isym = ', i
          ALLOCATE (CS(NANG(i), NANG(i)),CEL(NANG(i)))
          CALL SSUM3D(L(:NANG(i),i),RADA,NSPS(i),NANG(i),
     &        CK(:NSPS(i),i),CPHI(:NANG(i),:NSPS(i),i),DCMPLX(AK),CS)
          DO k = 1, NS
            TET = DBLE(k) * STEP
            DO nu = 1, NANG(i)
              CF = 0
              DO mu = 1, NANG(i) 
                CT = -CS(nu, mu)
                IF (nu.EQ.mu) CT = CT + 1.D0
                CT = CT * (0.D0, 2.D0) * PI / AK * 
     &                CDEXP(DCMPLX(0.D0,PI*DBLE(L(mu, i))/2.D0))
                CF = CF + CT * REHARM(L(mu, i), M(mu, i), TET, PHI)
              ENDDO
              SIGMAT(k) = SIGMAT(k) + CDABS(CF)**2
            ENDDO
          ENDDO

          DEALLOCATE(CS, CEL)
        ENDDO
        WRITE(SIGMANAME,"('sigmatet',I0,'_',I0,'_',I0,'_',I0)") 
     &          MODEL, INT(RADA), NDVR, LMAX
        OPEN(1, FILE=SIGMANAME)

        DO k = 1, NS
          WRITE(1, 76) DBLE(k) * STEP, SIGMAT(k)
        ENDDO
        CLOSE(1)
      CASE (2)
        TETK = 0.D0
        PHIK = 0.D0
        ALLOCATE(SIGMAT(0:NS), CSA(NANGMAX, NANGMAX, NSYM))
        SIGMAT = 0.D0
        STEP = PI / DBLE(NS)
        PHI = 0.D0
        AK = 0.10D0
        DO i = 1, NSYM
          ALLOCATE (CS(NANG(i), NANG(i)),CEL(NANG(i)))
          CALL SSUM3D(L(:NANG(i),i),RADA,NSPS(i),NANG(i),
     &        CK(:NSPS(i),i),CPHI(:NANG(i),:NSPS(i),i),DCMPLX(AK),CS)
          DO nu = 1, NANG(i)
            DO mu = 1, NANG(i)
              CSA(nu, mu, i) = CS(nu, mu)
            ENDDO
          ENDDO
          DEALLOCATE(CS, CEL)
        ENDDO
        DO k = 0, NS
          CF = 0
          DO i = 1, NSYM
            TET = DBLE(k) * STEP
            DO nu = 1, NANG(i)
              DO mu = 1, NANG(i) 
                CT = -CSA(nu, mu, i)
                IF (nu.EQ.mu) CT = CT + 1.D0
                CT = CT * (0.D0, 2.D0) * PI / AK * 
     &                CDEXP(DCMPLX(0.D0,PI*DBLE(L(mu,i)-L(nu,i))/2.D0))
                CF = CF + CT * REHARM(L(mu, i), M(mu, i), TETK, PHIK) *
     &                REHARM(L(nu, i), M(nu, i), TET, PHI)
              ENDDO
            ENDDO
          ENDDO
          SIGMAT(k) = CDABS(CF)**2
          IF (k .EQ. 0) PRINT *, 'CF = ', CF
        ENDDO
        DEALLOCATE(CSA)
        
        WRITE(SIGMANAME,"('dsigmaZ',I0,'_',I0,'_',I0,'_',I0)")
     &          MODEL, INT(RADA), NDVR, LMAX
        OPEN(1, FILE=SIGMANAME)

        DO k = 0, NS
          WRITE(1, 76) DBLE(k) * STEP, SIGMAT(k)
        ENDDO
        CLOSE(1)
      CASE (3)
        TETK = 0.D0
        PHIK = 0.D0
        ALLOCATE(SIGMAT(0:NS))
        SIGMAT = 0.D0
        STEP = PI / DBLE(NS)
        PHI = 0.D0
        AK = 0.2D0
        DO k = 0, NS
          CF = 0
          DO i = 1, NSYM
            ALLOCATE (CS(NANG(i), NANG(i)),CEL(NANG(i)))
            CALL SSUM3D(L(:NANG(i),i),RADA,NSPS(i),NANG(i),
     &        CK(:NSPS(i),i),CPHI(:NANG(i),:NSPS(i),i),DCMPLX(AK),CS)

            TET = DBLE(k) * STEP
            DO nu = 1, NANG(i)
              DO mu = 1, NANG(i) 
                CT = -CS(nu, mu)
                IF (nu.EQ.mu) CT = CT + 1.D0
                CT = CT * (0.D0, 2.D0) * PI / AK * 
     &                CDEXP(DCMPLX(0.D0,PI*DBLE(L(mu,i)-L(nu,i))/2.D0))
                CF = CF + CT * REHARM(L(mu, i), M(mu, i), TETK, PHIK) *
     &                REHARM(L(nu, i), M(nu, i), TET, PHI)
              ENDDO

            ENDDO
            DEALLOCATE(CS, CEL)
          ENDDO


          SIGMAT(k) = CDABS(CF)**2
        ENDDO
        WRITE(SIGMANAME,"('dsigmaZ',I0,'_',I0,'_',I0,'_',I0)")
     &          MODEL, INT(RADA), NDVR, LMAX
        OPEN(1, FILE=SIGMANAME)

        DO k = 0, NS
          WRITE(1, 76) DBLE(k) * STEP, SIGMAT(k)
        ENDDO
        CLOSE(1)
      CASE (4)
        IF (OX) THEN
          TETK = PI / 2.D0
          PHIK = 0.D0
          TMPNAME = 'X'
        ELSE IF (OY) THEN        
          TETK = PI / 2.D0
          PHIK = PI / 2.D0
          TMPNAME = 'Y'
        ELSE
          TETK = 0.D0
          PHIK = 0.D0
          TMPNAME = 'Z'
        ENDIF
        ALLOCATE(SIGMA2T(0:NS, 0:2*NS), CSA(NANGMAX, NANGMAX, NSYM))
        SIGMAT = 0.D0
        STEPT = PI / DBLE(NS)
        STEPP = PI / DBLE(NS)
        PHI = 0.D0
        AK = 0.10D0
        DO i = 1, NSYM
          ALLOCATE (CS(NANG(i), NANG(i)),CEL(NANG(i)))
          CALL SSUM3D(L(:NANG(i),i),RADA,NSPS(i),NANG(i),
     &        CK(:NSPS(i),i),CPHI(:NANG(i),:NSPS(i),i),DCMPLX(AK),CS)
          DO nu = 1, NANG(i)
            DO mu = 1, NANG(i)
              CSA(nu, mu, i) = CS(nu, mu)
            ENDDO
          ENDDO
          DEALLOCATE(CS, CEL)
        ENDDO
        ALLOCATE(XK(NANGMAX, NSYM))
        XK = 0
        DO i = 1, NSYM
          DO nu = 1, NANG(i)
            XK(nu, i) = REHARM(L(nu, i), M(nu, i), TETK, PHIK)
          ENDDO
        ENDDO
        DO k = 0, NS
          DO kp = 0, 2*NS
            CF = 0
            TET = DBLE(k) * STEPT
            PHI = DBLE(kp) * STEPP
            DO i = 1, NSYM
              DO nu = 1, NANG(i)
                XR = REHARM(L(nu, i), M(nu, i), TET, PHI)
                DO mu = 1, NANG(i) 
                  CT = -CSA(nu, mu, i)
                  IF (nu.EQ.mu) CT = CT + 1.D0
                  CT = CT * (0.D0, 2.D0) * PI / AK * 
     &                CDEXP(DCMPLX(0.D0,PI*DBLE(L(mu,i)-L(nu,i))/2.D0))
                  CF = CF+ CT*XK(mu, i)*XR        
                ENDDO
              ENDDO
            ENDDO
            SIGMA2T(k, kp) = CDABS(CF)**2            
            IF (k .EQ. 0 .AND. kp .EQ. 0) PRINT *, 'CF = ', CF
          ENDDO
        ENDDO
        DEALLOCATE(CSA, XK)
        
        WRITE(SIGMANAME,"('dsigma',A1,I0,'_',I0,'_',I0,'_',I0)")
     &          TMPNAME,MODEL, INT(RADA), NDVR, LMAX
        OPEN(1, FILE=SIGMANAME)

        DO k = 0, NS
          DO kp = 0, 2*NS
            WRITE(1, 76) DBLE(k) * STEPT, DBLE(kp) * STEPP,SIGMA2T(k,kp)
          ENDDO
        ENDDO
        CLOSE(1)
        WRITE(SIGMANAME,"('dsigma_matrix',A1,I0,'_',I0,'_',I0,'_',I0)")
     &          TMPNAME,MODEL, INT(RADA), NDVR, LMAX
        OPEN(1, FILE=SIGMANAME)
        DO k = 0, NS
          DO kp = 0, 2*NS
            WRITE(1, 76, ADVANCE = 'no') SIGMA2T(k,kp)
          ENDDO
          WRITE(1, *)
        ENDDO
        END SELECT

 76   FORMAT(3(E19.12,1X))
      
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
        
        FUNCTION CONSOLE_INPUT(NJ, NS, EIG_FROM_FILE, OX, OY)
        LOGICAL, INTENT(OUT) :: EIG_FROM_FILE
        INTEGER, INTENT(OUT) :: NJ
        LOGICAL CONSOLE_INPUT, OX, OY
        CHARACTER(LEN=20) :: STR

        CALL GETARG(IARGC(), STR)
        print *, str
        EIG_FROM_FILE = (STR .EQ. '-file')
        print *, EIG_FROM_FILE
        
        CONSOLE_INPUT = .FALSE.
        IF (IARGC() .LT. 1) RETURN
        CALL GETARG(1, STR)
        READ(STR,*) NJ

        CONSOLE_INPUT = .TRUE.
        
        IF (IARGC() .LT. 2) RETURN
        CALL GETARG(2, STR)
        READ(STR, *) NS
        IF (IARGC() .GT. 2) THEN
          CALL GETARG(3, STR)
          OX = (STR .EQ. 'ox')
          OY = (STR .EQ. 'oy')
        ENDIF
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
      END