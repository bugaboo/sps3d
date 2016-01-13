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
      ALLOCATABLE CZER(:,:)
      ALLOCATABLE L(:),M(:),PIR(:)
      NAMELIST /INF3D/MODEL,RADA,NDVR,LMAX,LAN,MAN,KSYM,NTET,NPHI,
     & 		KEYA,AMIN,AMAX,NUMA
      COMMON /POT_C/MODEL
      CHARACTER(LEN=30) :: SPSNAME, FILENAME, TMPNAME, SYMNAME
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
      
      
      NK = 1
      DO i = 1, NSPS
        IF (DABS(1 - 0.06D0 / DREAL(CK(i))) .LT. 0.25D0 .AND.
     &   DABS(DREAL(CK(i)) / DIMAG(CK(i))) .GT. 2.5D0) THEN
        NK = i
        ENDIF
      ENDDO
      GAMMA = 0.D0
      CVECS = 0.D0
      CVECSP = 0.D0
      VEC30 = 0.D0
      VEC3L = 0.D0
      ALLOCATE (CZER(LMAX,LMAX))
      CALL BESZERR(CZER,LMAX)
      DO i = 1, NANG
        TMP = DREAL(CK(NK)) * OMEGA(L(i), DCMPLX(0.D0, -RADA)
     &               * CK(NK)) * CDABS(CPHI(i, NK))**2
        GAMMA = GAMMA + TMP
        PRINT *, 'GAMMA', L(i), TMP
        DO k = 1, NDVR
          CVECS = CVECS + CVEC((i-1) * NDVR + k, NK)**2
        VEC3L = VEC3L + DIMAG(CK(NK)) * RADA * 
     &          CDABS(CVEC((i-1) * NDVR + k, NK))**2
        ENDDO
        CTMP=1.D0
        DO ip=1,L(i)
          CTMP=CTMP+CZER(L(i),ip)/((0.D0,1.D0)*CK(NK)*RADA+
     &      CZER(L(i),ip))**2
        ENDDO   
        CVECSP = CVECSP + (0.D0, 0.5D0) * CTMP * CPHI(i,NK)**2 / CK(NK)
        VEC30 = VEC30 + OMEGA(L(i),DCMPLX(0.D0, -RADA) * CK(NK)) * 
     &          CDABS(CPHI(i,NK))**2

      ENDDO
C      PRINT *, "CPHI", CPHI(1, NK)**2 / RADA**2
C      PRINT *, CK(NK), GAMMA / 2.D0 / DIMAG(CK(NK)) / DREAL(CK(NK))
C      PRINT *, 'CVECS', CVECS * RADA / 2.D0
C      PRINT *, '30', VEC30, VEC3L * 
      
      
      IF(NB.GT.0) WRITE(*,70) DREAL(CE(NB))
      IF (KSYM .NE. 2) THEN
        FILENAME = SYMNAME('L',MODEL,KSYM,NDVR,LMAX,INT(RADA),MAN)
        PRINT *, FILENAME
        OPEN(UNIT = 1, FILE=FILENAME, FORM = 'unformatted', 
     &          ACCESS = 'direct', RECL = NSPS * 16)
        WRITE(1, REC = 1) CK
        CLOSE(1)
        FILENAME = SYMNAME('C',MODEL,KSYM,NDVR,LMAX,INT(RADA),MAN)
        OPEN(UNIT = 1, FILE=FILENAME, FORM = 'unformatted', 
     &          ACCESS = 'direct', RECL = NSPS * NANG * 16)
        WRITE(1, REC = 1) CPHI
        CLOSE(1)
      ENDIF
      CAK=0.3D0
      CALL SSUM3D(L,RADA,NSPS,NANG,CK,CPHI,CAK,CS)
C      IF (.NOT. UCHECK(CS,NANG)) PRINT *, "UNITARITY PROBLEM"
      CALL CDETS(CS,CEL,CDET,NANG)
      PRINT *, "Sum det", CDET, AIMAG(CDLOG(CDET))/2.D0
      CALL SPRO(RADA,NSPS,CK,CAK,CDETP)
      CDETP = CDETP*CDEXP(-(0.D0,2.D0)*CAK*RADA*(NANG-1))
      PRINT *, "Prod det", CDETP, AIMAG(CDLOG(CDETP))/2.D0
      PRINT *, "Slength", SLENGTH(RADA,NSPS,NANG,CK)
      
C  Printing eigenvalues

      WRITE(SPSNAME,"(A7,I0,A1,I0,A1,I2.2,A1,I0,'_'I0)")
     & "spseig_",KSYM,'_',NDVR,'_',LMAX,'_',INT(RADA),MAN
      OPEN(1,FILE=SPSNAME)
      DO n=1,NSPS
        write(1,77) CK(n),CE(n)
      ENDDO
      CLOSE(1)
C  END
      DEALLOCATE(PIR,CK,CE,CVEC,CPHI,CS,CEL)
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
              IF (CDABS(CTMP) .GT. 1.D-5) THEN
                PRINT *, "NONUNITARY S MATRIX", i, j, CTMP
              ENDIF
            ENDDO
          ENDDO
          UCHECK = .TRUE.
          RETURN 
        END FUNCTION
        
        FUNCTION CONSOLE_INPUT(EMAX, NS)
          LOGICAL CONSOLE_INPUT
          CHARACTER(LEN=20) :: STR
          
          CONSOLE_INPUT = .FALSE.
          IF (IARGC() .EQ. 0) RETURN
          CALL GETARG(1, STR)
          READ(STR,*) EMAX
          CALL GETARG(2, STR)
          READ(STR,*) NS
          CONSOLE_INPUT = .TRUE.
          RETURN
        END FUNCTION
        
      END PROGRAM
