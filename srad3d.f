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
      CAK=0.3D0
      CALL SSUM3D(L,RADA,NSPS,NANG,CK,CPHI,CAK,CS)
      IF (.NOT. UCHECK(CS,NANG)) PRINT *, "UNITARITY PROBLEM"
      CALL CDETS(CS,CEL,CDET,NANG)
      PRINT *, "Sum det", CDET, AIMAG(CDLOG(CDET))/2.D0
      CALL SPRO(RADA,NSPS,CK,CAK,CDETP)
      CDETP = CDETP*CDEXP(-(0.D0,2.D0)*CAK*RADA*(NANG-1))
      PRINT *, "Prod det", CDETP, AIMAG(CDLOG(CDETP))/2.D0
      
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
