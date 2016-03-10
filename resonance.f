      PROGRAM RESONANCE
C
C  Finds and prints resonance in precomputed file 
C  
c  DELTARE - halfwidth of interval in which look for resonance with center in REK
c  COTANK  - minimum Re(k)/Im(k) relation for resonance
C  Other parameters are defined from command line arguments. 
C  Usage:
C  user@host# resonance model rada ndvr lmax lan man ksym REK
C
      USE ANGSYM
      PARAMETER(DELTARE = 0.1D0)          
      PARAMETER(COTANK = 1.5D0)
      REAL*8 :: RADA, REK, GAMMA, OMEGA
      INTEGER :: MODEL, NDVR, LMAX, MAN, LAN, KSYM, NANG, NSPS, NK,NUMK
      INTEGER, ALLOCATABLE :: L(:), M(:)
      COMPLEX*16, ALLOCATABLE :: CK(:), CPHI(:, :)
      CHARACTER(30) :: STR, FILENAME, SYMNAME
      
      IF (IARGC() .LT. 8) STOP 'AT LEAST 7 PARAMETERS REQUIRED'
      CALL GETARG(1, STR)
      READ(STR,*) MODEL
      CALL GETARG(2, STR)
      READ(STR,*) RADA
      CALL GETARG(3, STR)
      READ(STR,*) NDVR
      CALL GETARG(4, STR)
      READ(STR,*) LMAX
      CALL GETARG(5, STR)
      READ(STR,*) LAN
      CALL GETARG(6, STR)
      READ(STR,*) MAN
      CALL GETARG(7, STR)
      READ(STR,*) KSYM
      CALL GETARG(8, STR)
      READ(STR,*) REK
      
      CALL ANGBAS(KSYM,L,M,LMAX,NANG,LAN,MAN)
      NSPS = 0
      DO i= 1, NANG
        NSPS = NSPS + 2 * NDVR + L(i)
      ENDDO
     
      ALLOCATE(CK(NSPS), CPHI(NANG,NSPS))
      FILENAME = SYMNAME('L',MODEL,KSYM,NDVR,LMAX,INT(RADA),MAN)
      OPEN(UNIT=1,FILE=FILENAME,
     &     FORM = 'unformatted', ACCESS = 'direct', RECL = NSPS * 16)
      READ(1, REC = 1) CK
      CLOSE(1)
      OPEN(UNIT=1,FILE=SYMNAME('C',MODEL,KSYM,NDVR,LMAX,INT(RADA),MAN),
     &     FORM = 'unformatted', ACCESS='direct', RECL=NSPS*NANG*16)
      READ(1, REC = 1) CPHI
      CLOSE(1)
      
      NUMK = 0     
      DO i = 1, NSPS
        IF (DABS(1 - REK / DREAL(CK(i))) .LT. DELTARE .AND.
     &   DABS(DREAL(CK(i)) / DIMAG(CK(i))) .GT. COTANK) THEN
        NUMK = NUMK + 1
        NK = i
        ENDIF
      ENDDO
      IF (NUMK .EQ. 0) STOP 'STATE NOT FOUND'
      IF (NUMK .GT. 1) STOP 'MORE THAN ONE STATE FOUND'
      GAMMA = 0.D0
      DO i = 1, NANG
        GAMMA = GAMMA + DREAL(CK(NK)) * OMEGA(L(i), DCMPLX(0.D0, -RADA)
     &               * CK(NK)) * CDABS(CPHI(i, NK))**2
      ENDDO
      PRINT *, "Gamma", GAMMA, GAMMA/(2.D0*DREAL(CK(NK))*DIMAG(CK(NK)))
      DEALLOCATE(CK, CPHI)
      DEALLOCATE(L, M)
 
      END
