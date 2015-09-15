      FUNCTION REHARM(L,M,TET, PHI)
C=======================================================================
C  
C-----------------------------------------------------------------------
      IMPLICIT REAL*8(A,B,D-H,O-Z)     
      PARAMETER(PI=3.141592653589793238462643D0)
      INTEGER :: L, M, MP
      
      IF (M.EQ.0) THEN
	REHARM = POLJ(0,0,L,DCOS(TET))/DSQRT(2*PI)
	RETURN
      ENDIF
      IF (M.GT.0) THEN
	MP = M
	REHARM = DCOS(MP * PHI) / DSQRT(PI)
      ELSE
	MP = -M
	REHARM = DSIN(MP * PHI) / DSQRT(PI)
      ENDIF
      TMP = 1 - DCOS(TET)**2
      DO i = 1, MP / 2
	REHARM = REHARM * TMP
      ENDDO
      IF (MOD(MP, 2).EQ.1) REHARM = REHARM * DSQRT(TMP)
      REHARM = REHARM * POLJ(DBLE(MP), DBLE(MP), L - MP,DCOS(TET))
      RETURN
C     --------------------
      END
