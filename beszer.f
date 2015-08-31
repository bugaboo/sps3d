      SUBROUTINE BESZER(N,CZER)
C=======================================================================
C  KEY.EQ.0 - zeros of the Bessel polynomial y_N(z)
C  KEY.NE.0 - zeros of the reverse Bessel polynomial theta_N(z)
C
C  THIS ALGORITHM WORKS ONLY FOR N.LE.28!!!
C
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A,B,D-H,O-Z)
      IMPLICIT COMPLEX*16 (C)
      PARAMETER(NMAX=28,KTEST=0)
      DIMENSION CZER(*),P(0:NMAX),Q(NMAX),IND(NMAX)
C
      IF(N.LT.1) RETURN
      IF(N.EQ.1) THEN
        CZER(1)=-1.D0
        RETURN
      ENDIF
      IF(N.GT.NMAX) STOP ' *** BESZER: N > NMAX=28 ERROR'
C --- Skipping non-relevant zeroes
      OPEN(10, file='beszer')
      READ(10, *) X, Y
      IF (N.GT.2) THEN
	DO k=2, N-1
	  DO j=1, k
	    READ(10, *) X, Y
	  ENDDO
	ENDDO
      ENDIF
C --- Reading values
      DO k=1, N
	READ(10, *) X, Y
	CZER(k)=DCMPLX(X,Y)
      ENDDO
      CLOSE(10)
	
C --- Odering the zeros in decreasing order of Im(ZER)
      DO i=1,N
        P(i)=DREAL(CZER(i))
        Q(i)=-DIMAG(CZER(i))
      ENDDO
      CALL INDEXX(N,Q(1),IND)
      DO i=1,N
        ii=IND(i)
        CZER(i)=DCMPLX(P(ii),-Q(ii))
      ENDDO
C
 71   FORMAT('N=',I2,'  BESZER TEST: MAX|BES(ZER)|=',E9.3)
      RETURN 
      END
C=======================================================================