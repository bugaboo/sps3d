      SUBROUTINE BESZERR(CZER,LMAX)
      IMPLICIT REAL*8 (A,B,D-H,O-Z)
      IMPLICIT COMPLEX*16 (C)
      DIMENSION CZER(LMAX,LMAX)

C --- Reading values
      OPEN(10, file='beszer')
      DO k=1, LMAX
	DO j=1, k
	  READ(10, *) X, Y
	  CZER(k,j)=DCMPLX(X,Y)
	ENDDO
      ENDDO
      CLOSE(10)
      RETURN 
      END
C=======================================================================