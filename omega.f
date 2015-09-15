      FUNCTION OMEGA(L,CZ)
C=======================================================================
C  
C-----------------------------------------------------------------------
      IMPLICIT REAL*8(A,B,D-H,O-Z)
      IMPLICIT COMPLEX*16 (C)      
      INTEGER :: L
      ALLOCATABLE :: CZER(:,:)
	OMEGA = 1.D0
      IF (L.EQ.0) RETURN
      ALLOCATE(CZER(L, L))
      CALL BESZERR(CZER, L)
      COMEGA = 0.D0
      DO i = 1, L
	COMEGA = COMEGA+CZER(L,i)/(CZ-CZER(L,i))/(DCONJG(CZ)-CZER(L,i))
      ENDDO
      OMEGA = OMEGA + DREAL(COMEGA)
      DEALLOCATE(CZER)
      RETURN
C     --------------------
      END
