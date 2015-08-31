      PROGRAM TESTBES
      COMPLEX*16, DIMENSION(1,1) :: CZER
      
      CALL BESZERR(CZER,1)
      DO i=1,1
	DO j=1,i
	  PRINT *, CZER(i,j)
	ENDDO
	PRINT *,
      ENDDO
      END