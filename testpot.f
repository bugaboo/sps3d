      PROGRAM TEST
      USE ANGSYM
      REAL*8, ALLOCATABLE :: V(:,:)
      PARAMETER(LMAX=2)
      INTEGER, DIMENSION(:), ALLOCATABLE :: L, M

      COMMON /POT_C/ MODEL

      MODEL = 303
      CALL ANGBAS(0,L,M,LMAX,NANG)
      ALLOCATE (V(NANG,NANG))
      CALL POTMAT(200, 200, 5.D0, V, NANG, L, M, LMAX)
      DO i = 1,NANG
	DO j = 1,NANG
	  WRITE (*,77,ADVANCE='no') V(i,j)
	ENDDO
	PRINT *,
      ENDDO
      DEALLOCATE(V, L, M)
 77   FORMAT(99(E19.12,1X))
      END
