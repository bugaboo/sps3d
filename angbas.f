      MODULE ANGSYM
      CONTAINS
	SUBROUTINE ANGBAS(KSYM, L, M, LMAX, NANG, LAN, MAN)
	IMPLICIT REAL*8 (A-H,O-Z)
	INTEGER, INTENT(OUT) :: NANG
	INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: L, M

	n = 1
C
C     	KSYM = 0 no symmetry
C    	
	SELECT CASE(KSYM)
	CASE(0)
	  NANG = (LMAX + 1)**2
	  ALLOCATE(L(NANG), M(NANG))
	  DO i = 0, LMAX
              DO k = -i, i
                  L(n) = i
                  M(n) = k
                  n = n + 1
              ENDDO
          ENDDO  
C  KSYM=1 Angular symmetry          
	CASE(1)
	  NANG = LMAX+1
	  ALLOCATE(L(NANG), M(NANG))
	  DO i=1,NANG
	    L(i) = i-1
	    M(i) = MAN
	  ENDDO
C  KSYM=2 Spherical symmetry	  
	CASE(2)
	  NANG = 1
	  ALLOCATE(L(1),M(1))
	  L(1)=LAN
	  M(1)=0
	END SELECT
	RETURN
      END SUBROUTINE
C      
      END MODULE ANGSYM