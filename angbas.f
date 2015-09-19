      MODULE ANGSYM
      CONTAINS
	SUBROUTINE ANGBAS(KSYM, L, M, LMAX, NANG, LAN, MAN)
	IMPLICIT REAL*8 (A-H,O-Z)
	INTEGER, INTENT(OUT) :: NANG
	INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: L, M

	n = 1
C
C  KSYM = 0 no symmetry
C    	
	SELECT CASE(MOD(KSYM, 10))
	CASE(0)
	  SELECT CASE(KSYM / 10)
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
C  KSYM=20 Inversion symmetry even
	  CASE(2)
	  LTMP=LMAX/2
	  NANG=(LTMP+1)*(2*LTMP+1)
	  ALLOCATE(L(NANG), M(NANG))
	  DO i = 0, LMAX, 2
              DO k = -i, i
                  L(n) = i
                  M(n) = k
                  n = n + 1
              ENDDO
          ENDDO
C  KSYM=10 Inversion symmetry odd
	  CASE(1)
	    LTMP=LMAX/2
	    IF (MOD(LMAX,2).EQ.0) LTMP = LTMP - 1
	    NANG=(LTMP+1)*(2*LTMP+3)
	    ALLOCATE(L(NANG), M(NANG))
	    DO i = 1, LMAX, 2
		DO k = -i, i
		    L(n) = i
		    M(n) = k
		    n = n + 1
		ENDDO
	    ENDDO
	  END SELECT
C  KSYM=1 Angular symmetry          
	CASE(1)
	  NANG = LMAX+1-MAN
	  ALLOCATE(L(NANG), M(NANG))
	  DO i=1,NANG
	    L(i) = i-1+MAN
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