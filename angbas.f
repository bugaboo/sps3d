      MODULE ANGSYM
      CONTAINS
	SUBROUTINE ANGBAS(KSYM, L, M, LMAX, NANG, LAN, MAN)
	IMPLICIT REAL*8 (A-H,O-Z)
	LOGICAL :: SX, SY, SZ
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
C  KSYM=90 Inversion symmetry even
	  CASE(9)
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
C  KSYM=-90 Inversion symmetry odd
	  CASE(-9)
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
C KSYM = 10, 20, ... 80
	  CASE (1:8)
	    KT = KSYM / 10
	    SZ = (MOD(KT, 2).EQ.1)
	    KT = KT / 2
	    SY = (MOD(KT, 2).EQ.1)
	    KT = KT / 2
	    SX = (MOD(KT, 2).EQ.1)
	    NANG = 0
	    DO i = 0, LMAX
	      DO j = -i, i
		IF (SYMX(i,j,SX).AND.SYMY(i,j,SY).AND.SYMZ(i,j,SZ)) THEN
		  NANG = NANG + 1
		ENDIF
	      ENDDO
	    ENDDO
	    ALLOCATE(L(NANG), M(NANG))
	    k = 1
	    DO i = 0, LMAX
	      DO j = -i, i
		IF (SYMX(i,j,SX).AND.SYMY(i,j,SY).AND.SYMZ(i,j,SZ)) THEN
		  L(k) = i
		  M(k) = j
		  k = k + 1
		ENDIF
	      ENDDO
	    ENDDO
	  END SELECT
C  KSYM=1 Angular symmetry          
	CASE(1)
	  NANG = LMAX+1-ABS(MAN)
	  ALLOCATE(L(NANG), M(NANG))
	  DO i=1,NANG
	    L(i) = i-1+ABS(MAN)
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
      
      FUNCTION SYMX(L, M, B)
        LOGICAL :: SYMX, B
        INTEGER :: L, M                    
        SYMX = MOD(M, 2) .EQ. 0
        IF (M.GE.0 .NEQV. B) THEN
          SYMX = (.NOT.SYMX)
        ENDIF
        RETURN
      END FUNCTION
      
      FUNCTION SYMY(L, M, B)
        LOGICAL :: SYMY, B
        INTEGER :: L, M
        SYMY = (M.GE.0)
        IF (.NOT.B) SYMY = (.NOT.SYMY)
        RETURN
      END FUNCTION
C      
      FUNCTION SYMZ(L, M, B)
        LOGICAL :: SYMZ, B
        INTEGER :: L, M
        SYMZ = (MOD(L-M,2).EQ.0)
        IF (.NOT.B) SYMZ = (.NOT.SYMZ)
        RETURN
      END FUNCTION
      END MODULE ANGSYM