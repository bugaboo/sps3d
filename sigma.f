      FUNCTION SIGMA(L,M,NANG,CS,AK,TET,PHI)
      IMPLICIT REAL*8(A,B,D-H,O-Z)
      IMPLICIT COMPLEX*16(C)
      PARAMETER(PI=3.141592653589793238462643D0)
      DIMENSION L(*), M(*), CS(NANG, *)
      
      SIGMA = 0.D0
      CK = DCMPLX(0.D0, 2.D0*PI/AK)
      DO i = 1, NANG
        CF = 0.D0
        DO j = 1, NANG
          CT = -CS(i, j)
          IF (i.EQ.j) CT = CT + 1.D0
          CT = CT * CK
          CF=CF+CDEXP(DCMPLX(0.D0,PI*L(j)/2.D0))*CT*
     &		 REHARM(L(j), M(j),TET,PHI)
        ENDDO
        SIGMA = SIGMA + CDABS(CF)**2
      ENDDO
      RETURN
      END FUNCTION