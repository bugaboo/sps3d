      FUNCTION SYMNAME(MODE, MODEL, KSYM, N, L, R, MAN) RESULT(RES)
C  Filename for storing eigenvectors and eigenvalues
C  Modes: C - eigenvectors, L - eigenvalues
      CHARACTER(30) :: RES, TMP
      INTEGER, INTENT(IN) :: MODEL, KSYM, N, L, R, MAN
      CHARACTER, INTENT(IN) :: MODE
      SELECT CASE(KSYM)
      CASE (1)
        WRITE(TMP, 78) MODEL, N, L, R, MAN
      CASE (-90, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
        WRITE(TMP, 78) MODEL, KSYM, N, L, R
      END SELECT
      SELECT CASE(MODE)
      CASE ('C')
        RES = 'eigenvec' // TMP
      CASE ('L')
        RES = 'eigenval' // TMP
      END SELECT
      RETURN
 78   FORMAT(I0,'_',I0,'_',I0,'_',I0,'_',I0)
      END
