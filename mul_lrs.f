      SUBROUTINE MUL_LRS(NBAS,AS,VEC,WK,NBASD)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION AS(*),VEC(NBASD,*),WK(NBASD,*)
C
      CALL VMULSF(AS,NBAS,VEC,NBAS,NBASD,WK,NBASD)
      IJF=0
      DO JBF=1,NBAS
        DO IBF=1,JBF
          IJF=IJF+1
          SUM=0.D0
          DO k=1,NBAS
            SUM=SUM+VEC(k,IBF)*WK(k,JBF)
          ENDDO
          AS(IJF)=SUM
        ENDDO
      ENDDO
C
      RETURN
      END