      SUBROUTINE IKESTEP(KEYA1,AMIN1,AMAX,NUMA)
C=======================================================================
C  Momentum K and energy E step generator for spectral calculations
C-----------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      SAVE KEYA,AMIN,DA
C
      KEYA=KEYA1
      AMIN=AMIN1
      IF(KEYA.LT.0) THEN
        DA=DLOG10(AMAX/AMIN)/DBLE(NUMA)
      ELSE
        DA=(AMAX-AMIN)/DBLE(NUMA)
      ENDIF
      RETURN 
C     ----------------------
      ENTRY KESTEP(IA,AK,AE)
C     ----------------------
      SELECT CASE(KEYA)
C--------------------------
C  KEYA=1: equidistant in K
C--------------------------
      CASE(1)
      AK=AMIN+IA*DA
      AE=0.5D0*AK*AK
C--------------------------
C  KEYA=2: equidistant in E
C--------------------------
      CASE(2)
      AE=AMIN+IA*DA
      AK=DSQRT(2.D0*AE)
C-------------------------------
C  KEYA=-1: equidistant in lg(K)
C-------------------------------
      CASE(-1)
      AK=AMIN*10.D0**(IA*DA)
      AE=0.5D0*AK*AK
C-------------------------------
C  KEYA=-2: equidistant in lg(E)
C-------------------------------
      CASE(-2)
      AE=AMIN*10.D0**(IA*DA)
      AK=DSQRT(2.D0*AE)
      END SELECT
C
      RETURN 
      END
C=======================================================================