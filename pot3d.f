      FUNCTION POT3D(RAD,TET,PHI)
C=======================================================================
C  Potential energy V(r) for one-channel SPS EVP in the R problem
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /POT_C/MODEL
C
      SELECT CASE(MODEL)
C-----------------------------
C  model 1: Coulomb potential
C-----------------------------
      CASE(1)
      POT3D=-1.D0/RAD
C----------------------------
C  model 2: Yukawa potential
C----------------------------
      CASE(2)
      TMP=RAD*DSQRT((DSIN(TET)*DCOS(PHI)/1.D0)**2+
     & (DSIN(TET)*DSIN(PHI)/2.D0)**2+(DCOS(TET)/3.D0)**2)
      POT3D=-DEXP(-0.1D0*TMP)/TMP
C      
C----------------------------
C  model 300: Angle test
C----------------------------
      CASE(300)
      X = DCOS(TET)
c      POT3D = DSQRT(1-X**2)*POLJ(1.D0,1.D0,0.D0,X)*DSIN(PHI)/
c     &	 DSQRT(4*DATAN(1.D0))
      POT3D = DSIN(PHI)*DSIN(TET)*DSQRT(.75D0)
C----------------------------
C  model 301: Angle test
C----------------------------
      CASE(301)
      POT3D = DCOS(TET)*DSQRT(3.D0/DATAN(1.D0))/4.D0
      
C----------------------------
C  model 301: Angle test
C----------------------------
      CASE(302)
      POT3D = DCOS(PHI)*DSQRT(2.D0)
C----------------------------
C  model 301: Angle test
C----------------------------
      CASE(303)
      POT3D = DEXP(DCOS(TET) + DSIN(2.D0*PHI))
      END SELECT
      RETURN
      END
C=======================================================================