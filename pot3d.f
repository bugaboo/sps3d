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
      X=DSIN(TET)*DCOS(PHI)
      Y=DSIN(TET)*DSIN(PHI)
      Z=DCOS(TET)
      TX=0.5D0
      TY=2.D0
      R=RAD*DSQRT(X**2/TX+Y**2/TY+Z**2)
      POT3D=-DEXP(-0.1D0*R)/R
C----------------------------
C  model 3: Yukawa potential
C----------------------------
      CASE(3)
      X=8.D0*RAD*DSIN(TET)*DCOS(PHI)
      Y=9.5D0*RAD*DSIN(TET)*DSIN(PHI)
      Z=11.5D0*RAD*DCOS(TET)
      POT3D=0.5D0(X**2+Y**2+Z**2-2.29D2)
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