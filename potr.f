      FUNCTION POTR(RAD)
C=======================================================================
C  Potential energy V(r) for one-channel SPS EVP in the R problem
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /POTR_C/MODELR
C
      SELECT CASE(MODELR)
C--------------------------------
C  model 1: rectangular potential
C--------------------------------
      CASE(1)
      POTR=-10.D0
C--------------------------------
C  model 2: exponential potential
C--------------------------------
      CASE(2)
      POTR=-DEXP(-RAD)
C-----------------------------
C  model 3: Gaussian potential
C-----------------------------
      CASE(3)
      POTR=-1.D0*DEXP(-RAD*RAD)
C---------------------------
C  model 4: Eckart potential
C---------------------------
      CASE(4)
      POTR=-21.D0/DCOSH(RAD)**2
C--------------------------
C  model 5: Morse potential
C--------------------------
      CASE(5)
      TMP=DEXP(1.5D0-RAD)
      POTR=12.5D0*TMP*(TMP-2.D0)
C-----------------------------------
C  model 6: potential with a barrier
C-----------------------------------
      CASE(6)
      POTR=7.5D0*RAD*RAD*DEXP(-RAD)
C---------------------------------------
C  model 7: Model 1 potential for H^{-}
C  this:  E0=-0.027750994 and SL=5.96505
C  exact: E0=-0.0277510   and SL=5.965
C---------------------------------------
      CASE(7)
      R=RAD/2.5026D0
      POTR=-0.3831087D0*DEXP(-R*R)
C---------------------------------------
C  model 8: Model 2 potential for H^{-}
C  this:  E0=-0.027751007 and SL=5.96500
C  exact: E0=-0.0277510   and SL=5.965
C---------------------------------------
      CASE(8)
      TMP=DSQRT(RAD*RAD+2.76D0**2)
      POTR=-19.74974D0*DEXP(-TMP)/TMP
C---------------------------------------
C  model 9: Model 3 potential for H^{-}
C  this:  E0=-0.02775103  and SL=4.957
C  this:  E0=-0.02775100  and SL=5.305
C  this:  E0=-0.02775101  and SL=5.642
C  this:  E0=-0.02775101  and SL=5.963
C  exact: E0=-0.0277510   and SL=5.965
C---------------------------------------
      CASE(9)
c      R=RAD/1.0D0
c      POTR=-1.723380D0*DEXP(-R*R)
c      R=RAD/1.5D0
c      POTR=-0.8598014D0*DEXP(-R*R)
      R=RAD/2.0D0
      POTR=-0.5398646D0*DEXP(-R*R)
c      R=RAD/2.5D0
c      POTR=-0.3837017D0*DEXP(-R*R)
C-----------------------------
C  model 10: Coulomb potential
C-----------------------------
      CASE(10)
      POTR=-1.D0/RAD
C----------------------------
C  model 11: Yukawa potential
C----------------------------
      CASE(11)
      POTR=-DEXP(-0.01D0*RAD)/RAD
C---------------------------------------
C  model 12: Gauss-Coulomb potential
C  this: E0=-0.4854833635, SL=32.2704512
C---------------------------------------
      CASE(12)
      POTR=-DEXP(-(RAD/10.D0)**2)/RAD
C----------------------------------
C  model 13: Model potential for He
C  this:  E(1s)=-0.90372448
C----------------------------------
      CASE(13)
      POTR=-(1.D0+DEXP(-2.132405D0*RAD))/RAD
C----------------------------------
C  model 14: Model potential for Li
C  this:  E(2s)=-0.198146915844
C----------------------------------
      CASE(14)
      BET=1.655705D0
      POTR=-(1.D0+2.D0*(1.D0+BET*RAD)*DEXP(-2.D0*BET*RAD))/RAD
C-----------------------------------
C  model 15: Model potential for Ne*
C-----------------------------------
      CASE(15)
      ZEFF=1.D0+9.D0*(
     &    (1.29942636D0*DEXP(-3.06518063D0*RAD)-
     &     1.12569309D0*DEXP(-4.51696279D0*RAD))*RAD*RAD+
     &    (0.47271993D0*DEXP(-2.90207510D0*RAD)+
     &     0.42068967D0*DEXP(-5.80903874D0*RAD))*RAD+
     &    (0.53912121D0*DEXP(-2.41322960D0*RAD)+
     &     0.46087879D0*DEXP(-4.68014471D0*RAD)))
      POTR=-ZEFF/RAD
C---------------------
C  model 0: free model
C---------------------
      CASE(0) 
      P1=4.418d0
      P2=1.351d0
      Z0=36.d0
      TMP=1.D0/((P1/P2)*(DEXP(P2*RAD)-1.D0)+1.D0)
      ZEFF=Z0-(Z0-1.D0)*(1.D0-TMP)
      POTR=-ZEFF/RAD
C
      END SELECT
      RETURN
      END
C=======================================================================