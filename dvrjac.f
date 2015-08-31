      SUBROUTINE DVRJAC(A,B,NDVR,X,W,T,KHAM,H,NDIM)
C=======================================================================
C
C  JACOBI DISCRETE VARIABLE REPRESENTATION
C
C  A,B       - PARAMETERS OF THE JACOBI POLYNOMIALS
C  NDVR      - NUMBER OF DVR QUADRATURE POINTS 
C  X(i),W(i) - GAUSS-JACOBI ABSCISSAS AND WEIGHTS
C  T(n,i)    - FBR-DVR TRANSFORMATION MATRIX; n=POLYNOMIAL, i=POINT
C  KHAM      - 0 (1) - DON'T (DO) CALCULATE H(ij)
C  H(ij)     - KINETIC PART OF THE DVR HAMILTONIAN IN SSM
C  NDIM      - DIMENSION PARAMETER, NDIM.GE.NDVR+1
C
C-----------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (EPS=1.D-14,ITRMAX=30)
      DIMENSION X(*),W(*),T(NDIM,*),H(*)
C
C  Calculation of X(i), W(i), and T(n,i)
C
      AN=DBLE(NDVR)
      AB=A+B
      ANAB=2.D0*AN+AB
      P1=(A-B)*AN/ANAB
      P2=-AN
      TMP=DSQRT(ANAB+1.D0)
      TMP1=2.D0*DSQRT(AN*(AN+A)*(AN+B)*(AN+AB)/(ANAB-1.D0))/ANAB
      P3=TMP*TMP1
      WW=TMP/TMP1
      DO 1 i=NDVR,1,-1
        IF(i.EQ.NDVR) THEN
          AAN=A/AN
          BAN=B/AN
          R1=(1.D0+A)*(2.78D0/(4.D0+AN*AN)+0.768D0*AAN/AN)
          R2=1.D0+1.48D0*AAN+0.96D0*BAN+0.452D0*AAN*AAN+0.83D0*AAN*BAN
          Z=1.D0-R1/R2
        ELSE IF(i.EQ.NDVR-1) THEN
          R1=(4.1D0+A)/((1.D0+A)*(1.D0+0.156D0*A))
          R2=1.D0+0.06D0*(AN-8.D0)*(1.D0+0.12D0*A)/AN
          R3=1.D0+0.012D0*B*(1.D0+0.25D0*DABS(A))/AN
          Z=Z-(1.D0-Z)*R1*R2*R3
        ELSE IF(i.EQ.NDVR-2) THEN
          R1=(1.67D0+0.28D0*A)/(1.D0+0.37D0*A)
          R2=1.D0+0.22D0*(AN-8.D0)/AN
          R3=1.D0+8.D0*B/((6.28D0+B)*AN*AN)
          Z=Z-(X(NDVR)-Z)*R1*R2*R3
        ELSE IF(i.EQ.2) THEN
          R1=(1.D0+0.235D0*B)/(0.766D0+0.119D0*B)
          R2=1.D0/(1.D0+0.639D0*(AN-4.D0)/(1.D0+0.71D0*(AN-4.D0)))
          R3=1.D0/(1.D0+20.D0*A/((7.5D0+A)*AN*AN))
          Z=Z+(Z-X(4))*R1*R2*R3
        ELSE IF(i.EQ.1) THEN
          R1=(1.D0+0.37D0*B)/(1.67D0+0.28D0*B)
          R2=1.D0/(1.D0+0.22D0*(AN-8.D0)/AN)
          R3=1.D0/(1.D0+8.D0*A/((6.28D0+A)*AN*AN))
          Z=Z+(Z-X(3))*R1*R2*R3
        ELSE
          Z=3.D0*(X(i+1)-X(i+2))+X(i+3)
        ENDIF
        itr=0
        DZ=0.D0
 99     itr=itr+1
        IF(itr.GT.ITRMAX) STOP ' *** DVRJAC ITERATIONS ERROR'
        Z1=Z
        Z=Z1-DZ
        CALL JACPOL(Z,T(1,i),NDVR)
        PP=((P1+P2*Z)*T(NDVR+1,i)+P3*T(NDVR,i))/(1.D0-Z*Z)
        DZ=T(NDVR+1,i)/PP
        IF(DABS(DZ).GT.EPS) GOTO 99
        X(i)=Z
        W(i)=WW/(T(NDVR,i)*PP)
        TMP=DSQRT(W(i))
        DO n=1,NDVR
          T(n,i)=TMP*T(n,i)
        ENDDO
 1    CONTINUE
      IF(KHAM.EQ.0) RETURN
C
C  Calculation of H(ij)
C
      E0=0.25D0*AB*(AB+2.D0)
      ij=0
      DO 2 j=1,NDVR
        DO 3 i=1,j
          ij=ij+1
          H(ij)=0.D0
          DO n=1,NDVR
            AN=DBLE(n)
            En=E0+(AN-1.D0)*(AN+AB)
            H(ij)=H(ij)+T(n,i)*En*T(n,j)
          ENDDO
 3      CONTINUE
 2    CONTINUE
C
      RETURN 
      END
C=======================================================================