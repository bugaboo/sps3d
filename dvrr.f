      SUBROUTINE DVRR(KPOL,LAN,NDVR,X,W,T,DVR,PIR)
C=======================================================================
C  DVR arrays for the SPS EVP in the R problem 
C  [see Appendix C in PRA 58, 2077 (1998)]
C  KPOL.EQ.0 - Legendre polynomials (LAN is not used in this case)
C  KPOL.NE.0 - Jacobi(0,2*LAN) polynomials
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(*),W(*),T(NDVR+1,*),DVR(*),PIR(*)
      ALLOCATABLE VK(:),WK(:,:)
C
C  DVR arrays X, W, T, and DVR and right-end polynomial amplitudes VK
C
      ALLOCATE(VK(NDVR),WK(NDVR+1,NDVR+1))
      IF(KPOL.EQ.0) THEN
        CALL LEGPOM(NDVR)
        CALL DVRLEG(NDVR,X,W,T,0,0,ERR,NDVR+1)
        DO n=1,NDVR
          VK(n)=DSQRT(DBLE(n)-0.5D0)
        ENDDO
      ELSE
        BET=2.D0*DBLE(LAN)
        CALL JACPOM(0.D0,BET,NDVR)
        CALL DVRJAC(0.D0,BET,NDVR,X,W,T,[0],ERR,NDVR+1)
        CALL JACPOL(1.D0,VK,NDVR-1)
        TMP=2.D0**LAN
        DO n=1,NDVR
          VK(n)=TMP*VK(n)
        ENDDO
      ENDIF
C
C  Right-end DVR amplitudes PIR
C
      DO i=1,NDVR
        SUM=0.D0
        DO n=1,NDVR
          SUM=SUM+T(n,i)*VK(n)
        ENDDO
        PIR(i)=SUM
      ENDDO
C
C  Kinetic matrix DVR
C
      SUM=0.D0
      DO n=1,NDVR
        FF=VK(n)**2
        nm=(n*(n+1))/2
        DVR(nm)=2.D0*FF*SUM+0.5D0*(FF-0.5D0)**2
        TMP=VK(n)*(2.D0*SUM+FF-0.5D0)
        DO m=n+1,NDVR
          nm=nm+m-1
          DVR(nm)=VK(m)*TMP
        ENDDO
        SUM=SUM+FF
      ENDDO
      CALL MUL_LRS(NDVR,DVR,T,WK,NDVR+1)
      DEALLOCATE(VK,WK)
C
      RETURN
      END
C=======================================================================