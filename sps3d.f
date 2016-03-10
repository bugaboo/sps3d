      SUBROUTINE SPS3D(LMAX,NANG,L,M,RADA,NDVR,NTET,NPHI,PIR,
     &                 CK,CE,CVEC,CPHI,NB,NA,NOI)
C=======================================================================
C  SPS basis for the R problem with L.GE.0 [PRA 75, 062704 (2007)]
C  'product' formula for surface amplitudes
C-----------------------------------------------------------------------
      USE ANGSYM
      IMPLICIT REAL*8 (A,B,D-H,O-Z)
      IMPLICIT COMPLEX*16 (C)
      PARAMETER(KTEST=0)
      ALLOCATABLE CZER(:,:)
      ALLOCATABLE X(:),W(:),T(:,:),DVR(:)
      DIMENSION CK(*),CE(*),CVEC(NDVR*NANG,*),CPHI(NANG,*),PIR(*)
      DIMENSION L(*),M(*)
      ALLOCATABLE HAM(:,:),VKR(:),VKI(:),VK(:),WK(:,:),CWK(:,:),RAD(:)
      ALLOCATABLE V(:,:,:), VT(:,:)
      
C
C  Allocation
C      
      
      NDVRS=(NDVR*(NDVR+1))/2
      ALLOCATE(X(NDVR),W(NDVR),T(NDVR+1,NDVR+1),DVR(NDVRS),RAD(NDVR))
      
      NSPS = 0
      DO i= 1, NANG
        NSPS = NSPS + 2 * NDVR + L(i)
      ENDDO
      
C
C  Setting up DVR arrays
C
      SR=0.5D0*RADA
      NDVR2=2*NDVR
      CALL DVRR(0,0,NDVR,X,W,T,DVR,PIR)
      ij=0
      DO j=1,NDVR
        RAD(j)=SR*(1.D0+X(j))
        PIR(j)=PIR(j)/RAD(j)
        DO i=1,j
          ij=ij+1
          DVR(ij)=DVR(ij)/(RAD(i)*RAD(j))+PIR(i)*PIR(j)
        ENDDO
      ENDDO
C
C  Zeros of the reverse Bessel polynomial
C
      ALLOCATE (CZER(LMAX,LMAX))
      CALL BESZERR(CZER,LMAX)
C
C  Calculating potential
C
      ALLOCATE(VT(NANG,NANG),V(NANG,NANG,NDVR))
      DO k=1,NDVR
      CALL POTMAT(NTET,NPHI,RAD(k),VT,NANG,L,M,LMAX)
      DO i=1,NANG
        DO j=1,NANG
        V(i,j,k) = VT(i,j)
        ENDDO 
      ENDDO
      ENDDO
C
C  Setting up SPS Hamiltonian
C
      ALLOCATE(HAM(NSPS,NSPS))
      DO j=1,NSPS
        DO i=1,NSPS
          HAM(i,j)=0.D0
        ENDDO
      ENDDO
      nnu = 0
      DO nu = 1, NANG
        DO i=1,NDVR
        HAM(nnu+i, nnu+NDVR+i)=1.D0
      ENDDO

      ij=0
      DO j=1,NDVR
        DO i=1,j
          ij=ij+1
          TMP=-2.D0*DVR(ij)
          HAM(nnu+NDVR+i, nnu+j)=TMP
          HAM(nnu+NDVR+j, nnu+i)=TMP
          TMP=2.D0*RADA*PIR(i)*PIR(j)
          HAM(nnu+NDVR+i, nnu+NDVR+j)=TMP
          HAM(nnu+NDVR+j, nnu+NDVR+i)=TMP
        ENDDO
        HAM(nnu+NDVR+j,nnu+j)=HAM(nnu+NDVR+j,nnu+j)-
     &            DBLE(L(nu)*(L(nu)+1))/RAD(j)/RAD(j)
      ENDDO
      TMP=DSQRT(2.D0/RADA)
      DO ip = 1, L(nu)/2
        ip1=NDVR2+2*ip-1
        ip2=ip1+1
        ZER1=DREAL(CZER(L(nu),ip))
        ZER2=DIMAG(CZER(L(nu),ip))
        DO i=1,NDVR
          TMPi=TMP*PIR(i)
          HAM(nnu+NDVR+i,nnu+ip1)=2.D0*TMPi
          HAM(nnu+ip1,nnu+i)=-ZER1*TMPi
          HAM(nnu+ip2,nnu+i)=-ZER2*TMPi
        ENDDO
        HAM(nnu+ip1,nnu+ip1)=-ZER1/RADA
        HAM(nnu+ip1,nnu+ip2)=+ZER2/RADA
        HAM(nnu+ip2,nnu+ip1)=-ZER2/RADA
        HAM(nnu+ip2,nnu+ip2)=-ZER1/RADA
      ENDDO
      IF(MOD(L(nu),2).EQ.1) THEN
        ZER=CZER(L(nu),L(nu)/2+1)
        DO i=1,NDVR
          TMPi=TMP*PIR(i)
          HAM(nnu+NDVR+i,nnu+2*NDVR+L(nu))=TMPi
          HAM(nnu+2*NDVR+L(nu),nnu+i)=-ZER*TMPi
        ENDDO
        HAM(nnu+2*NDVR+L(nu),nnu+2*NDVR+L(nu))=-ZER/RADA
      ENDIF
      
      nmu = 0
      DO mu = 1, NANG
        DO i = 1,NDVR
          HAM(nnu+NDVR+i,nmu+i)=HAM(nnu+NDVR+i,nmu+i)-2.D0*V(nu,mu,i)
        ENDDO
      nmu = nmu + 2*NDVR + L(mu)
      ENDDO
      nnu = nnu + 2*NDVR + L(nu)
      ENDDO
C
C  Solution of the SPS EVP
C
C --- IMSL
c     ALLOCATE(CWK(NSPS,NSPS))
c     CALL DEVCRG(NSPS,HAM,NSPS,CK,CWK,NSPS)
c      DO n=1,NSPS
c        CK(n)=-(0.D0,1.D0)*CK(n)
c        DO i=1,NDVR
c          CVEC(i,n)=CWK(i,n)
c        ENDDO
c      ENDDO
c      DEALLOCATE(CWK)
C --- LAPACK
      NVK=10*NSPS
      ALLOCATE(VKR(NSPS),VKI(NSPS),VK(NVK),WK(NSPS,NSPS))
      PRINT *, "Call DGEEV"
      CALL DGEEV('N','V',NSPS,HAM,NSPS,VKR,VKI,ERR,1,WK,NSPS,
     &           VK,NVK,info)
      PRINT *, "DGEEV done"
      DO n=1,NSPS
        CK(n)=-(0.D0,1.D0)*DCMPLX(VKR(n),VKI(n))
        IF(VKI(n).EQ.0.D0) THEN
        nnu = 0
        DO nu = 0, NANG-1
          DO i=1,NDVR
            CVEC(nu*NDVR+i,n)=DCMPLX(WK(nnu+i,n),0.D0)
          ENDDO
          nnu = nnu + 2*NDVR + L(nu+1)
        ENDDO
        ELSE IF(VKI(n).GT.0.D0) THEN
        nnu = 0
        DO nu = 0, NANG-1
          DO i=1,NDVR
            CVEC(nu*NDVR+i,n)  =DCMPLX(WK(nnu+i,n),+WK(nnu+i,n+1))
            CVEC(nu*NDVR+i,n+1)=DCMPLX(WK(nnu+i,n),-WK(nnu+i,n+1))
          ENDDO
          nnu = nnu + 2*NDVR + L(nu+1)
        ENDDO
        ENDIF
      ENDDO
      DEALLOCATE(VKR,VKI,VK,WK)
      DEALLOCATE(HAM)

C
C  Ordering the SPS
C
       CALL SPSODR(NDVR*NANG,NSPS,CK,CE,CVEC,NB,NA,NOI)
      
C
C  Normalization
C
      DO n=1,NSPS
        CTK = (0.D0,1.D0)*RADA/CK(n)
        CSUMC = 0.D0
        DO nu=1,NANG
          CTMP=1.D0
          DO ip=1,L(nu)
            CTMP=CTMP+CZER(L(nu),ip)/((0.D0,1.D0)*CK(n)*RADA+
     &            CZER(L(nu),ip))**2
          ENDDO        
          CSUMF = 0.D0
          DO i=1,NDVR
            DO j=1,NDVR
              CSUMF=CSUMF+CVEC((nu-1)*NDVR+i,n)*PIR(i)*PIR(j)
     &          * CVEC((nu-1)*NDVR+j,n)
            ENDDO
            CSUMC = CSUMC + CVEC((nu-1)*NDVR+i,n)**2
          ENDDO
          CSUMC = CSUMC + CTK*CTMP*CSUMF
        ENDDO

        CTMP = CDSQRT(1/CSUMC/SR)
        DO i=1,NDVR*NANG
          CVEC(i,n)=CVEC(i,n)*CTMP
        ENDDO
        
        DO nu=1,NANG
        CPHI(nu,n)=0.D0
        DO i=1,NDVR
          CPHI(nu,n)=CPHI(nu,n)+CVEC((nu-1)*NDVR+i,n)*PIR(i)
        ENDDO
        CPHI(nu,n)=CPHI(nu,n)*RADA
        ENDDO
      ENDDO
      
      IF (KTEST.EQ.1) THEN
C
C  Checking orthogonality
C
      DO n=1,10
        DO k=n,20
          CTK = (0.D0,2.D0)*RADA/(CK(n)+CK(k))
          CSUMC = 0.D0
          DO nu=1,NANG
            CTMP=1.D0
            DO ip=1,L(nu)
              CTMP=CTMP+CZER(L(nu),ip)/((0.D0,1.D0)*CK(n)*RADA+
     &          CZER(L(nu),ip))/((0.D0,1.D0)*CK(k)*RADA+CZER(L(nu),ip))
            ENDDO         
            CSUMF = 0.D0
            DO i=1,NDVR
              DO j=1,NDVR
                CSUMF=CSUMF+CVEC((nu-1)*NDVR+i,n)*PIR(i)*PIR(j)
     &                  *CVEC((nu-1)*NDVR+j,k)
              ENDDO
              CSUMC = CSUMC+CVEC((nu-1)*NDVR+i,n)*CVEC((nu-1)*NDVR+i,k)
            ENDDO
            CSUMC = CSUMC + CTK*CTMP*CSUMF
          ENDDO
          TMP = CDABS(CSUMC)
          PRINT *, CSUMC, n, k
        ENDDO
      ENDDO
      ENDIF
      
C
 77   FORMAT(4(E19.12,1X))
      DEALLOCATE(X,W,T,DVR,RAD,V,VT)
      RETURN
      END
C=======================================================================