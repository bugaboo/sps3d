      SUBROUTINE POTMAT(NTET,NPHI,RAD,V,NANG,L,M,LMAX)
      USE OMP_LIB
      IMPLICIT REAL*8 (A,B,D-H,O-Z)
      IMPLICIT COMPLEX*16 (C)
      INTEGER, DIMENSION(NANG) :: L,M
      DIMENSION V(NANG,NANG)
      ALLOCATABLE X(:,:), W(:,:), T(:), XT(:), WT(:)
c
      ALLOCATE(X(NTET,0:2*LMAX),W(NTET,0:2*LMAX), T(2*NTET))
      ALLOCATE(XT(NTET), WT(NTET))
      PI=2*DACOS(0.D0)
      DO i=0,2*LMAX
	CALL JACPOM(DBLE(i)/2.D0,DBLE(i)/2.D0,NTET)
	CALL GAUJAC(DBLE(i)/2.D0,DBLE(i)/2.D0,NTET,XT,WT,T)
	DO j=1,NTET
	  X(j,i) = XT(j)
	  W(j,i) = WT(j)
	ENDDO
      ENDDO
      
!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(NANG,NPHI,NTET,L,M,V,X,W,PI,RAD)
!$OMP DO
      DO i = 1,NANG
        DO j = 1,NANG
	  SUMP = 0.D0
	  SUMP1 = 0.D0
	  DO k = 1,NPHI
	    PHI = DBLE(2*k-1)*PI/DBLE(NPHI)
	    M1 = M(i)
	    IF (M(i).EQ.0) THEN
	      TMP = 1.D0
	    ELSE 
	      IF (M(i).GT.0) THEN
		TMP = DSQRT(2.D0)*DCOS(DBLE(M(i))*PHI)
	      ELSE
		M1 = -M1
		TMP = DSQRT(2.D0)*DSIN(-DBLE(M(i))*PHI)
	      ENDIF
	    ENDIF
	    M2 = M(j)
	    IF (M(j).EQ.0) THEN
	    ELSE 
	      IF (M(j).GT.0) THEN
		TMP = TMP*DSQRT(2.D0)*DCOS(DBLE(M(j))*PHI)
	      ELSE 
		M2 = -M2
		TMP = TMP*DSQRT(2.D0)*DSIN(-DBLE(M(j))*PHI)
	      ENDIF
	    ENDIF
	    SUMT = 0.D0
	    DO n = 1,NTET
      SUMT=SUMT+W(n,M1+M2)*POLJ(DBLE(M1),DBLE(M1),L(i)-M1,X(n,M1+M2))
     & 		*POLJ(DBLE(M2),DBLE(M2),L(j)-M2,X(n,M1+M2)) *
     & 		POT3D(RAD,DACOS(X(n,M1+M2)),PHI)
	    ENDDO
	    SUMP = SUMP + TMP*SUMT
	  ENDDO
	  V(i,j) = SUMP/DBLE(NPHI)
	ENDDO
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
      DEALLOCATE(X,W,XT,WT,T)
      RETURN
      END