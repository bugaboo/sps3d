      SUBROUTINE CDETS(CS,CW,CDET,N)
*  CS is S-matrix
*  CW are eigenphase shifts
*  CDET is Det[S]
      IMPLICIT REAL*8 (A,B,D-H,O-Z)
      IMPLICIT COMPLEX*16 (C)
      DIMENSION CW(*),CS(N,*)
      ALLOCATABLE :: CWK(:), WK(:), CWL(:,:),CWR(:,:)
      
      ALLOCATE (WK(2*N),CWK(10*N),CWL(N,N),CWR(N,N))
      CALL ZGEEV('N','N',N,CS,N,CW,CWL,N,CWR,N,CWK,10*N,WK,INFO)
      
      DEALLOCATE (WK,CWK,CWL,CWR)
      CDET = 1.D0
      DO i=1,N
	CDET = CDET*CW(i)
	CW(i)=(0.D0,-0.5D0)*CDLOG(CW(i))
      ENDDO
      END SUBROUTINE