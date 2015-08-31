      SUBROUTINE SPSODR(NDVR,NSPS,CK,CE,CVEC,NB,NA,NOI)
C=======================================================================
C  Two-level ordering of SPS (for SPSR0, SPSRL, and SPSX)
C  (1) in three groups: bound, antibound, and outgoing/incoming pairs
C  (2) within each group: in increasing order of |k|
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A,B,D-H,O-Z)
      IMPLICIT COMPLEX*16 (C)
      DIMENSION CK(*),CE(*),CVEC(NDVR,*)
      ALLOCATABLE IND(:),VK(:),VKR(:),VKI(:),CVK(:)
C
C  Ordering in |k|
C
      ALLOCATE(IND(NSPS),VK(NSPS),VKR(NSPS),VKI(NSPS),CVK(NSPS))
      DO n=1,NSPS
        VK(n)=CDABS(CK(n))
      ENDDO
      CALL INDEXX(NSPS,VK,IND)
      DO n=1,NSPS
        CVK(n)=CK(n)
      ENDDO
      DO n=1,NSPS
        CK(n)=CVK(IND(n))
      ENDDO
      DO i=1,NDVR
        DO n=1,NSPS
          CVK(n)=CVEC(i,n)
        ENDDO
        DO n=1,NSPS
          CVEC(i,n)=CVK(IND(n))
        ENDDO
      ENDDO
C
C  Ordering in groups
C
      DO n=1,NSPS
        VKR(n)=DREAL(CK(n))
        VKI(n)=DIMAG(CK(n))
      ENDDO
C --- bound states
      NB=0
      DO n=1,NSPS
        IF(VKR(n).EQ.0.D0.AND.VKI(n).GT.0.D0) THEN
          NB=NB+1
          IND(NB)=n
        ENDIF
      ENDDO
C --- antibound states
      NA=0
      DO n=1,NSPS
        IF(VKR(n).EQ.0.D0.AND.VKI(n).LE.0.D0) THEN
          NA=NA+1
          IND(NB+NA)=n
        ENDIF
      ENDDO
      NBA=NB+NA
C --- outgoing states
      NO=0
      DO n=1,NSPS
        IF(VKR(n).GT.0.D0.AND.VKI(n).LE.0.D0) THEN
          NO=NO+1
          IND(NBA+2*NO-1)=n
        ENDIF
      ENDDO
C --- incoming states
      NI=0
      DO n=1,NSPS
        IF(VKR(n).LT.0.D0.AND.VKI(n).LE.0.D0) THEN
          NI=NI+1
          IND(NBA+2*NI)=n
        ENDIF
      ENDDO
      PRINT *, "Bound", NO, NI, NBA, NSPS, NA, NB
C --- test
C      IF(NO.NE.NI.OR.NBA+NO+NI.NE.NSPS) STOP ' *** SPSODR: ERROR'
      NOI=NO
      WRITE(*,70) NB,NA,NOI
C --- ordering
      DO n=1,NSPS
        CVK(n)=CK(n)
      ENDDO
      DO n=1,NSPS
        CK(n)=CVK(IND(n))
        CE(n)=0.5D0*CK(n)*CK(n)
      ENDDO
      DO i=1,NDVR
        DO n=1,NSPS
          CVK(n)=CVEC(i,n)
        ENDDO
        DO n=1,NSPS
          CVEC(i,n)=CVK(IND(n))
        ENDDO
      ENDDO
      DEALLOCATE(IND,VK,VKR,VKI,CVK)
C
 70   FORMAT(' SPSODR: NB=',I3,'  NA=',I3,'  NO=NI=',I4)
      RETURN
      END
C=======================================================================