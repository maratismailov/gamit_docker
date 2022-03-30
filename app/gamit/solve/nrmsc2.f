      SUBROUTINE NRMSC2(A,B,SCALE,NCOEF,IFLAG,IFLAG2)
C
      character*256 message
      
      integer*4 ncoef,iflag,iflag2,idiag,ij,i,j

      REAL*8 A,B,SCALE(NCOEF),DIGT
      DIMENSION A((NCOEF*(NCOEF+1))/2),B(NCOEF)
C
C.....
      GO TO (10,20), IFLAG+1
C
   10 CONTINUE
CD      WRITE(6,750)
CD   750 FORMAT(' NORMALIZATION IN PROGRESS')
C
C NORMALIZE NORMAL EQUATIONS
C.....COMPUTE SCALE FACTORS FOR THE NORMAL EQUATIONS
      DO 100 I=1,NCOEF
      IDIAG=I*(I+1)/2
      DIGT=A(IDIAG)
      IF (DIGT.GT.0.0D0) GO TO 50
      WRITE (message,'(a,i4,a,f15.4,a)') 'Diagonal term ',i
     . ,'of normal equations not positive (=',digt,'); set = 1.0'
      call report_stat('WARNING','SOLVE','nrmsc2',' ',message,0)
      SCALE(I)=1.D0
* MOD TAH 940415 Set diagonal to 1 as well to stop inversion problem
      a(idiag) = 1.d0
      GO TO 100
   50 SCALE(I)=1.D0/DSQRT(DIGT)
  100 CONTINUE
C.....SCALE THE COEFFICIENT MATRIX AND RIGHT SIDES OF THE
C.....NORMAL EQUATIONS
   20 CONTINUE
CD      WRITE(6,*) IFLAG
      IJ=0
      DO 600 I=1,NCOEF
      DO 500 J=1,I
      IJ=IJ+1
      A(IJ)=A(IJ)*SCALE(I)*SCALE(J)
CD      WRITE (6,1000) SCALE(I),A(IJ)
cd 1000 FORMAT (1X,2(E22.15,1X))
  500 CONTINUE
      IF(IFLAG2.EQ.1) B(I)=B(I)*SCALE(I)
Cd      WRITE (6,1000) SCALE(I),B(I)
  600 CONTINUE
C
      RETURN
      END
