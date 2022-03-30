Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.

      Subroutine SOLVE1( nded,ierinv )

      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'        
      include 'parameters.h'

      character*256 message

      integer nb12,nb22,ij22,indx,ij11,nded
     .      , ierinv,nctemp,il,id,jd,jl,ii,jj,i,j,k,ioerr
        
c*** DEBUG
c      integer*4 jel ,ifr(10)

      real*8 scale(maxprm),rcond

C APRIL 1984 REWRITE SOLVE1 FOR ORIGINAL AMATRIX
C WRITTEN ON DISK
C
C      STATEMENT OF PROBLEM:
C   GIVEN A COEFFICIENT MATRIX A
C   A RIGHT HAND VECTOR B
C   NUMBER OF VARIABLES NCOEF (ntpart)
C   A VECTOR WHICH FLAGS FIXED OR FREE PARAMS (FREE)
C
C     WE WISH TO
C   A. SOLVE FOR THE FREE PARAMS FOR FIXED VALUES OF THE
C      FIXED PARAMS
C
C   B.   CALCULATE CHI2 FOR VALUES OF THE FIXED PARAMS,
C        CHI2 BEING MINIMIZED W.R.T FREE PARAMS
C        (NOT NECESSARILY IN THAT ORDER)
C
C
C     THE EQUATIONS CAN BE WRITTEN SIMPLY IN <BRA|KET> NOTATION
C
C
C      IF:
C       |X> IS THE PARAMTER VECTOR
C       |B> IS THE R.H.S VECTOR
C        A IS THE COEFFICIENT MATRIX
C       R2SUM IS THE WEIGHTED SUM OF THE SQUARED RESIDUALS
C
C
C WE USE THE SUBSCRIPT 1 TO MARK THE SUBSPACE
C      OF THE FREE PARAMTERS, 2 FOR THE FIXED
C
C       THUS WE WRITE |X>=|X1>*|X2>    WHERE * IS TENSOR PRODUCT
C LIKEWISE    |B>=|B1>*|B2>
C
C WHILE A IS BROKEN DOWN:
C
C
C     I  A11       | A12 I
C     I            |     I
C A=  I............|.....I
C     I  A21       | A22 I
C
C
C     THE ONLY MATRIX WHICH HAS TO BE INVERTED IS A11
C
C THE SOLUTIONS:
C
C LET |B1'> = |B1>-A12|X2>
C
C THEN |X1> = INVERSE(A11)|B1'> = INV(A11)|B1> - INV(A11)A12|X2>
C
C AND CHI2 = R2SUM - 2*<B2|X2> + <X2|A22|X2> - <X1|A11|X1>
C          = R2SUM - 2<B|X> + <X|A|X>
C
C          IF WE WISH TO SOLVE FOR |X1> FIRST
C
C AND CHI2 = R2SUM - 2*<B2|X2> + <X2|A22|X2> - <B1'|INV(A11):B1'>
C
C       IF WE DON'T
C
C
C    THE FUNCTION OF SOLVE1 IS TO FORM THE MATRIX A11 AND INVERT
C
C  THE FUCTION OF SOLVE2 IS TO SOLVE FOR :X1> AND CHI2
C   FOR THE PARTICULAR VALUES OF :X2> WITH WHICH
C   SOLVE2 IS CALLED
C
C
C    THE ROUTINES WORK IF ALL FREE = 1 OR ALL FREE = 0
C
C--------------------------------------------------------------
C
CD     WRITE(6,334)
CD 334 FORMAT(' IN SOLVE1')
      REWIND 27
      READ(27,iostat=ioerr) NCTEMP
      if( ioerr.ne.0 )  call report_stat('FATAL','SOLVE','solve1',' '
     .  ,'Error reading first record of temp file 27',ioerr)
 
      IF(NCTEMP.NE.ntpart) then
        WRITE(message,2) ntpart,NCTEMP
2       FORMAT(' Problem reading normal eqns: ntpart,nctemp=',2i5)
        call report_stat('FATAL','SOLVE','solve1',' ',message,0)
      endif

C-----------------------------------------------
C
      IJ11=0
      NB22=NLIVE*(NLIVE+1)/2
      NB12=NB22+NDED*(NDED+1)/2
      IJ22=NB22
C
      IL=0
      ID=0

      DO 100 I=1,ntpart
C
      READ(27,iostat=ioerr) (B(J),J=1,I) 
      if( ioerr.ne.0 )  call report_stat('FATAL','SOLVE','solve1',' '
     .  ,'Error reading second record of temp file 27',ioerr)
 
CD     WRITE(6,671) (B(J),J=1,I)
cd  671 FORMAT(6D15.8)
C
      IF(FREE(I).EQ.0) ID=ID+1
      IF(FREE(I).NE.0) IL=IL+1
C
      JD=0
      JL=0
      DO 100 J=1,I
C
      IF(FREE(J).EQ.0) JD=JD+1
      IF(FREE(J).NE.0) JL=JL+1
C
C
      IF((FREE(I).EQ.0).OR.(FREE(J).EQ.0)) GO TO 10
C A11
      IJ11=IJ11+1
      A(IJ11)=B(J)
CD     WRITE(6,*) I,J,IJ11,A(IJ11)
      GO TO 100
C
10    CONTINUE
      IF((FREE(I).NE.0).OR.(FREE(J).NE.0)) GO TO 20
C A22
      IJ22=IJ22+1
      A(IJ22)=B(J)
CD     WRITE(6,*) I,J,IJ22,A(IJ22)
      GO TO 100
C
20    CONTINUE
C A12
C  FOR A12 We store the whole matrix. It is Nlive * Nded
C  in dimension.It could be stored in a linear manner using
C
      IF(FREE(I).EQ.0) GO TO 25
C  I IS THE LIVE ONE
      II=IL
      JJ=JD
      GO TO 30
25    CONTINUE
      II=JL
      JJ=ID
30    CONTINUE
C
C**   FILL A12 PORTION - IT'S STORED AFTER A11 AND A22
      INDX=NB12+NLIVE*(JJ-1)+II
      A(INDX)=B(J)
CD     WRITE(6,*) I,J,II,JJ,INDX,A(INDX)
C**   A12(II,JJ)=B(J)
100   CONTINUE
C
C  FINALLY READ IN B (RIGHT HAND SIDE)
      READ(27,iostat=ioerr) (B(J),J=1,ntpart) 
      if( ioerr.ne.0 )  call report_stat('FATAL','SOLVE','solve1',' '
     .  ,'Error reading third record of temp file 27',ioerr)
 
C
       DO 115 I=1,ntpart
115     BORG(I)=B(I)
C-------------------------------------------------------------
C
       IL=0
       ID=0
       DO 120 I=1,ntpart
CD      WRITE(6,651) I,FREE(I)
cd  651  FORMAT(2I5)
       IF(FREE(I).EQ.0) GO TO 118
       IL=IL+1
       K=IL
       GO TO 119
118    CONTINUE
       ID=ID+1
       K=ID+NLIVE
119    CONTINUE
       B(K)=BORG(I)
CD      WRITE(6,335) K,B(K)
cd 335    FORMAT(2X,I5,2X,F25.8)
120    CONTINUE
C
C B NOW CONTAINS NLIVE B'S FOLLOWED BY NDED
C-------------------------------------------------------------


c  *** debug print the ratio of offdiagonal to diagonal
c      and scale the diagonal
c      print *,'SOLVE1 ntpart nlive nded ',ntpart,nlive,nded
c      call printa( ntpart ) 
c      j= 0
c      do i=1,ntpart 
c        j = j +1
c        if( j.gt.10 ) stop 10   
c        ifr(j) = free(i)
c        if( mod(i,10).eq.0)  then     
c          write(*,'(i6,10i2)') i-9,(ifr(jj),jj=1,10)
c          j=0 
c        endif
c      enddo     
c      print *,'SOLVE1 checking aij/aii ratio '
c      do i=1,nlive    
c        do j=1,i-1
c          if( dabs(a(jel(i,j)))/dsqrt(a(jel(i,i))*a(jel(j,j))).gt.1.d0 )
c     .      then
c             print *,'i,j, aii aij ',i,j,a(jel(i,i)),a(jel(i,j))
c          endif
c        enddo
c        a(jel(i,i)) = a(jel(i,i))*(1.d0 + 1.d-6) 
c      enddo  
c      print *,'SOLVE1 after diagonal '
c      call printa( 1440 )    


C
C NORMALIZE NORMAL EQUATIONS
CD     WRITE(6,*) NLIVE
CD     WRITE(6,*) A(1),B(1)
      IF(NLIVE.GT.0) CALL NRMSCL(SCALE,NLIVE,0,1)
CD     WRITE(6,*) A(1),B(1),SCALE(1)
      if( logprt ) WRITE(6,888)
  888 FORMAT(' Solving Normal Equations')
C      IF(NLIVE.GT.0)
C     1  CALL VINV2(3,NLIVE,IERINV)
      IF(NLIVE.GT.0) CALL INVER2(A,B,3,NLIVE,rcond,IERINV)

CD     WRITE(6,*) A(1),B(1),SCALE(1)
      IF(NLIVE.GT.0) CALL NRMSCL(SCALE,NLIVE,1,1)
CD     WRITE(6,*) A(1),B(1),SCALE(1)
      IF(IERINV.NE.0)
     . call report_stat('WARNING','SOLVE','solve1',' '
     .      , 'Bad inversion--too many parameters adjusted',0)
C
C      NOTE:  VECTOR B1  NOW ACTUALLY CONTAINS INV(A11):B1>
C
      RETURN
C
      END
