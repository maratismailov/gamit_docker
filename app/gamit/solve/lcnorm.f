Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
C
      Subroutine LCNORM( ierinv,ijob )
C
C       IJOB = 1 : Copy LC mode normal matrix
C       IJOB = 2 : Solve LC mode normal equation ( = SOLVE1)

      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'            
      include 'parameters.h'
             
      integer ierinv,ijob,lrow,lrow0,lrow1,lb,il,id
     .      , i4,ij4,ij,i,j,k,ioerr

c     This added for bias-paramter debug
c      integer*4 jel

      real*8 scale(maxprm),r2lc,rcond

      character*256 message

      GOTO (10,200), IJOB
c             
c      print *,'LCNORM 1 free(211-270 '
c      write(*,'(10i7)') (free(i),i=211,270)
 10   CONTINUE
      REWIND 28
      R2SUM = 0.0D0
C---- Copy N11 to a  
      ij4 = lpart*(lpart+1)/2 
c      print *,'LCNORM lpart ij4 alc ',lpart,ij4,(alc(i),i=1,ij4) 
      call copy1d(1,ij4,0,alc,a)
      LROW = LPART
      LROW0 = LROW
      LROW1 = LROW
C---- Read N21 and N22 from file 28.
C---- Read r2sum
         READ (28,iostat=ioerr) R2LC 
        if( ioerr.ne.0 )  call report_stat('FATAL','SOLVE','lcnorm',' '
     .    ,'Error reading first record of temp file 28',ioerr)  
c        print *,'LCNORM read 28 r2lc r2sum ',r2lc,r2sum
         R2SUM = R2SUM + R2LC
         READ (28,iostat=ioerr) LB  
         if( ioerr.ne.0 )  call report_stat('FATAL','SOLVE','lcnorm',' '
     .     ,'Error reading second record of temp file 28',ioerr)
c         print *,'r2lc r2sum lb ',r2lc,r2sum,lb
C---- Read N21
         DO 110 I = 1,LB
            LROW = LROW + 1
            IJ4 = (LROW-1)*LROW/2
            READ (28,iostat=ioerr) (a(IJ4 + I4),I4 = 1,LPART)
            if( ioerr.ne.0 )  call report_stat('FATAL','SOLVE'
     .        ,'lcnorm',' ','Error reading third record of temp file 28'
     .        ,ioerr)
c           print *,'ij4 a ',ij4,(a(ij4+i4),i4=1,lpart)
 110     CONTINUE         
c      print *,'LCNORM after reading 28'
c      call printa(lpart)
C---- Read N22
         DO 120 I = 1,LB
            LROW0 = LROW0 + 1
            IJ4 = (LROW0 - 1)*LROW0/2 + LROW1
            READ (28,iostat=ioerr) (a(IJ4 + I4),I4 = 1,I)   
            if( ioerr.ne.0 )  call report_stat('FATAL','SOLVE','lcnorm'
     .        ,' ','Error reading fourth record of temp file 28',ioerr)
c            print *,'lrow0 a ',(a(IJ4 + I4),I4 = 1,I)
 120     CONTINUE
         LROW = LROW + LB
         LROW0 = LROW
         LROW1 = LROW0
c       call printan22(l1bias(1))
C---- Copy BLC to b
      call copy1d(1,ntpart,0,blc,b)
      GOTO 400
c
C-------------------------------------------------------------
 200  CONTINUE
      call copy1d(1,ntpart,0,b,borg)
c      print *,'LCNORM 2 free(211-270 '
c      write(*,'(10i7)') (free(i),i=211,270)
       IL = 0
       ID = 0
       DO 220 I = 1,NTPART
          IF(FREE(I).EQ.0) GO TO 218
          IL = IL + 1
          K = IL
          I4 = K*(K - 1)/2
          IJ4 = I*(I - 1)/2
          DO 240 J = 1,I
             IF (FREE(J).EQ.0) GOTO 240
             I4 = I4 + 1
             IJ = IJ4 + J
             a(I4) = a(IJ)
 240      CONTINUE
          GO TO 219
 218      CONTINUE
          ID = ID + 1
          K = ID + NLIVE
 219      CONTINUE
          b(K) = BORG(I)
 220   CONTINUE
C
C B NOW CONTAINS NLIVE B'S FOLLOWED BY NDED
C-------------------------------------------------------------
      IF (NLIVE.EQ.0) GOTO 400   
                
c       print *,'LCNORM ntpart nlive ',ntpart,nlive
c       print *,'LCNORM before  reweight '
c       call printa( nlive )    
c      print *,'LCNORM checking aij/aii ratio '
c      do i=1,nlive    
c        do j=1,i-1
c          if( dabs(a(jel(i,j)))/dsqrt(a(jel(i,i))*a(jel(j,j))).gt.1.d0 )
c     .      then
c            print *,'i,j, aii aij ',i,j,a(jel(i,i)),a(jel(i,j))
c          endif
c        enddo
c         a(jel(i,i)) = a(jel(i,i))*(1.d0 + 1.d-6) 
c      enddo  
           
c      print *,'LCNORM after reweight '
c      call printa( 1440 ) 
c      call printan22( 490 )   

C
C NORMALIZE NORMAL EQUATIONS
      CALL NRMSCL(SCALE,NLIVE,0,1)
      call report_stat('STATUS','SOLVE','lcnorm',' '
     .          , 'Solving normal equations in LC mode',0)
c      CALL VINV2(3,NLIVE,IERINV)
      CALL INVER2(A,B,3,NLIVE,rcond,IERINV)
      CALL NRMSCL(SCALE,NLIVE,1,1)
      IF(IERINV.NE.0) then
        WRITE(message,2342)
2342    FORMAT(' Bad matrix inversion in LCNORM'
     .        ,' One of the parameters should not have been adjusted')
        call report_stat('WARNING','SOLVE','lcnorm',' ',message,0)
      endif
C
C      NOTE:  VECTOR B1  NOW ACTUALLY CONTAINS INV(A11):B1>
C
 400  CONTINUE
      RETURN
      END


