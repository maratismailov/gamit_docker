Copyright (c) Massachusetts Institute of Technology and the University of
California, San Diego. 1995.   All rights reserved.

      SUBROUTINE CHDRED ( IOBS,IPRNT,NEPOCH,INTER,NTIME,IY,IM,ID
     1                  , IHR,MIN,SEC,NCHAN,ISCHAN,NDAT,DATTYP,LAMBDA
     2                  , offarp,OFFSL1,OFFSL2,SITNAM,rcvrsw,swver
     3                  , rctype,rcvnum,anttyp,antnum
     4                  , ircint,isessn,NTEXT,TEXT 
     5                  , NSLIP, ISLIP, ISLPST)

C     Read a C-file header.  R. King 24 June 1987

      implicit none

      include '../includes/dimpar.h'

      character*1  skd,gnss
      character*3  rcvrsw
      character*4 antmod,svantmod,dryzen,wetzen,drymap,wetmap,ionsrc
      character*5 frame,precmod,nutmod,gravmod,sprmod 
     .,                  otidemod,atmtide,hydrolmod
      character*6        etidemod,magfield
      character*8        atmlmod,cextra(maxext)  
      character*10 antmod_snx,svantmod_snx
      CHARACTER*32 SITNAM
      CHARACTER*16 OBFILN,TFILN,jfiln
      CHARACTER*20 RLABEL(MAXLAB),rctype,rcvnum,anttyp,antnum
     .            , svantbody(maxsat)
      CHARACTER*80 TEXT(MAXTXT)
      character*16 latstr,lonstr

      INTEGER*4 NTEXT,NDAT,LAMBDA(MAXSAT,MAXDAT),DATTYP(MAXDAT)
     1      , ISCHAN(MAXSAT),IDMS(MAXLAB),ISLOT(MAXLAB)
     2      , NTIME,INTER,NPART,IPRNT,IHR,MIN,NLABEL,NEPOCH,ID,IM,IY
     3      , MAXPAR,NCHAN,I,J,NPARAM,IOBS,avlmet,nslip
     4      , ircint,isessn,nclock,ietide,isptide,norb,iblk
     5      , niextra,iextra(maxext),nrextra,ncextra     
      integer*2 islip(maxcsb),islpst(maxcsb)
      integer*4 jde,jdr

      PARAMETER (MAXPAR=7+MAXSAT*MAXORB)
C
      REAL*8 offarp(3),OFFSL1(3),OFFSL2(3),PREVAL(MAXPAR),rextra(maxext)
     1     , SEC
     2     , te,tr,ut1,ut1dot,xp,xpdot,yp,ypdot
     3     , psi,psidot,eps,epsdot
     4     , elvcut,clock(4)
     .     , atmlavg,hydrolavg,fL1(maxsat),fL2(maxsat)

c**rwk 190826: Temporary to avoid changing the c-file format for 10.71, we've 
c              written only a a single-frequency value for the SV antenna PCO 
* MOD TAH 1909120: Copied from model version of routine with same name.  
*     Added to readc2 call here
*     Dummy name here since svantdx is declared in includes/model.h
C TAH: Problem here that needs fixing.  In  model.h, svantdx(3,2,maxsat) is
C     declared but only the first a 3,maxsat array is written to the c-file.
C MOD TAH 200126: Made svantdx(3,2,maxsat) through out so problem goes away.
C     Commemted line below (Defined locally here so needed).
      real*8 svantdx(3,2,maxsat) 

      real*4 swver
C
C
C------------------------------------------------------------------------
C
C
C        Read the Header Records
C
      call readc1 (iobs,ntext,text )

      call readc2 ( iobs
     .            , sitnam,rctype,rcvnum,rcvrsw,swver,anttyp,antnum
     .            , npart,norb,gnss,nchan,ischan,fL1,fL2 
     .            , ndat,dattyp,lambda
     .            , skd,nepoch,inter,ircint,ntime,isessn
     .            , iy,im,id,ihr,min,sec
     .            , offarp,offsl1,offsl2, svantdx 
     .            , obfiln,tfiln,jfiln
     .            , frame,precmod,nutmod,gravmod,sprmod  
     .            , ietide,isptide
     .            , etidemod,otidemod,atmtide,atmlmod,hydrolmod
     .            , atmlavg,hydrolavg,dryzen,wetzen,drymap,wetmap
     .            , ionsrc,magfield
     .            , antmod_snx,antmod,svantbody,svantmod_snx,svantmod
     .            , elvcut,nclock,clock
     .            , jde,te,jdr,tr
     .            , ut1,ut1dot,xp,xpdot,yp,ypdot
     .            , psi,psidot,eps,epsdot
     .            , avlmet
     .            , nslip,islip,islpst     
     .            , niextra,iextra,nrextra,rextra,ncextra,cextra)

      call readc3 ( iobs,nlabel,nparam,islot,idms,rlabel,preval )
C
C
C        Print the C-File header information
C
      WRITE(IPRNT,10) (TEXT(I),I=1,NTEXT)
10    FORMAT(/,1X,'Input C-file Header Information:',//
     1        ,3X,' Comments:',/
     2        ,50(A80/) )
      call wdms(1,preval(1),latstr)
      call wdms(2,preval(2),lonstr)
      WRITE(IPRNT,15) SITNAM,PREVAL(1),latstr,PREVAL(2),lonstr,PREVAL(3)
     .              , rcvrsw,swver
15    FORMAT(/,'--------------------------',/
     1        ,1X,'Long site name: ',A32,/
     2        ,1X,'latitude  = ',f12.8,' radians or ',a16,/
     2        ,1X,'longitude = ',f12.8,' radians or ',a16,/
     3        ,1x,'radius    = ',f12.6,' kilometers',//
     4        ,1x,'Receiver, software: ',a3,1x,f5.2)

      write(iprnt,20) offarp,OFFSL1,OFFSL2
     1              , IY,ID,IM,IHR,MIN,SEC,isessn,skd,NCHAN,NEPOCH
     2              , INTER,ircint,NPART,NLABEL,NDAT
     3              , OBFILN,TFILN,jfiln,jde,te,jdr,tr,ut1,ut1dot
     4              , xp,xpdot,yp,ypdot,psi,psidot,eps,epsdot
     5              , avlmet,nslip,(islip(i),islpst(i),i=1,nslip)
20    format (1X,'Antenna offsets (Up North East) ARP:',3f7.4
     .        ,'  L1:',3F7.4,'  L2:',3F7.4,'  meters',/
     6        ,1X,'Start time:',I5,2I3,2X,2I3,F7.3,3x,'Session:',i2,//
     7        ,1X,'SKD = ',a,3x,'NCHAN = ',I2,//
     X        ,1X,'NEPOCH = ',I4,3X,'INTER=',I4,3x,'IRCINT=',i4,//
     8        ,1X,'NPART = ',I2,3X,'NORB=',I2,3X,'NLABEL=',I2
     9        ,3X,'NDAT = ',i3,/
     B        ,1X,'Raw data file: ',A16,/
     .        ,1x,'T-file:        ',A16,/
     .        ,1x,'J-file:        ',A16,/
     .        ,1x,'IC epoch:      ',i8,f12.2,/
     .        ,1x,'E-Rot. epoch:  ',i8,f12.2,/
     .        ,1x,'UT1, rate:     ',f9.6,d12.3,/
     .        ,1x,'x-pole, rate:  ',f10.5,d12.3,/
     .        ,1x,'y-pole, rate:  ',f10.5,d12.3,/
     .        ,1x,'Delta-psii,rate:',f10.5,d12.3,/
     .        ,1x,'Delta-eps,rate:',f10.5,d12.3,/
     .        ,1x,'AVLMET:        ',i4,/
     .        ,1x,'Bias flags (',i4,') : ',14i4,/,24x,14i4 )  
      WRITE(IPRNT,30)
30    FORMAT(/,2x,'PRN   LAMBDAS')

      DO 35 I=1,NCHAN
         WRITE(IPRNT,31) gnss,ISCHAN(I),(LAMBDA(I,J),J=1,NDAT)
31       FORMAT(1X,a1,i2,4I4)
35    CONTINUE

      write (iprnt,*) ' '
      write (iprnt,*) 'RLABEL    ISLOT IDMS '
      WRITE(IPRNT,40) (RLABEL(I),ISLOT(I),IDMS(I),I=1,NLABEL)
40    FORMAT (20(1X,A20,I5,I3,/))

      WRITE(IPRNT,*) 'PREVAL'
      WRITE(IPRNT,50) (PREVAL(I),I=1,NPARAM)
50    FORMAT(10(12X,3D22.14,/))

      write(iprnt,60) (islip(i),islpst(i),i=1,nslip)
60    format(1x,'Extra biases:  ',i3,/,17x,12i5,/17x,12i5)
               
      write(iprnt,'(5x,a,1x,20i4)') 
     .    'niextra iextra',(iextra(i),i=1,niextra)
      write(iprnt,'(5x,a,1x,20d12.3)') 
     .    'nrextra rextra',(rextra(i),i=1,nrextra)
      write(iprnt,'(5x,a,1x,20(1x,a8))') 
     .    'ncextra cextra',(cextra(i),i=1,ncextra)

      write (iprnt,*) ' '
      write (iprnt,*) ' '
C
      RETURN
      END
