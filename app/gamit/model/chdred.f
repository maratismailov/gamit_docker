Copyright (c) Massachusetts Institute of Technology and the University of
California, San Diego. 1994.   All rights reserved.
  
       Subroutine CHDRED ( iy,im,id,ihr,min,sec ) 

C     Read a C-file header.  R. King 24 June 1987
c     Last modified R King 28 November 2014

      implicit none

      include '../includes/dimpar.h'
      include '../includes/units.h'    
      include '../includes/global.h'
      include '../includes/model.h'
         
      integer*4 iy,im,id,ihr,min,norb 
     .        , idms(maxlab),islot(maxlab)
     .        , nclock,jde,jdr,i,j

c     function
      integer*4 julday
   
      real*8  tr,clock(4),sec
                     
c     -the SKD codes now local and not used
      character*1 skd
      character*16 xfiln 
                                                
c**rwk 190826: Temporary to avoid changing the c-file format for 10.71, we've 
c              written only a a single-frequency value for the SV antenna PCO 
C TAH: Problem here that needs fixing.  In  model.h, svantdx(3,2,maxsat) is
C     declared but only the first a 3,maxsat array is written to the c-file.
C MOD TAH 200126: Made svantdx(3,2,maxsat) through out so problem goes away.
C     Commemted line below.
!     real*8 svantdxx(3,maxsat) 
 

                   
C-----------------------------------------------------------------------

c     ***read the c-file headers***

      call readc1 (iobs
     .,            ntext,text)

      call readc2 (iobs
     .,            sitnam,rctype,rcvrsn,rcvers,swver,anttyp,antsn
     .,            npart,norb,gnss,nchan,ischan,fL1,fL2
     .,            ndat,dattyp,lambda
     .,            skd,nepoch,inter,ircint,mtime,isessn
     .,            iy,im,id,ihr,min,sec
     .,            offarp,offl1,offl2, antdaz, svantdx
     .,            xfiln,tfiln,jfiln
     .,            frame,precmod,nutmod,gravmod,srpmod                 
     .,            eradmod,antradmod
     .,            ietide,isptide,speopmod
     .,            etidemod,otidemod,atidemod,atmlmod,hydrolmod 
     .,            atmlavg,hydrolavg,dryzen,wetzen,drymap,wetmap  
     .,            ionsrc,magfield
     .,            antmod_snx,antmod,svantbody,svantmod_snx,svantmod
     .,            elvcut,nclock,clock
     .,            jdet,tet,jdr,tr
     .,            ut1,ut1dot,xp,xpdot,yp,ypdot
     .,            psi,psidot,eps,epsdot
     .,            avlmet
     .,            nslip,islip,islpst
     .,            niextra,iextra,nrextra,rextra,ncextra,cextra)

c** rwk 190826; see above 
* MOD TAH 200126: Removed code below becuase not needed anymore with
*     L1/L2 PCO svantdx saved in c-file.   
C     do j=1,nchan
C       do i=1,3
C         svantdx(i,1,j) = svantdxx(i,j)   
C         svantdx(i,2,j) = svantdxx(i,j)
C       enddo
C    enddo 

      call check_y2k(iy)

C     Read the third header, although none of the information is needed

      call readc3 (iobs,nlabel,nparam,
     .             islot,idms,rlabel,preval)


C        Print the C-File header information
C
            write(iprnt,25)
25    format('----------------------------------------------------',/,
     1       1x,'** INPUT C-FILE HEADER INFORMATION **',/)
      WRITE(IPRNT,'(A80)') (TEXT(I),I=1,NTEXT)
      WRITE(IPRNT,'(//,3X,A,//)') 'Record 2 :'
      WRITE(IPRNT,'(1X,2A,/)') 'Site: ',SITNAM
      WRITE(IPRNT,'(1X,A,2F12.8,A,F14.6,/)')
     1     'Coordinates (lat,lon,rad):',PREVAL(1),PREVAL(2)
     2,    ' radians',PREVAL(3)
      write(iprnt,'(1x,a,a20,a,a20)')
     .     'Receiver type: ',rctype,'   Serial #: ',rcvrsn
      write(iprnt,'(1x,a,a3,a,f5.2)')
     .     'Receiver software: ',rcvrswx,'   Version #: ',swver
      write(iprnt,'(1x,a,a20,a,a20)')
     .     'Antenna type: ',anttyp,'   Serial #: ',antsn
      WRITE(IPRNT,'(1X,A,3F7.4,A,3F7.4,A,3F7.4,A,/)')
     1     'Antenna offsets (Up North East) ARP:',offarp,'  L1:',offl1
     2,    '  L2:',offl2,'  meters'                             
      WRITE(IPRNT,'(1x,a,a,3x,a,i2,3x,a,I4,3x,a,i4)')
     1   'GNSS=',gnss,'mtime=',mtime,'nepoch=',nepoch,'inter=',inter
     2 , 'ircint=',ircint
      WRITE(IPRNT,'(1X,A,I5,2I3,2X,2I3,F7.3,a,i2,/)')
     1     'Start time:',IY,ID,IM,IHR,MIN,SEC ,'  Session:',isessn
      WRITE(IPRNT,'(1X,A,I2,3X,A,I2)') 'nchan=',nchan,'ndat=',ndat
      WRITE(IPRNT,'(1X,A,10I4,/)') 'dattyp=',(dattyp(i),i=1,ndat)
      write(iprnt,'(a)') 'Chan PRN     lambda         fL1           fL2'
      do i=1,nchan
         write(iprnt,'(2i4,4i4,2f12.8)') 
     .       i,ischan(i),(lambda(i,j),j=1,ndat),fl1(i),fl2(i)
      enddo
      WRITE(IPRNT,'(/,1X,4A,/)')
     1   'Raw data file: ',xfiln,'  T-file: ',TFILN,'  J-file: ',jfiln
      write(iprnt,'(1x,a,7(1x,a5))')
     .   'Frame and models: ',frame,precmod,nutmod,gravmod,srpmod
     .         ,eradmod,antradmod
      write(iprnt,'(1x,a,2i4)') 'ietide isptide: ',ietide,isptide
      write(iprnt,'(1x,a,3(1x,a8))') 
     .    'speopmod etidemod otidemod: ',speopmod,etidemod,otidemod
      write(iprnt,'(1x,a,3(1x,a8))') 
     .    'atidemod atmlmod hydrolmod: ',atidemod,atmlmod,hydrolmod  
      write(iprnt,'(1x,2(a,1x,a4,1x,a4))') 'Zenith-delay models: '
     .   ,dryzen,wetzen,'  Mapping functions: ',drymap,wetmap
      write(iprnt,'(1x,a4,1x,a6)') 'Ionosphere models ',ionsrc,magfield
      write(iprnt,'(1x,a,a10,1x,a4,f7.2)')
     .      'Ground antenna model, elev-cut: ',antmod_snx,antmod,elvcut
      write(iprnt,'(a)') 'Satellites PN  svantbody  ant model ',
     .                   ' L1 offsets     L2 offsets'
      do i=1,nchan
* MOD TAH 200126: Updated to output L1 and L2 offsets.
         write(iprnt,
     .      '(2x,i2,3x,a1,i2,1x,a20,1x, a10,1x,a3,2(1x,3f7.3))') 
     .      i,gnss,ischan(i),svantbody(i)
     .      ,svantmod_snx(i),svantmod(i),(svantdx(:,j,i),j=1,2)
      enddo
      write(iprnt,'(1x,a,i2,4d12.3)')    
     .      'nclock clock: ', nclock,(clock(i),i=1,nclock)
      WRITE(IPRNT,110) nslip,(islip(i),islpst(i),i=1,nslip)
110   FORMAT(/,3X,'Bias flags (',i4,') : ',14i4,/,24x,14i4)
      write(iprnt,120) avlmet
120   format(/,5x,'Met data (series binary flag) :',i3)
      WRITE(IPRNT,'(5X,A,I2,3X,A,I2,3X,A,I2)')
     1   'NPART=',NPART,'NLABEL=',NLABEL
      write(iprnt,'(5x,a,1x,20i4)') 
     .    'niextra iextra',(iextra(i),i=1,niextra)
      write(iprnt,'(5x,a,1x,20d12.3)') 
     .    'nrextra rextra',(rextra(i),i=1,nrextra)
      write(iprnt,'(5x,a,1x,20(1x,a8))') 
     .    'ncextra cextra',(cextra(i),i=1,ncextra)
      WRITE(IPRNT,'(//,3X,A,//)') 'Record 3 :'
      WRITE(IPRNT,40) (RLABEL(I),ISLOT(I),IDMS(I),I=1,NLABEL)
40    FORMAT(5X,'RLABEL, ISLOT, IDMS =',A20,I5,I3,14(/,26X,A20,I5,I3))
      WRITE(IPRNT,50) (PREVAL(I),I=1,NPARAM)
50    FORMAT(5X,'PREVAL=',3D22.14,1000(/,12X,3D22.14))
C         

c   Convert the start time to JD, seconds-of-day for storing in model.h

      jd0 = julday(im,id,iy)      
      t0  = ihr*3600.d0 + min*60.d0 + sec

      RETURN
      END
