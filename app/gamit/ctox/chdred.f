Copyright (c) Massachusetts Institute of Technology and the University of
California, San Diego, 1995. All rights reserved.

      Subroutine CHDRED( dump,iy,im,id,ihr,min,sec )      

C     Read a C-file header.  R. King 24 June 1987/last revised 29 November 2014
C
      implicit none

      include '../includes/dimpar.h'    
      include '../includes/units.h'
      include '../includes/global.h'
      include '../includes/model.h'

      LOGICAL dump

      character*1 skd 
      character*4  c_time_type
      character*16 latstr,lonstr

      character*80 SV_format  ! Format string for sv output due to 
                              ! possible difference number of observables

      INTEGER*4 IDMS(MAXLAB),ISLOT(MAXLAB),iy,im,id,ihr,min
     3      , I,J
     4      , nclock,norb
      integer*4 jdr
c     function
      integer*4 julday 

C
      REAL*8 SEC,tr,clock(4)
c**rwk 190826: Temporary to avoid changing the c-file format for 10.71, we've 
c              written only a a single-frequency value for the SV antenna PCO 
* MOD TAH 1909120: Copied from model version of routine with same name.  
*     Dummy name here since svantdx is declared in includes/model.h
C TAH: Problem here that needs fixing.  In  model.h, svantdx(3,2,maxsat) is
C     declared but only the first a 3,maxsat array is written to the c-file.
C MOD TAH 200126: Made svantdx(3,2,maxsat) through out so problem goes away.
C     Commemted line below.
!     real*8 svantdxx(3,maxsat) 

C
C
C------------------------------------------------------------------------
C
C
C        Read the Header Records
C          
      call readc1 (iuc,ntext,text )

      call readc2 ( iuc                    
     .            , sitnam,rctype,rcvrsn,rcvrswx,swver,anttyp,antsn
     .            , npart,norb,gnss,nchan,ischan,fL1,fL2
     .            , ndat,dattyp,lambda
     .            , skd,nepoch,inter,ircint,mtime,isessn
     .            , iy,im,id,ihr,min,sec
     .            , offarp,offl1,offl2, antdaz, svantdx
     .            , obfiln,tfiln,jfiln
     .            , frame,precmod,nutmod,gravmod,srpmod     
     .            , eradmod,antradmod
     .            , ietide,isptide,speopmod
     .            , etidemod,otidemod,atidemod,atmlmod,hydrolmod    
     .            , atmlavg,hydrolavg,dryzen,wetzen,drymap,wetmap
     .            , ionsrc,magfield
     .            , antmod_snx,antmod,svantbody,svantmod_snx,svantmod
     .            , elvcut,nclock,clock
     .            , jdet,tet,jdr,tr
     .            , ut1,ut1dot,xp,xpdot,yp,ypdot
     .            , psi,psidot,eps,epsdot
     .            , avlmet
     .            , nslip,islip,islpst   
     .            , niextra,iextra,nrextra,rextra,ncextra,cextra)          
cd      print *,'CHDRED anttyp swver rctype irctint '
cd     .     ,anttyp,swver,rctype,ircint
                        
      call readc3 ( iuc,nlabel,nparam,islot,idms,rlabel,preval )   


c        Print the C-File header information

      if (dump) write(iprnt,10) ntext,(text(i),i=1,ntext)
      write(iscrn,10) ntext,(text(i),i=1,ntext)
10    FORMAT(/,1X,'Input C-file Header Information ',//
     1        ,3X,' Comments (ntext=',i3,'): ',/
     2        ,50(A80/) )
      call wdms(1,preval(1),latstr)
      call wdms(2,preval(2),lonstr)
      if (dump)       
     .WRITE(IPRNT,15) SITNAM,PREVAL(1),latstr,PREVAL(2),lonstr,PREVAL(3)
     .              , rctype,rcvrsn,anttyp,antsn,rcvrswx,swver
      WRITE(iscrn,15) SITNAM,PREVAL(1),latstr,PREVAL(2),lonstr,PREVAL(3)
     .              , rctype,rcvrsn,anttyp,antsn,rcvers,swver
15    FORMAT(/,'--------------------------',/
     1        ,1X,'Long site name: ',A32,/
     2        ,1X,'latitude  = ',f12.8,' radians or ',a16,/
     2        ,1X,'longitude = ',f12.8,' radians or ',a16,/
     3        ,1x,'radius    = ',f12.6,' kilometers',//
     4        ,1x,'rctype=',a20,' rcvrsn=',a20
     5        ,1x,'anttyp=',a20,' antsn=',a20,/
     4        ,1x,'Receiver, software: ',a3,1x,f5.2)

      if( mtime.eq.1 ) then
        c_time_type = 'UTC '                          
      elseif (mtime.eq.2 ) then
        c_time_type = 'GPST'
      else
        call report_stat('FATAL','CTOX','chdred',mtime,
     .  'Bad C-file time type ',0)
      endif
      write(iscrn,'(/,a,a4)') 'C-file time is ',c_time_type
      if( dump ) then
        write(iprnt,'(a,3(3f7.4))') 
     .    ' Antenna offsets (U N E) ARP:',offarp,offl1,offl2
        write(iprnt,'(a,i5,4i3,2x,f7.3)')
     .    ' Start time:,',iy,id,im,ihr,min,sec
        write(iprnt,'(a,3(1x,a16))') 
     .    ' X, J, T files:',obfiln,tfiln,jfiln
        write(iprnt,'(a,a1,1x,a1,i5,i5,i5)')
     .    ' gnss, skd,nepoch inter ircint ',gnss,skd,nepoch,inter,ircint
        write(iprnt,'(a,5i4)')
     .    ' mtime ndat chan npart ',mtime,ndat,nchan,npart
        write(iprnt,'(a,7(1x,a8))')
     .    ' frame precmod nutmod gravmod srpmod eradmod antradmod: '
     .      ,frame,precmod,nutmod,gravmod,srpmod,eradmod,antradmod
        write(iprnt,'(a,2i4,2(1x,a8))')
     .    ' ietide,isptide,speopmod,etidemod:'
     .     ,ietide,isptide,speopmod,etidemod    
        write(iprnt,'(a,4(1x,a8))')
     .     ' otidemod atidemod,atmlmod,hydrolmod:'
     .      ,otidemod,atidemod,atmlmod,hydrolmod
        write(iprnt,'(a,4(1x,a4))')
     .    ' dryzen wetzen drymap wetmap: '
     .     ,dryzen,wetzen,drymap,wetmap
        write(iprnt,'(a,a4,1x,a6)')
     .    ' ionsrc magfield: ',ionsrc,magfield
        write(iprnt,'(a,a10,1x,a4,f6.2,i4)')
     .    ' antmod_snx antmod elvcut avlmet: '
     .     ,antmod_snx,antmod,elvcut,avlmet
        write(iprnt,'(a,i2,4d12.4)')
     .    ' nclock clock:   ',nclock,(clock(i),i=1,nclock)
        write(iprnt,'(a,i8,f12.2)')
     .    ' IC epoch:       ',jdet,tet
        write(iprnt,'(a,i8,f12.2)')
     .    ' E-Rot epoch:    ',jdr,tr
        write(iprnt,'(a,f10.6,d12.3)')
     .    ' UT1, rate:      ',ut1,ut1dot
        write(iprnt,'(a,f10.6,d12.3)')
     .    ' x-pole, rate:   ',xp,xpdot
        write(iprnt,'(a,f10.6,d12.3)')
     .    ' y-pole, rate:   ',yp,ypdot
        write(iprnt,'(a,f10.6,d12.3)')
     .    ' Delta-psi,rate: ',psi,psidot
        write(iprnt,'(a,f10.6,d12.3)')
     .    ' Delta-eps,rate: ',eps,epsdot 
* MOD TAH 101111: Added back output antenna offsets               
        write(iprnt,'(2a)')
     .    ' Chan PRN  svantbody            fL1       fL2      ',
     .    ' svantmod_snx svantmod  Lambdas       Offsets'
*       Create format that is needed
        write(SV_format,21) ndat
cc21    format('(i3,2i5,4x,a10,1x,a4,',I2.2,'i3,1x,3f7.3))')
  21    format('(i3,2x,a1,i2,1x,a20,1x,2f10.3,2x,a10,4x,a4,3x,',I2.2
     .        ,'i3,1x,3f7.3,1x,3f7.3))')
* MOD TAH 190920: Changed svantdx to svantdxx
* MOD TAH 200126: Updated for full 2,3,nchan svantdx array
        do i=1,nchan
          write(iprnt,SV_format)  i,gnss,ischan(i),svantbody(i)
     .,    fl1(i)*1.d-6,fl2(i)*1.d-6
     .,    svantmod_snx(i),svantmod(i),(lambda(i,j),j=1,ndat)
     .,    (svantdx(:,j,i),j=1,2)
        enddo   
      endif
      if ( dump ) write(iprnt,23) nslip,(islip(i),islpst(i),i=1,nslip)    
      if( dump ) then
         write(iprnt,'(1x,a,i3,20i6)') 'niextra iextra: '
     .            ,niextra,(iextra(i),i=1,niextra)   
         write(iprnt,'(1x,a,i3,20d12.3)') 'nrextra rextra: '
     .            ,nrextra,(rextra(i),i=1,nrextra)   
         write(iprnt,'(1x,a,i3,20a8)') 'ncextra cextra: '
     .            ,ncextra,(cextra(i),i=1,ncextra)   
      endif
                                               
  23   format(1x,'Bias flags (',i4,') : ',100(/,2i8))

      if ( dump ) then

        write (iprnt,*) ' '
        write (iprnt,*) 'rlabel    islot idms '
        WRITE(IPRNT,40) (RLABEL(I),ISLOT(I),IDMS(I),I=1,NLABEL)
40      FORMAT (29(1X,A20,I5,I3,/))

        WRITE(IPRNT,*) 'preval'
        WRITE(IPRNT,50) (PREVAL(I),I=1,NPARAM)
50      FORMAT(10(12X,3D22.14,/))

        write(iprnt,60) (islip(i),islpst(i),i=1,nslip)
60      format(1x,'Extra biases:  ',i3,/,17x,12i5,/17x,12i5)
                        
        write(iprnt,*) ' '
        write(iprnt,*) 'niextra= ',niextra
        write(iprnt,*) 'iextra = ',(iextra(i),i=1,niextra)     
        write(iprnt,*) 'nrextra= ',nrextra
        write(iprnt,*) 'rextra = ',(rextra(i),i=1,nrextra)
        write(iprnt,*) 'ncextra= ',ncextra
        write(iprnt,*) 'cextra = ',(cextra(i),i=1,ncextra)
        write(iprnt,*) ' '
        write(iprnt,*) ' ' 

      endif

c   Save the start time for storing in model.h and writing the x-file 
      jd0 = julday(im,id,iy)      
      t0  = ihr*3600.d0 + min*60.d0 + sec
      
c RWK 160419 This remmoved in favor of preval, stored in model.h
c     TMPLAT=PREVAL(1)
c     TMPLON=PREVAL(2)
c     TMPRAD=PREVAL(3)

cd    print *,'end of chdred '
      RETURN
      END







