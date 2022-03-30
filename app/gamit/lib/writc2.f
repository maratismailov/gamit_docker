      subroutine writc2 (lunit
     .,                sitnam,rctype,rcvnum,rcvrsw,swver,anttyp,antnum
     .,                npart,norb,gnss,nsat,isprn,fl1,fl2 
     .,                ndat,dattyp,lambda
     .,                skd,nepoch,inter,ircint,mtime,isessn
     .,                iy,im,id,ihr,min,sec
     .,                offarp,offsl1,offsl2, antdaz, svantdx
     .,                obfiln,tfiln,jfiln
     .,                frame,precmod,nutmod,gravmod,srpmod
     .,                eradmod,antradmod 
     .,                ietide,isptide,speopmod
     .,                etidemod,otidemod,atmtide,atmlmod,hydrolmod      
     .,                atmlavg,hydrolavg,dryzen,wetzen,drymap,wetmap
     .,                ionsrc,magfield
     .,                antmod_snx,antmod,svantbody,svantmod_snx,svantmod
     .,                elvcut,nclock,clock
     .,                jdet,tet,jdr,tr
     .,                ut1,ut1dot,xp,xpdot,yp,ypdot
     .,                psi,psidot,eps,epsdot
     .,                avlmet
     .,                nslip,islip,islpst 
     .,                niextra,iextra,nrextra,rextra,ncextra,cextra)
     
c     Write the second part of a C-file which is a header -- see readc2 for list 

c     K. Feigl/R. King  July 1989    
c     R. King Sep 1998:  Modfied to add 'norb' parameter   
c     R. King Jan 2003  Modified to add more model names, avg loading ocean/atm loading
c                       values, and more extra slots            
c     R. King Aug 2010  Modified for release 10.4 to lengthen the names of the
c                       antenna models 
c     R. King Jan 2013  Modified for release 10.41 to add dryzen and wetzen 
c     R. King Mar 2014  Added eradmod, antradmod, ionsrc, magfield
c     T. Herring Jan 2020: Converted svantdx into (2,3,NS) to save L1 and L2
c                       values (needed for Galileo).TAH 200126 (See readc2.f)
c                       Added antdaz for antenna azimuth (repr03/itrf2020)
      implicit none

      include '../includes/dimpar.h'

c     input variables:
      integer*4          nsat,norb,npart,ndat,dattyp(maxdat)    
     .,                  nepoch,inter,mtime,isprn(maxsat)
     .,                  iy,im,id,ihr,min,lambda(maxsat,maxdat)
     .,                  avlmet,nslip,ircint,isessn,ietide,isptide   
     .,                  nclock,niextra,iextra(maxext),nrextra,ncextra
      integer*2          islip(maxcsb),islpst(maxcsb)
      integer*4          jdet,jdr
      real*8             sec,offarp(3),offsl1(3),offsl2(3),tet,tr
     .,                  ut1,ut1dot,xp,xpdot,yp,ypdot
     .,                  psi,psidot,eps,epsdot
     .,                  rextra(maxext),elvcut,clock(4) 
C MOD TAH 200126: Changed to svantdx(2,3,maxsat) (added L1/L2 index)
     .,                  atmlavg(3),hydrolavg(3),svantdx(2,3,maxsat)
     .,                  fL1(maxsat),fL2(maxsat)
      real*8 antdaz       ! Antenna azimith (deg) (model.h TAH 200205))
      character*1        skd,gnss
      character*3        rcvrsw
      character*4        antmod,svantmod(maxsat),dryzen,wetzen
     .,                  drymap,wetmap,ionsrc
      character*5        frame,precmod,nutmod,gravmod,srpmod
     .,                  eradmod,antradmod   
      character*6        magfield
      character*8        etidemod,atmlmod,otidemod,atmtide,hydrolmod
     .,                  speopmod,cextra(maxext) 
      character*10       antmod_snx,svantmod_snx(maxsat)
      character*16       tfiln,obfiln,jfiln
      character*20       rctype,rcvnum,anttyp,antnum,svantbody(maxsat)
      character*32       sitnam
      real*4             swver

      integer            iflag,lunit,i,j
      integer*4          ioerr,len,rcpar
      character*16       cfname
      character*80       prog_name

      logical            debug
      data debug/.false./


      iflag = 2        
      if( debug ) then
       print *,'WRITC2: sitnam= ',sitnam
       print *,'WRITC2: rctype,rcvnum= ',rctype,rcvnum
       print *,'WRITC2: rcvrsw,swver= ',rcvrsw,swver
       print *,'WRITC2: anttyp,antnum= ',anttyp,antnum
       print *,'WRITC2: npart = ',npart
       print *,'WRITC2: norb  = ',norb
       print *,'WRITC2: gnss  = ',gnss
       print *,'WRITC2: nsat  = ',nsat
       print *,'WRITC2: isprn = ',(isprn(i),i=1,nsat)
       print *,'WRITC2: ndat  = ',ndat
       print *,'WRITC2: dattyp= ',(dattyp(i),i=1,ndat)
       do i=1,nsat      
          print *,'WRITC2: fl1 fl2 ',fl1(i),fl2(i)
          print *,'WRITC2: lambda ',(lambda(i,j),j=1,ndat)
       enddo
       print *,'WRITC2: skd   = ',skd
       print *,'WRITC2: nepoch= ',nepoch
       print *,'WRITC2: inter = ',inter
       print *,'WRITC2: ircint= ',ircint
       print *,'WRITC2: mtime = ',mtime
       print *,'WRITC2: isessn = ',isessn
       print *,'WRITC2: iy,im,id,ihr,min,sec ',iy,im,id,ihr,min,sec
       print *,'WRITC2: offapr= ',(offarp(i),i=1,3)
       print *,'WRITC2: offsl1= ',(offsl1(i),i=1,3)
       print *,'WRITC2: offsl2= ',(offsl2(i),i=1,3)
       print *,'WRITC2: svantdx= ',(svantdx(:,:,j),j=1,nsat)
       print *,'WRITC2: obfiln= ',obfiln
       print *,'WRITC2: tfiln = ',tfiln
       print *,'WRITC2: jfiln = ',jfiln
       print *,'WRITC2: frame,precmod,nutmod= ',frame,precmod,nutmod
       print *,'WRITC2: gravmod,srpmod,eradmod,antradmod= '
     .                  ,gravmod,srpmod,eradmod,antradmod
       print *,'WRITC2: ietide,isptide= ',ietide,isptide  
       print *,'WRITC2: speopmod,etidemod,otidemod,atmtide= '
     .                 ,speopmod,etidemod,otidemod,atmtide
       print *,'WRITC2: atmlmod,hydrolmod= ',atmlmod,hydrolmod       
       print *,'WRITC2: dryzen,wetzen,drymap,wetmap= '
     .                 ,dryzen,wetzen,drymap,wetmap  
       print *,'WRITC2: ionsrc,magfield= ',ionsrc,magfield
       print *,'WRITC2: antmod_snx,antmod,elvcut= '
     .                 ,antmod_snx,antmod,elvcut
       do i=1,nsat
         print *,'WRITC2: svantbody,svantmod_snx, svantmod= '
     .             ,svantbody(i),svantmod_snx(i),svantmod(i) 
       enddo        
       print *,'WRITC2: nclock,clock= ',nclock,(clock(i),i=1,nclock)
       print *,'WRITC2: jdet,tet = ',jdet,tet
       print *,'WRITC2: jdr,tr = ',jdr,tr
       print *,'WRITC2: ut1,ut1dot= ',ut1,ut1dot
       print *,'WRITC2: xp,xpdot= ',xp,xpdot
       print *,'WRITC2: yp,ypdot= ',yp,ypdot
       print *,'WRITC2: psi,psidot,eps,epsdot= ',psi,psidot,eps,epsdot
       print *,'WRITC2: avlmet= ',avlmet
       print *,'WRITC2: nslip= ',nslip
       print *,'WRITC2: islip,islpst= ',(islip(i),islpst(i),i=1,nslip)
       print *,'WRITC2: niextra= ',niextra
       print *,'WRITC2: iextra = ',(iextra(i),i=1,niextra)     
       print *,'WRITC2: nrextra= ',nrextra
       print *,'WRITC2: rextra = ',(rextra(i),i=1,nrextra)
       print *,'WRITC2: ncextra= ',ncextra
       print *,'WRITC2: cextra = ',(cextra(i),i=1,ncextra)
      endif
c
c     get calling program name and C-file name for report_stat
      len = rcpar(0,prog_name)
      inquire( unit=lunit, name=cfname, iostat=ioerr )
      if( ioerr.ne.0 ) goto 1000

      write (unit    =  lunit
     .,     iostat  =  ioerr
     .,     err     =  1000)
     .      iflag
     .,     sitnam,rctype,rcvnum,rcvrsw,swver
     .,     anttyp,antnum
     .,     npart,norb,gnss,nsat,(isprn(i),i=1,nsat)
     .,     (fL1(i),fL2(i),i=1,nsat)
     .,     ndat,(dattyp(i),i=1,ndat)
     .,     ((lambda(i,j),i=1,nsat),j=1,ndat)
     .,     skd,nepoch,inter,ircint,mtime,isessn
     .,     iy,im,id,ihr,min,sec
     .,     (offarp(i),i=1,3),(offsl1(i),i=1,3),(offsl2(i),i=1,3)
     .,     antdaz
     .,     (svantdx(:,:,j),j=1,nsat)
     .,     obfiln,tfiln,jfiln
     .,     frame,precmod,nutmod,gravmod,srpmod,eradmod,antradmod
     .,     ietide,isptide 
     .,     speopmod,etidemod,otidemod,atmtide,atmlmod,hydrolmod  
     .,     (atmlavg(i),i=1,3),(hydrolavg(i),i=1,3)
     .,     dryzen,wetzen,drymap,wetmap,ionsrc,magfield
     .,     antmod_snx,antmod
     .,     (svantbody(i),i=1,nsat),(svantmod_snx(i),i=1,nsat)
     .,     (svantmod(i),i=1,nsat)
     .,     elvcut,nclock,(clock(i),i=1,nclock)
     .,     jdet,tet,jdr,tr
     .,     ut1,ut1dot,xp,xpdot,yp,ypdot,psi,psidot,eps,epsdot
     .,     avlmet
     .,     nslip,(islip(i),i=1,nslip),(islpst(i),i=1,nslip)
     .,     niextra,(iextra(i),i=1,niextra)
     .,     nrextra,(rextra(i),i=1,nrextra)
     .,     ncextra,(cextra(i),i=1,ncextra)   

 1000 if (ioerr .ne. 0) then
         call ferror (ioerr,6)
         call report_stat('FATAL',prog_name,'lib/writc2',cfname
     .                   ,'Error writing second C-file record',ioerr)
      endif

      return
      end
