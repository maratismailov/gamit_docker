Copyright (c) Massachusetts Institute of Technology and the University of
California at San Diego, 1995.  All rights reserved.

      Subroutine getcfs(iprnt,iuc,filnam,ioerr,
     $                  mdual,nparam,nsat,npart,nepoch,isprn,
     $                  obfiln,tfiln,sitnam,elvcut,
     $                  preval,nlabel,clabel,islot,
     $                  norb,inter,it0,t00)

c Purpose: 

c     Open c-files and retrieve necessary header information

      implicit none

      include '../includes/dimpar.h'

      integer dattyp(maxdat)
      integer*4 i
      integer*4 id
      integer*4 ihr
      integer*4 im
      integer*4 inter
      integer*4 ioerr
      integer*4 iphase
      integer*4 iprnt
      integer*4 iuc
      integer*4 it0(3)
      integer*4 iy
      integer*4 idms(maxlab)
      integer*4 isprn(maxsat)
      integer*4 islot(maxlab)
      integer*4 lambda(maxsat,maxdat)
      integer*4 min
      integer*4 mtime
      integer*4 ndat
      integer*4 nepoch
      integer*4 iblk(maxsat)
      integer*4 niextra,iextra(maxext),nrextra,ncextra
      integer*4 nlabel
      integer*4 norb
      integer*4 nparam
      integer*4 nsat
      integer*4 npart
      integer*4 ntext
      integer*4 avlmet,nslip,ircint,isessn
      integer*2 islip(maxcsb),islpst(maxcsb)
      integer*4 jde,jdr
      integer*4 ietide,isptide,nclock

      logical mdual

      real*4 swver

      real*8  rextra(maxext)
* MOD TAH 200127: Addeed L1/L2 frequencies to svantdx
      real*8  offarp(3),offsl1(3),offsl2(3),svantdx(3,2,maxsat)
      real*8  antdaz  ! Antenna aligment from True N (deg).
      real*8  preval(maxprm)
      real*8  sec
      real*8  t00(3)
      real*8  te,tr
      real*8  ut1,ut1dot,xp,xpdot,yp,ypdot
      real*8  psi,psidot,eps,epsdot
      real*8  elvcut,clock(4)   
      real*8  atmlavg(3),hydrolavg(3)                
      real*8  fL1(maxsat),fL2(maxsat)

      character*1  skd,gnss
      character*3  rcvrsw
      character*4  antmod,svantmod(maxsat)
     .,            dryzen,wetzen,drymap,wetmap,ionsrc
      character*5 frame,precmod,nutmod,gravmod,srpmod,eradmod,antradmod
      character*6 magfield
      character*8 etidemod,atmlmod,otidemod,atmtide,hydrolmod
     .,           speopmod,cextra(maxext) 
      character*10 antmod_snx,svantmod_snx(maxsat)
      character*12 sitnam
      character*32 sitnam32
      character*16 filnam
      character*16 lowerc
      character*16 obfiln
      character*16 tfiln,jfiln
      character*20 clabel(maxlab),rctype,rcvnum,anttyp,antnum
     .           , svantbody(maxsat)
      character*80 text(maxtxt)
      character*256 message
                   
c      print *,'GETCFS filnam ',filnam

      call copens (lowerc(filnam), 'old', iuc, ioerr)
      if (ioerr.ne.0) return
c
      call readc1 (iuc,ntext,text)
      call readc2 (iuc
     .,            sitnam32,rctype,rcvnum,rcvrsw,swver,anttyp,antnum
     .,            npart,norb,gnss,nsat,isprn,fL1,fL2
     .,            ndat,dattyp,lambda
     .,            skd,nepoch,inter,ircint,mtime,isessn
     .,            iy,im,id,ihr,min,sec
     .,            offarp,offsl1,offsl2,antdaz, svantdx
     .,            obfiln,tfiln,jfiln
     .,            frame,precmod,nutmod,gravmod,srpmod,eradmod,antradmod
     .,            ietide,isptide,speopmod
     .,            etidemod,otidemod,atmtide,atmlmod,hydrolmod  
     .,            atmlavg,hydrolavg,dryzen,wetzen,drymap,wetmap    
     .,            ionsrc,magfield
     .,            antmod_snx,antmod,svantbody,svantmod_snx,svantmod
     .,            elvcut,nclock,clock
     .,            jde,te,jdr,tr
     .,            ut1,ut1dot,xp,xpdot,yp,ypdot
     .,            psi,psidot,eps,epsdot
     .,            avlmet
     .,            nslip,islip,islpst  
     .,            niextra,iextra,nrextra,rextra,ncextra,cextra)
      call readc3 (iuc
     .,            nlabel,nparam
     .,            islot,idms,clabel,preval)
                          
      close (unit=iuc)

c     Get the # of orbital parameters from radiation pressure model 
c     ** added to C-file 980911 -- no longer needed
c      norb = 0 
c      if( npart.ge.6 ) then
c        if( srpmod.eq.'BERNE' .or. srpmod.eq.'BERN1' ) then
c          norb = 15
c        elseif( srpmod.eq.'BERN2 ' ) then
c          norb = 12
c       else
c          norb = 9
c        endif
c      endif

c     currently only c*12 site descriptions are actually created
      sitnam = sitnam32(1:12)
      it0(1) = im
      it0(2) = id
      it0(3) = iy
      t00(1) = ihr
      t00(2) = min
      t00(3) = sec
c
c     determine if L1 and L2 exist
c
      iphase = 0
      do 10 i = 1, ndat
         if (dattyp(i).eq.1 .or. dattyp(i).eq.2) iphase=iphase+1
   10 continue
      if (iphase.eq.1) then
         mdual = .false.
      elseif (iphase.eq.2) then
         mdual = .true.
      else
         write(message,'(a,i2)')
     1      'Error, bad number of phase observables: ',iphase
         call report_stat('FATAL','CFMRG','getcfs',' ',message,0)
      endif

      write(message,'(2a)') ' Site: ',sitnam
c      call report_stat('STATUS','CFMRG','getcfs',' ',message,0)
c     write (*,'(1x,2a)')    'Site: ',site
      write (iprnt,'(1x,2a)')    'Site: ',sitnam
c      write (*,'(1x,a,32i4)') 'Satellites: ',(isprn(i),i=1,nsat)
* MOD TAH 200618: Updated 32I to 50I to allow for 35 Beidou satellites
      write (iprnt,'(1x,a,50i4)') 'Satellites: ',(isprn(i),i=1,nsat)
c
      if (mdual) then
c         write(*,'(1x,a)') 'Dual-frequency C-file'
         write(iprnt,'(1x,a)') 'Dual-frequency C-file'
      else
c         write(*,'(1x,a)') 'Single-frequency C-file'  
         write(iprnt,'(1x,a)') 'Single-frequency C-file'
      endif

      return
      end
