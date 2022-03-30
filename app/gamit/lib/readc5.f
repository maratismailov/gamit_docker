      subroutine readc5 (lunit
     .,                  iprn
     .,                  elev,elvdot,azim,azmdot,nadang
     .,                  atmdel,svcepc,svcl1
     .,                  tau,drate
     .,                  ierfl,data_flag
     .,                  ndats,obsv,omcs,isnr,ampl1,ampl2
     .,                  nspare,spare
     .,                  nparts,tmpart)

c     read the fifth part of a C-file which has the data and partials
c     for one sat at one station at one epoch
c
c     notes: Several integer variables are stored as integer*2 in the C-file,
c            but passed around in the GAMIT software as integer*4.
c
c            Similarly, the met variables elev, etc. are stored as real*4 in
c            the C-files, but real*8 everywhere in the GAMIT package.
c
c            This is to save space, but avoid confusion.
c
      implicit none

      include '../includes/dimpar.h'

c     ***passed***
      integer*4           lunit
      integer*4           iprn
      real*8              elev,elvdot,azim,azmdot,nadang
      real*8              atmdel
      real*8              svcepc,svcl1
      real*8              tau,drate
      integer*4           ierfl,data_flag
      integer*4           ndats
      real*8              obsv(maxdat),omcs(maxdat)
      real*4              ampl1,ampl2
      integer*4           isnr(maxdat)
      integer*4           nspare
      real*8              spare(maxspr)
      integer*4           nparts
      real*8              tmpart(maxlab)

c     ***C-file and internal to this routine***
      integer             i,iflag,len,rcpar
      integer*2           iprn2
      integer*2           ierfl2
      integer*2           ndats2
      integer*2           isnr2(maxdat)
      integer*2           nparts2
      integer*2           nspare2
      integer*4           ioerr
      real*4              elev4,elvdot4,azim4,azmdot4,nadang4
      real*4              atmdel4
      character*80        prog_name

c RWK 150205: ndats and nparts seem to be included on both Record 2 and Record 5,
c            with no logical distinction since they cannot change during the session.
               
c     get the calling module name for report_stat
      len = rcpar(0,prog_name)

      read (unit    =  lunit
     .,     iostat  =  ioerr
     .,     end     =  1000
     .,     err     =  1000)
     .      iflag
     .,     iprn2
     .,     elev4,elvdot4,azim4,azmdot4,nadang4
     .,     atmdel4
     .,     svcepc,svcl1
     .,     tau,drate
     .,     ierfl2,data_flag
     .,     ndats2
     .,     (obsv  (i),i=1,ndats2)
     .,     (omcs  (i),i=1,ndats2)
     .,     (isnr2 (i),i=1,ndats2)
     .,     ampl1,ampl2
     .,     nspare2,(spare(i),i=1,nspare2)
     .,     nparts2,(tmpart(i),i=1,nparts2)

 1000 if (ioerr .ne. 0) then    
         call report_stat('FATAL',prog_name,'lib/readc5',' '
     .                   ,'Error reading C-file',ioerr)
      endif

      if (iflag .ne. 5) then  
         call report_stat('FATAL',prog_name,'lib/readc5',' '
     .                   ,'Wrong iflag',0)
      endif

      do 10 i = 1, ndats2
         isnr(i) = isnr2(i)
   10 continue

c     convert i*2 variables to i*4
      iprn     = iprn2
      ierfl    = ierfl2
      ndats    = ndats2
      nspare   = nspare2
      nparts   = nparts2

c     convert real*4 to real*8
      elev   = elev4
      elvdot = elvdot4
      azim   = azim4
      azmdot = azmdot4  
      nadang = nadang4
      atmdel = atmdel4

c      print *,'READC5: iprn ',           iprn
c      print *,'READC5: elev, elvdot ',   elev,elvdot
c      print *,'READC5: azim, azmdot ',   azim,azmdot
c      print *,'READC5: nadang '      ,   nadang
c      print *,'READC5: atmdel ', atmdel
c      print *,'READC5: svcepc,svcl1 ',   svcepc,svcl1
c      print *,'READC5: tau,drate ',      tau,drate
c      print *,'READC5: ierfl data_flag ',ierfl,data_flag
c      print *,'READC5: ndats ',          ndats
c      print *,'READC5: obsv ',           (obsv  (i),i=1,ndats)
c      print *,'READC5: omcs ',           (omcs  (i),i=1,ndats)
c      print *,'READC5: isnr ',           (isnr  (i),i=1,ndats)
c      print *,'READC5: ampl1, ampl2 ',   ampl1,ampl2
c      print *,'READC5: nspare ',         nspare
c      print *,'READC5: spare ',          (spare(i),i=1,nspare)
c      print *,'READC5: nparts ',         nparts
c      print *,'READC5: tmpart ',         (tmpart(i),i=1,nparts)

      return
      end
