      subroutine writc5 (lunit
     .,                  iprn
     .,                  elev,elvdot,azim,azmdot,nadang
     .,                  atmdel
     .,                  svcepc,svcl1
     .,                  tau,drate
     .,                  ierfl,data_flag
     .,                  ndats,obsv,omcs,isnr,ampl1,ampl2
     .,                  nspare,spare
     .,                  nparts,tmpart)


c     write the fifth part of a C-file which has the data and partials
c     for one sat at one station at one epoch
c
c     notes: Several integer variables are stored as integer*2 in the C-file,
c            but passed around in the GAMIT software as integer*4.
c
c            Similarly, the variables elev and azim are stored as real*4 in the C-files,
c            but real*8 everywhere in the GAMIT package.
c
c            This is to save space, but avoid confusion.
c

      implicit none

      include '../includes/dimpar.h'
      include '../../libraries/includes/const_param.h'  ! Temporary for debug

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
      integer             i
      integer             iflag
      integer*2           iprn2
      integer*2           ierfl2
      integer*2           ndats2
      integer*2           isnr2(maxdat)
      integer*2           nparts2
      integer*2           nspare2
      integer*4           ioerr
      integer*4 len,rcpar
      real*4              elev4,elvdot4,azim4,azmdot4,nadang4
      real*4              atmdel4
      character*16        cfname
      character*80        prog_name


      iflag = 5

c     convert i*4 variables to i*2

      iprn2    = iprn
      ierfl2   = ierfl
      ndats2   = ndats
      nspare2  = nspare
      nparts2  = nparts

      do 10 i = 1, ndats
         isnr2(i) = isnr(i)
   10 continue

c     convert real*8 to real*4
      elev4   = elev
      elvdot4 = elvdot
      azim4   = azim
      azmdot4 = azmdot
      nadang4 = nadang
      atmdel4 = atmdel


      write(unit    =  lunit
     .,     iostat  =  ioerr
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
c       get calling program name and m-file name for report_stat
        len = rcpar(0,prog_name)
        inquire( unit=lunit, name=cfname, iostat=ioerr )
        call report_stat('FATAL',prog_name,'lib/writc4',cfname
     .                      ,'Error writing C-file data record',ioerr)
      endif

!     write(*,200) iprn, elev*180/pi, azim*180/pi, omcs(1:ndats)
!200  format('WRITC5 ',I3,2F8.2,1x,10(F20.4,1x))
c      print *,'WRITC5: iprn2 ',              iprn2
c      print *,'WRITC5: elev4, elvdot4 ',     elev4,elvdot4
c      print *,'WRITC5: azim4, azmdot4 ',     azim4,azmdot4
c      print *,'WRITC5: nadang4         '     nadang4
c      print *,'WRITC5: atmdel4 ',            atmdel4 
c      print *,'WRITC5: svcepc,svcl1',        svcepc,svcl1
c      print *,'WRITC5: tau,drate ',          tau,drate
c      print *,'WRITC5: ierfl2,data_flag ',   ierfl2,data_flag
c      print *,'WRITC5: ndats2 ',             ndats2
c      print *,'WRITC5: obsv ',               (obsv  (i),i=1,ndats)
c      print *,'WRITC5: omcs ',               (omcs  (i),i=1,ndats)
c      print *,'WRITC5: isnr2 ',              (isnr2  (i),i=1,ndats)
c      print *,'WRITC5  ampl1, ampl2 ',       ampl1,ampl2
c      print *,'WRITC5: nspare2 ',            nspare2
c      print *,'WRITC5: spare ',              (spare(i),i=1,nspare2)
c      print *,'WRITC5: nparts2 ',            nparts2
c      print *,'WRITC5: tmpart ',             (tmpart(i),i=1,nparts)
      return
      end
