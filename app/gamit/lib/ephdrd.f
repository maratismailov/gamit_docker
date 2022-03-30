Copyright (c) Massachusetts Institute of Technology and the University of
California at San Diego, 1995.  All rights reserved.

      Subroutine EPHDRD(fjd)

c  a) Read lunar and solar ephemeris file headers,
c  b) be sure that their spans cover the times of interest, and
c  c) make initializations for ephemerides reads

c Rick Abbot - November 1984
c R. King - March 2018  Minor mods to make compatible with new ARC code with
c                       variables in common

 
c     Input: fjd used to test boundaries of ephemeris files
c            frame (in arc.h) used to check frame of the ephemerides

c     Output (all in arc.h): 
c       fjdbmn ffjdemn fjdbsn fjdesn - start/stop of ephemeris files
c       ilun isun - logical unit numbers for files 
c       sdeltm sdelts - tabular intervals of the values on the files
c       nintrsm nintrss - number of coordinate components (set to 3)
c       initialization of interpolation pointers

      implicit none
                         
      include '../includes/dimpar.h'  
      include '../includes/units.h'
      include '../includes/global.h' 
      include '../includes/arc.h'


c     local
      real*8 fjd,fjdemn,fjdesn
      character*5 sun_frame,moon_frame

               
c     variable for status reporting
      integer*4 len,rcpar,ioerr
      character*80 prog_name
      character*256 message

cd      print *,'EPHDRD frame ',frame 
             
c     get the calling program name for report_stat
      len = rcpar(0,prog_name)

c      open the files

cd       print *,'EPHDRD units ',iterm,iscrn,iprnt,ibody,isun,ilun
cd     .                       ,inut,iut1,ipole                  
c.....open lunar ephemeris file
      ilun = 33
      open(unit=ilun,file='luntab.',status='old',iostat=ioerr)
      if (ioerr .ne. 0) then
        call report_stat('FATAL',prog_name,'lib/ephdrd','luntab.',
     .  'Error opening luntab. file',ioerr) 
      else
        call report_stat('STATUS',prog_name,'lib/ephdrd','luntab.',
     .  'Open lunar ephemeris file',ioerr)
      endif
c
c.....open solar ephemeris file
      isun = 34
      open(unit=isun,file='soltab.',status='old',iostat=ioerr)
      if (ioerr .ne. 0) then
        call report_stat('FATAL',prog_name,'lib/ephdrd','soltab.',
     .  'Error opening soltab. file',ioerr)
      else
        call report_stat('STATUS',prog_name,'lib/ephdrd','soltab.',
     .  'Open solar ephemeris file',ioerr)
      endif


c.....read sun and moon ephemeris frames
      moon_frame = ' '            
      sun_frame = ' '  
c     skip the first header record and read the dates and frame
      read(ilun,'(1x)')
      read(ilun,'(31x,f7.0,1x,f7.0,23x,a5)') fjdbmn,fjdemn,moon_frame 
      read(isun,'(1x)') 
      read(isun,'(31x,f7.0,1x,f7.0,23x,a5)') fjdbsn,fjdesn,sun_frame
cd      print *,'EPHDRD frame moon_frame sun_frame ',frame,moon_frame
      if( moon_frame.ne.'J2000' ) moon_frame = 'B1950'
      if( sun_frame.ne.'J2000' )  sun_frame  = 'B1950'
      if( moon_frame.ne.sun_frame ) then
        write(message,'(a,2a)')
     .       'Frames of soltab. and luntab. are different: '
     .       ,sun_frame,moon_frame
        call report_stat('FATAL',prog_name,'lib/ephdrd',' ',message,0)
      elseif( sun_frame.ne.frame ) then
         write(message,'(a,1x,a5,1x,a5)')
     .      'Input frame and lunar/solar frame are different:'
     .      , frame,sun_frame
        call report_stat('FATAL',prog_name,'lib/ephdrd',' ',message,0)
      endif                           

      sdeltm = 0.5d0
      sdelts = 4.0d0
      nintrsm = 3
      nintrss = 3
      iendfm=0
      iendfs=0

c     See if the lunar ephemeris bounds the time of interest
      if( (fjd - (fjdbmn+0.5d0*sdeltm)). lt.0.d0 ) then
        write(iscrn,110) fjdbmn,fjdemn
  110   format(//,' Lunar tabular ephemeris parameters',
     .          /,f9.0,2x,'Ephemeris start jd' ,
     .         /,f9.0,2x,'Ephemeris end jd')
        write(message,111) fjd,fjdbmn+5.d0*sdeltm
  111   format ('Lunar tabular ephemeris starts too late '
     .,'Observation start jd = ',f14.5,' Lunar ephem start jd = ',f14.5)
        call report_stat('FATAL',prog_name,'lib/ephdrd',' ',message,0)
      endif

      if( ((fjdemn-5.d0*sdeltm) - fjd).lt.0.d0 ) then
        write(iscrn,110) fjdbmn,fjdemn
        write(message,112) fjd,fjdemn-5.d0*sdeltm
  112   format ('Lunar tabular ephemeris ends too early. ',
     .  'Observation end jd = ',f14.5,' Lunar ephem end jd = ',f14.5)
        call report_stat('FATAL',prog_name,'lib/ephdrd',' ',message,0)
      endif
               
c     See if the solar ephemeris bounds the time of interest
      if( (fjd - (fjdbsn+0.5d0*sdelts)). lt.0.d0 ) then
        write(iscrn,113) fjdbsn,fjdesn
  113   format(//,' Solar tabular ephemeris parameters',
     .          /,f9.0,2x,'Ephemeris start jd' ,
     .         /,f9.0,2x,'Ephemeris end jd')
        write(message,114) fjd,fjdbsn+5.d0*sdelts
  114   format ('Solar tabular ephemeris starts too late '
     .,'Observation start jd = ',f14.5,' Lunar ephem start jd = ',f14.5)
        call report_stat('FATAL',prog_name,'lib/ephdrd',' ',message,0)
      endif
      if( ((fjdesn-5.d0*sdelts) - fjd).lt.0.d0 ) then
        write(iscrn,115) fjdbsn,fjdesn
        write(message,115) fjd,fjdesn-5.d0*sdelts
  115   format ('Solar tabular ephemeris ends too early. ',
     .  'Observation end jd = ',f14.5,' Solarr ephem end jd = ',f14.5)
        call report_stat('FATAL',prog_name,'lib/ephdrd',' ',message,0)
      endif

      ji0m=0
      jilm=0
      iy1m=0
      iy2m=0
      jlastm=0
      jnowm=0                      
      sdeltm=1.d0/sdeltm      
      ji0s=0
      jils=0
      iy1s=0
      iy2s=0
      jlasts=0
      jnows=0
      sdelts=1.d0/sdelts

      return
      end
