Copyright (c) Massachusetts Institute of Technology and the University of
California at San Diego, 1994.  All rights reserved.
      subroutine ngsopn( spfile,tfile,tfilef,iungs,iux,iut )
C
C Open the files for NGSTOT and/or TTONGS
C R.W. King    18 Nov 1991
C
C Enter NGS file and T-file names
C
      implicit none
                   
      integer*4 iungs,iux,iut,iper,iflag,gpsweek,gpsday,idoy,iy,iy1
     .         ,len,rcpar,ioerr,nblen

      character*16 spfile,tfile,tfilef
      character*11 tfname
      character*80 tmpname,prog_name
c
      len = rcpar(0,prog_name)
c
c......Open the NGS fromat ephemeris file
      
      open(unit=iungs,file=spfile,status='old',form='formatted',
     .    iostat=ioerr)
      if (ioerr .ne. 0 ) then
          call report_stat('FATAL',prog_name,'orbits/ngsopn',spfile,
     .    'Error opening NGS Ephemeris file: ',ioerr)
      endif
c
c.....Extract GPS week and day for t-file name  

c* rwk 080321 this code doesn't work for ultra-rapids, which have the form iguWWWWD_HH.sp3
c*      iper=index(spfile,'.')
c*      read(spfile(iper-5:iper-2),'(i4)',iostat=ioerr) gpsweek  
c*    this ought to work, assume the form tttWWWW*
      read(spfile(4:7),'(i4)',iostat=ioerr) gpsweek
      if(ioerr .ne. 0) then
        call report_stat('WARNING',prog_name,'orbits/ngsopn',spfile,
     .  'Error reading gpsweek from sp3 file:',ioerr)
      endif
c*   ditto
c*     read(spfile(iper-1:iper-1),'(i1)',iostat=ioerr) gpsday
      read(spfile(8:8),'(i1)',iostat=ioerr) gpsday
      if(ioerr .ne. 0) then
        call report_stat('WARNING',prog_name,'orbits/ngsopn',spfile,
     .  'Error reading gpsday from sp3 file: ',ioerr)
      endif
      idoy=0
      iy=0
      call doygwk(idoy,iy,gpsweek,gpsday)
      tfname='           '
      tfname(1:5)='txxxx'
      iy1 = mod(iy,10)
      write(tfname(6:6),'(i1)') iy1
      tfname(7:7)='.'
      write(tfname(8:10),'(i3)') idoy
      if(tfname(8:8).eq.' ') tfname(8:8)='0'
      if(tfname(9:9).eq.' ') tfname(9:9)='0'
c
c.....Define earth-fixed T-file name

      tfilef=tfile
      iper=index(tfile,'.')

c.....Correct format (need a period in name) ?
      if(iper.eq.0) call report_stat('FATAL','NGSTOT','orbits/ngsopn'
     .             ,tfile,'No period in T-file name',ioerr)
      tfilef(iper-1:iper-1)='e'
c
c.....Open inertial T-file
      open(unit=iut,file=tfile,status='new',form='unformatted'
     .     ,iostat=ioerr)
      if (ioerr .ne. 0) then
        open(unit=iut,file=tfile,status='unknown',form='unformatted'
     .           ,iostat=ioerr)
        if (ioerr .ne. 0) then
          call report_stat('FATAL',prog_name,'orbits/ngsopn',tfile,
     .             'Error opening inertial T-file: ',ioerr)
        endif
      endif
      close(unit=iut)
c
c.....Open earth-fixed T-file
      open(unit=iut,file=tfilef,status='unknown',form='unformatted'
     .    ,iostat=ioerr)
        if (ioerr .ne. 0) then
          call report_stat('FATAL',prog_name,'orbits/ngsopn',tfilef,
     .    'Error opening earth fixed T-file: ',ioerr)
        endif

      return
      end
