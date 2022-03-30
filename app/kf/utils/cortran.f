      program cortran

      implicit none 

*     General coordinate system transformation program with
*     datum shifts
*
*     Runstring:
*     % cortran <infile> <outfile> <intype> <insys> <outtype> <outsys>
*     where <infile> is a globk style coordinate file with the input 
*           coordiates of type <intype>. Minimum entries are 
*            SiteName 3-coordinates in free format unless the intype is
*           UTM in which case the zone and hemisphere (N/S) must follow
*           the coordinates
*           <outfile> is the output coordinate file of <outtype>
*           <intype> is the type of input coordintes.  The choices are
*              XYZ  -- Cartesian coordinates
*              GEOD -- Geodetic latitude, longitude (+ve East) and Height
*              UTM  -- UTM coordinates with geodtic height as the 3rd enty
*                      (Zone and Hemisphere are included in this output)
*           <insys> is the input system.  The allowable systems are defined
*              in kf/gen_util/datum_def.f.  They currently are
*              WGS84 -- WGS84
*              NAD27 -- North American Datum 1927
*              OMAN  -- Oman geodetic system
*              (Note: the system includes both the ellipsoid and any 
*                     translation with respect to ITRF2000)
*           <outtype> is the output coordinate type
*           <outsys> is the output system.

      include '../includes/const_param.h'

* PROGRAM VARIBALES

      real*8 incrd(3), outcrd(3)  ! Input and output coordinates of site
      
      character*128 infile, outfile ! Names of input and output files
      character*4 intype, outtype ! Input and output types (XYZ, GEOD, UTM)
      character*8 insys,  outsys  ! Input and output systems (WGS84,NAD27,OMAN)
*                                   (See kf/gen_util/datum_def.f)

      integer*4 zone  ! UTM Zone number
      character*1 hemi   ! Hemisphere for UTM (N/S)
      character*8 site_name  ! Name of site (no spaces)

      integer*4 ierr, jerr  ! IOSTAT errors
      integer*4 len_run, len_fi, len_fo, rcpar, indx ! Lengths of runstring entries
      integer*4 trimlen  ! Length of string
      integer*4 ounit    ! Output unit number, allows 6 to be used

      character*256 line  ! Line read from file


***** Read the entries from the runstring and print help if not complete
      len_fi = rcpar(1,infile)
      if( len_fi.eq.0 ) call proper_runstring('cortran.hlp',
     .                    'cortran/infile',1)
      open(50,file=infile,status='old',iostat=ierr)
      call report_error('IOSTAT',ierr,'open',infile,1,'cortran') 
      len_fo = rcpar(2,outfile)
      if( len_fo.eq.0 ) call proper_runstring('cortran.hlp',
     .                    'cortran/outfile',1)
      if( outfile(1:2).eq.'6 ' ) then
          ounit = 6
      else
          ounit = 60
          open(ounit,file=outfile,status='unknown',iostat=ierr)
          call report_error('IOSTAT',ierr,'creat',outfile,1,'cortran')
      endif 

****  Now get the input system characeristics
      len_run = rcpar(3,intype)
      if( len_run.eq.0 ) call proper_runstring('cortran.hlp',
     .                    'cortran/intype',1)
      call casefold(intype)
      len_run = rcpar(4,insys)
      if( len_run.eq.0 ) call proper_runstring('cortran.hlp',
     .                    'cortran/insys',1)
      call casefold(insys)

****  Get the output system characteristics
      len_run = rcpar(5,outtype)
      if( len_run.eq.0 ) call proper_runstring('cortran.hlp',
     .                    'cortran/outtype',1)
      call casefold(outtype)
      len_run = rcpar(6,outsys)
      if( len_run.eq.0 ) call proper_runstring('cortran.hlp',
     .                    'cortran/outsys',1)
      call casefold(outsys)

***** OK: Write some new header records to the output file (other comment
*     line are copied accross)
      write(*,110) infile(1:len_fi), intype, insys,
     .             outfile(1:len_fo), outtype, outsys
 110  format('* CORTRAN: Transforming coordinates in ',a,1x,a,1x,a,/,
     .       '*                                   to ',a,1x,a,1x,a)
      if( ounit.ne.6 )
     .write(ounit,110,iostat=ierr) infile(1:len_fi), intype, insys,
     .             outfile(1:len_fo), outtype, outsys
      call report_error('IOSTAT',ierr,'writ',outfile,1,'cortran')
      if( outtype(1:3).eq.'XYZ' ) then
         write(ounit,120) 
 120     format('*  Site          X (m)          Y (m)          Z (m)')
      else if( outtype(1:4).eq.'GEOD' ) then
         write(ounit,140) 
 140     format('*  Site      Latitude (deg)  Longitude (deg)',
     .          '  Height (m)')
      else if( outtype(1:3).eq.'UTM' ) then
         write(ounit,160) 
 160     format('  Site   Northing (m)    Easting (m)    Height (m)',
     .          ' Zone Hemi')
      else
         call report_stat('FATAL','CORTRAN','main',outtype,
     .           'Unknown output type',0)
      endif

****  Loop over input, writing comments to output and converting the 
*     coordinates
      ierr = 0
      do while( ierr.eq.0 )
          read(50,'(a)',iostat=ierr) line
*         See if EXTENDED entry:
          indx = 1 
          call GetWord(line, site_name,indx)
          if( ierr.eq.0 .and. site_name.ne.'EXTENDED' ) then
*             See if comment
              if ( line(1:1).ne.' ' .or. trimlen(line).eq.0 ) then
                  write(ounit,'(a)') line(1:max(1,trimlen(line)))
              else
*                 We seem to have a none comment line
                  if( intype(1:3).eq.'UTM' ) then
                     read(line,*,iostat=jerr) site_name, incrd, 
     .                                        zone, hemi
                  else
                     read(line,*,iostat=jerr) site_name, incrd
                  endif
                  call report_error('IOSTAT',jerr,'decod',line,1,
     .                              'cortran')
                  if( jerr.eq.0 ) then

****                  If geod input convert from lat long in degrees
*                     to co-lat and long in radians
                      if( intype(1:4).eq.'GEOD' ) then
                          incrd(1) = (90-incrd(1))*pi/180
                          incrd(2) = incrd(2)*pi/180
                      end if

                      call geod_to_geod(incrd, outcrd, intype, outtype,
     .                                  insys, outsys, zone, hemi)

*                     Now write out the results
                      if( outtype.eq.'XYZ' ) then
                         write(ounit,220) site_name, outcrd
 220                     format(1x,a8,3F15.4)
                      else if( outtype.eq.'GEOD' ) then
                         outcrd(1) = (pi/2-outcrd(1))*180/pi
                         outcrd(2) = outcrd(2)*180/pi
                         write(ounit,240) site_name, outcrd
 240                     format(1x,a8,2F15.9,F15.4)
                      else if( outtype.eq.'UTM' ) then
                         write(ounit,260) site_name, outcrd, zone, hemi
 260                     format(1x,a8,3F15.4,1x,i2.2,1x,a1)
                      endif
                  endif
              endif 
          endif
      end do

***** Thats all
      close(50)
      if( ounit.ne.6 ) close(ounit)
      end


                     




