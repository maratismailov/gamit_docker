CTITLE PHSMOD_to_CODE
      subroutine phsmod_to_code( ant_mod, ant_cod, dir)
      
      implicit none

*     Rouitne to encode from 8 characters to 4 characters 
*     known antenna model.  If not known, then simply copied.

* ant_mod - full model name
* ant_cod - CODE name
* MOD TAH 960926: Changed the ant model to be encoded for the
*           igs type models.
* dir     - Direction to go:  ENCODE for mod -> cod
*                             DECODE for cod -> mod


      character*(*) ant_mod, ant_cod, dir      
      
*     See which way we are going

      if( dir(1:1).eq.'E' ) then
C         call casefold( ant_mod )
C         if( index(ant_mod,'IGS_01').gt.0 ) then
C             ant_cod = 'IGS1'
C         else if ( ant_mod.eq.'IGS_TEST05' .or. 
C    .              ant_mod.eq.'IGS_T05'.or. 
C    .              ant_mod(1:10).eq.'IGS05_1402')  then
C             ant_cod = 'ELI5'
C         else if ( ant_mod(1:10).eq.'IGS05_1451') then
C             ant_cod = 'ELI2'
C         else if ( ant_mod(1:5).eq.'IGS05' ) then
C             ant_cod = 'ELI5'
C         else
C             ant_cod = ant_mod
C         end if
C
* MOD TAH 101015: Code not needed with extended phase center name:
*         Simply copy model name
          ant_cod = ant_mod         
      else

*         Convert code back to model (The GS_0 line covers an old bug).
          if( ant_cod(1:4).eq.'IGS1' .or. ant_cod(1:4).eq.'GS_0' ) then
              ant_mod = 'IGS_01'
* MOD TAH 970724: Added more name conversions to IGS_01
          else if( ant_cod(1:4).eq.'ELIG' .or. 
     .             ant_cod(1:4).eq.'ELIx' ) then
              ant_mod = 'IGS_01'
          else if( ant_cod(1:4).eq.'ELI5' .or.
     .             ant_cod(1:4).eq.'ELI2' .or.
     .             ant_cod(1:4).eq.'AZI2'  )then
              ant_mod = 'IGS05'
          elseif( ant_cod(1:5).eq.'ELI2' .or.
     .            ant_cod(1:4).eq.'AZI2'  ) then
              ant_mod = 'IGS05'
          else
              ant_mod = ant_cod
          end if
      end if
      
*     Thats all
      return 
      end

CTITLE SMOD_TO_FULL

      subroutine SMod_to_Full( AntMod, name, type)

      implicit none

*     Routine to convert 4-character code from gamit into a
*     model nanme and type.

      character*(*) AntMod  ! GAMIT Ant mode e.g., ELI5
     .,             name    ! Name of model.  This will need to be
                            ! kept unuique
     .,             type    ! Sinex definition of type 
*                             (E Elev or F full)

***** Extract the information.  As models are used this code
*     will need to be updated.
      type(3:4) = ' '
*     Get the model name
      if( AntMod(3:4).eq.'I5' .or. 
     .    AntMod(3:4).eq.'I7'     ) then
          name = 'IGS05' 
          type(3:3) = 'A'    ! This model is absolute
      else  if( AntMod(3:4).eq.'I2' ) then
          name = 'IGS05' 
          type(3:3) = 'A'    ! This model is absolute
      elseif (AntMod(1:4).eq.'NONE' ) then
          name = 'NONE'
          type(3:3) = 'R'
      elseif( antmod(1:3).eq.'IGS' .and. antmod(6:6).eq.'_' ) then
          name = antmod
          type(3:4) = 'AE'
      else
          name = 'UNKNOWN'
          type(3:3) = '?'     ! leave as ? since we dont know
      endif

*     See if Elevation or full
      if( AntMod(1:2).eq.'EL' ) then
           type(4:4) = 'E'
      elseif( AntMod(1:2).eq.'AZ' ) then
           type(4:4) = 'F'
      elseif( AntMod(1:2).eq.'NO' ) then
           type(4:4) = '-'
      else
           type(4:4) = '?'
      endif

*

***** Thats all
      return
      end

CTITLE GET_SVNUM

      subroutine Get_SVnum( epoch )

      implicit none

*      Routine to read svnav.dat and get the SV number for each
*      PRN in the solution.

      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include 'htoglb_comm.h'

* PASSED 
      real*8 epoch   ! Start epoch of the solution.  All PRNs
*                      before this epoch will be used.

* LOCAL
      integer*4 ie   ! Counter to find directoty name
     .,      trimlen ! Length string
     .,      ierr, jerr    ! IOSTAT error
     .,      pn, sv  ! PRN and SV number read from file
     .,      date(5) ! Date of start of interval
     .,      endd(5) ! Date end of interval (for igs_snx YDS)
     .,      i       ! Loop counters
     .,      indx    ! Position of Version string

      real*4 svnver  ! svnav.dat version

      real*8 sectag  ! Second tag
     .,      jde, jds      ! Julian date of satellite entry start and stop

      character*256 svnav_file  ! Name of svnav.dat file.
     .,       line    ! Line read from file

      integer*4 code   ! Block types
      integer*4 conoff       ! Funxtion to return PRN offet from
                             ! constellation type 
      integer*4 off    ! Constellation offset value
      integer*4 sod     ! Second of day from igs_metadata

      character*22 antbody   ! Antenna/Satellite Body type
      character*1 sys  ! System type G == GPS
      character*4 svname       ! sys + satellite number as string.
      character*3 prname       ! sys + PRN number

      logical igs_meta  ! Set true if file linked to svnav.dat is igs_metadata.snx
      logical found     ! Set true when block and valid line found,

****  Get the name of svnav.dat.  Get the directory where the 
*     hfile is located.  Get the last / in the name
      ie = trimlen(hfile)
      do while ( ie.ge.1 .and. hfile(ie:ie).ne.'/' )
         ie = ie - 1
      end do
      if( ie.eq.0 ) then
          svnav_file = 'svnav.dat'
      else
          svnav_file = hfile(1:ie) // 'svnav.dat'
      endif
*     Try to open file
      open(103,file=svnav_file, status='old', iostat=ierr)
      call report_error('IOSTAT',ierr,'open',svnav_file,0,'Get_SVnum')

*     If there is an error opening; try the version on tables
      if( ierr.ne. 0 ) then
         call getenv('HOME',svnav_file)
         svnav_file(trimlen(svnav_file)+1:) = '/gg/tables/svnav.dat'
         open(103,file=svnav_file, status='old', iostat=ierr)
         call report_error('IOSTAT',ierr,'open',svnav_file,0,
     .                     'Get_SVnum')
*        if there is still an error then we can't extract.  Values
*        are initialized to zero so leave it at that.
         if( ierr.ne.0 ) RETURN
      end if

      write(*,120) svnav_file(1:trimlen(svnav_file))
 120  format('Opened svnav.dat file ',a)

*     OK: Mow read the lfile.  Skip first two lines
      read(103,'(a)' ) line
* MOD TAH 190705: See of igs_metadata.snx
      if( line(1:5).eq.'%=SNX' ) then
          igs_meta = .true.
          read(line(6:),*,iostat=ierr) svnver
          call report_error('IOSTAT',ierr,'read',line,0,
     .         'svnav.dat->igs_metadata.snx version')
      else           ! Original svnav.dat
         igs_meta = .false. 

* MOD TAH 050115: See what version of the file.
         svnver = 2.0   ! Default to new version
         indx = index(line,'Version') 
         if( indx.gt.0 ) then
             read(line(indx+9:),*,iostat=ierr) svnver
             call report_error('IOSTAT',ierr,'read',line,0,
     .                        'svnav.dat version line')
         else
            svnver = 1.0
         endif
      endif

* MOD TAH 190705: Drop support for vers 1.0
      if( svnver.lt. 2.0 ) RETURN

*     If this is igs_metadata.snx, skip down to +SATELLITE/IDENTIFIER block
      if( igs_meta ) then
         found = .false.
         do while ( .not.found )
            read(103,'(a)',iostat=ierr) line
            if( line(1:14).eq.'+SATELLITE/PRN' ) then
                found = .true.
            elseif( line(1:14).eq.'-SATELLITE/PRN' .or. 
     .          ierr.ne.0 ) then
                call report_stat('WARNING','HTOGLB','get_ginfo', 
     .             'svnav.dat/igs_metadata.snx',
     .             'Failed to find SATELLITE/PRN block',ierr)
                RETURN
            endif
         enddo      ! File not positioned at start of block
      else

*        Skip header line          
         read(103,'(a)' ) line
      endif

      sectag = 0.d0

      date(2) = 1   ! Set jan since 3 argument is doy.
      endd(2) = 1   ! Set jan since 3 argument is doy.

      do while ( ierr.eq.0 )
          read(103,'(a)',iostat=ierr) line
          if ( trimlen(line).eq.0 ) ierr = -1
          if ( ierr.eq.0 ) then
             if( .not. igs_meta ) then
* MOD TAH 191126: Read start and stop time.
* G    34  18   0  BLOCK IIA               972900.     P      0.1233  2018  24  0  0  2100   1  0  0   #resumed transmitting sometime after 23 January
                read(line,145,iostat=jerr) sys, sv, pn, antbody, 
     .                date(1), date(3), date(4), date(5), 
     .                endd(1), endd(3), endd(4), endd(5) 

 145            format(1x,a1,2x,I4,I4,4x,2x,a22,29x,I4,1x,
     .                 i3,1x,I2,1x,I2,2x,I4,1x,i3,1x,I2,1x,I2 )  
C MOD TAH 190705: This does not seem to be meeded?
C               call name_to_blk(+1, antbody, code)
                found = .true.
             elseif( line(1:1).eq.' ' ) then
                read(line,150, iostat=jerr) svname, date(1),date(3),sod, 
     .               endd(1),endd(3), endd(5), prname
 150            format(1x,a4,2(1x,i4,1x,i3,1x,I5),1x,a3)
                date(4) = sod/3600
                date(5) = (sod-date(4)*3600)/60 
                endd(4) = endd(5)/3600
                endd(5) = (endd(5)-endd(4)*3600)/60 
                if( endd(1).eq.0 ) endd(1) = 2100

                sys = svname(1:1)
                read(svname(2:),*) sv
                read(prname(2:),*) pn
                found = .true.
             else    ! Comment line
                 if( line(1:14).eq.'-SATELLITE/PRN') ierr = -1
                 found = .false.
             end if

* MOD TAH 180401: Update the PRN number for the sys type
             if( found ) then 
                pn = pn + conoff(sys)

!SYS SVN  PRN CHAN ANT/BODY               MASS(G) YAW BIAS  YAW RATE  START           STOP             COMMENTS
! G     1   4   0  BLOCK I                 453800.     U      0.1999  1978  53  0  0  1985 199  0  0                                           
!PRN SV BLK  MASS(G)  BIASED  YAW RATE  YR MO DY HR MN   DX      DY      DZ     (see key at bottom)
! 4,  1  1   453800.     U     0.1999 1978 02 22 00 00   0.210   0.0     0.854 

*               Convert the date to JD
                call ymdhms_to_jd(date,sectag,jds)
                call ymdhms_to_jd(endd,sectag,jde)
* MOD TAH 190705: Skip code test, does not seem to be needed.
C               if( jd.lt. epoch .and. code.gt.0 ) then
C               if( jd.lt. epoch ) then
* MOD TAH 191126: Make sure epoch falls in side range (SVs are turned off and on
*               and can have different PRNs.
                if( epoch.ge.jds .and. epoch.le.jde ) then
                    do i = 1,qnum_svs
                       if( qsvi_prn(i).eq.pn ) then
                           qsvi_svn(i) = sv
                           qsvi_launch(i) = jds
                       endif
                    end do
                end if
             endif       ! Found 
          end if
      end do 

****  Thats all
      close(103)
      return
      end


