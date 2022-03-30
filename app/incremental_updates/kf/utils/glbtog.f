      program glbtog
 
*     This program will read either globk output files or
*     glbak output files and generate GAMIT G-files for all
*     the epheremis elements that it finds. 
*
*     The runstring is:
*     % glbtog <globk/glbak file> [4-5 character code] ['Comment']
*     where <globk/glbak file> is the name of the globk output or
*         glbak file and
*     [4-5 character code] is an optional 4 or 5 character code
*         to be used in the G-file name.  If four characters are given
*         the fifth ill be the last digit of the year.  If no name is
*         given then 'glbkY' will be used where Y is the last digit
*         of the year.
*     ['Comment']  is an optional comment for the g-file

* MOD TAH 200609: Modified with addition of optional GNSS argument.
*     Naming scheme made consistent with GAMIT scheme:
*       <orbt>    (replaces -cent and -pre, which will still work for now)
*                                                         3-char sp3 name 
*        igsf   AC-combined GPS-only final                     -> igs  
*        mitf   MIT  GPS-only final                            -> mit
*        mitr   MIT  GLONASS-only                              -> mir
*        mitm   MIT  Multi-GNSS                                -> mim
*       If no orbt option given, all GNSS output, file name not changed
*       If G  orbt option given, only GPS file name not changed


 
*     If multiple G-files will be produced then an end numerical value
*     will be added to the name.
 
      include 'glbtog.h'
 
* LOCAL VARIABLES
 
*   len_run     - Length of runstring element
*   rcpar       - Routine to read runstring
*   ierr,jerr   - IOSTAT error
*   i,j         - Loop counters
*   indx_date   - Index of date string on line read
*   indx_prn    - Index of prn  string on line read
*   indx_run    - Index of run  string on line read
*   indx_svs    - Index of IC epoch string on line read
*   indx_sys    - Index to GPS System line
*   indx        - Generic index return for read_line
*   n           - Index to radiation parameter
*   code        - Radiation code number
*   pin         - Start of radiation names (e.g. "PRN_02   Direct" points to D 
      integer*4 len_run, rcpar, ierr, jerr, i, indx_date,
     .    indx_prn, indx_run, indx_svs, date(5), indx, n,
     .    indx_sys, code, pin

 
*   sectag      - seconds tage in jd conversions
*   rdum        - Dummy real*8 value for skipping over adjustments
*                 during reads
 
      real*8 sectag, rdum

*   cdum        - Dummy character string for read_line
*   arel        - Arc element name

      character*4 cdum, arel
      
*   gkel        - Globk element name
*   svel        - GAMIT radiation name

      character*14 gkel
      character*16 svel
 
*   line        - Line read from input file
 
      character*256 line

      character*1 ssys  ! GNSS system 

 
****  START: Get the runstring elements
 
      write(*,100)
 100  format(/' glbtog: Make g-files from GLOBK/GLBAK output',/,
     .            ' --------------------------------------------',/)
 
      len_run = rcpar(1, input_file)
      if( len_run.le.0 ) then
          call proper_runstring('glbtog.hlp','glbtog',1)
      else
 
*         Try to open the input file
          open(100,file=input_file, iostat=ierr, status='old')
          call report_error('IOSTAT',ierr,'open',input_file,1,
     .                'glbtog/Open input')
      end if

****  Initialize the time and frame variables
      gtime   = 'UNKN'
      gframe  = 'UNKN'
      gprec   = 'UNKN'
      gsrpmod = 'UNKN'
* MOD TAH 200602: Changed default values
      gnut    = 'IAU0A'
      ggrav   = 'EGR08'
      geradmod = 'NONE'
      gantradmod = 'NONE'
      
****  See if g-file code passed:
      len_run = rcpar(2,code_gfile)
      if( len_run.eq.0 ) then
          code_gfile = 'glbk'
      end if

***** See if comment passed
      len_run = rcpar(3,comment_line)
      if( len_run.eq.0 ) comment_line = ' '

* MOD TAH 200611: See if option to select GNSS passed (to make
*     gfiles compatible with 10.71
      len_run = rcpar(4,sel_gnss)
      if( len_run.eq.0 ) sel_gnss = ' '
      call casefold( sel_gnss ) 
 
**** Now start reading the input file:
 
      num_svs = 0
      prev_gfile = ' '
      gindx = 0
 
      do while ( ierr.eq.0 )
          read(100,'(a)', iostat=ierr ) line
 
*         If no error, see what we have:
          if( ierr.eq.0 ) then
 
*             See if we have found new experiment date:
              indx_date = 0
              indx_date = index(line,'Solution refers to     :')
              if( indx_date.eq.0 ) then 
                  indx_date= index(line,'EXPERIMENT date :')
              end if
 
*             If we found a new experiment date, then output current
*             gfile and decode new date.
              if( indx_date.ne.0 ) then
                  if( num_svs.gt.0 ) then
                      call write_gfile
*                     Reset the number of satellites
                      num_svs = 0
                      gtime   = 'UNKN'
                      gframe  = 'UNKN'
                      gprec   = 'UNKN'
                      gsrpmod = 'UNKN'
                  end if
 
*                 Now decode current date
                  indx_date = index(line,':')
                  read(line(indx_date:),200) date
 200              format(2x,I4,1x,I2,1x,I2,1x,i2,1x,i2,1x,f5.2)
                  sectag = 0.d0
                  call ymdhms_to_jd(date, sectag, exp_jd)
              end if
 
*             See if we found runtime
              indx_run = index(line,'Run time')
              if( indx_run.gt.0 ) then
                  indx_run = index(line,':')
                  read(line(indx_run:),200) date
                  sectag = 0.d0
                  call ymdhms_to_jd(date, sectag, run_jd)
              end if
 
*             See if we found satellite IC
              indx_svs = index(line,'Satellite IC epoch')
              if( indx_svs .gt. 0 ) then
                  indx_svs = index(line,':')
                  read(line(indx_svs:),200) date, sectag
                  call ymdhms_to_jd(date, sectag, svs_jd)
              end if
              
*             See if we GPS System information found
              indx_sys = index(line,'GPS System Information')
              if( indx_sys .gt. 0 ) then
                  indx_sys = index(line,':')
* MOD TAH 080103: New entryies for nutation and gravity
* MOD TAH 140327: Added ggeradmod and ggantradmod to read      
                  read(line(indx_sys:),220) gtime, gframe, 
     .                   gprec, gsrpmod, gnut, ggrav,
     .                   geradmod, gantradmod
 220              format(2x,5x,a4,7x,a5,
     .                   12x,a5,17x,a5,10x, a5,9x,a5,
     .                   10x,a5, 11x, a5 )
c220              format(2x,'Time ',a4,' Frame ',a5,
c    .                   ' Precession ',a5,
c    .                   ' Radiation model ',a5)
c  165     format(' GPS System Information : Time ',a4,' Frame ',a5,
c     .           ' Precession ',a5,' Radiation model ',a5,
c     .           ' Nutation ',a5,' Gravity ',a5,
c     .           ' EarthRad ',a5,' AntThrust ',a5)
               end if
 
*             See if we found a satellite
              indx_prn = index(line,'Inert.  X')
              if( indx_prn.gt.0 ) then
 
*                 OK found start of satellite, Get this elment
*                 and the rest.
** WARNING**  Current code assumes all 9 elements are present.
*
*                 Get the PRN Number
                  num_svs = num_svs + 1
                  indx_prn = index(line,'PRN_')
* MOD TAH 180401: Update for GNSS satellite name format.
                  if( indx_prn.gt.0 ) then 
                      new_gnss = .false.
* MOD TAH 180608: Update to read SVN as well
                      read(line(indx_prn:),250) prns(num_svs),
     .                       svns(num_svs)
* MOD TAH 150825: Read all 4-digits; later check that SV matches PN
 250                  format(4x,2I2)
                      prns(num_svs) = prns(num_svs)*100+svns(num_svs)
                      gnss_sys(num_svs) = 'G'
c                 read(line,260) orb_svs(1,num_svs)
c260              format(39x,F20.4)
                  else
*                     New format e.g., R796_01, E210_01, G063_01
                      new_gnss = .true.
                      read(line(8:),'(a,I3,1x,I2)') ssys, 
     .                     svns(num_svs), prns(num_svs)
                      gnss_sys(num_svs) = ssys
                  endif

                  indx = 40
                  call read_line(line,indx,'R8',jerr,orb_svs(1,num_svs),
     .                           cdum)
                  call read_line(line,indx,'R8',jerr,rdum, cdum)
                  call read_line(line,indx,'R8',jerr,orb_sig(1,num_svs),
     .                           cdum)
 
*                 Read the 5 remaing position and velocity values
                  do i = 1,5
                      read(100,'(a)') line
                      indx = 40
                      call read_line(line,indx,'R8',jerr,
     .                     orb_svs(i+1,num_svs), cdum)
                      call read_line(line,indx,'R8',jerr,rdum, cdum)
                      call read_line(line,indx,'R8',jerr,
     .                     orb_sig(i+1,num_svs), cdum)

*                     See if error on this read to handle case when
*                     apriori value of radiation parameters is zero
*                     on only adjustment is output.  In this case
*                     rdum contains the sigma
                      if( jerr.eq.-1 ) then
                          orb_sig(i+1,num_svs) = rdum
                      end if
     
                  end do
                  
*                 Now get the radiation model
                  num_svs_elem = 0
                 
                  do while (index(line,'Semimajor axis').eq.0 )
                      read(100,'(a)') line
* MOD TAH 150615: Replaced 16:28 with ns:ns+12 based on PRN location.
                      pin = index(line,'PRN') + 9
* MOD TAH 180401: Check for new format
                      if( pin.eq.9 ) then   ! PRN not found
                          pin = 17
                      endif
                       
                      gkel = line(pin:pin+12)
                      call svel_to_code( svel, gkel, arel, code,
     .                                   'GTOC')
* MOD TAH 981103: Only use the radiation parameters.  Ingore the last
*                 three entries which are antenna offsets.
                      if( code.ge.7 .and. code.le.max_svs_elem-3 ) then
                          num_svs_elem = num_svs_elem + 1
                          n = num_svs_elem + 6
                          arels(n) = arel                    
                          indx = 40
                          call read_line(line,indx,'R8',jerr,
     .                         orb_svs(n,num_svs), cdum)
                          call read_line(line,indx,'R8',jerr,rdum,cdum)
                          call read_line(line,indx,'R8',jerr,
     .                         orb_sig(n,num_svs), cdum)

*                         See if error on this read to handle case when
*                         apriori value of radiation parameters is zero
*                         on only adjustment is output.  In this case
*                         rdum contains the sigma
                          if( jerr.eq.-1 ) then
                              orb_sig(n,num_svs) = rdum
                          end if
                      end if
                  end do
                  
******            Convert the units
                  do i = 1,3
                      orb_svs(i,num_svs) = orb_svs(i,num_svs)/1000.d0
                  end do
                  do i = 4,6
                      orb_svs(i,num_svs) = orb_svs(i,num_svs)/1000.d3
                  end do

* MOD TAH 201228: Check sigma on position.  If larger than 1 meter do 
*                 not output.  (Happens with back solution and deleted
*                 satellite in GLX file by use apr_svs <SV> 0 0 0 0 0 0 0R 0A
                  if( orb_sig(1,num_svs)+orb_sig(2,num_svs)+
     .                orb_sig(3,num_svs).ge.1.0 ) then
                      write(*,310) gnss_sys(num_svs),svns(num_svs),
     .                   prns(num_svs),orb_sig(1:3,num_svs) 
 310                  format('DELETING ',a1,I3.3,'_',I2.2,' due to ', 
     .                       'large XYZ sigma ',3F10.3,' m')
                       num_svs = num_svs - 1 
                  endif

              end if
          end if
      end do
 
***** See if we have one last g-file to output
      if( num_svs.gt.0 ) then
          call write_gfile
      end if
 
****  Thats all
      end
 
CTITLE WRITE_GFILE
 
      subroutine write_gfile
 
*     Routine to write out the g-file.  If first creates the name
*     the file, compares it the previous one written and modifies
*     it if it is the same.  It then writes the new file.
 
      include 'glbtog.h'
 
* LOCAL VARIABLES
 
*   i,j     - Loop counters
*   trimlen - Length of string
*   date(5) - Generic date string
*   doy     - Day-of-year
*   ierr    - IOSTAT error
*   date_exp(5) - Experiment date (ymdhm)
*   date_svs(5) - IC date (ymdhm)
*   doy_svs     = Day-of_year for orbits
 
      integer*4 i,j, trimlen, date(5), doy, ierr, date_exp(5),
     .    date_svs(5), doy_svs
 
*   sectag  - Generic seconds tag
*   secexp  - experiment seconds tag
*   secsvs  - IC seconds tag
 
 
      real*8 sectag, secexp, secsvs

      integer*4 pnl, svl  ! PRN and SVN number split from 4-digit PRN SVN
                          ! combination.
      integer*4 num_out   ! Number of satelllites output.  Used to check to
                          ! see if gfile should be removed due to no
                          ! satellites.  TAH 200611.

      logical full   ! Set true when svnav.dat ver 2.0 is found and 
                     ! full satellite information is output
      logical cont   ! Set true to output PRN if SVN matches expected value

      logical out_sv  ! Set true if multi-gnss satellite to be output
  
***** First convert the date information back to the forms
*     needed.
*     The experiment date is used to get the gfile day number.
*     To get the day number generate the jd on the first of the year
*     for the experiment date that we have:
 
      call jd_to_ymdhms(exp_jd, date_exp, secexp)
 
*     Now get jd at Jan 1 of year.
      call ymd_to_doy( date_exp, doy )
 
*     If needed get the last digit of year
      yr_last = mod(date_exp(1),10)
 
*     Generate the g-file name
      full_code = code_gfile
* MOD TAH 200611: See if we are change name based on gnss
      if( trimlen(code_gfile).eq.4 ) then
          write(full_code(5:5),'(I1)') yr_last
      end if
* MOD TAH 200611: See if we are change name based on gnss
      if ( sel_gnss.eq.'R' ) full_code(4:4) = 'r'
      if ( sel_gnss.eq.'E' ) full_code(4:4) = 'e'
      if ( sel_gnss.eq.'C' ) full_code(4:4) = 'c'
      if ( sel_gnss.eq.'M' ) full_code(4:4) = 'm'

 
*     Complete the gfile name
      curr_gfile = 'g' // full_code(1:5)
      write(curr_gfile(7:),200) doy
 200  format('.',I3.3)
 
*     Now see if we have used this name before.
      if( curr_gfile.eq.prev_gfile ) then
 
*         We have already generated a file of this name, so increment
*         index number and append to name.  (Assume here that no more
*         100 gfiles of the same name)
          gindx = gindx + 1
          write(curr_gfile(11:),210) gindx
 210  format('.',i2.2)
      else
 
*         Reset the index number and previous name
          gindx = 0
          prev_gfile = curr_gfile
      end if
 
****  OPen the gfile
      open(200,file=curr_gfile, iostat=ierr,status='unknown')
      call report_error('IOSTAT',ierr,'open',curr_gfile,0,'write_gfile')
      if( ierr.eq.0 ) then
 
*         Start write gfile
          call jd_to_ymdhms(svs_jd, date_svs, secsvs)
          call ymd_to_doy(date_svs,doy_svs)
c**          write(200,300) date_svs(1)-1900, doy_svs, date_svs(4),   
          write(200,300) mod(date_svs(1),100), doy_svs, date_svs(4),
     .                   date_svs(5), nint(secsvs), 
     .                   gtime, gframe, gprec, gsrpmod, gnut, ggrav,
     .                   geradmod, gantradmod

 300      format(I2,1x,i3,1x,i2,1x,i2,1x,i2,20x,a4,1x,a5,1x,
     .           a5,1x,a5,1x,a5,1x,a5,1x,a5,1x,a5)
 
*         Write out elements (agiain we have set for 9 elements')
          write(200,310) num_svs_elem+6, (arels(i),i=7,num_svs_elem+6)
 310      format(i2,' X    Y    Z    XDOT YDOT ZDOT',13a5)
 
*         Now write out the comment line:
          call jd_to_ymdhms(run_jd, date, sectag)
          write(200,320) input_file(1:trimlen(input_file)), date,sectag
 320      format('G-file generated from analysis ',a,/,
     .        'Run-time of analysis: ',i4,2('/',i2),1x,i2,':',i2,1x,
     .          F5.2)

*         If comment has been passed then write out now
          if( trimlen(comment_line).gt.0 ) then
              write(200,'(a)') comment_line(1:trimlen(comment_line))
          end if
 
*         Write the end message:
          write(200,340)
 340      format('END')
 
*         Now start writing out satellites:  Tell user file name
*         being written:
* MOD TAH 200602: Switch svn and PRN order to match orginal format
          if( new_gnss ) then
              write(*,400) trim(curr_gfile), num_svs, 
     .              (gnss_sys(i), svns(i), prns(i),  i = 1,num_svs)
 400          format('Writing ',A,' with ',i3,' satellites.',
     .            100(1x,a1,I3.3,'_',I2.2))
          else 
              write(*,410) trim(curr_gfile), num_svs, 
     .               (prns(i), i= 1,num_svs)
 410           format('Writing ',A,' with ',i3,' satellites. PRN',
     .            100I5.4)
          endif
*         Start the writing:  WARN the user if indirect pressures
*         greater than 0.3)
          call Get_ginfo( svs_jd , full )
          num_out = 0
          do i = 1,num_svs
* MOD TAH 180401: See if new gnss formated names for satellites
              if( .not.new_gnss ) then 
                 pnl = prns(i)/100
                 svl = prns(i)-100*pnl
                 svns(i) = svl
* MOD TAH 160815: See if no SV number (old version of .org file).
                 if( pnl.eq.0 ) then
                     pnl = svl
                     svl = 0
                 endif
              else
* MOD TAH 180401: Values already set so just save
                 pnl = prns(i)
                 svl = svns(i)
              endif

              cont = .false.
* MOD TAH 200611: Only write satellite if correct GNSS selection
              out_sv = .false.
              if(sel_gnss.eq.'G'.and.gnss_sys(i).eq.'G') out_sv = .true. 
              if(sel_gnss.eq.'R'.and.gnss_sys(i).eq.'R') out_sv = .true. 
              if(sel_gnss.eq.'E'.and.gnss_sys(i).eq.'E') out_sv = .true. 
              if(sel_gnss.eq.'C'.and.gnss_sys(i).eq.'C') out_sv = .true. 
              if(sel_gnss.eq.'M'.or. sel_gnss.eq.' ' ) out_sv = .true.
              if( out_sv ) then
                 num_out = num_out + 1  ! Increment so we know how many SVs out.
                 if( .not.full ) then 
                     write(200,500) pnl
 500                 format('PRN ',i2)
                     cont = .true.
                 elseif( svns(i).eq.svl .and. svl.gt.0 ) then   ! Write new line type 
* MOD TAH 180401: Used the have GNSS type.
                     write(200,505) gnss_sys(i), pnl, svns(i), btype(i) 
 505                 format(a1,I2.2,1x,I4,1x,a)
!G 2   61 IIR-B  
                     cont = .true.
                 elseif( svl.eq.0 ) then  ! Old .org file with no SV number
                     write(200,505) gnss_sys(i), pnl, svns(i), btype(i) 
                     cont = .true.
                 else
                     write(*,510) i, prns(i), svns(i), svl
 510                 format('Mis-match SV number. Satellite ',i2,
     .                      ' Full PRN ',i4,' SV expected ',I2.2, 
     .                      ' Found ',i2.2)
                 endif 
                 if( cont ) then 
                    do j = 1,num_svs_elem+6
                        write(200,520) orb_svs(j,i), orb_sig(j,i)
 520                    format(1x,D19.13,5x,F10.4)
                    end do

*                   Look at direct pressure deviation from unity.
                    if( abs(orb_svs(7,i)-1.d0).gt. 0.3d0 ) then
                        write(*,540) prns(i),arels(7)
                    end if 
 
*                   Check the indirect values
                    do j = 8,num_svs_elem+6
                        if( abs(orb_svs(j,i)).gt.0.3d0) then
                            write(*,540) prns(i),arels(j)
 540                    format(' ** WARNING ** For PRN/SVN ',i4,
     .                    ' Radiation Pressure ',a4,' deviates > 0.3')
                        end if
                    end do
                 endif
              endif       ! Correct GNSS or multi-gnss
          end do
 
****      Close the file
          write(200,600)
 600      format('END')
* MOD TAH 200611: If no satellites written delete file and reset
*         previous name so new name not generated
          if( num_out .gt.0 ) then
             close(200)
          else
             write(*,620) trim(curr_gfile), sel_gnss
 620         format('Deleting ',a,' due no ',a,' satellites')
             close(200,status='delete')
             prev_gfile = ' '   ! so new name not generated
          endif
      end if
 
****  Thats all
      return
      end
 
CTITLE GET_GINFO

      subroutine Get_ginfo( epoch , full )

      implicit none

*     Routine to read svnav.dat and information about
*     satellites: e.g., G 3   68 IIF   

      include 'glbtog.h'

* PASSED 
      real*8 epoch   ! Start epoch of the solution.  All PRNs
*                      before this epoch will be used.
      logical full   ! Returns true if all information set (svnav.dat V 2.0)

* LOCAL
      integer*4 ie   ! Counter to find directoty name
     .,      trimlen ! Length string
     .,      ierr, jerr    ! IOSTAT error
     .,      pn, sv  ! PRN and SV number read from file
     .,      pnl, svl    ! Extracted PRN number and SV number
     .,      date(5) ! Date of start of interval
     .,      endd(5) ! Date end of interval (for igs_snx YDS)
     .,      i       ! Loop counters
     .,      indx    ! Position of Version string

      real*4 svnver  ! svnav.dat version

      real*8 sectag  ! Second tag
     .,      jds, jde     ! Julian date of satellite entry

      character*256 svnav_file  ! Name of svnav.dat file.
     .,       line    ! Line read from file

      integer*4 code   ! Block types
      character*22 antbody   ! Antenna/Satellite Body type
      character*22 antbody_in  ! Name read from igs_metdata (needs to be
                               ! changed to match GAMIT names.
* MOD TAH 200602:  Added gamit_bodytyp to save the trailing string used
*     by GAMIT.  Old code only worked for GPS.
      character*8 gamit_bodytyp   ! Body type for GAMIT
      character*4 svname       ! sys + satellite number as string.

      character*1 sys  ! System type G == GPS
      character*64 prog_name
      character*10 cospar_id   ! ID not used in GAMIT
      character*6  SatCat      ! Catalog number; not used in GAMIT.

      logical igs_meta  ! Set true if file linked to svnav.dat is igs_metadata.snx
      logical found     ! Set true when block and valid line found,


****  Get the name of svnav.dat.  Try local version of 
*     file name
      prog_name = 'GLBTOG' 
      svnav_file = 'svnav.dat'
      full = .false.

*     Try to open file
      open(103,file=svnav_file, status='old', iostat=ierr)
!     call report_error('IOSTAT',ierr,'open',svnav_file,0,'Get_SVnum')

*     If there is an error opening; try the version on tables
      if( ierr.ne. 0 ) then

*        Try in ../tables/
         svnav_file = '../tables/svnav.dat'
         open(103,file=svnav_file, status='old', iostat=ierr)
         if( ierr.ne.0 ) then 
*           Try gg location
            call getenv('HOME',svnav_file)
            svnav_file(trimlen(svnav_file)+1:) = '/gg/tables/svnav.dat'
            open(103,file=svnav_file, status='old', iostat=ierr)
            call report_error('IOSTAT',ierr,'open',svnav_file,0,
     .                     'Get_SVnum')
         endif
*        if there is still an error then we can't extract.  Values
*        are initialized to zero so leave it at that.
         if( ierr.ne.0 ) RETURN
      end if

!     write(*,120) svnav_file(1:trimlen(svnav_file))
!120  format('Opened svnav.dat file ',a)

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

      if( svnver.lt. 2.0 ) RETURN
      full = .true.

*     If this is igs_metadata.snx, skip down to +SATELLITE/IDENTIFIER block
      if( igs_meta ) then
         found = .false.
         do while ( .not.found )
            read(103,'(a)',iostat=ierr) line
            if( line(1:13).eq.'+SATELLITE/ID' ) then
                found = .true.
            elseif( line(1:13).eq.'-SATELLITE/ID' .or. 
     .          ierr.ne.0 ) then
                call report_stat('WARNING',prog_name,'get_ginfo', 
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
             if( .not.igs_meta ) then 
* MOD TAH 191126: Read start and stop time.
* G    34  18   0  BLOCK IIA               972900.     P      0.1233  2018  24  0  0  2100   1  0  0   #resumed transmitting sometime after 23 January
                read(line,145,iostat=jerr) sys, sv, pn, antbody, 
     .                date(1), date(3), date(4), date(5), 
     .                endd(1), endd(3), endd(4), endd(5) 

 145            format(1x,a1,2x,I4,I4,4x,2x,a22,29x,I4,1x,
     .                 i3,1x,I2,1x,I2,2x,I4,1x,i3,1x,I2,1x,I2 )  
                call ymdhms_to_jd(date,sectag,jds)
                call ymdhms_to_jd(endd,sectag,jde)
                call name_to_blk(+1, antbody, code)
* MOD TAH 200625: Save the gamit_bodytyp based on system type.
                gamit_bodytyp = '??'
                if( sys.eq.'G' ) gamit_bodytyp = antbody(7:)
                if( sys.eq.'R' ) gamit_bodytyp = antbody(9:)
                if( sys.eq.'E' ) gamit_bodytyp = antbody(9:)
                if( sys.eq.'C' ) gamit_bodytyp = antbody(8:)

                found = .true.
             elseif( line(1:1).eq.' ' ) then   ! Decode line
                read(line ,150, iostat=ierr) svname, cospar_id, SatCat, 
     .                   antbody_in
 150            format(1x,a4,1x,a9,1x,A6,1x,a15)
                if( ierr.ne.0 )  
     .              call report_stat('FATAL',prog_name,'get_ginfo',
     .             'lib/read_svsinex','svnav.dat',
     .              'Error decoding ID block of SINEX file',ierr)
!*SVN_ COSPAR ID SatCat Block__________ Comment__________________________________
!*                                                                               
! G001 1978-020A  10684 GPS-I           Launched 1978-02-22; NAVSTAR 1
*                Due to sub_char allowing multiple replacements
*                we can't directly replace'GLO' with 'GLONASS'
*                (sub_char knows and won't replace, so we need
*                to do in two steps
* MOD TAH 200602: Extract the shortened body type (i.e., IIIA,M,K1A etc
*                needed by arc g-files.
                 gamit_bodytyp = antbody_in(5:)
                 call sub_char( antbody_in,'GPS-','BLOCK_')
                 call sub_char( antbody_in,'_',' ')
                 call sub_char( antbody_in,'GLO','XXX')
                 call sub_char( antbody_in,'XXX','GLONASS')     
                 call sub_char( antbody_in,'GAL','XXX')
                 call sub_char( antbody_in,'XXX','GALILEO')
                 call sub_char( antbody_in,'BDS','BEIDOU')
                 antbody = antbody_in
                 call name_to_blk(+1, antbody, code)
                 jds = 2415020.5   ! SVnumber to antbody is unique so no date (1900)
                 jde = 2488069.5   ! SVnumber to antbody is unique so no date (2100)
                 sys = svname(1:1)
                 read(svname(2:),*) sv 
                 found = .true.
             else
                if( line(1:13).eq.'-SATELLITE/ID' ) ierr = -1
                found = .false.
             end if

!SYS SVN  PRN CHAN ANT/BODY               MASS(G) YAW BIAS  YAW RATE  START           STOP             COMMENTS
! G     1   4   0  BLOCK I                 453800.     U      0.1999  1978  53  0  0  1985 199  0  0                                           
!PRN SV BLK  MASS(G)  BIASED  YAW RATE  YR MO DY HR MN   DX      DY      DZ     (see key at bottom)
! 4,  1  1   453800.     U     0.1999 1978 02 22 00 00   0.210   0.0     0.854 
*            Convert the date to JD
             if( found ) then 
* MOD TAH 191126: Make sure epoch falls in side range (SVs are turned off and on
*               and can have different PRNs. (For SVNAV.dat)
C               if( jd.lt. epoch .and. code.gt.0 ) then
                if( epoch.ge.jds .and. epoch.le.jde 
     .                           .and. code.gt.0 ) then
                    do i = 1,num_svs
* MOD TAH 150825: Get high order 2-digits of PRN
* MOD TAH 160815: Account for old .org with just PRN number (no SV number)
* MOD TAH 180401: See if new GNSS satellite names being used
* MOD TAH 180609: .not.new_gnss only works with GPS
                       if( .not.new_gnss .and. sys.eq.'G' ) then 
                          pnl = prns(i)/100
                          svl = prns(i) - pnl*100
* MOD TAH 190705: SVnumber is unique, use only this.
                          if( sv.eq.svl .and. svl.ne.0 ) then
* MOD TAH 200602: Use the antbody_in string where the gamit body type is
*                 in the 5th position (GPS-IIIA -> IIIA; GLO-M+ -> M+ etc)
                              btype(i) = gamit_bodytyp
* These for old solutions and svnav.dat proper is likely needed.
                          elseif ( pnl.eq.pn .and. svl.eq.0 ) then
                              svns(i) = sv   ! Save SV number since this is unknown
                              prns(i) = pnl*100 + sv
                              btype(i) = gamit_bodytyp
                          endif
                       else
* MOD TAH 190705: SVnumber is unique, use only this.
!                         if( pn.eq.prns(i) .and. sv.eq.svns(i) ) then
                          if( sv.eq.svns(i) ) then
                             btype(i) = gamit_bodytyp
                             exit
                          endif
                       end if
                    end do
                end if
              endif            ! Found
          endif
      end do 

****  Thats all
      close(103)
      return
      end
