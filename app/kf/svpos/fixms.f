 
      program fixms

      implicit none 

*     Program to fix millisec offsets between CA code range
*     measurements in trimble receiverers.
 
*     The runstring of the program is
*     % fixms [navfile] [data file]
*     where nav_file is nave file name
*         data_file is the name of the Rinex data file.
* MOD TAH 971219: Added feature to fix phase if there is a jump
*     in these values as well.

*     The output rinex file is called [data file].fix
 
      include 'fixms.h'
      include '../../libraries/includes/const_freq.h'

 
* Main pogram variables

      integer*4 i
 
*   eof     - Indicates end of file.
 
      logical eof
 
***** Get the runstring runstring
 
      call get_svsrun
 
*     Read in the ephemeris file
      call read_nav
 
*     Now loop over the data and get an estimate of the site position.
*     This rouitne will also echo the header records of the rinexfile
*     into the output.
 
      call read_data_head
 
      write(*,100) 
 100  format(' FIXMS: Program to fix erroneous CA code ranges') 

      eof = .false.
      i = 0
      
      do while ( .not.eof)
          call read_range(eof)
          i = i + 1 
          if( .not.eof ) call fix_CArange(i)   
 
      end do
 
****  Thats all
      end
 
CTITLE fix_CArange
 
      subroutine fix_CArange(ep)

      implicit none 
 
*     Routine to remove anomalous millisecond offsets in the 
*     CA code ranges.
 
      include 'fixms.h'
      include '../includes/const_param.h'
      include '../../libraries/includes/const_freq.h'    

* sec10 - minutes and seconds in tenths of seconds units.

      integer*4 i,j, num_av, date(5), ep,
     .          ctop(max_sat)
      real*8  average, svclk(max_sat),  sectag

      real*8  clock_epoch

* dms, dms_L1, dms_L2 -- Millisecond offsets during calcualtion
*     and saved for L1 and L2.

      real*8  dms, rms, dms_L1, dms_L2
     
      character*80 outline

****  compute the epheremis position at the measurement
*     time
 
c     write(*,*) 'For epoch ',data_epoch
*     See if debug wanted:
C     if( ep.ge.debug_start .and. ep.le.debug_end ) then
C         write(*,110) 'Satellite Position at t-66.66 us', ep, 0.0
C110      format('+',1x,a,' Epoch ',i5,' Site clock ',f17.2,' (m)',/,
C    .           '+PRN',5x,'Xe (m)',7x,'Ye (m)',7x,'Ze (m)',
C    .           5x,'SV Clock (m)')
C     end if

*     prev_site_clk is initialized in get_svs_run (can be passed by
*     user or set to zero).
      average = prev_site_clk 
      if( ep.eq.1 ) then
         init_dms = data_sec - nint(data_sec)
         write(*,115) init_dms
 115     format('Setting initial ms offset to ',f6.4,' ms')
      end if

      do i = 1, num_sat
 
          call eph_to_xyz( data_epoch-0.066666/86400.d0, i, 'E')
          svclk(i) = (af0(i)+
     .                af1(i)*(data_epoch-toe_jd(i))*86400.d0)*
     .               vel_light

*         See if debug wanted:
C         if( ep.ge.debug_start .and. ep.le.debug_end ) then
C             write(*,120) prn(i), (svs_xyz(j,i),j=1,3), svclk(i)
C120          format('+',i3, 3F13.2, F13.2)
C         end if
 
      end do
 
      do i = 1, num_chan
 
          omc_OK(i) = .false.
          do j = 1, num_sat
              if( chan(i).eq.prn(j) ) then
 
*                 Compute a rough range
                  p1c(i) = sqrt( (site_xyz(1)-svs_xyz(1,j))**2+
     .                           (site_xyz(2)-svs_xyz(2,j))**2+
     .                           (site_xyz(3)-svs_xyz(3,j))**2)
                  ctop(i) = prn(j)
               end if
           end do
      end do 

      do i = 1, num_chan
          do j = 1, num_sat
              if( chan(i).eq.prn(j) ) then
 
*                 Compute the range accounting for the propagation
*                 delay.
                  clock_epoch = data_epoch -
     .                   (average+p1c(i))/vel_light/86400.d0
                  svclk(i) = (af0(i)+
     .                        af1(i)*(clock_epoch -
     .                                toe_jd(i))*86400.d0)*
     .                        vel_light
                  call eph_to_xyz( data_epoch-
     .                      (average+p1c(i))/vel_light/86400.d0,
     .                       j, 'E')
                  p1c(i) = sqrt( (site_xyz(1)-svs_xyz(1,j))**2+
     .                           (site_xyz(2)-svs_xyz(2,j))**2+
     .                           (site_xyz(3)-svs_xyz(3,j))**2)

*                 See if debug wanted:
C                 if( ep.ge.debug_start .and. ep.le.debug_end ) then
C                     write(*,120) prn(j), (svs_xyz(k,j),k=1,3), 
C    .                             svclk(j)
C                 end if
 
****              compute O-C
                  omc_cl1(i) = 0.d0
                  omc_cl2(i) = 0.d0
                  omc_pl1(i) = 0.d0
                  omc_pl2(i) = 0.d0
                  if( c1f(i).ne.99 )
     .            omc_cl1(i) = C1o(i) -p1c(i)+svclk(j)-prev_site_clk 
                  if( c2f(i).ne.99 )
     .            omc_cl2(i) = C2o(i) -p1c(i)+svclk(j)-prev_site_clk 
                  if( p1f(i).ne.99 )
     .            omc_pl1(i) = P1o(i) -p1c(i)+svclk(j)-prev_site_clk 
                  if( p2f(i).ne.99 )
     .            omc_pl2(i) = P2o(i) -p1c(i)+svclk(j)-prev_site_clk 
 
              end if
          end do
      end do

*     See if debug for o-minus-c wanted
      if( ep.ge.debug_start .and. ep.le.debug_end ) then
          write(*,140) 'Range differences', ep
 140      format('+ ',a,' Epoch ',i4,/,
     .           '+CHAN PRN',3x,'Obs (m)',6x,'Range (m)',3x,
     .           'Obs-Calc L1 ',3x,'Obs-Calc L2')
       
          do i = 1, num_chan
              write(*,150) i,ctop(i), c1o(i), p1c(i), omc_cl1(i),
     .                     omc_pl2(i) 
 150          format('+',i3,1x,i4, 4F13.2)
          end do
      end if
      
*    
*     Loop over all the channels and update the clock as we go
      average = 0.d0
      num_av = 0
      rms = 0

*     Now check the consistency of ranges
*     Loop over all the channels and update the clock as we go

      do i = 1, num_chan
      
*        Clear the millisec offset values
         dms_L1 = 0
         dms_L2 = 0      

*        Check each of range values
         if( p1f(i).ne.99 .and. p1o(i).ne.0 ) then
             if( abs(omc_pl1(i)).gt.0.25d-3*vel_light ) then
*                The range seems to be on wrong millisecond.
*                correct the value, and then let it con contribute
*                to estimate of clock
                 dms = nint(omc_pl1(i)/(1.d-3*vel_light))*
     .                                (1.d-3*vel_light)
                 dms_L1 = dms
                 omc_pl1(i) = omc_pl1(i) - dms
                 p1o(i) = p1o(i) - dms
                 write(*,300) ep, i, ctop(i), 'P1 Range',
     .                        dms/(1.d-3*vel_light)
 300             format('Epoch ',i5,' CH ',i2,' PRN ',i2,1x,a,
     .                   1x,F8.1,' ms')
             end if
             average = average + omc_pl1(i)
             rms = rms + omc_pl1(i)**2
             num_av = num_av + 1
         end if
         if( c1f(i).ne.99 .and. c1o(i).ne.0 ) then
             if( abs(omc_cl1(i)).gt.0.25d-3*vel_light ) then
*                The range seems to be on wrong millisecond.
*                correct the value, and then let it con contribute
*                to estimate of clock
                 dms = nint(omc_cl1(i)/(1.d-3*vel_light))*
     .                                (1.d-3*vel_light)
                 dms_L1 = dms
                 omc_cl1(i) = omc_cl1(i) - dms
                 c1o(i) = c1o(i) - dms
                 write(*,300) ep, i, ctop(i), 'C1 Range',
     .                     dms/(1.d-3*vel_light)
             end if
             average = average + omc_cl1(i)
             rms = rms + omc_cl1(i)**2
             num_av = num_av + 1
         end if
         if( p2f(i).ne.99 .and. p2o(i).ne.0 ) then
             if( abs(omc_pl2(i)).gt.0.25d-3*vel_light ) then
*                The range seems to be on wrong millisecond.
*                correct the value, and then let it con contribute
*                to estimate of clock
                 dms = nint(omc_pl2(i)/(1.d-3*vel_light))*
     .                                (1.d-3*vel_light)
                 dms_L2 = dms
                 omc_pl2(i) = omc_pl2(i) - dms
                 p2o(i) = p2o(i) - dms
                 write(*,300) ep, i, ctop(i), 'P2 Range',
     .                        dms/(1.d-3*vel_light)
             end if
             average = average + omc_pl2(i)
             rms = rms + omc_pl2(i)**2
             num_av = num_av + 1
         end if

* MOD TAH 971219: See if we should adjust the phase as well.
         if( fix_phase ) then
             if( L1o(i).ne.0 .and. dms_L1.ne.0 ) then
                 L1o(i) = L1o(i) - dms_L2/vel_light*fL1
                 write(*,300) ep, i,  ctop(i), 'L1 Phase',
     .                        dms_L1/(1.d-3*vel_light)
             end if
             if( L2o(i).ne.0 .and. dms_L2.ne.0 ) then
                 L2o(i) = L2o(i) - dms_L2/vel_light*fL2
                 write(*,300) ep, i,  ctop(i), 'L2 Phase',
     .                        dms_L2/(1.d-3*vel_light)
             end if
         end if
      end do
      

****  See if we need to adjust time tag (Removed -0.5 from calc)
      if( num_av.gt.0 ) then
          average = average / num_av
          prev_site_clk = prev_site_clk + average
      end if

      dms = nint(prev_site_clk/(1.d-3*vel_light)-0.5d0)*1.d-3
      data_sec = nint(data_sec) + dms + init_dms

      call jd_to_ymdhms(data_epoch+1.d-4/86400.d0, date, sectag)     
      write(*,600) ep, num_av, prev_site_clk/(1.d-3*vel_light), 
     .             sectag+1.d-6, data_sec, average
 600  format(' Epoch ',i5,' clk from ',i4,' Value ',F8.3,' ms,',
     .       ' Old sectag ',F7.3,' New ',f7.3,' sec, Adj ',f9.2, ' m')

****  Now set up to write out the rinex file.  See if we need
*     to adjust time tag
      write(200, 710) date_obs(1), date_obs(2), date_obs(3), 
     .                date_obs(4), date_obs(5),
     .                data_sec, plf, num_chan, (ctop(i), i=1,num_chan)
 710  format(5i3,f11.7,20i3)
      do i = 1, num_chan
         outline = ' '
         do j = 1, num_data_types
            if( data_types(j).eq.'C1' ) then
                if( c1o(i).ne.0.0d0 ) then
                    write(outline((j-1)*16+1:),'(F14.3)') c1o(i)
                end if
                if( c1i(i).ne.0 ) then
                    write(outline(j*16-1:),'(I1)') c1i(i)
                end if
                if( c1f(i).ne.0 .and.c1f(i).ne.99 ) then
                    write(outline(j*16  :),'(I1)') c1f(i)
                end if
            end if
            if( data_types(j).eq.'P1' ) then
                if( p1o(i).ne.0.0d0 ) then
                    write(outline((j-1)*16+1:),'(F14.3)') p1o(i)
                end if
                if( p1i(i).ne.0 ) then
                    write(outline(j*16-1:),'(I1)') p1i(i)
                end if
                if( p1f(i).ne.0 .and.p1f(i).ne.99 ) then
                    write(outline(j*16  :),'(I1)') p1f(i)
                end if
            end if
            if( data_types(j).eq.'L1' ) then
                if( L1o(i).ne.0.0d0 ) then
                    write(outline((j-1)*16+1:),'(F14.3)') L1o(i)
                end if
                if( L1i(i).ne.0 ) then
                    write(outline(j*16-1:),'(I1)') L1i(i)
                end if
                if( L1f(i).ne.0 .and.L1f(i).ne.99 ) then
                    write(outline(j*16  :),'(I1)') L1f(i)
                end if
            end if
            if( data_types(j).eq.'P2' ) then
                if( p2o(i).ne.0.0d0 ) then
                    write(outline((j-1)*16+1:),'(F14.3)') p2o(i)
                end if
                if( p2i(i).ne.0 ) then
                    write(outline(j*16-1:),'(I1)') p2i(i)
                end if
                if( p2f(i).ne.0 .and.p2f(i).ne.99 ) then
                    write(outline(j*16  :),'(I1)') p2f(i)
                end if
            end if
            if( data_types(j).eq.'L2' ) then
                if( L2o(i).ne.0.0d0 ) then
                    write(outline((j-1)*16+1:),'(F14.3)') L2o(i)
                end if
                if( L2i(i).ne.0 ) then
                    write(outline(j*16-1:),'(I1)') L2i(i)
                end if
                if( L2f(i).ne.0 .and.L2f(i).ne.99 ) then
                    write(outline(j*16  :),'(I1)') L2f(i)
                end if
            end if
         end do

*        Now write out the outline
         write(200,'(a80)') outline(1:80)
      end do



****  Thats all
      return
      end

CTITLE READ_RANGE
 
      subroutine read_range(eof)

      implicit none 
 
*     This routine will read the next group of ranges from the rinex file
*     eof if returned true if we run out of data.
 
      include '../includes/const_param.h'
      include 'fixms.h'
 
*   date(5)     - Ymdhm of observation
*   flags(5)    - Flags read from file.
*   ierr        - IOSTAT Error
*   i,j         - Loop counter
*   id          - Dummy entry for power fail flag
 
      integer*4  flags(5), ierr, i,j,  trimlen, ili(5)
 
*   sectag      - Seconds tag
*   vals(5)     - Range values read
 
      real*8 sectag, vals(5)
 
*   eof     - Indicates end of file.
 
      logical eof
 
      character*256 line
 
****  Read in the next line
      read(101,'(a)', iostat=ierr) line
      if( ierr.ne.0 .or. trimlen(line).eq.0 ) then
          eof = .true.
          RETURN
      end if
 
*     Decode the line
* MOD TAH 980919: Updated read for ne rinex standard with G for
*     GPS satellites.
c     read(line,*,iostat=ierr) date_obs, sectag, plf, num_chan,
c    .            (chan(i),i=1,num_chan)
      read(line,100,iostat=ierr) date_obs, sectag, plf, num_chan,
     .            (chan(i),i=1,num_chan)
 100  format(5i3,f11.7,i3,i3,12(1x,i2))

      call ymdhms_to_jd( date_obs, sectag, data_epoch)
      data_sec = sectag
 
*     Now loop over the data records.
      do i = 1, num_chan
          read(101,120, iostat=ierr) (vals(j), ili(j), flags(j), j =1,5)
 120      format(10(f14.3,i1,i1))
 
*         Now assign the phase and range measurements
          P1f(i) = 99
          P2f(i) = 99
          C1f(i) = 99
          C2f(i) = 99
          do j = 1,5
              if( data_types(j).eq.'L1' ) L1o(i) = vals(j)
              if( data_types(j).eq.'L2' ) L2o(i) = vals(j)
              if( data_types(j).eq.'P1' ) P1o(i) = vals(j)
              if( data_types(j).eq.'C1' ) C1o(i) = vals(j)
              if( data_types(j).eq.'P2' ) P2o(i) = vals(j)
              if( data_types(j).eq.'C2' ) C2o(i) = vals(j)
*             Save the SNR flags
              if( data_types(j).eq.'L1' ) L1f(i) = flags(j)
              if( data_types(j).eq.'L2' ) L2f(i) = flags(j)
              if( data_types(j).eq.'P1' ) P1f(i) = flags(j)
              if( data_types(j).eq.'C1' ) C1f(i) = flags(j)
              if( data_types(j).eq.'P2' ) P2f(i) = flags(j)
              if( data_types(j).eq.'C2' ) C2f(i) = flags(j)
*             Loss of lock indicators
              if( data_types(j).eq.'L1' ) L1i(i) = ili(j)
              if( data_types(j).eq.'L2' ) L2i(i) = ili(j)
              if( data_types(j).eq.'P1' ) P1i(i) = ili(j)
              if( data_types(j).eq.'C1' ) C1i(i) = ili(j)
              if( data_types(j).eq.'P2' ) P2i(i) = ili(j)
              if( data_types(j).eq.'C2' ) C2i(i) = ili(j)
          end do

*         Should repeat the above code but wont be needed at the
*         moment.
          if( num_data_types.gt.5 ) 
     .        read(101,120, iostat=ierr) (vals(j), flags(j), 
     .                   j =1,num_data_types-5)
      end do

      if( ierr.ne.0 ) eof = .true.
 
****  Thats all
      return
      end
 
 
CTITLE READ_DATA_HEAD
 
      subroutine read_data_head

      implicit none 
 
*     Routine to read the header from the data files.
*     Returns will be:
*     Approximate site position
*     station name

      include '../includes/const_param.h' 
      include 'fixms.h'
 
* Local variables
 
*   ierr    - IOSTAT error
*  trimlen  - Length of string
*   indx    - Position of substring in a string
 
      integer*4 ierr, trimlen, indx, i, date(5)
      real*8 sectag
 
*   eoh     - End of header flags
*   rinex   - True if rinex file
 
      logical eoh, rinex
 
*   line    - Line read from file
 
 
      character*256 line, rxout_name
      character*1 cr
 
****  Open the data file
      cr = char(13) 
      open(101,file=data_file, iostat=ierr, status='old')
      call report_error('IOSTAT',ierr,'open',data_file, 1,
     .        'read_data_file/fixms')

*     OPen the output file
      rxout_name = data_file(1:trimlen(data_file)) // '.fix'
      open(200,file=rxout_name , iostat=ierr, status='new')
      call report_error('IOSTAT',ierr,'creat',rxout_name, 1,
     .        'read_data_file/fixms')

 
*     Start looping over the header
      eoh = .false.
      rinex = .false.
      do while ( .not.eoh )
          read(101,'(a)', iostat=ierr) line
          call sub_char(line,cr,' ')
          call casefold(line)
          call report_error('IOSTAT',ierr,'read', data_file,1,
     .                'End of hedader not found')
          if( trimlen(line).lt.10) eoh = .true.
          if( index(line,'END OF HEAD').gt.0 ) eoh = .true.
 
*         See if we find Rinex file
          indx = index(line,'OBSERVATION DATA')
          if( indx.gt.0 ) rinex = .true.

*         See if we find Rinex file
          indx = index(line,'COMMENT')
          if( indx.gt.0 ) rinex = .true.
 
*         See if marker name
          indx = index(line,'MARKER NAME')
          if( indx.gt.0 ) site_name = line(1:8)
 
*         See if position
          indx = index(line,'APPROX POSITION')
          if( indx.gt.0 .and. site_xyz(1).eq.0 ) then
              read(line,*,iostat=ierr) site_xyz
          end if
 
*         See if data types
          indx = index(line,'TYPES OF OBS')
          if( indx .gt.0 ) then
              read(line,*, iostat=ierr) num_data_types, 
     .               (data_types(i), i=1,num_data_types)
          end if

*         Write the line out. If we are at the end of headesr
*         then write comment out about fix:
          if( eoh ) then
              call systime( date, sectag )
              write(200,120) (date(i),i=1,3), 
     .                        prev_site_clk/(1.d-3*vel_light)
 120          format('Ranges fixed with FIXMS ',i5,'/',i2,
     .               '/',i2,' Init clk ',f7.3,' ms     ',
     .               'COMMENT             ')             
              if( fix_phase ) write(200,125)
 125          format('Phases adjusted by millisec jumps in FIXMS',
     .               '                  COMMENT             ')             

          end if
          write(200,'(a80)') line(1:80)
      end do
 
*     See if we found rinex file
      if( rinex ) then
          write(*,150) data_file(1:trimlen(data_file)),
     .        site_name, site_xyz
 150      format('* For RINEX data file ',a,/,
     .        '* Site ',a,' Aprrox. position ',3F15.3)
      else
          write(*,170) data_file(1:trimlen(data_file))
 170      format(a,' Does not appear to be RINEX file')
          stop 'fixms: Wrong type data file'
      end if
 
****  Thats all
      return
      end
 
 
 
CTITLE GET_SVSRUN
 
      subroutine get_svsrun
 
      implicit none 

*     Routine to get the cunstring.
 
      include '../includes/const_param.h'
      include 'fixms.h'
 
* LOCAL VARIABLES
 
*  len_run     - Length of the runstring element
*  trimlen     - Length of string
* date_start(5), date_stop(5) - Start and stop dates
*               for generating results
*  i           - Loop counter
*  rcpar       - Gets runstring entry
 
 
      integer*4 len_run, i, rcpar
 
*  runstring   - Elements of runstring
 
      character*256 runstring
 
*     Get the first runstring parameter
      len_run = rcpar(1,nav_file)
      if( len_run.le.0 ) then
          call proper_runstring('fixms.hlp','fixms/nav file',1)
      end if
 
*     See if latitude and longitude passed
      len_run = rcpar(2,data_file)
      if( len_run.le.0 ) then
          call proper_runstring('fixms.hlp','fixms/data file',1)
      end if

*     See if phase shold be fixed
      len_run = rcpar(3,runstring)
      fix_phase = .false.
      if( len_run.gt.0 ) then
          if( runstring(1:1).eq.'p' .or. runstring(1:1).eq.'P') then
              fix_phase = .true.
          end if
      end if

*     Now see if positions passed 
      do i = 1,3
         len_run = rcpar(3+i,runstring)
         if(  len_run.gt.0 ) then
             read(runstring,*) site_xyz(i)
         else
             site_xyz(i) = 0.d0
         end if
      end do

*     See if initial clock offset passed. (in millisec)
      len_run= rcpar( 7, runstring)
      if( len_run.gt.0 ) then
          read(runstring,*) prev_site_clk  
*         Convert values to meters
          prev_site_clk = prev_site_clk*vel_light*1.d-3
      else
          prev_site_clk = 0.d0
      end if
 
*     See if debug wanted
      len_run= rcpar( 8, runstring)
      if( len_run.gt.0 ) then
          read(runstring,*) debug_start
      else
          debug_start = -1   
      end if
 
*     See if debug wanted
      len_run= rcpar( 9, runstring)
      if( len_run.gt.0 ) then
          read(runstring,*) debug_end   
      else
          if( debug_start.ne.-1   ) then
              debug_end = debug_start
          else
              debug_end = 0
          end if 
      end if
 
****  Thats all
      return
      end
 
CTITLE READ_NAV
 
      subroutine read_nav
 
      implicit none 

*     Thuis routine will read the navigation file.  Only the first
*     occurence of a satellite ephemeris entry will be used.  (This
*     is set by the PRN number still being zero)
 
      include 'fixms.h'
 
* LOCAL VARIABLES
 
*   ierr        - IOSTAT error
*   i           - Loop counter
*   rp          - Read Prn number
*   trimlen     - length of string
*   ll          - Length of line string
*   date(5)     - ymdhm
*   rinex_version - Version of rinex file
 
 
      integer*4 ierr, i, rp, trimlen, ll, date(5) 
      real*4 rinex_version
 
*   sectag      - Seconds tag in date
 
       real*8 sectag
 
*   still_header    - Indicates that we are sill in the header
*                   - part of file.
*   have            - Set true if we already have this satellite
 
       logical still_header, have
 
*   line            - line read from file
 
       character*256 line

*   cr - Carriage return (for handling dos files)
      character*1 cr

      cr = char(13)
 
****  Open the nav_file
 
      open(100, file=nav_file, status='old', iostat=ierr)
      call report_error('IOSTAT',ierr,'open',nav_file,1,
     .            'fixms/read_nav')
 
*     Loop over the file reading the ephemeris entries.  First clear
*     all of the PRN numbers so we know when a PRN has been read
 
      do i = 1, max_sat
          prn(i) = 0
      end do
      rinex_version = 1
 
      still_header = .true.
      do while ( still_header )
          read(100,'(a)', iostat=ierr) line
          call sub_char(line,cr,' ')
          call casefold(line)
          ll = trimlen(line)
          if( ll.eq.0 .or. index(line,'END OF HEAD').gt.0 ) then
              still_header = .false.
          end if
          call report_error('IOSTAT',ierr,'read','NAV FILE HEADER',
     .                      1,'read_nav')
          if( index(line,'RINEX VERSION').gt.0 ) then
              read(line,*, iostat=ierr) rinex_version
              write(*,120) rinex_version
 120          format('* Rinex version ',f5.1,' Nav file found')
          end if
      end do
 
      num_sat = 0
 
*     Now start reading the entries
      do while ( ierr.eq.0 )
          read(100,'(a)',iostat=ierr) line
          call sub_char(line,cr,' ')
          if( trimlen(line).eq.0 ) ierr = -1
*                                     ! See which PRN
          if( ierr.eq.0 ) then
              read(line,*) rp
 
*             See if we already have this prn
              have = .false.
              do i = 1, num_sat
                  if( rp.eq.prn(i) ) then
                      have = .true.
                  end if
              end do
 
*                                         ! This is first
              if( .not.have ) then
*                                             ! on this prn so read
                  num_sat = num_sat + 1
                  rp = num_sat
                  read(line,200) prn(rp), date, sectag,
     .                        af0(rp),af1(rp), af2(rp)
 200              format(i2,5i3,f5.1,3d19.8)
 
                  call ymdhms_to_jd( date, sectag, toe_jd(rp))
 
*                 Read the rest of the entires
                  read(100,210) aode(rp), crs(rp), dn(rp), m0(rp)
                  read(100,210) cuc(rp), ecc(rp), cus(rp), art(rp)
                  read(100,210) toe(rp), cic(rp), om0(rp), cis(rp)
                  read(100,210) i0(rp) , crc(rp), w(rp)  , omd(rp)
                  read(100,210) idt(rp), cflg12(rp), weekno, pflg12(rp)
                  read(100,210) svacc, svhealth, tgd, aodc(rp)
 210              format(3x,4d19.8)
                  if( rinex_version.gt.1 ) read(100,'(a)') line
              else
*                 Skip 6 lines in file
                  do i = 1,6
                      read(100,'(a)') line
                  end do
                  if( rinex_version.gt.1 ) read(100,'(a)') line
              end if
          end if
      end do
 
      write(*,300) num_sat, nav_file(1:trimlen(nav_file))
 300  format('* ', i5,' satellites found in ',a)
 
      if( num_sat.eq.0 ) stop ' fixms: No satellites found'
 
      return
      end
CTITLE EPH_TO_XYZ
 
      subroutine eph_to_xyz(t, i, sys)

      implicit none 
 
*     Routine to compute XYZ coordinates of satellite at time t for
*     sattellite number i.  SYS is 'I' or 'E' for inertial or Earth
*     fixed.
 
      include 'fixms.h'
 
*   t       - Time for computation (day number)
 
 
      real*8 t
 
*   i       - Satellite number
 
 
      integer*4 i
 
*   sys     - System for results.
 
 
      character*(*) sys
 
* LOCAL VARIABLES
 
*   j           - Loop counter in Kepler's equation
 
 
      integer*4 j
 
*   gm      - GM
*   eom     - Earth rotation rate (rads/sec)
*   a       - Semimajor axis
*   n0      - Mean motion
*   tk      - Time of epoch from toe (seconds)
*   n       - COrrected mean motion
*   mk      - Mean anomaly
*   ek      - Eccentric anomaly
*   vk      - true anomaly
*   sinvk, cosvk    - Sin and cos of true anomaly
*   pk      - argument of latitude
*   duk, drk, dik   - Coorections to arg of lat, radius and
*           - inclinations
*   uk      - Argumenr of latitude
*   rk      - radius at time tk
*   ik      - Inclination at tk
*   xpk, ypk    - Inplane coordiantes
*   omk     - Longitude of the asecdning node
*   rot_mat(3,3)    - Rotation matrix from XYZ to NEU
 
 
 
 
      real*8 gm, eom, a, n0, tk, n, mk, ek, vk, sinvk, cosvk, pk,
     .    duk, drk, dik, uk, rk, ik, xpk, ypk, omk, rot_mat(3,3)
 
      gm = 3.986005d14
      eom = 7.2921151467d-5
 
      a = art(i)*art(i)
      n0 = sqrt(gm/a**3)
 
      tk = (t-toe_jd(i))*86400.0d0
      n = n0 + dn(i)
      mk = m0(i) + n*tk
 
****  Solve Keplers equation
      ek = mk
      do j = 1, 10
          ek = mk + ecc(i)*sin(ek)
      end do
 
***** Get the true anomaly
      sinvk = sqrt(1-ecc(i)**2)*sin(ek)/(1 - ecc(i)*cos(ek))
      cosvk = (cos(ek)-ecc(i))/(1-ecc(i)*cos(ek))
 
      vk = atan2(sinvk, cosvk)
 
*     Argument of latitude
      pk = vk + w(i)
 
*     Correction terms
      duk = cus(i)*sin(2*pk) +cuc(i)*cos(2*pk)
      drk = crs(i)*sin(2*pk) +crc(i)*cos(2*pk)
      dik = cis(i)*sin(2*pk) +cic(i)*cos(2*pk)
 
      uk = pk + duk
      rk = a*(1-ecc(i)*cos(ek)) + drk
      ik = i0(i) + dik + idt(i)*tk
 
*     Get inplane coordinates
      xpk = rk*cos(uk)
      ypk = rk*sin(uk)
 
*     Compute the longitude of the ascending node
      omk = om0(i) + omd(i)*tk
 
*     If we are in Earth fixed frame account for rotation of Earth
      if( sys(1:1).eq.'E' .or. sys(1:1).eq.'e') then
          omk = omk - eom*(tk+toe(i))
      end if
 
*     Get three_d coordinates
      svs_xyz(1,i) = xpk*cos(omk) - ypk*sin(omk)*cos(ik)
      svs_xyz(2,i) = xpk*sin(omk) + ypk*cos(omk)*cos(ik)
      svs_xyz(3,i) = ypk*sin(ik)
 
*     Now convert to latitude and longitude
      call xyz_to_geod( rot_mat, svs_xyz(1,i), svs_loc(1,i))
 
****  Thats all
      return
      end
 
 
 
