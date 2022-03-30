CTITLE READ_ALL_RINEX

      subroutine read_all_rinex(ns)

      implicit none

*     Routine to read all the rinex data into memory.  The indexing is
*     based on the ref_start mjd.
      
      include '../includes/const_param.h'   
      include '../includes/xfile_def.h'
      include 'track_com.h'

* PASSED VARIABLES
* ns  -- Station number to read (used to connect to the logical unit also)

      integer*4 ns,marker

* LOCAL VARIABLES
* i,j,k  -- Loop counters
* ep     -- Epoch number for the current measurements

      integer*4 i,j,k, ep, trimlen

* frep -- Fractional epoch.  Used to check to see if we should take
*     this data epoch
* ad_mjd -- Adjusted time tag based on psuedorange values
* min_psr -- Minimum psuedorange used to get the time tag adjustment
* max_eperr, max_epusd -- Maximum epoch error seen and used

      real*8 frep, ad_mjd, min_psr, max_eperr, max_epusd

* eof    -- EOF (set true at end of file)
* used   -- Logical to indicate that a PRN has already been added to list
* SP3_OK -- Set true if satellite appears in SP3 file
* reported_sp3(max_sat) -- Set true once a satellite has been reported as
*           not being in SP3 file
* curr_stopgo -- Current status stopgo mode (false when reciever is static
*     and true when kinematic (ID 2 record changes to true, ID 3 changes 
*     to false)
      logical eof
      logical used, sp3_OK, reported_sp3(max_sat), rx2_read,
     .        curr_stopgo

* buffer -- Character buffer used to read rinex/xfiles
      character*256 buffer
* MOD AZ 190305: additional variable for issue explained in Line 282
* marker -- defined to solve some unconventional issue
      marker = 0
      eof = .false.
      curr_stopgo = .false.
      max_eperr = 0.d0
      max_epusd = 0.d0

*     Restore the header information about the types of data in
*     this file

      call get_headinfor(ns)

      obs_lu = 50 + ns  ! rd_rinex_obs assumes ns from obs_lu
      
*     Clear the pointer arrays to the data so that we know where
*     we actually have data
      do i = 1, max_epochs
         num_chan_se(ns,i) = 0
      end do
      do i = 1, max_sat 
         reported_sp3(i) = .false.
      end do 

*     Loop over the data file reading the rinex file and saving the
*     results

c      eof = .false.
C      curr_stopgo = .false.
c      max_eperr = 0.d0
c      max_epusd = 0.d0

      do while ( .not.eof  )

*         Get line with list of satellites from either the xfile or
*         rinex file.  (The read code 0 indicates to read directly
*         from the file)

          rx2_read = .false.
          if(obs_file_type(ns).eq."R")  call read_rinex_obs_1(eof,
     .             buffer,0, curr_stopgo)
          if(obs_file_type(ns).eq."X") call read_xf_obs_1(eof,
     .             buffer, 0)

*         Based on the time, compute the epoch number for this measurement
          frep = (xf_mjd - ref_start)*86400.d0/usr_interval + 1
          ep = nint((xf_mjd - ref_start)*86400.d0/usr_interval) + 1
          if( ep.gt.usr_nepochs .and. usr_nepochs.gt.0  ) eof = .true.
c          write(*,*) ep,' ',xf_mjd,' ',eof
* MOD TAH 000823: Changed the epoch tolerance calculation to account for
*         possible offsets in the receiver clock (based on psuedorange).
*         If the epock looks OK, read the data record and then use
*         the psuedorange to get it basically correct.
          if( ep.gt.0 .and. ep.le.max_epochs .and.  .not.eof ) then

*             OK, data is past first epoch time and we have not reached
*             EOF so decode the data.
c          write(*,*) eof
              if(obs_file_type(ns).eq."R") 
     .                     call read_rinex_obs_2(eof, ep)
              rx2_read = .true.
c           write(*,*) '101 ',eof
              if( ep.ge.debug_start .and. ep.le.debug_end ) then
                 do j = 1, xf_msat
!                  print *,'RAW OBS ', ep, ns, ' CH ',j, xf_obsv(1:6,j)
                 end do 
              endif   

****          Now check the "fine tuning" of the epoch.  Here we use the
*             L1 psuedo-range to approximately adjust the epoch to make sure
*             that clock has not drifted to much.  Use the minimum value of
*             the pseudorange
              min_psr = 1.d20
c              write(*,*) '113 ',eof
              do j = 1, xf_msat
* MOD TAH 100717: Check P1 and C1 separate
C                if ( xf_obsv(3,j).lt. min_psr ) min_psr = xf_obsv(3,j)
                 if ( xf_obsv(3,j).ne.0 .and.
     .                xf_obsv(3,j).lt. min_psr ) min_psr = xf_obsv(3,j)
                 if ( xf_obsv(5,j).ne.0 .and.
     .                xf_obsv(5,j).lt. min_psr ) min_psr = xf_obsv(5,j)

              end do

* MOD TAH 050624: Check that the min_psr is not the inital value (indicating
*             pseudornage data).  If not changed set back to nominal
              if( min_psr.eq. 1.d20 ) min_psr  = 20.d6

*             Adjust the time tag by accounting for psuedorange
              ad_mjd = xf_mjd - (min_psr-20.d6)/vel_light/86400.d0

*             Now re-compute the epochs with adjusted time tags
              frep = (ad_mjd - ref_start)*86400.d0/usr_interval + 1
              if ( ep.ge.debug_start .and. ep.le.debug_end ) 
     .            print *,'EP TEST ',ep, frep ,
     .            nint((ad_mjd - ref_start)*86400.d0/usr_interval) + 1  
              ep = nint((ad_mjd - ref_start)*86400.d0/usr_interval) + 1
c             write(*,215) 'TT ',ep, frep, (ad_mjd-xf_mjd)*86400.d0,
c     .              (xf_mjd - 
c     .              ((ep-1)*usr_interval/86400.d0+ref_start))*86400.d0,
c     .              (frep-ep)*usr_interval
c 215          format(a,1x,i8,F12.3, 3F12.4)
              if( abs((frep-ep)*usr_interval).gt.max_eperr .and.
     .            (abs(frep-ep)*usr_interval).lt.0.01d0  )
     .            max_eperr = abs((frep-ep)*usr_interval)
           end if

c           write(*,*) '146 ',eof
*          Continue only if the epoch is with 0.01 seconds of the correct
*          value.
           if( ep.gt.0 .and. ep.le.max_epochs .and.  .not.eof .and.
     .        (abs(frep-ep)*usr_interval).lt.0.1d0 ) then

*             Save the largest epoch error actually used 
              if( abs((frep-ep)*usr_interval).gt.max_epusd )
     .            max_epusd = abs((frep-ep)*usr_interval)

c              write(*,*) '155 ',eof
*             Get the offset from the nominal time tag for this 
*             measurement
              sec_offset(ns,ep) = (xf_mjd - 
     .              ((ep-1)*usr_interval/86400.d0+ref_start))*86400.d0
              if( ep-(ep/10000)*10000.eq.0 ) then
                  write(*,220) ns, ep, xf_mjd, sec_offset(ns,ep)
 220              format('Site ',i2,' Epoch ',i8,' MJD ',F16.8,
     .                  ' Time offset (sec) ',F8.5)
              end if

*             Now save the data in the memory arrays
              if ( ep.gt.debug_start .and. ep.le.debug_end ) 
     .        print *,'DEBUG: EP, NS ',ep, ns, ' XF_MSAT ', xf_msat
              num_chan_se(ns,ep) = xf_msat
              if ( xf_msat.gt.max_chan ) then
                  write(*,225) xf_msat, max_chan
 225              format('** DISASTER ** Number of channels ',I4,
     .                   ' exceeds max_chan ',i4)
                  stop 'TRACK: Too many channels needed'
              endif
 
c              write(*,*) '176 ',eof
              do j = 1, xf_msat

*                Pull of the different types of data that we need
                 L1o_all_cse(j,ns,ep) = xf_obsv(1,j)
                 L2o_all_cse(j,ns,ep) = xf_obsv(2,j)
                 P1o_all_cse(j,ns,ep) = xf_obsv(3,j)
* MOD TAH 100710: If P1 range is zero; use the C1 range
                 if( xf_obsv(3,j).eq.0 ) P1o_all_cse(j,ns,ep) = 
     .                                  xf_obsv(5,j)
                 P2o_all_cse(j,ns,ep) = xf_obsv(4,j)
* MOD TAH 130425: Use C2 if P2 is not available.
                 if( xf_obsv(4,j).eq.0 ) P2o_all_cse(j,ns,ep) = 
     .                                  xf_obsv(6,j)

*                Save some additional information about this epoch
                 data_flag_cse(j,ns,ep) = 1
* MOD TAH 040316: See if in kinematic mode
                 if( curr_stopgo ) then
                     call sbit(data_flag_cse(j,ns,ep),6,1)
                 endif

c                 write(*,*) '197 ',eof
                 if( debug_start.eq.-2 ) then
                     write(*,230) ep, site_names(ns), xf_iprn(j), 
     .                     xf_obsv(1:5,j),
     .                     xf_obsv(1,j) - (fR1/fR2)*xf_obsv(2,j),
     .                     xf_obsv(1,j) - xf_obsv(2,j) -
     .                     dfsf*(xf_obsv(3,j)*fR1/vel_light+
     .                           xf_obsv(4,j)*fR2/vel_light),
     .                     xf_obsv(1,j) - xf_obsv(2,j) -
     .                     dfsf*(xf_obsv(5,j)*fR1/vel_light+
     .                           xf_obsv(4,j)*fR2/vel_light)

 230                 format('XOBS ',i5,1x,a4,1x,' PRN ',i3,1x,
     .                    5(F15.3,1x),' WL ',3(F15.3,1x))
                 end if

                 amb_point_cse(j,ns,ep) = 0
                 ctop_cse(j,ns,ep) = xf_iprn(j)

c                 write(*,*) '215 ',eof
****             Keep track of PRN's used
                 used = .false.
                 do k = 1, num_prn
* MOD TAH 150411: Restrict selection to GPS with PRN <= 32.
                    if( xf_iprn(j).eq.prn_used(k)  ) used = .true.
                 end do

                 if( .not.used  ) then
                    num_prn = num_prn + 1
                    prn_used(num_prn) = xf_iprn(j)
                 end if

c                 write(*,*) '227 ',eof
*****            See if we have an edit for this data: 
* MOD TAH 101113: Allow -1 as PRN number to edit all.
                 do i = 1, num_edits
                    if( ss_edit(1,i).eq.ns .and. 
     .                  (ss_edit(2,i).eq.xf_iprn(j) .or.
     .                   ss_edit(2,i).eq. -1)  .and.
     .                  xf_mjd.ge.tt_edit(1,i) .and.
     .                  xf_mjd.le.tt_edit(2,i)      ) then
*                       User has said not to use this data
                        call sbit(data_flag_cse(j,ns,ep),2,1)
                    end if
                 end do

c                 write(*,*) '240 ',eof
*****            See if usr_bflag to be applied
                 do i = 1,num_abf
                    if( ss_abf(1,i).eq.ns .and. 
     .                  ss_abf(2,i).eq.xf_iprn(j) .and.
     .                  xf_mjd.ge.tt_abf(i) ) then
*                       Add bias flag
                        call sbit(data_flag_cse(j,ns,ep),3,1)
                        call sbit(data_flag_cse(j,ns,ep),4,1)
                        call sbit(data_flag_cse(j,ns,ep),7,1)  ! Set bit to 
*                                            | show user added.
                        ss_abf(1,i) = -ns   ! Stops being reapplied
                        write(*,240) i, ep, 
     .                     site_names(abs(ss_abf(1,i))), ss_abf(2,i) 
 240                    format('Adding user bias ',i3,' at Epoch ',i5,
     .                         ' site ',a4,' PRN ',i3.2) 
                    endif
                 end do
                    
c                 write(*,*) '258 ',eof
*****            See if the satellite is in the SP3 file 
                 sp3_OK = .false.
                 do i = 1, num_sat
                    if( prn_sp3(i).eq.xf_iprn(j) ) sp3_OK = .true.
                 end do
c                 write(*,*) '272 ',eof
                 if( .not.sp3_OK ) then
                     data_flag_cse(j,ns,ep) = 3
c                     write(*,*) '275 ',eof
* MOD AZ 190305: in some Linux system, the "eof" will turn ".true."
*                automatically, investigation around this seems to indicate
*                some unconventional definition of bool variable "eof"
*                The "marker" here is to solve this issue.
                     if (.not.eof) marker = 1
                     if( .not.reported_sp3(xf_iprn(j)) ) then
                         reported_sp3(xf_iprn(j)) = .true.
c                         write(*,*) '278 ',eof
                         write(*,280) xf_iprn(j)
 280                     format('**WARNING** PRN ',i3,' Not in',
     .                          ' SP3 file')
*                        Add to excluded 
c                         write(*,*) '281 ',eof
                         num_exclude = num_exclude + 1
                         if( num_exclude.gt. max_prn ) then
                            write(*,290)  num_exclude
  290                       format('** ERROR ** Too many excluded ',
     .                             'satellites ',I4)
                            stop 'TRACK: Too many exclude satellites'
                         endif
                         ss_exclude(num_exclude) = xf_iprn(j)
                     end if
                     if ( marker.eq.1 ) eof = .false.
                 end if
              end do
          else
*             Skip over these measurements
c              write(*,*) '285 ',eof
              if(obs_file_type(ns).eq."R" .and. .not.rx2_read
     .           .and. .not.eof ) then
                     call read_rinex_obs_2(eof, ep) 
              endif
c              write(*,*) '290 ',eof
   
          end if
c         write(*,*) eof 
         if( ep.gt.max_epochs ) then
              write(*,320) ns, max_epochs 
320           format('File ',i2,' has more epochs than allowed:',
     .               ' max_epochs is ',i6)
              eof = .true.
              ep = ep - 1
          end if
c          write(*,*) eof
      end do

* MOD TAH 160812: Remove any prn's that are in the ss_exclude list
*    (needed when excluded satellites are added for not being in SP3
*     file)
      do j = 1, num_exclude
         i = 0
         do while ( i.lt. num_prn )
            i = i + 1
            if( prn_used(i).eq.ss_exclude(j) ) then
*               Remove this PRN
                do k = i, num_prn-1
                   prn_used(k) = prn_used(k+1)
                enddo
                i = i - 1
                num_prn = num_prn - 1
             end if
         end do
      end do


****  Tell user about the file:
      write(*,410) site_names(ns), ep, num_prn, 
     .            (prn_used(i),i=1,num_prn)
 410  format('For site ',a4,' last epoch found ',i6,/,
     .       'Found ',i3,' PRNS -> ',100(I3.2,1x))
      num_epochs = max(num_epochs,ep)
      data_end(ns) = (ep-1)*usr_interval/86400.d0 + ref_start
      write(*,430) max_eperr, max_epusd
 430  format('Maxiumum epoch errors: Seen ',F9.4,' sec, and ',
     .       'Used ',F8.4,' sec')

      if( trimlen(rwl_root).gt.0 ) then
          print *,'Outputing Raw Widelanes with ',trim(rwl_root)
          call output_wls(ns,'RWL')
      end if

****  Thats all 
      return
      end

CTITLE FLAG_GAPS

      subroutine flag_gaps(ns)

      implicit none

*     Routine to scan the data run from the data files and flag
*     all the gaps greater than min_gap.  If the amount of good data 
*     between gaps is less than min_dur then the data is marked bad.
*     This routine now actually uses the bias flags pre-inserted by
*     scan_ddlg subroutine.
                       
      include 'track_com.h'

* PASSED VARIABLES
* ns  -- Station number to read (used to connect to the logical unit also)

      integer*4 ns

* LOCAL VARIABLES
* i, j, k -- Loop counters
* ep      -- epoch number loop counter
* ch      -- Channel number corresponding to PRN number
* num_in_gap -- Number of missing data in gap
* last_good_ep -- Last epoch at which we had good data
* last_bf_ep   -- Epoch of start of last bias flag (used to make sure
*            we have enough data to determine bias flag)

      integer*4 i, k, ep, ch, num_in_gap, last_good_ep,
     .          last_bf_ep, num_good 


* first_found -- Logical to indicate that we found the first good
*     observation for a particular site and satellite
* marked      -- Set true when first bias flag for site/PRN marked
* data_good   -- Set true for good data point
* data_OK     -- Logical function that returns true for good data
* kbit        -- Checks status of bit

      logical first_found, marked, data_good, data_OK,  
     .        kbit


****  Start by looping over all satellites.  At each gap we assign
*     a bias flag number

      do i = 1, num_prn
         marked = .false.
         first_found = .false.
         num_in_gap  = 0
         num_good = 0
         do ep = 1, num_epochs
           
*           See if we have the prn at this epoch
            data_good = .false.
            do k = 1, num_chan_se(ns,ep)
               if( ctop_cse(k,ns,ep).eq.prn_used(i) ) then

*                  We have a measurement to the satellite at this
*                  epoch.  See if measurement marked bad
                   ch = k
                   if( data_OK(data_flag_cse(k,ns,ep),data_mask) ) then
                       data_good = .true.
                       num_good = num_good + 1
                       if( .not.first_found ) first_found = .true.
                   end if
               end if
            end do

*           Now see if we have good data 
            if( data_good  ) then
*               See if have missed data before this and/or we are
*               marking the first good observation
C               if( (first_found .and. .not. marked) 
C    .              .or. kbit(data_flag_cse(ch,ns,ep),4)) then
                if( kbit(data_flag_cse(ch,ns,ep),4) ) then
*                   Save the last_good_ep of current bias flag.
*                   (Because we may have removed and ambiquity due
*                   to short data, the end epoch of prevous bias flag
*                   may already have been set, so do not overwrite.
                  
                    num_ambs = num_ambs + 1
                    bf_ents(1,num_ambs) = ns
                    bf_ents(2,num_ambs) = prn_used(i)
                    bf_ents(3,num_ambs) = ep
                    bf_ents(4,num_ambs) = 0
                    bf_ents(5,num_ambs) = 0

*                   Save the information about this bf_flag
                    
                    if( num_ambs.gt.max_ambs ) then
                        write(*,220) max_ambs
 220                    format('**DISASTER** Too many ambiquities are',
     .                         ' needed.  Max allowed is ',i5)
                        stop 'TRACK: Abort too many ambiquities needed'
                    end if
                    marked = .true.
 
*                   Mark the L1 and L2 bias flags
                    call sbit(data_flag_cse(ch,ns,ep),3,1)
                    call sbit(data_flag_cse(ch,ns,ep),4,1)
                    last_bf_ep = ep

                end if
 
*               Reset the number in gap and save the ambiquity information
C               num_in_gap = 0
                amb_point_cse(ch,ns,ep) = num_ambs
                last_good_ep = ep
*               Save last data point epoch of bias flags (updated with
*               each good measurements)
                bf_ents(4,num_ambs) = ep
            end if
         end do      ! Looping over epochs

      end do         ! Looping over satellites

****  Thats all
      return
      end

CTITLE DATA_OK

      logical function data_ok(flag,mask)

      implicit none

*     Routine which returns true is none of the bits in mask are
*     set in flag.

* PASSED VARIABLES
* flag  -- data flag being checked
* mask  -- Masking bits to check only those for bad data

      integer*4 flag, mask

* cand  -- Common library version of the and function

      integer*4 cand

* kbit  -- Logical to test but setting
      logical kbit

****  OK see if any bits match
      data_OK = .true.
      if( cand(flag,mask).gt.0 ) data_OK = .false.
*     Bit 1 needs to be set for there to be data, so check
      if( .not.kbit(flag,1) )  data_OK = .false.
      

****  Thats all
      return
      end

CTITLE SELECT_ARB_BF

      subroutine select_arb_bf

      implicit none

*     Routine to select the bias flags to be arbitarily set
*     as known.
      
      include 'track_com.h'

* LOCAL VARIABLES
* i,j,k -- Loop counter
* ns    -- Site loop counter
* last_main_bf -- Last BF number for the first site
* max_over_lap -- Maxiumum duration of the overlap between bias
*          flags
* max_ol_prn   -- PRN with maximum overlap
* max_ol_bf    -- Number of the bias flag with maximum overlap
* ol_start     -- First epoch of overlap
* ol_end       -- Last epoch of overlap
* over_lap     -- Duration of the overlap
* mxx_ent      -- Entry in the bf_ents that correspond to the longest
*                 overlap.  Used for the initial selection of which
*                 ambiquities will be set to 0.
* mxx_ovl      -- Overlap associated with maximum overlap.  This was used
*                 in test code to see if we should fix only some of reference
*                 station bias flags.

      integer*4 i,j,k, ns, last_main_bf, max_over_lap, max_ol_prn,
     .          max_ol_bf, ol_start, ol_end, over_lap, mxx_ent,
     .          mxx_ovl


****  Start by finding the last bias flag at the primary site.
      do i = 1, num_ambs
*        Check for the first site and mark bias flag as resolved
         if( bf_ents(1,i).eq.1 ) then
             last_main_bf = i
C            bf_ents(5,i) = 1 
         end if
      end do

****  Start:   Loop over all stations except the first.
      mxx_ovl = 0
      do ns = 2, num_site
*        Check the bias flags at this site.  Find the overlap
*        with the primary site
         max_ol_prn = 0
         max_over_lap = 0
         do j = last_main_bf+1, num_ambs
*           Find the over lapping bias flags with main site
            if( bf_ents(1,j).eq.ns ) then
*               BF entry for this station.  Find match to reference
*               site
                do k = 1, last_main_bf

                  if( bf_ents(2,k).eq.bf_ents(2,j) ) then
*                     Satellites are the same.  Compute the length
*                     of common data
                      ol_start = max(bf_ents(3,k),bf_ents(3,j))
                      ol_end   = min(bf_ents(4,k),bf_ents(4,j))
                      over_lap = ol_end - ol_start
                      if( over_lap.gt.max_over_lap ) then
                          max_over_lap = over_lap
                          max_ol_prn = bf_ents(2,j)
                          max_ol_bf  = j
                      end if
                   end if
                end do
             end if
         end do

*        Tell user what is happening.
         if( max_ol_prn.gt.0 ) then
*            We found an overlapping PRN, tell the user which we
*            have selected
             if ( 1.eq.1 )
     .       write(*,220) max_ol_prn, site_names(ns), 
     .            max_over_lap, max_ol_bf
 220         format('BIAS Set for PRN',i3.2,' Site ',a4,
     .              ' Duration ',i8,' epochs, BF # ',i4)
             bf_ents(5,max_ol_bf) = 1
             wls_ref(1,max_ol_bf) = -1

*            See if this overlap is greater than all others
             if( max_over_lap.gt.mxx_ovl ) then
                 mxx_ent = max_ol_bf
             end if
         else
             write(*,240) site_names(ns)
 240         format('**DISASTER** No overlapping data at site ',a4)
             stop
         end if
      end do

****  Now set the bias flags to be assumed known at the reference
*     site.  These will be all bias flags that overlap by at least
*     min_good with the longest overlap entry.
      do i = 1, last_main_bf
*        Check for the first site and mark bias flag as resolved
         if( bf_ents(1,i).eq.1 ) then

*            Check the over lap between this bias flag and the 
*            maximum duration one
             ol_start = max(bf_ents(3,i),bf_ents(3,mxx_ent))
             ol_end   = min(bf_ents(4,i),bf_ents(4,mxx_ent))
             over_lap = ol_end - ol_start
             if( over_lap.gt.min_good ) then
                 bf_ents(5,i) = 1
             end if
         end if
         bf_ents(5,i) = 1
      end do

****  Thats all
      return
      end 

CTITLE REPORT_BF

      subroutine report_bf(un,type)

      implicit none

*     Routine to report the status of the bias flags
*     in the solution
      
      include 'track_com.h'

* PASSED VARIABLES
* type  -- Character string with type of report

      character*(*) type
      integer*4 un     ! Unit number for output

* LOCAL VARIABLES
* i, j -- Loop counters

      integer*4 i,j

*     Write out the values
      write(un,120) type, type
 120  format('BIAS FLAG REPORT: Type ',1x,a,/, 
     .       '.  # Site  S  PRN      Start    Stop Fixd      ',
     .       'L1 cycles      L2 cycles      DD Bias Refs        ',
     .       'MW-WL      +-       EX-WL      +-    ',a)

      do i = 1, num_ambs
         write(un,140) i, site_names(bf_ents(1,i)),
     .                  (bf_ents(j,i),j=1,5), 
     .                  (ambiq_all(j,i),j=1,2),
     .                  (wls_ref(j,i),j=1,3),
     .                  (wls_all(j,i),wls_sig(j,i), j=1,2), type
 140     format(i4,1x,a4,1x,i2,1x,' PRN ',i3.2, 2I8,3x,o2.2,
     .         2F15.1,3i6,' WL ',4(F9.3,1x),1x,a)
      end do

****  Thats all
      return 
      end

CITLE RESOLVE_WL

      subroutine resolve_wl(ns, iter, tot_resolved)

      implicit none

*     Routine to resolve the MW-WL and ionospheric constraint
*     wide lane ambiquities.
      
      include 'track_com.h'

* PASSED VARIABLES
* ns -- Station number being resolved
* iter -- Iteration number, with each iteration we make it easier
*         to resolve the ambiquity to ensure that we all have
*         estimates by the time we are done.
* tot_resolved -- Total number of resolved ambiquities

      integer*4 ns, iter, tot_resolved

* LOCAL VARIABLES
* num_new_resolved -- Number of new ambiquities resolved with 
*     each iteration. Iteration stops if = 0
* i  -- Loop counter over ambiquities
* tar_b1  -- Target site ambiquity numbers to match with one being
*            checked
* ref_b1, ref_b2  -- Reference site ambiquity numbers

      integer*4 num_new_resolved, i, tar_b1,
     .          ref_b1, ref_b2
                                                     
* all_resolved -- Logical set true when all ambiquities are
*     resolved
* resolved     -- Set true when particular ambiquity resolved
* possible     -- Logical that stays true while ever there are
*     more combinations to check

      logical all_resolved, possible, resolved

      integer*4 count, cnt2

      integer*4 bf5_copy(max_ambs)


***** Start looping over the bias flags computing the MW-WL and
*     ionospheric constraint.

      all_resolved = .false.
      tot_resolved = 0
      num_new_resolved = 1
      count = 0
      do i = 1,num_ambs
         bf5_copy(i) = -1
      enddo

      do while ( .not.all_resolved .and.num_new_resolved.gt.0 
     .           .and. count.lt.500 )

*        Loop over the bias flags at this site 
         num_new_resolved = 0
         all_resolved = .true.
         count = count + 1
 
         do i = 1, num_ambs
            if( bf_ents(1,i).eq.ns .and. 
     .          bf_ents(5,i).eq.0        ) then
*               OK, bf matches this site and it is not resolved
*               Now find maximum overlap between this bias flag
*               and a resolved bias at this site
                all_resolved = .false.

*               Find the maximum overlap between a resolved and 
*               non-resolved bias flag.
                resolved = .false.
                possible = .true.
                ref_b1 = 0
                ref_b2 = 0
                tar_b1 = 0
                cnt2 = 0
                
                do while ( .not.resolved .and. possible 
     .                     .and. cnt2.lt. 500 )
*                  The tar_b1 value is incremented in side the
*                  resolution loop.
                   cnt2 = cnt2 + 1
                   call res_one_wl( i, tar_b1, ref_b1, ref_b2,
     .                              resolved, possible, iter )
                   if( i.eq.-1 ) then
                      print *,'RES_OW ',i, tar_b1, ref_b1, ref_b2,
     .                        count, cnt2, resolved, possible,iter
                   endif
                end do
                if( resolved ) num_new_resolved = num_new_resolved + 1
                
            end if
         end do
         tot_resolved = tot_resolved + num_new_resolved
*        Now see if anything has changed
         if ( .not.all_resolved ) then
             all_resolved = .true.
             do i = 1,num_ambs
                if( bf_ents(5,i).ne.bf5_copy(i) ) 
     .                            all_resolved = .false.
             enddo
         end if
*        Now save copy of bf5 for checking next iteration
         do i = 1,num_ambs
            bf5_copy(i) = bf_ents(5,i)
         enddo

      end do

****  Thats all 
      return 
      end

CITLE RES_ONE_WL

      subroutine res_one_wl(tar_b2, tar_b1, ref_b1, ref_b2,
     .                      resolved, possible, iter )

      implicit none

*     Routine to search through the solved bias flags to find
*     the next one possible of resolve the WL ambiquity and to
*     resolve it.

      include 'track_com.h'

* PASSED VARIABLES
* tar_b2 -- Target site bias flag being resolved
* tar_b1 -- Target site bias flag used as reference one for double
*           differences
* ref_b1, ref_b2 -- Corresponding BF numbers for the reference site
*           (usually  site 1).
* iter   -- Iteration number

       integer*4 tar_b1, tar_b2, ref_b1, ref_b2, iter 

* resolved -- Set true if we are able resolve WL bias
* possible -- Set true if there is still more combinations of data
*             that we may able to use.

       logical resolved, possible

* LOCAL VARIABLES
* ns  -- Site number
* lv  -- Satellite PRN trying to be resolved
* kv  -- Reference PRN for double differences
* fv  -- Final PRN we are looking for (depends if lv or kv found first)
* over_lap -- Overlap of data for combination found
* ol_start -- First epoch of overlap
* ol_end   -- Last epoch of overlap 

       integer*4 ns, lv, kv, fv, over_lap, ol_start, ol_end 


* found  -- Logical to indicate we have found combination
* found_2nd -- Logical for second PRN at first site

       logical found, found_2nd

****   Set the logicals
       resolved = .false.
       possible = .true.

****   Get the station and PRN for target site
       ns = bf_ents(1,tar_b2)
       lv = bf_ents(2,tar_b2)

****   Now start scanning up.  Find the next overlapping
*      satellite with the target PRN.  Move to the next bf flag
*      at the target site
       found = .false.
       do while ( .not.found )
          tar_b1 = tar_b1 + 1
          if( tar_b1.gt.num_ambs ) then
*             We have run out of bias flags to try, this bias is
*             not possible yet
              possible = .false.
              if( 1.eq.2 )
     .        write(*,120) lv, ns, bf_ents(3,tar_b2), bf_ents(4,tar_b2)
 120          format('BIAS FLAG PRN',i3.2,' site ',i2,' Epochs ',
     .               2I8,' Not resolvable yet')
              RETURN
          end if
*         See if this bias is correct station and resolved
          if( bf_ents(1,tar_b1).eq.ns .and. 
     .        bf_ents(5,tar_b1).ge.1  .and. 
     .        bf_ents(2,tar_b1).ne.lv        ) then
              found = .true.
              kv = bf_ents(2,tar_b1)
              ref_b1 = 0
              ref_b2 = 0
          end if
      end do

****  OK, we have found a possible satellite to use.  Now find the
*     corresponding entries for site 1.
      found = .false.
      do while ( .not.found )
         ref_b1 = ref_b1 + 1

*        See if matches the PRN we need and their is overlap of the
*        data.   Check that the site is different, that the bias
*        have been resolved and the satellite is one of the pair
*        we are looking for.
         if( bf_ents(1,ref_b1).ne.ns .and. bf_ents(5,ref_b1).ge.1 .and.
     .      (bf_ents(2,ref_b1).eq.lv.or.bf_ents(2,ref_b1).eq.kv) ) then

*            Save the PRN of the satellite we are now looking
             if( bf_ents(2,ref_b1).eq.lv ) fv = kv
             if( bf_ents(2,ref_b1).eq.kv ) fv = lv

*            OK, found one satellite; now try the next
             ref_b2 = ref_b1
             found_2nd = .false.
             do while ( .not.found_2nd )
                ref_b2 = ref_b2 + 1
*               Check if this the same reference site, that the PRN
*               matchs the final one we are looking for and the
*               ambiquity has been resolved.
                if( bf_ents(1,ref_b2).eq.bf_ents(1,ref_b1) .and.
     .              (bf_ents(2,ref_b2).eq.fv .and.
     .               bf_ents(5,ref_b2).ge.1) .or.
     .              (bf_ents(2,ref_b1).eq.fv .and.
     .              bf_ents(5,ref_b1).ge.1)  ) then

*                   Check the amount of overlap
                    ol_start = max(bf_ents(3,tar_b1),bf_ents(3,tar_b2),
     .                             bf_ents(3,ref_b1),bf_ents(3,ref_b2))
                    ol_end   = min(bf_ents(4,tar_b1),bf_ents(4,tar_b2),
     .                             bf_ents(4,ref_b1),bf_ents(4,ref_b2))
                    over_lap = ol_end - ol_start

                    if( over_lap.gt.0 ) then
                        call get_wl(tar_b2, tar_b1, ref_b2, ref_b1,
     .                              ol_start, ol_end, resolved, iter)
                        if( tar_b2.eq.-1 ) then
                            print *,'GET_WL ',tar_b2, tar_b1, ref_b2, 
     .                         ref_b1,ol_start, ol_end, resolved, iter
                        endif
* MOD TAH 070109: If We have found a double difference and the next
*                 ambiquity is a different satellite; say we found_2nd
                        if ( ref_b2.lt. num_ambs .and. 
     .                       bf_ents(2,ref_b2+1).ne.fv ) 
     .                                             found_2nd = .true.
                        
                    end if
                    if( resolved ) found_2nd = .true.
                end if

*               See we have run out of site 1 data, set the found_2nd
*               to true, so we will exit loop.  If not resolved then
*               we will continue searching.
                if( ref_b2.ge.num_ambs ) found_2nd = .true.
             end do

*            If we have resolved the bias, then set found to exit.  
*            otherwise check if we still have data
             if( resolved ) found = .true.
         end if
         if( ref_b1.ge.num_ambs ) found = .true.
 
      end do

****  Thats all
      return 
      end

CTITLE GET_WL

      subroutine get_wl( tar_b2, tar_b1, ref_b2, ref_b1, 
     .                   ol_start, ol_end, 
     .                   resolved, iter)

      implicit none 

*     Routine to determine the MW-WL and ionospheric WL and see if
*     L1-L2 ambiquities can be resolved.
      
      include '../includes/const_param.h'
      include 'track_com.h'

* PASSED VARIABLES
* tar_b2 -- Target site bias flag being resolved
* tar_b1 -- Target site bias flag used as reference one for double
*           differences
* ref_b2 -- Reference site bf_ent being used (to get site)
* ref_b1 -- Final bias to which things are done relative
* ol_start -- First epoch of overlap
* ol_end   -- Last epoch of overlap
* iter     -- Iteration number.  If iteration is negative then
*             ambiquities are not updated.

      integer*4 tar_b1, tar_b2, ref_b2, ref_b1, 
     .          ol_start, ol_end, iter 

* resolved -- Set true if we are able resolve WL bias

      logical resolved

* LOCAL VARIABLES
* i   -- Loop counter
* ns  -- Site number of target site
* ks  -- Site number of reference site
* lv  -- Satellite PRN trying to be resolved
* kv  -- Reference PRN for double differences
* ep  -- Epoch counter
* nwl -- Number of wide lane estimates
* wl_eps(max_epochs) -- Epoch number for each WL value
* min_mwl -- Minimum number of WLS estimates (starts at min_good
*     and decreases with iterations)

      integer*4 i, j, ns, ks, lv, kv, ep, nwl, wl_eps(max_epochs),
     .    min_mwl

      integer*4 PtoL  ! PRN to list location
     .,         sv    ! Satellite number

* D11(4), D12(4), D21(4), D22(4) -- The four
*     one-ways sets of values that make up the double differences
* data_dd(4) -- Double differenced data (L1,L2,P1,P2)
* mw_wl, ex_wl -- MW and Extra-wide lanes estimates
* mw_wl_all(max_epochs), ex_wl_all(max_epochs) -- Saved values
*       of the wide lanes that we can check for slips
* sum_wl, var_wl, wgh_wl  -- Sum, sum**2 and weights of MW-WL
* sum_ex, var_ex, wgh_wl  -- Sum, sum**2 and wieghts of Extra-wide lane
* weight  -- Weight given to widelane estimates.  Based on inverse of
*            sin(elev)
* mean_wl, rms_wl -- Mean and RMS of MW-WL
* mean_ex, rms_ex -- Mean and RMS of Extra-Wide lane
* NL1, NL2  -- Estimates of the number of cycles at L1 and L2.
* mw_wls_sig -- Sigma estimate for the WM-WL
* ex_wls_sig -- Sigma estimate for the Extra widelane
* dL1n       -- Difference of L1 from an integer value
* wl_tol     -- WL tolerance (set depending on iteration).
* mt_tol     -- Tolerance on how to match previous value

      real*8 D11(4), D12(4), D21(4), D22(4), weight,
     .       mw_wl_all(max_epochs), ex_wl_all(max_epochs),
     .       sum_wl, var_wl, wgh_wl, sum_ex, var_ex, wgh_ex,
     .       mean_wl, rms_wl, mean_ex, rms_ex, NL1, NL2,
     .       mw_wls_sig, ex_wls_sig, dL1n, wl_tol, mt_tol

      real*8 mw_wl(4), ex_wl(4), mw_wldd, ex_wldd  ! MW EX WLs by site,
             ! satellite and then double difference version

* elev_11, elev_12, elev_21, elev_22 -- Elevation angles (deg) of the
*     four one-way measurements

      real*8 elev_11, elev_12, elev_21, elev_22

* new_met, save_met -- Metric for judging if wide-lane double difference
*     should be saved.  Current is min(diff,0.1)/sigma
c       real*8 new_met, save_met


* OK  -- Logical to indicate that data is OK
* save_wls_all -- Logical set true when the current widelane looks
*        better than previous ones and should be saved.

      logical OK, save_wls_all, L1, L2, L3

* df  -- Data flag returned by get_ow
      integer*4 df

***** OK, loop over all of the data for this combination
*     computing the mean and rms of the MW-WL and ION-WL
      nwl = 0
      sum_wl = 0
      sum_ex = 0
      var_wl = 0
      var_ex = 0
      wgh_wl = 0
      wgh_ex = 0

      ns  = bf_ents(1,tar_b2)
      ks  = bf_ents(1,ref_b2)
      kv  = bf_ents(2,tar_b2)
      lv  = bf_ents(2,tar_b1)

      do ep = ol_start, ol_end

*        Get each of the oneway values that we need.
         call get_ow(ep,ns,lv, D11, elev_11, OK, df)
         if( OK ) call comp_wl(ep, ns,lv, D11,mw_wl(1), ex_wl(1))

         if( OK ) call get_ow(ep,  ns,kv, D12, elev_12, OK, df)
         if( OK ) call comp_wl(ep, ns,kv, D12,mw_wl(2), ex_wl(2))

         if( OK ) call get_ow(ep,  ks,lv, D21, elev_21, OK, df)
         if( OK ) call comp_wl(ep, ks,lv, D21, mw_wl(3), ex_wl(3))

         if( OK ) call get_ow(ep,  ks,kv, D22, elev_22, OK, df)
         if( OK ) call comp_wl(ep, ks,kv, D22, mw_wl(4), ex_wl(4))

*        If all is still OK, compute the MW-WL and ION
         if( OK ) then

*            Compute the weight based on the elevation angles
             weight =  1.d0/(1.d0/sin(elev_11*pi/180)**2 +
     .                    1.d0/sin(elev_12*pi/180)**2 +
     .                    1.d0/sin(elev_21*pi/180)**2 +
     .                    1.d0/sin(elev_22*pi/180)**2 )

*            Now compute the ML-WL and ion_wl (extra-wide lane)
             mw_wldd = (mw_wl(1)-mw_wl(2))-(mw_wl(3)-mw_wl(4))
             ex_wldd = (ex_wl(1)-ex_wl(2))-(ex_wl(3)-ex_wl(4))
             if( ep.ge.debug_start .and. ep.le. debug_end )
     .       write(*,180) 'GETWLDD EP ',ep, ns, ks, lv, kv,
     .                    mw_wldd, mw_wl, ex_wldd, ex_wl 
 180         format(a12,I5,' S12 ',2i3,' P12 ',2I3, 
     .             ' MW ',5F15.4,' EX ',5F14.4)
 

*            Save the WL's in an array for further processing
             nwl = nwl + 1
             mw_wl_all(nwl) = mw_wldd
             ex_wl_all(nwl) = ex_wldd
             wl_eps(nwl)    = ep

*            Sum for mean and RMS
             sum_wl = sum_wl + mw_wldd*weight
             var_wl = var_wl + mw_wldd**2*weight
             wgh_wl = wgh_wl + weight
             sum_ex = sum_ex + ex_wldd*weight
             var_ex = var_ex + ex_wldd**2*weight 
             wgh_ex = wgh_ex + weight
         end if
      end do

* MOD TAH 090114: Changed limit here to have at least min_good data
*     rather than just 1 or more
*     Set the tolerance say WL set.  Keep increasing so it will be set
      min_mwl = min_good
      wl_tol = 10
      mt_tol = 100
      if( iter.le.25 ) then
         wl_tol = 0.1*float(iter)**2
         mt_tol = 1.d-3
      elseif ( iter.lt.40 ) then
         wl_tol = 0.004*iter**2
         mt_tol = 0.0004*iter
      endif
      if( iter.gt.5 ) then
          min_mwl = max(5,min_good/iter**2)
      endif

****  Compute mean and RMS
      if( nwl.gt.0 ) then
          mean_wl = sum_wl/wgh_wl
          mean_ex = sum_ex/wgh_ex
          if( nwl.gt.1 ) then
             rms_wl  = sqrt(abs((var_wl-mean_wl**2*wgh_wl)/wgh_wl))
             rms_ex  = sqrt(abs((var_ex-mean_ex**2*wgh_wl)/wgh_wl))
          else
*            If only one value, then set large sigma
             rms_wl  = 99.d0
             rms_ex  = 99.d0
          end if

*         Compute the estimates of the number of amgiquities
*         at L1 and L2.  (Note: use of nint(mean_wl) this should
*         yield better estimates of the cycles)
C         NL1 = 77.d0/17.d0*nint(mean_wl) - 60.d0/17.d0*mean_ex
C         NL2 = 60.d0/17.d0*nint(mean_wl) - 60.d0/17.d0*mean_ex
          NL1 = fR1/(fR1-fR2)*nint(mean_wl) - fR2/(fR1-fR2)*mean_ex
          NL2 = fR2/(fR1-fR2)*nint(mean_wl) - fR2/(fR1-fR2)*mean_ex
         
          if( tar_b2.eq.-1 )
     .    write(*,220) tar_b2, kv, ns, lv, ks, nwl, mean_wl, rms_wl,  
     .                 mean_ex, rms_ex, NL1, NL2, 
     .                 wl_tol, mt_tol, 
     .                 ol_start, ol_end, min_good
 220      format('MWL',i4,4i3,i7,2(F15.3,1x,F8.3,1x),2F15.3,2f8.2,2i8,
     .           1x,i4)

      end if

*     Now see if we should try resolve the ambiquity.
C     if( nwl.gt.0 .and. iter.gt.0 ) then
      if( tar_b2.eq.-1 )
     .    print *,'NW: No Test ', tar_b2,nwl, min_mwl, min_good, 
     .                            rms_wl, iter
      if( nwl.ge.min_mwl .and. iter.gt.0 ) then

          if( tar_b2.eq.-1 )
     .    print *,'NW: Test ', tar_b2,nwl, min_mwl, min_good, rms_wl

*         Now see if can mark as reliable.  Compute the sigma of the
*         Mean.  Here we assume a 10 minute correlation time.
          mw_wls_sig = rms_wl/max(1.d0,sqrt(nwl*usr_interval/wl_tau))
          ex_wls_sig = rms_ex/max(1.d0,sqrt(nwl*usr_interval/wl_tau))
          dL1n = abs(NL1-nint(NL1))

*         Now compute if we will accept this value and fix the number
*         of cycles
          resolved = .true.
* MOD TAH 090115: Check both this widelane and the saved one since
*         the saved value from a previous iteration may not have
*         passed the resolved test.
          if( mw_wls_sig.gt.wl_tol .or.
     .        wls_sig(1,tar_b2).gt.wl_tol ) 
     .        resolved = .false.
c         if( resolved .and. 
c    .        dL1n .gt.0.1d0*(float(iter)) )    resolved = .false.
* MOD TAH 070109: Do not resolve any values on first iteration
*         unless it is to reference site.
          if ( iter.eq.1 .and. wls_ref(1,tar_b1).ne.-1 )
     .         resolved = .false.
                   
*         See if this estimate is better than any previous ones
*         we have.
          save_wls_all = .false.
          if( wls_ref(1,tar_b2).gt.0 ) then 
*             OK, we have a previous estimate.  Now see if this is
*             better
* MOD TAH 090114: Added test to make sure that number of a reasonable
*             size as well as sigma.
C             if( mw_wls_sig.le.wls_sig(1,tar_b2) .and. 
C             if( mw_wls_sig-wls_sig(1,tar_b2).le.0.d0 .and. 
C    .            nwl.ge. min_good ) then
              if( 1.eq.2 ) 
     .        print *,'MW_WLS ',mw_wls_sig,wls_sig(1,tar_b2),
     .                mw_wls_sig-wls_sig(1,tar_b2), iter

*             Compute metric: Not implemented at this time
              if( mw_wls_sig-wls_sig(1,tar_b2).le.mt_tol .and.  
     .            nwl.ge. min_good ) then
*                 OK, this one has a better sigma estimate. Save values
                  save_wls_all = .true.
                  if( 1.eq.2 )
     .            write(*,*) ' Replacing BF# ',tar_b2,' from ',
     .               wls_ref(2,tar_b2), wls_all(1,tar_b2),mw_wls_sig
              end if

c             if( tar_b2.eq.54 ) then
              if( .not.save_wls_all  .and. iter.gt.1 ) then
c                print *,'SAVE_WLS_ALL ',mw_wls_sig, wls_sig(1,tar_b2),
c     .              mw_wls_sig-wls_sig(1,tar_b2), nwl, min_good, 
c     .              save_wls_all, iter
                 if( mw_wls_sig-wls_sig(1,tar_b2).le.mt_tol ) then
c                     print *,'Setting save_wls_all ',nwl, min_good, iter
                     if( nwl.ge.min_good ) then
                        print *,'save_wls_all test failed ',iter,
     .                           tar_b2, nwl, min_good, 
     .                           mw_wls_sig-wls_sig(1,tar_b2)
                        save_wls_all = .true.
                     endif
                 endif
              endif
          else
              save_wls_all = .true.
          end if

          if( tar_b2.eq.-1 ) then
              print *,'Amb ',tar_b2,iter, save_wls_all, nwl, 
     .            resolved, tar_b1, ref_b2, ref_b1, 
     .            mw_wls_sig, wls_sig(1,tar_b2),
     .            mw_wls_sig - wls_sig(1,tar_b2)

          end if

          if( save_wls_all ) then
*             Save the residuals to wls and sigmas
              wls_all(1,tar_b2) = mean_wl - (nint(NL1)-nint(NL2))
              wls_all(2,tar_b2) = mean_ex - 
     .                           (nint(NL1)-(fR1/fR2)*nint(NL2))
              wls_sig(1,tar_b2) = mw_wls_sig
              wls_sig(2,tar_b2) = ex_wls_sig
*             Now save the double difference combination that we
*             used.
              wls_ref(1,tar_b2) = tar_b1
              wls_ref(2,tar_b2) = ref_b2
              wls_ref(3,tar_b2) = ref_b1

*             Save estimates of ambiquities: kv is target ambiquity
              sv = PtoL( kv ) 
C             ambiq_all(1,tar_b2) = ambiq_all(1,tar_b2) + nint(NL1)*
C    .                                    fR1/fL1(sv)
C             ambiq_all(2,tar_b2) = ambiq_all(2,tar_b2) + nint(NL2)*
C    .                                    fR2/fL2(sv)
              ambiq_all(1,tar_b2) = ambiq_all(1,tar_b2) + nint(NL1)
              ambiq_all(2,tar_b2) = ambiq_all(2,tar_b2) + nint(NL2)
              if( 1.eq.2 )
     .        write(*,235) kv, ns, lv, tar_b2, 
     .                        (ambiq_all(j,tar_b2),j=1,2),
     .                        (wls_all(j,tar_b2),j=1,2)
 235          format('SAVING PRN',i3.2,' Site ',i2,' with PRN',i3.2,
     .                   ' BF# ',i4,1x,2F15.1,' Resid ',2F9.3)


*             OK, we saved information.  Now see if we can resolve it
              if( resolved ) then  
                  bf_ents(5,tar_b2) = 1
                  if( 1.eq.2 )
     .            write(*,240) kv, ns, lv, tar_b2, 
     .                        (ambiq_all(j,tar_b2),j=1,2),
     .                        (wls_all(j,tar_b2),j=1,2)
 240              format('FIXING PRN',i3.2,' Site ',i2,' with PRN',i3.2,
     .                   ' BF# ',i4,1x,2F15.1,' Resid ',2F9.3)
              end if
          else
* MOD TAH 090115: If we have not saved the values; do not return resolved as 
*         true.
              resolved = .false.
          end if
      end if
 

****  OK For moment; see how it goes
      return
      end

CTITLE GET_OW

      subroutine get_ow(ep, ns, kv, data, elev, OK, df)

      implicit none 

*     Routine to return the 4 observable at epoch ep, for site ns
*     and PRN kv.  If data is OK, returned in array data


      include 'track_com.h'

* PASSED VARIABLES
* ep  -- Epoch number
* ns  -- Site number
* kv  -- PRN number 
* df  -- Data flag itself.  Initialize to zero in case no data found
 
      integer*4 ep, ns, kv, df

* data(4) -- Values of L1, L2, P1, P2 (cycles and meters)
* elev    -- Elevation angle (deg) 

      real*8 data(4), elev

* OK -- Set true if data is OK
* data_OK -- Logical function to test data

      logical OK, data_OK

* LOCAL VARIABLES
* na  -- Ambiquity number of this one-way
* i   -- Loop counters
* ch  -- Channel number 

      integer*4 na, i, ch

* data_OK    -- Logical function; returns true if data is OK

****  OK, see if have measurement and if it is OK
      OK = .false.
      df = 0
      do i = 1,4
         data(i) = 0.d0
      end do
      elev = 90.d0

      do i = 1, num_chan_se(ns,ep)
         if( ctop_cse(i,ns,ep).eq. kv ) then
             ch = ctop_cse(i,ns,ep)
*            OK, found channel see if data OK
             OK = data_OK(data_flag_cse(i,ns,ep),data_mask)
             df = data_flag_cse(i,ns,ep)

*            Get ambiquity number for these measurements
             na = amb_point_cse(i,ns,ep)
             if( na.gt.0 ) then 
                 data(1) = L1o_all_cse(i,ns,ep) + ambiq_all(1,na)
                 data(2) = L2o_all_cse(i,ns,ep) + ambiq_all(2,na)
                 data(3) = P1o_all_cse(i,ns,ep)
                 data(4) = P2o_all_cse(i,ns,ep)
                 elev    = elev_cse(i,ns,ep)
             else if( OK ) then
C                write(*,120) ep, ns, kv
C120             format('**ERROR** At epoch ',i8,' Site ',i2,
C    .                  ' PRN',i2.2,' has no ambquity pointer')
C                OK = .false. 
C                write(*,240) ns, ep, num_chan_se(ns,ep),
C    .                (ctop_cse(j,ns,ep),amb_point_cse(j,ns,ep),
C    .                 j=1,num_chan_se(ns,ep))
C240             format('SITE: ',I2,i8,i3,' Chans ',32(I2.2,1x,i3,1x))
*                Probably an initial call before ambiquities have been
*                set.  Just return the raw data
                 data(1) = L1o_all_cse(i,ns,ep)
                 data(2) = L2o_all_cse(i,ns,ep)
                 data(3) = P1o_all_cse(i,ns,ep)
                 data(4) = P2o_all_cse(i,ns,ep)
                 elev    = elev_cse(i,ns,ep)
             end if
         end if
      end do

****  Thats all
      return
      end

CTITLE COMP_WL
            
      subroutine comp_wl( ep, ns, pn, data, mw_wl, ex_wl)

      implicit none 

*     Routine to compute the MW and extra-wide lane.  We apply
*     corrections here for DCB's and iondel if use_ionex is true.
                           
      include '../includes/const_param.h'
      include 'track_com.h'

* PASSED VARIABLES

      integer*4 ns   ! Site number
     .,         pn   ! Satellite PRN number
     .,         ep   ! Epoch number

* data(4)  -- L1, L2, P1, P2 data corrected for known ambiquities
* mw_wl, ex_wl -- MW and Extra-wide lane values

      real*8 data(4), mw_wl, ex_wl

* LOCAL Variables
      real*8 dcbap(2)   ! DCB coorections (could be zero if not needed).

      real*8 data_upd(4)  ! Data with DCB and ION accounted for,
      real*8 rcv_mjd, rcv_sec  ! Receiver time MJD and seconds of day
                          ! These are adquate for the ion delay calc.
      real*8 elev, az, dnadir    ! Elevation and azimuth (deg)
      real*8 tec_los      ! LOS TEC value
      real*8 dsitL1, dsitL2, dsvsL1, dsvsL2
      real*8 IntAntMod    ! Function to compute antenna model

      real*8 range(2), phase(2)  ! Theorectical range and phase at L1/L2
      real*8 dry_map, wet_map    ! Dry and wet mapping functions
      real*8 rcv_time     ! Time of reception approximately.
      real*8 trn_time, trn_sec  ! Time at transmitter.
      real*8 svs_ear(3)   ! SVS position in earth fixed frame.


      logical ob(2)       ! True is el/nadir out of bounds.

      integer*4 ch, j     ! channel number

****  Get channel
      ch = -1
      do j = 1,num_chan_se(ns,ep)
         if( ctop_cse(j,ns,ep).eq.pn ) ch = j
      end do
!      print *,'COMP_WL ep, ns, pn, ch ',ep, ns, pn, ch
      if( ch.eq.-1 ) then
         print *,'BAD Channel in comp_wl ns,pn,ep ',ns,pn,ep
         stop 'BAD Channel in comp_wl'
      endif


* MOD TAH 130621: Just call theory to get the calc'd values and remove
*     from data to use o-minus-c
      call get_kine_apr(ep)
      rcv_time = ((ep-1)*usr_interval+
     .           sec_offset(ns,ep))/86400.d0+ref_start
      trn_time = rcv_time - data(3)/vel_light/86400.d0
      trn_sec  = ((ep-1)*usr_interval+sec_offset(ns,ep)) +
     .                 ref_sec - data(3)/vel_light
      
      call theory(trn_time, trn_sec, ns, ch,
     .            curr_site_xyz(1,ns), 
     .            svs_ear,  range, phase, elev, az, dry_map, 
     .            wet_map, 'F', ep) 

      data_upd(1) = data(1) - phase(1)
      data_upd(2) = data(2) - phase(2)
      data_upd(3) = data(3) - range(1)
      data_upd(4) = data(4) - range(2)

****  Get DCB corrections
      call get_dcb(ns, pn, dcbap) 
      data_upd(3) = data_upd(3) + dcbap(1)
      data_upd(4) = data_upd(4) + dcbap(2)

***** Get antenna phase center model
C     call get_kine_apr(ep)
C     rcv_mjd = ((ep-1)*usr_interval+sec_offset(ns,ep))/86400.d0+
C    .            ref_start
C     rcv_sec = ((ep-1)*usr_interval+sec_offset(ns,ep)) +
C    .            ref_sec

C     call eph_to_xyz( rcv_mjd, rcv_sec, pn, 'E', ep )
C     call get_elev(curr_site_xyz(1,ns), svs_xyz(1,pn), elev, az)
C     call get_nadir(curr_site_xyz(1,ns), svs_xyz(1,pn), dnadir)

C     dsitL1 = IntAntMod(num_sit_dph(1,ns), sit_dzn(1,1,ns),
C    .                   sit_dphs(1,1,1,ns), 90-elev, az, 
C    .                   max_zen, max_az, ob(1))
C     dsitL2 = IntAntMod(num_sit_dph(1,ns), sit_dzn(1,1,ns),
C    .                   sit_dphs(1,1,2,ns), 90-elev, az, 
C    .                   max_zen, max_az, ob(1))

*     Compute the satellite nadir angle
C     dsvsL1 = IntAntMod(num_svs_dph(1,pn), svs_dna(1,1,pn),
C    .                   svs_dphs(1,1,1,pn), dnadir, az,
C    .                   max_dna, 1, ob(2))
C     dsvsL2 = IntAntMod(num_svs_dph(1,pn), svs_dna(1,1,pn),
C    .                   svs_dphs(1,1,2,pn), dnadir, az,
C    .                   max_dna, 1, ob(2))

C     data_upd(1) = data_upd(1) - (dsitL1 + dsvsL1)*fR1/vel_light
C     data_upd(2) = data_upd(2) - (dsitL2 + dsvsL2)*fR2/vel_light
!      data_upd(3) = data_upd(3) - (dsitL1 + dsvsL1)
!      data_upd(4) = data_upd(4) - (dsitL2 + dsvsL2)

****  If IONEX file read; Add correction to model
C     if( use_ionex ) then
C        call interp_ionex( ep, ns, ch, curr_site_xyz(1,ns), az, elev,  
C    .                     rcv_mjd, tec_los )
*        Remove ion-delay contribution (Opposite sign to theory where
*        values are added to theory model).
C        data_upd(1) = data_upd(1) + l1tecu*tec_los*fR1/vel_light ! L1 phase
C        data_upd(2) = data_upd(2) + l2tecu*tec_los*fR2/vel_light ! L2 phase
C        data_upd(3) = data_upd(3) - l1tecu*tec_los ! L1 range
C        data_upd(4) = data_upd(4) - l2tecu*tec_los ! L2 range
C     end if

****  OK, compute the values
      mw_wl = data_upd(1) - data_upd(2) - 
     .  dfsf*(data_upd(3)*fR1/vel_light+data_upd(4)*fR2/vel_light)
      ex_wl = data_upd(1) - (fR1/fR2)*data_upd(2)

****  Thats all
      return
      end

CTITLE REPORT_WLS

      subroutine report_wls

      implicit none 

*     Routine to report the WLS status 
      
      include 'track_com.h'

* LOCAL VARIABLES
* i, j -- Loop counters

      integer*4 i,j

*     Write out the values
      write(*,120) 
 120  format('WIDELANE STATUS REPORT:',/, 
     .       '   #    Site   S  PRN  Fixd  Ref BF#s   WM-Residual',
     .       '     +-  EX-Residual  +- ')
      do i = 1, num_ambs
         write(*,140) i, site_names(bf_ents(1,i)),
     .          (bf_ents(j,i),j=1,2), bf_ents(5,i),
     .                  (wls_ref(j,i),j=1,2),
     .                  (wls_all(j,i),wls_sig(j,i), j=1,2)
 140     format(i4,1x,a4,1x,i2,1x,' PRN ',i3.2,5x,o2.2,2x,2i6, 
     .         4(F9.3,1x))
      end do

****  Thats all
      return 
      end
             
CTITLE ALL_WLS_RES

      logical function all_wls_res(ns) 

      implicit none

*     Routine to see if have resolved all the widelanes at site ns 
      
      include 'track_com.h'

* PASSED VARIABLES
* ns -- Site number

       integer*4 ns

* LOCAL VARIABLES
* i  -- Loop counters

      integer*4 i

****  Loop over all ambiquities and see if resolved
      all_wls_res = .true.
      do i = 1, num_ambs
         if( bf_ents(1,i).eq.ns .and. bf_ents(5,i).eq.0 )
     .         all_wls_res = .false.
      end do

****  Thats all
      return 
      end
             
CTITLE CHECK_ARB_BF

      subroutine check_arb_bf(ns, iter, tot_iter) 

      implicit none

*     Routine to see if we have to set a new arbitary bias flag
*     because of a gap in the data.
      
      include 'track_com.h'

* PASSED VARIABLES
* ns -- Site number being checked
* iter -- Iteration in the cleaning.
* tot_iter -- Total number of iterations

      integer*4 ns, iter, tot_iter

* LOCAL VARIABLES
* i, j  -- Loop counters
* nr    -- Number of resolved biases (needs to be more than 1 for anything
*          to happen here.  Routine is only called after not new biases
*          have been resolved.
* max_dur -- Longest duration non-overlapping bias flag
* max_ent -- Bias flag entry that is non-overlapping and has max_duration
* ol_start, ol_end -- Start and end of overlap interval

      integer*4 i,j, max_dur, max_ent, ol_start, ol_end, nr 


* over_lap -- Set true if there is overlap (in which we do not want
*     fix arbirarily 

      logical over_lap

****  OK, start by making sure more than one bias is set (ie., this
*     not early in the iteration)
      nr = 0
      do i = 1, num_ambs
         if( bf_ents(1,i).eq.ns .and. bf_ents(5,i).ge.1 ) nr = nr + 1
      end do

*     If nr (number resolved is not greater than one, just return)
      if( nr. le. 1 ) RETURN

*     OK, bias have been resolved so check the ones that have now entry.
      max_dur = 0
      do i = 1, num_ambs
         if( bf_ents(1,i).eq.ns .and. wls_ref(1,i).eq.0 ) then 

*            OK, we have a bias flags which there has been no attempt
*            to resolve.  See if there is a group which are all at the
*            same time.
             over_lap = .false.
             do j = 1, num_ambs
                if( bf_ents(1,j).eq.ns .and. wls_ref(1,i).ne.0 ) then
*                   OK, this bias flag has been estimated, see if it
*                   overlaps with the one we are checking
                    ol_start = max(bf_ents(3,i),bf_ents(3,j))
                    ol_end   = min(bf_ents(4,i),bf_ents(4,j))
                    if( ol_end.gt.ol_start ) over_lap = .true.
                end if
             end do

*            If no overlap then save the duration of the bias we just 
*            tested
             if( .not.over_lap .and. 
     .            max_dur.lt.bf_ents(4,i) - bf_ents(3,i) ) then
                 max_dur = bf_ents(4,i) - bf_ents(3,i)
                 max_ent = i
             end if
         end if
      end do

****  OK, if max_dur is greater than 0, then we have a disconnected set
*     of bias flags, so fix the one with the longest duration
      if( max_dur.gt.0 .and. iter.gt. 0 ) then 
          if( 1.eq.1 )
     .    write(*,220) max_ent, max_dur, iter, tot_iter
 220      format('NON-OVERLAPPING BIAS SET: Fixing BF# ',i4,
     .           ' Duration ', i8,' at iterations ',2(i4,1x))
          bf_ents(5,max_ent) = 1
          wls_ref(1,max_ent) = -1
          iter = 0
      end if

****  Thats all
      return
      end

CTITLE SCAN_DDLG

      subroutine scan_ddlg

      implicit none

*     Routine to scan the double difference ex-wl or LG and flag
*     data with jumps or missing data
                           
      include '../includes/const_param.h'
      include 'track_com.h'

* PASSED VARIABLES
* None

* LOCAL VARIABLES
* i,j,k,l, ep  -- Loop counters
* kv, lv       -- PRN's corresponding to sequential list of satellites
* prev_ep      -- Previous epoch number 
* num_good     -- Counter for number of good measurements
* last_good    -- Last good data point
* last_good_bf -- Last good data point with a bias flag

      integer*4 i,j,k, kv, ep,  num_good, last_good,
     .          last_good_bf

* exwl  -- Estimate of the extra wide-lane (also known in gamit as LG)
* prev_exwl -- Estimate from previous epoch
* data(4) -- Estimates of the phase and range for each of the one-ways
*              in the double difference.  The first index is over phase
*              and range (L1,L2,P1,P2); the second index is by site
* dexwl     -- Change in widelane
* elev      -- elevation angle

      real*8 exwl, prev_exwl, data(4),  elev
      real*8 mwwl, prev_mwwl

* dmw_wl -- Absolute value of change in MW-WL (may be corrected for 
*           missed milli-second jumps
* dmw_wl_ms -- Change in widelane expected due to missed milli-second 
*           time jumps
* dmv_test -- Test offset

      real*8 dmw_wl, dex_wl, dmw_wl_ms, dmw_test

* av_mwwl, av_exwl -- average value, +- 5points
* dmw, dex  -- Change in MW-wl and EX-wl

      real*8 av_mwwl, av_exwl 

* first     -- Logical to indicate first measurement from start or
*              after a gap.
* dd_OK     -- Set true is double difference should be OK
* OK        -- The status of each of the oneways.

      logical first, first_ex, first_mw,  OK, kbit


* df      -- Data flag returned by get_ow
      integer*4 df

      integer*4 date(5)
      real*8 sectag
      real*8 dcberr(2)  ! DCB contribution for this receiver/prn

      integer*4 PtoL   ! Generate the satellite list number from the PRN
     .,         lv     ! Satellite list number 

****  First flag and remove any short segments from the first
*     site 
C      i = 1
      do i = 1, num_site
         do k = 1, num_prn
            kv = prn_used(k)
            num_good = 0
            first = .true.
            do ep = 1, num_epochs
*              Check for bad data and low elevation
               data_mask = 18
               call get_ow(ep,i, kv, data,  elev, OK, df)
               if( OK ) then

*                  Check Data flag to see if bias parameter
C                  data_mask = 8
C                  call get_ow(ep,i,kv,data,elev,OK, df)
                   OK = .not.kbit(df,4)
                   if( OK ) then
*                      No bias flag on this point.  If the
*                      number of good data points is 0, then
*                      add a bias flag
                       if( num_good.eq.0 ) then
                           call rep_addbf(ep,ep,i,kv,'First BF')
                           call mark_slip(ep, i, kv, 4, 0.d0)
                           last_good_bf = ep 
                       end if
                       num_good = num_good + 1
                       last_good = ep
                       first = .false.
                   else
*                      Bias flag is set, see how many good
*                      data we have.  If we don't have enough
*                      then delete the previous data
                       if( num_good.lt.min_good .and. 
     .                     num_good.gt.0             ) then
                           call rep_addbf(last_good_bf,ep-1,i,kv,
     .                             'Editing not enough good data')
                           do j = last_good_bf, ep-1
                              call mark_slip(j,i,kv,2,-1.d0)
                           end do
                       end if

                       last_good_bf = ep
                       last_good    = ep
                       num_good = 1
                       first = .false.
                   end if
               else
*                  We have encounted a bad data point
*                  Reset the number good data
* MOD TAH 030725: Added check on gaps since last epoch.
                   if( num_good.lt.min_good .and. 
     .                 num_good.gt.0 .and. 
     .                 ep-last_good.ge.max_gap  ) then
                       call rep_addbf(last_good_bf,ep-1,i,kv,
     .                        'Editing not enough good data reset') 
                       do j = last_good_bf, ep-1
                          call mark_slip(j,i,kv,2,-1.d0)
                       end do
                   end if
                   if( ep-last_good.ge.max_gap ) then
                      num_good = 0
                   end if
               end if
            end do

*           Clean up any tail on the data
            if( num_good.lt.min_good .and. 
     .          num_good.gt.0 ) then
                call rep_addbf(last_good_bf,last_good,i,kv,
     .               'Editing not enough good data tail')
                do j = last_good_bf, last_good
                   call mark_slip(j,i,kv,2,-1.3d0)
                end do
            end if
         end do
      end do

****  Check the one-way ionospheric delay (called the exwl).
*     We check again in the double differences with a finer
*     tolerance condition. We don't want to make the tolerance 
*     to tight here because antenna rotations effect the EXWL.
* MOD TAH 030727: Added check of MWWL as well.
      data_mask = 18

* MOD TAH 080625: Set MW-WL jump due to milli-second offset
      dmw_wl_ms = dfsf*(fR1+fR2)*1.d-3
      write(*,205) 
 205  format('Checking OW Widelanes')
      do i = 1, num_site
         do k = 1, num_prn
            kv = prn_used(k)
            first_ex = .true.
            first_mw = .true.
            
            do ep = 1, num_epochs
*              Check for bad data and low elevation
               call get_ow(ep,i, kv, data,  elev, OK, df)
               call get_dcb(i, kv, dcberr) 
               if( OK ) then
                   exwl =  data(1) - (fR1/fR2)*data(2)
                   mwwl =  data(1) - data(2) -
     .                     dfsf*((data(3)+dcberr(1))*fR1/vel_light+
     .                           (data(4)+dcberr(2))*fR2/vel_light)


*                  Compute size of mw-wl jump and see if due to
*                  ms clock reset
                   dmw_wl = abs(prev_mwwl-mwwl)

                   if( abs(dmw_wl/dmw_wl_ms).gt.1.d-2 ) then
                       if( abs(dmw_wl-nint(dmw_wl/dmw_wl_ms)*dmw_wl_ms)
     .                                          .lt.mwwl_jmp ) then
* MOD TAH 090111: Do not update so that bias flag will be added
                           dmw_test = abs(dmw_wl-
     .                              nint(dmw_wl/dmw_wl_ms)*dmw_wl_ms)
* MOD TAH 091017: Re-added option so that bias will not be added.  Seems
*                 to be some cases where adding the bias flag is needed.
*                 Added new option that will set bias flag
                           if( .not. set_msec_bf ) then
                               dmw_wl = dmw_test
                               write(*,207) ep, site_names(i), kv, 
     .                            dmw_test, abs(prev_mwwl-mwwl),
     .                            nint(dmw_wl/dmw_wl_ms)*dmw_wl_ms 
 207                           format('MilliSec WL Jump Epoch ',i6,
     .                             ' Site ',a4,' SVS ',i2,' Res ',F6.2,
     .                            ' Total ',2F10.2,' cyc')
                            end if
                       endif
                   endif

*                  Now check the widelane continuity
                   if( first_ex ) then
                       first_ex = .false.
                   else if( abs(prev_exwl-exwl).gt.
     .                 2.d0*max_ion_jmp )  then
                       dex_wl = exwl-prev_exwl
                       call mjd_to_ymdhms((ep-1)*usr_interval/86400.d0
     .                      +ref_start,date, sectag)
                       write(*,210) ep, date,sectag, site_names(i), kv,
     .                              dex_wl
 210                   format('Marking OW EX-WL slip at epoch ',
     .                      i6,1x,i4,4(1x,i2.2),1x,F4.1,
     .                      ' Site ',a4,' PRN ',I3.2,' Size ',
     .                      F15.2,' Cycles')
* MOD TAH 090120: Add check of average value
                       if ( abs(dex_wl).lt.
     .                     10*max_ion_jmp )  then
                           call get_avwl( ep, i, kv, av_mwwl, av_exwl)
                           dex_wl = av_exwl
                           write(*,215) ep, date,sectag, site_names(i),
     .                             kv,  dex_wl, 2*max_ion_jmp
 215                       format('Marking AV EX-WL slip at epoch ', 
     .                          i6,1x,i4,4(1x,i2.2),1x,F4.1,
     .                         ' Site ',a4,' PRN ',I3.2,' Size ',
     .                         F15.2,' Cycles, Tol ',F6.2)
                       end if
                       if( abs(dex_wl).gt.2.d0*max_ion_jmp)    
     .                 call mark_slip(ep,i,kv,4, dex_wl)
                       first_ex = .true.   ! Reset so that value will be re-initialized
                   end if
                   if( first_mw ) then
                       first_mw = .false.
                   else if( dmw_wl.gt. mwwl_jmp )  then
                       call mjd_to_ymdhms((ep-1)*usr_interval/86400.d0
     .                      +ref_start,date, sectag)
                       write(*,230) ep, date,sectag, site_names(i), kv,
     .                              mwwl-prev_mwwl
 230                   format('Marking OW MW-WL slip at epoch ',
     .                      i6,1x,i4,4(1x,i2.2),1x,F4.1,
     .                      ' Site ',a4,' PRN ',I3.2,' Size ',
     .                      F15.2,' Cycles')

* MOD TAH 090120: Add check of average value
                       if ( dmw_wl.lt. 10*mwwl_jmp ) then
                           call get_avwl( ep, i, kv, av_mwwl, av_exwl)
                           dmw_wl = av_mwwl
                           write(*,235) ep, date,sectag, site_names(i),
     .                              kv, dmw_wl,mwwl_jmp 
 235                       format('Marking AV MW-WL slip at epoch ', 
     .                          i6,1x,i4,4(1x,i2.2),1x,F4.1,
     .                         ' Site ',a4,' PRN ',I3.2,' Size ',
     .                         F15.2,' Cycles, Tol ',F6.2)
                       end if
                       if( abs(dmw_wl).gt.mwwl_jmp)    
     .                 call mark_slip(ep,i,kv,4,dmw_wl)
                       first_mw = .true.    ! Reset to start nexy segment
                   end if
                   prev_exwl = exwl
                   prev_mwwl = mwwl
               end if
            end do
         end do
      end do


****  OK, Loop over all possible double differences.  Here we
*     are also check on bias flags.
      data_mask = 26

****  Now remove short segments of data.  We make sure that there
*     is sufficient data between bias flags.
      write(*,'(a)') 'Removing short segments of data'
      do i = 1, num_site
         do k = 1, num_prn
            kv = prn_used(k)
            num_good = 0
            last_good = -max_gap
            first = .true.
            do ep = 1, num_epochs
               data_mask = 18
               call get_ow(ep,i, kv, data,  elev, OK, df)
               if( OK ) then
C                  data_mask = 8
C                  call get_ow(ep,i,kv,data,elev,OK,df)
                   OK = .not.kbit(df,4)
                   if( OK ) then
*                      No bias flag on this point.  If the
*                      number of good data points is 0, then
*                      add a bias flag
                       if( num_good.eq.0 ) then
                           call rep_addbf(ep,ep,i,kv,'First BF2')
                           call mark_slip(ep, i, kv, 4, 0.1d0)
                           last_good_bf = ep 
                       end if
                       num_good = num_good + 1
                       last_good = ep
                       first = .false.
                   else
*                      Bias flag is set, see how many good
*                      data we have.  If we don't have enough
*                      then delete the previous data
                       if( num_good.lt.min_good .and. 
     .                     num_good.gt.0  ) then
                           call rep_addbf(last_good_bf,ep-1,i,kv,
     .                              'Editing not enough good data2')
                           do j = last_good_bf, ep-1
                              call mark_slip(j,i,kv,2,-1.1d0)
                           end do
                       end if
                       last_good_bf = ep
                       last_good    = ep
                       num_good = 1
                       first = .false.
                   end if
               else
*                  We have encounted a bad data point
*                  Reset the number good data 
* MOD TAH 030725: Added check on gaps size.
                   if( num_good.lt.min_good .and.         
     .                 num_good.gt.0  .and.
     .                 ep-last_good.ge.max_gap ) then              
                        call rep_addbf(last_good_bf,ep-1,i,kv,
     .                         'Editing not enough good data reset2')
                       do j = last_good_bf, ep-1          
                          call mark_slip(j,i,kv,2,-1.1d0) 
                       end do
                       num_good = 0                             
                   end if
                   if( ep-last_good.ge.max_gap ) then                                 
                      num_good = 0
                   end if
               end if
            end do

****        Now clean up any tails on the data
            if( num_good.lt.min_good .and. 
     .          num_good.gt.0               ) then
                call rep_addbf(last_good_bf,last_good,i,kv,
     .               'Editing not enough good data tail2')
 
                do j = last_good_bf, last_good
                   call mark_slip(j,i,kv,2,-1.2d0)
                end do
            end if

         end do
      end do  

****  Thats all
      return 
      end

CTITLE MARK_SLIP

      subroutine mark_slip( ep, ns, lv, bit, jump )

      implicit none

*     Routine to mark a cycle slip or show gap in data

      include 'track_com.h'

* PASSED VARIABLES
* ep  -- Epoch number
* ns  -- site number
* lv  -- PRN number
* bit -- Bit to set in data_flag

      integer*4 ep, ns, lv, bit

* jump  -- Size of cycle slip
      real*8 jump

* LOCAL VARIABLES 
* i   -- Loop counter

      integer*4 i

* found  -- Logical to indicate that we found the channel

      logical found

***** OK, Loop over channels and find the one to mark
      found = .false.
      do i = 1, num_chan_se(ns,ep)
         if( ctop_cse(i,ns,ep).eq.lv ) then
             if( ep.le.0 )
     .       write(*,220) ep, ns, lv, site_names(i), jump, 
     .                    data_flag_cse(i,ns,ep)
 220         format('FLAGGING GAP/CYCLE SLIP: Epoch ',i6,
     .              ' Site ',a4,' PRN ',i3.2,' Chan ',i2,
     .              ' DEX_WL ',F20.2,' cycles, DF ',i4)
             call sbit(data_flag_cse(i,ns,ep),bit,1)
             found = .true.
         end if
      end do

***** Make sure we found it
      if( .not.found ) then
C         write(*,320) ep,ns,lv
C320      format('**Error** Channel not found for Epoch ',
C    .           I6,' Site ',i2,' PRN ',I2.2)
      end if

****  Thats all
      return
      end 

CTITLE READ_BF

      subroutine read_bf

      implicit none

*     Routine to read the bias flag values in.  The layout of the file 
*     must exactly match that which is output by the program.

      include 'track_com.h'

* PASSED VARIABLES
* None

* LOCAL VARIABLES
* i, j -- Loop counters
* ibf_ents(5) -- Input bias flag reference
* iwls_ref(2) -- Input widelane references

      integer*4 i,j, ibf_ents(5), iwls_ref(3), ierr, trimlen, ii,
     .          jerr

* iambig_all(2) -- Input ambiquitity values
* iwls_all(2)   -- Input WLS values
* iwls_sig(2)   -- Input WLS sigma

      real*8 iambig_all(2), iwls_all(2), iwls_sig(2)

* match  -- Check true if input matches expectation

      logical match

* line   -- Line read from file
      character*256 line


*     Open the ambin file
      if( trimlen(ambin_file).eq.0 ) RETURN

      open(99, file=ambin_file, status='old',iostat=ierr)
      call report_error('IOSTAT',ierr,'open',ambin_file,1,'READ_BF')

*     Now loop over the file reading it
      i = 0
      do while( ierr.eq.0 )
         read(99,'(a)', iostat=ierr) line
         if( ierr.eq.0 .and. line(1:1).eq.' ' ) then

*            Increment the ambiquiuty number, and make sure it matches
*            what is expected
             i = i + 1
             read(line,140,iostat=jerr) ii,(ibf_ents(j),j=1,5), 
     .                   (iambig_all(j),j=1,2),
     .                   (iwls_ref(j),j=1,3),
     .                   (iwls_all(j),iwls_sig(j), j=1,2)

 140         format(i4,1x,4x,1x,i2,1x,5x,i3.2, 2I8,3x, o2,2F15.1,3i6,
     .              4x,4(F9.3,1x),1x,a)
c140     format(i4,1x,a4,1x,i2,1x,' PRN ',i2.2, 2I8,3x,o2.2,
c    .         2F15.1,3i6,' WL ',4(F9.3,1x),1x,a)

             call report_error('IOSTAT',jerr,'read',line,0,'READ_BF')

             match = .true.
             do j = 1, 4
                if( ibf_ents(j).ne.bf_ents(j,i)) match = .false.
             end do
             if( .not.match ) then
                write(*,*) 'Entry ',i,' of ambin_file does not match',
     .               'Input BFS ',(ibf_ents(j),j=1,5),
     .               'Expected  ',(bf_ents(j,i),j=1,5)
                stop 'Track: Ambin file does not match'
             end if

*            OK, Copy over the new values
             do j = 1,5
                bf_ents(j,i) = ibf_ents(j)
             end do

             do j = 1,2
                ambiq_all(j,i) = iambig_all(j)
                wls_all(j,i) = iwls_all(j)
                wls_sig(j,i) = iwls_sig(j)
             end do
* MOD TAH 130730: Correctly copy 3 wl_ref values
             do j = 1,3
                wls_ref(j,i) = iwls_ref(j)
             end do

         end if
      end do 

***   Recompute the wide-lanes
      call recomp_wls

      write(*,220) ambin_file(1:trimlen(ambin_file))
 220  format('AMBIGUITIES READ FROM ',a)
      call report_bf(6,'INPUT')
      if( lus.ne.6 ) then
          write(lus,220) ambin_file(1:trimlen(ambin_file))
          call report_bf(lus,'INPUT')
      end if
      close(99) 

****  Thats all
      return 
      end

CTITLE REP_ADDF

      subroutine rep_addbf(ep_start,ep_end,ns,lv, descp)

      implicit none

*     Routine to report when mark slip is called either for
*     editing or for adding a bias flag

      include 'track_com.h'

* PASSED VARIABLES
* ep_start, ep_end -- Start and stop epoch numbers
* ns  -- Site number
* lv  -- PRN number

      integer*4 ep_start, ep_end, ns, lv

* descp  -- Description of why edited or added

      character*(*) descp

****  Tell user what and why added
c     write(*,120) ep_start, ep_end, site_names(ns), lv, descp
c120  format('Marking from Ep ',i6,' to ',i6,' Site ',a4,' PRN ',
c    .       I2.2,' Reason: ',a)

      return
      end


CITLE CHECK_WL

      subroutine check_wl

      implicit none

*     Routine to check the MW-WL and ionospheric constraint
*     wide lane ambiquities after all the amibquities have
*     been resolves

      include 'track_com.h'

* PASSED VARIABLES
* ns -- Station number being resolved
* iter -- Iteration number, with each iteration we make it easier
*         to resolve the ambiquity to ensure that we all have
*         estimates by the time we are done.

      integer*4 ns, ks, iter

* LOCAL VARIABLES
* i  -- Loop counter over ambiquities
* tar_b1  -- Target site ambiquity numbers to match with one being
*            checked
* ref_b1, ref_b2  -- Reference site ambiquity numbers

      integer*4 i, j,k,l, tar_b1,  tar_b2, ref_b1, ref_b2

* lv  -- Satellite PRN trying to be resolved
* kv  -- Reference PRN for double differences
* fv  -- Final PRN we are looking for (depends if lv or kv found first)
* over_lap -- Overlap of data for combination found
* ol_start -- First epoch of overlap
* ol_end   -- Last epoch of overlap 

       integer*4 lv, kv, kp, over_lap, ol_start, ol_end 


                                                     
* resolved     -- Set true when particular ambiquity resolved
      logical resolved, possible


***** Start looping over the bias flags computing the MW-WL and
*     ionospheric constraint.
      iter = -1

      do i = 1, num_ambs
         tar_b1 = i
         if( wls_ref(1,i).gt.0        ) then

*            This is a non-arbitary fixed bias             
             resolved = .false.
             possible = .true.
             ns = bf_ents(1,i)
             lv = bf_ents(2,i)
*            Find another station and satellites
             do ks = 1, num_site
                if( ks.ne.ns ) then
                do kp = 1, num_prn
                   kv = prn_used(kp)
                   if( kv.ne.lv ) then
*                      Find ambiquities for the other
*                      satellite at this site and
*                      other site and both satellites.
                       do j = 1,num_ambs
                          if( bf_ents(1,j).eq.ns .and.
     .                        bf_ents(2,j).eq.kv ) then
                              tar_b2 = j
                              do k = 1,num_ambs
                              if( bf_ents(1,k).eq.ks .and.
     .                            bf_ents(2,k).eq.lv ) then
                                  ref_b1 = k
   
                                  do l = 1,num_ambs
                                     if( bf_ents(1,l).eq.ks .and.
     .                                   bf_ents(2,l).eq.kv ) then
                                         ref_b2 = l
             
*                   Check the amount of overlap
                    ol_start = max(bf_ents(3,tar_b1),bf_ents(3,tar_b2),
     .                             bf_ents(3,ref_b1),bf_ents(3,ref_b2))
                    ol_end   = min(bf_ents(4,tar_b1),bf_ents(4,tar_b2),
     .                             bf_ents(4,ref_b1),bf_ents(4,ref_b2))
                    over_lap = ol_end - ol_start

             
                    if( over_lap.gt.0 ) then
                        call get_wl(tar_b1, tar_b2, ref_b1, ref_b2,
     .                              ol_start, ol_end, resolved, iter)
                    end if
                                     endif
                                  end do
                              end if
                              end do
                          end if
                       end do
                   end if
                end do
                end if
            end do
         end if
      end do

****  Thats all 
      return 
      end

CTITLE SCAN_SDLG

      subroutine scan_sdlg

      implicit none

*     Routine to scan the single differences (site-site1) ex-wl or LG and flag
*     data with jumps or missing data

      include '../includes/const_param.h'
      include 'track_com.h'

* PASSED VARIABLES
* None

* LOCAL VARIABLES
* i,j,k,l, ep  -- Loop counters
* kv       -- PRN's corresponding to sequential list of satellites
* prev_ep      -- Previous epoch number 
* num_good     -- Counter for number of good measurements
* last_good    -- Last good data point
* last_good_bf -- Last good data point with a bias flag

      integer*4 i,j,k, kv, ep, num_good, last_good,
     .          last_good_bf

* exwl  -- Estimate of the extra wide-lane (also known in gamit as LG)
* prev_exwl -- Estimate from previous epoch
* data(4,2) -- Estimates of the phase and range for each of the one-ways
*              in the double difference.  The first index is over phase
*              and range (L1,L2,P1,P2); the second index is by site
* dexwl     -- Change in widelane
* elev(2)   -- Four elevation angles

      real*8 exwl, prev_exwl, data(4,2),  elev(2)


* first     -- Logical to indicate first measurement from start or
*              after a gap.
* dd_OK     -- Set true is double difference should be OK
* OK(2)     -- The status of each of the oneways.

      logical first, OK(2), kbit

* df(2)      -- Data flags returned by get_ow
      integer*4 df(2), date(5)
      real*8 sectag

****  First flag and remove any short segments from all sites
      do i = 1, num_site
         do k = 1, num_prn
            kv = prn_used(k)
            num_good = 0
            first = .true.
            do ep = 1, num_epochs
*              Check for bad data and low elevation
               data_mask = 18
               call get_ow(ep,i, kv, data,  elev(1), OK(1), df(1))
               if( OK(1) ) then

*                  Check Data flag to see if bias parameter
C                  data_mask = 8
C                  call get_ow(ep,i,kv,data,elev(1),OK(1), df(1))
                   OK(1) = .not.kbit(df(1),4)
                   if( OK(1) ) then
*                      No bias flag on this point.  If the
*                      number of good data points is 0, then
*                      add a bias flag
                       if( num_good.eq.0 ) then
                           call rep_addbf(ep,ep,i,kv,'First BF')
                           call mark_slip(ep, i, kv, 4, 0.d0)
                           last_good_bf = ep 
                       end if
                       num_good = num_good + 1
                       last_good = ep
                       first = .false.
                   else
*                      Bias flag is set, see how many good
*                      data we have.  If we don't have enough
*                      then delete the previous data
                       if( num_good.lt.min_good .and. 
     .                     num_good.gt.0             ) then
                           call rep_addbf(last_good_bf,ep-1,i,kv,
     .                             'Editing not enough good data')
                           do j = last_good_bf, ep-1
                              call mark_slip(j,i,kv,2,-1.d0)
                           end do
                       end if

                       last_good_bf = ep
                       last_good    = ep
                       num_good = 1
                       first = .false.
                   end if
               else
*                  We have encounted a bad data point
*                  Reset the number good data
* MOD TAH 030725: Added check on gap size
                   if( num_good.lt.min_good .and. 
     .                 num_good.gt.0 .and.
     .                 ep-last_good.ge.max_gap ) then
                       call rep_addbf(last_good_bf,ep-1,i,kv,
     .                        'Editing not enough good data reset') 
                       do j = last_good_bf, ep-1
                          call mark_slip(j,i,kv,2,-1.d0)
                       end do
                   end if
                   if( ep-last_good.ge.max_gap ) then
                      num_good = 0
                   end if
               end if
            end do

*           Clean up any tail on the data
            if( num_good.lt.min_good .and. 
     .          num_good.gt.0 ) then
                call rep_addbf(last_good_bf,last_good,i,kv,
     .               'Editing not enough good data tail')
                do j = last_good_bf, last_good
                   call mark_slip(j,i,kv,2,-1.3d0)
                end do
            end if
         end do
      end do

****  Check the single difference ionospheric delay (called the exwl).
*     Any bias flags on site 1 are transferred to the other sites. 
*     So this loop is done from site 2 on (referenced to site 1)
      data_mask = 18   ! Only checking bad data and elev cutoff
      do i = 2, num_site
         do k = 1, num_prn
            kv = prn_used(k)
            first = .true.
            
            do ep = 1, num_epochs
*              Check for bad data and low elevation
               call get_ow(ep,i, kv, data(1,1), elev(1), OK(1),df(1))
               call get_ow(ep,1, kv, data(1,2), elev(2), OK(2),df(2))
               if( OK(1) .and. OK(2) ) then
*                  Both one ways are good.  Now see if we should compare
*                  i.e., is there a bias flag on one of the oneways.
                   exwl =  (data(1,1)-data(1,2)) - 
     .                      (fR1/fR2)*(data(2,1) - data(2,2))

*                  If there is a bias flag on site 1, mark it also on 
*                  site i
                   if( kbit(df(2),4) ) then
                       call mark_slip(ep,i,kv,4,2.0d0)
                       call sbit(df(1),4,1)
                   end if 
                   if( first .or. kbit(df(1),4) ) then
                       first = .false.
                   else if( abs(prev_exwl-exwl).gt.
     .                      max_ion_jmp )  then
                       call mark_slip(ep,i,kv,4,
     .                              exwl-prev_exwl)
c                      write(*,210) ep, i, kv,
c    .                              exwl-prev_exwl
c210                   format('Marking SD EX-WL slip at epoch ',
c    .                      i6,' Site ',i1,' PRN ',I2.2,' Size ',
c    .                      F15.2,' Cycles')
                       call mjd_to_ymdhms((ep-1)*usr_interval/86400.d0
     .                      +ref_start,date, sectag)

                       write(*,210) ep, date,sectag, site_names(i), kv,
     .                              exwl-prev_exwl
 210                   format('Marking SD EX-WL slip at epoch ',
     .                      i6,1x,i4,4(1x,i2.2),1x,F4.1,
     .                      ' Site ',a4,' PRN ',I3.2,' Size ',
     .                      F15.2,' Cycles')

                       first = .true.
                   end if
                   prev_exwl = exwl
               end if
            end do
         end do
      end do

****  Now remove short segments of data.  We make sure that there
*     is sufficient data between bias flags.
      do i = 1, num_site
         do k = 1, num_prn
            kv = prn_used(k)
            num_good = 0
            last_good = -max_gap
            first = .true.
            do ep = 1, num_epochs
               data_mask = 18
               call get_ow(ep,i, kv, data,  elev(1), OK(1), df(1))
               if( OK(1) ) then
C                  data_mask = 8
C                  call get_ow(ep,i,kv,data,elev,OK,df)
                   OK(1) = .not.kbit(df(1),4)
                   if( OK(1) ) then
*                      No bias flag on this point.  If the
*                      number of good data points is 0, then
*                      add a bias flag
                       if( num_good.eq.0 ) then
                           call rep_addbf(ep,ep,i,kv,'First BF2')
                           call mark_slip(ep, i, kv, 4, 0.1d0)
                           last_good_bf = ep 
                       end if
                       num_good = num_good + 1
                       last_good = ep
                       first = .false.
                   else
*                      Bias flag is set, see how many good
*                      data we have.  If we don't have enough
*                      then delete the previous data
                       if( num_good.lt.min_good .and. 
     .                     num_good.gt.0  ) then
                           call rep_addbf(last_good_bf,ep-1,i,kv,
     .                              'Editing not enough good data2')
                           do j = last_good_bf, ep-1
                              call mark_slip(j,i,kv,2,-1.1d0)
                           end do
                       end if
                       last_good_bf = ep
                       last_good    = ep
                       num_good = 1
                       first = .false.
                   end if
               else
*                  We have encounted a bad data point
*                  Reset the number good data 
* MOD TAH 030725: Added check on gap size
                   if( num_good.lt.min_good .and.         
     .                 num_good.gt.0 .and.
     .                 ep-last_good.ge.max_gap ) then              
                        call rep_addbf(last_good_bf,ep-1,i,kv,
     .                         'Editing not enough good data reset2')
                       do j = last_good_bf, ep-1          
                          call mark_slip(j,i,kv,2,-1.1d0) 
                       end do
                       num_good = 0                             
                   end if
                   if( ep-last_good.ge.max_gap ) then                                 
                      num_good = 0
                   end if
               end if
            end do

****        Now clean up any tails on the data
            if( num_good.lt.min_good .and. 
     .          num_good.gt.0               ) then
                call rep_addbf(last_good_bf,last_good,i,kv,
     .               'Editing not enough good data tail2')
 
                do j = last_good_bf, last_good
                   call mark_slip(j,i,kv,2,-1.2d0)
                end do
            end if

         end do
      end do  

****  Thats all
      return 
      end

CTITLE OUTPUT_WLS

      subroutine output_wls(ns,type)

      implicit none

*     Routine to report the WL values for data checking

      include '../includes/const_param.h'
      include '../includes/xfile_def.h'
      include 'track_com.h'

* PASSED VARIABLES
      integer*4 ns   ! Station number
      character*(*) type

* LOCAL VARIABLES
      integer*4 pn, ep, nc, j ! PRN, Epoch, num channels and chan num
     .,         ierr          ! IOSTAT error
     .,         k             ! Counter over satellites
      integer*4 trimlen, lenr, na
      integer*4 df   ! Data flag from OW
      integer*4 PtoL     ! Function to get list number from PRN number

      logical file_not_open   ! True when the output file has not been
                              ! opened
      logical OK              ! Set true if oneway OK

      real*8 time   ! Time of value
     .,      mwwl, exwl, drng   ! MW and EX Widelines and range diff (cycles)
     .,      dcberr(2)   ! DCB contribution 
     .,      ow_data(4)  ! Phase and range one-way data
     .,      ow_elev     ! Elevation angle

      character*256 fn     ! Name of output file
      character*256 root
      logical raw

****  Loop over each satellite at this site
      if( type.eq.'RWL' ) then
          lenr = trimlen(rwl_root)
          root = rwl_root
          raw = .true.
      else
          lenr = trimlen(wls_root)
          root = wls_root
          raw = .false.
      end if

      print *,'PRN USE # ',num_prn, ' PRN ', prn_used(1:num_prn)
      print *,'PRN SP3 # ',num_sat, ' PRN ', prn_sp3(1:num_sat)


      do k = 1,num_prn
         pn = prn_used(k)   ! Get the PRN number from list

         file_not_open = .true.
         do ep = 1,num_epochs

*           See if we have data at this epoch
            nc =  num_chan_se(ns,ep)
* MOD TAH 200618: Updated 32I to 50I to allow for 35 Beidou satellites
            if( ep.gt.debug_start .and. ep.le.debug_end ) then
               write(*,110) ep, ns, nc, k, pn, ctop_cse(1:nc,ns,ep)
 110           format('OUTPUT_WLS EP ',2I3,' NC ',I3,' K PN ',2i4, 
     .                ' CTOP ',50I4)
            end if 
            do j = 1, nc
                if( ctop_cse(j,ns,ep).eq.pn ) then
*                  OK, we an value
                   if( file_not_open ) then
                      write(fn,120) root(1:lenr),site_names(ns),
     .                       pn
 120                  format(a,'_',a4,'_PRN',i3.3)
                      open(110,file=fn,iostat=ierr,status='unknown')
                      call report_error('IOSTAT',ierr,'open',fn,1,
     .                     'report_wls')
                      file_not_open = .false.
                      write(110,140)
 140                  format('*  EPOCH PRN          MW-WL (cyc)',
     .                       '         EX-WL (cyc)     L1-L2 Range (m)',
     .                       ' DF      MJD               Elev (deg)')
                   endif

*                  Compute the EX and MW widelanes
                   call get_dcb(ns, pn, dcberr) 

                   if( raw ) then
                       mwwl = L1o_all_cse(j,ns,ep)-L2o_all_cse(j,ns,ep)-
     .                    dfsf*((P1o_all_cse(j,ns,ep)+dcberr(1))*
     .                                                    fR1/vel_light+
     .                          (P2o_all_cse(j,ns,ep)+dcberr(2))*
     .                                                    fR2/vel_light)
                       exwl = L1o_all_cse(j,ns,ep)-
     .                         (fR1/fR2)*L2o_all_cse(j,ns,ep)
                       ow_elev = 0.0
                       na = 1
                       if( ep.gt.debug_start.and.ep.le.debug_end ) then
                          print *,'PRN ',pn, ' EP', ep, ns, j,   
     .                       data_flag_cse(j,ns,ep), ' DATA ',
     .                       L1o_all_cse(j,ns,ep),
     .                       L2o_all_cse(j,ns,ep),P1o_all_cse(j,ns,ep),
     .                       P2o_all_cse(j,ns,ep)
                       end if  
                   else
* MOD TAH 130621: Make output based on o-minus-c calcualtion of widelane

                       na = amb_point_cse(j,ns,ep)
                       if( na.gt.0 ) then 
C                          mwwl = L1o_all_cse(j,ns,ep)+ambiq_all(1,na)-
C    .                           (L2o_all_cse(j,ns,ep)+ambiq_all(2,na))-
C    .                      dfsf*((P1o_all_cse(j,ns,ep)+dcberr(1))*
C    .                                                    fR1/vel_light+
C    .                            (P2o_all_cse(j,ns,ep)+dcberr(2))*
C    .                                                    fR2/vel_light)
C                           exwl = L1o_all_cse(j,ns,ep)+ambiq_all(1,na)-
C    .                           (fR1/fR2)*
C    .                            (L2o_all_cse(j,ns,ep)+ambiq_all(2,na))
                            call get_ow(ep,ns,pn,ow_data,ow_elev,OK, df)
                            call comp_wl(ep,ns,pn,ow_data, mwwl,exwl)
                       end if  
                   endif
            
                   drng = (P1o_all_cse(j,ns,ep)-P2o_all_cse(j,ns,ep))*
     .                    (fR1/vel_light)
                   time = ref_start+(ep-1)*usr_interval/86400.d0
* Only output if data ambiguity associated with widelanes (e.g., < elev cut
*                  will not be output on finals wls but will be there for
*                  raw widelanes.
                   if ( na.gt.0 )
     .             write(110,220) ep, pn, mwwl, exwl, drng, 
     .                            data_flag_cse(j,ns,ep), time, ow_elev
 220               format(i8,1x,i3.2,1x,3F20.3,o4,F20.8,1x, F10.3)
                endif
            end do
         end do
         close(110)
       end do

****   Thats all
       return
       end


CTITLE RECOMP_WLS

      subroutine recomp_wls

      implicit none

*     Routine to loop over all widelanes in the ambiguity table and recompute

      include '../includes/const_param.h'
      include 'track_com.h'

* LOCAL
* na       -- Ambiquiuty number
* ol_start -- First epoch of overlap
* ol_end   -- Last epoch of overlap

      integer*4 na, ol_start, ol_end 


* LOCAL VARIABLES
* i   -- Loop counter
* ns  -- Site number of target site
* ks  -- Site number of reference site
* lv  -- Satellite PRN trying to be resolved
* kv  -- Reference PRN for double differences
* ep  -- Epoch counter
* nwl -- Number of wide lane estimates
* wl_eps(max_epochs) -- Epoch number for each WL value

      integer*4 i, ns, ks, lv, kv, ep, nwl, wl_eps(max_epochs)

* D11(4), D12(4), D21(4), D22(4) -- The four
*     one-ways sets of values that make up the double differences
* data_dd(4) -- Double differenced data (L1,L2,P1,P2)
* mw_wl, ex_wl -- MW and Extra-wide lanes estimates
* mw_wl_all(max_epochs), ex_wl_all(max_epochs) -- Saved values
*       of the wide lanes that we can check for slips
* sum_wl, var_wl, wgh_wl  -- Sum, sum**2 and weights of MW-WL
* sum_ex, var_ex, wgh_wl  -- Sum, sum**2 and wieghts of Extra-wide lane
* weight  -- Weight given to widelane estimates.  Based on inverse of
*            sin(elev)
* mean_wl, rms_wl -- Mean and RMS of MW-WL
* mean_ex, rms_ex -- Mean and RMS of Extra-Wide lane
* NL1, NL2  -- Estimates of the number of cycles at L1 and L2.
* mw_wls_sig -- Sigma estimate for the WM-WL
* ex_wls_sig -- Sigma estimate for the Extra widelane
* dL1n       -- Difference of L1 from an integer value

      real*8 D11(4), D12(4), D21(4), D22(4), weight,
     .       mw_wl_all(max_epochs), ex_wl_all(max_epochs),
     .       sum_wl, var_wl, wgh_wl, sum_ex, var_ex, wgh_ex,
     .       mean_wl, rms_wl, mean_ex, rms_ex, NL1, NL2,
     .       mw_wls_sig, ex_wls_sig, dL1n

      real*8 mw_wl(4), ex_wl(4), mw_wldd, ex_wldd  ! MW EX WLs by site,
             ! satellite and then double difference version

!      real*8 dcb_11(4), dcb_12(4), dcb_21(4), dcb_22(4)  ! DCB 
!              ! corrections (1:2 are zero for phase always)

* elev_11, elev_12, elev_21, elev_22 -- Elevation angles (deg) of the
*     four one-way measurements

      real*8 elev_11, elev_12, elev_21, elev_22


* OK  -- Logical to indicate that data is OK
* save_wls_all -- Logical set true when the current widelane looks
*        better than previous ones and should be saved.

      logical OK

* df  -- Data flag returned by get_ow
      integer*4 df

***** OK, loop over all ambiguities
      do na = 1, num_ambs

*        See if we meed to compute
         if( wls_ref(1,na).gt.0 ) then
*
*           Start the calculation   
            nwl = 0
            sum_wl = 0
            sum_ex = 0
            var_wl = 0
            var_ex = 0
            wgh_wl = 0
            wgh_ex = 0

            ns  = bf_ents(1,na)             ! Site 1
            ks  = bf_ents(1,wls_ref(2,na))  ! Site 2
            lv  = bf_ents(2,na)             ! PRN 1
            kv  = bf_ents(2,wls_ref(1,na))  ! PRN 2

            ol_start = max(bf_ents(3,na),bf_ents(3,wls_ref(1,na)),
     .          bf_ents(3,wls_ref(2,na)),bf_ents(3,wls_ref(3,na)))
            ol_end   = min(bf_ents(4,na),bf_ents(4,wls_ref(1,na)),
     .          bf_ents(4,wls_ref(2,na)),bf_ents(4,wls_ref(3,na)))


            do ep = ol_start, ol_end

*              Get each of the oneway values that we need.
               call get_ow(ep,ns,lv, D11, elev_11, OK, df)
               if( OK ) call comp_wl(ep, ns,lv, D11,mw_wl(1), ex_wl(1))

               if( OK ) call get_ow(ep,  ns,kv, D12, elev_12, OK, df)
               if( OK ) call comp_wl(ep, ns,kv, D12,mw_wl(2), ex_wl(2))

               if( OK ) call get_ow(ep,  ks,lv, D21, elev_21, OK, df)
               if( OK ) call comp_wl(ep, ks,lv, D21, mw_wl(3), ex_wl(3))

               if( OK ) call get_ow(ep,  ks,kv, D22, elev_22, OK, df)
               if( OK ) call comp_wl(ep, ks,kv, D22, mw_wl(4), ex_wl(4))

*              If all is still OK, compute the MW-WL and ION
               if( OK ) then

*                  For the double difference of the data
!                  do i = 1, 4
!                     data_dd(i) = (data_11(i)-data_12(i))-
!    .                             (data_21(i)-data_22(i)) +
!    .                             (dcb_11(i)-dcb_12(i))-
!    .                             (dcb_21(i)-dcb_22(i))
!                  end do

*                  Compute the weight based on the elevation angles
                   weight =  1.d0/(1.d0/sin(elev_11*pi/180)**2 +
     .                          1.d0/sin(elev_12*pi/180)**2 +
     .                          1.d0/sin(elev_21*pi/180)**2 +
     .                          1.d0/sin(elev_22*pi/180)**2 )

*                  Now compute the ML-WL and ion_wl (extra-wide lane)
                   mw_wldd = (mw_wl(1)-mw_wl(2))-(mw_wl(3)-mw_wl(4))
                   ex_wldd = (ex_wl(1)-ex_wl(2))-(ex_wl(3)-ex_wl(4))

                   if( ep.ge.debug_start .and. ep.le. debug_end )
     .             write(*,180) 'RECOMPWL EP ', ep, ns, ks, lv, kv, 
     .                           mw_wldd, mw_wl, ex_wldd, ex_wl 
 180               format(a12,I5,' S12 ',2i3,' P12 ',2I3, 
     .                   ' MW ',5F18.3,' EX ',5F18.3)

*                  Save the WL's in an array for further processing
                   nwl = nwl + 1
                   mw_wl_all(nwl) = mw_wldd
                   ex_wl_all(nwl) = ex_wldd
                   wl_eps(nwl)    = ep

*                  Sum for mean and RMS
                   sum_wl = sum_wl + mw_wldd*weight
                   var_wl = var_wl + mw_wldd**2*weight
                   wgh_wl = wgh_wl + weight
                   sum_ex = sum_ex + ex_wldd*weight
                   var_ex = var_ex + ex_wldd**2*weight 
                   wgh_ex = wgh_ex + weight
               end if
            end do

****        Compute mean and RMS
            if( nwl.gt.0 ) then
                mean_wl = sum_wl/wgh_wl
                mean_ex = sum_ex/wgh_ex
                if( nwl.gt.1 ) then
                   rms_wl  = sqrt(abs((var_wl-mean_wl**2*wgh_wl)/
     .                                                  wgh_wl))
                   rms_ex  = sqrt(abs((var_ex-mean_ex**2*wgh_wl)/
     .                                                  wgh_wl))
                else
*                  If only one value, then set large sigma
                   rms_wl  = 99.d0
                   rms_ex  = 99.d0
                end if

*               Compute the estimates of the number of amgiquities
*               at L1 and L2.  (Note: use of nint(mean_wl) this should
*               yield better estimates of the cycles)
C               NL1 = 77.d0/17.d0*nint(mean_wl) - 60.d0/17.d0*mean_ex
C               NL2 = 60.d0/17.d0*nint(mean_wl) - 60.d0/17.d0*mean_ex
                NL1 = fR1/(fR1-fR2)*nint(mean_wl)-fR2/(fR1-fR2)*mean_ex
                NL2 = fR2/(fR1-fR2)*nint(mean_wl)-fR2/(fR1-fR2)*mean_ex
                
                if( 1.eq.2 )
     .          write(*,220) na, kv, ns, lv, ks, nwl, mean_wl, rms_wl,  
     .                    mean_ex, rms_ex, NL1, NL2, 
     .                    ol_start, ol_end
 220            format('MWL',i4,4i3,i7,2(F15.3,1x,F8.3,1x),2F15.3,2i8)

            end if

*           Now see if we should try resolve the ambiquity.
            if( nwl.gt.0 ) then

*               Now see if can mark as reliable.  Compute the sigma of the
*               Mean.  Here we assume a 10 minute correlation time.
                mw_wls_sig = rms_wl/
     .                       max(1.d0,sqrt(nwl*usr_interval/wl_tau))
                ex_wls_sig = rms_ex/
     .                       max(1.d0,sqrt(nwl*usr_interval/wl_tau))
                dL1n = abs(NL1-nint(NL1))


*               Save the residuals to wls and sigmas
                wls_all(1,na) = -mean_wl 
                wls_all(2,na) = -mean_ex 
                wls_sig(1,na) = mw_wls_sig
                wls_sig(2,na) = ex_wls_sig

            end if
         end if
      end do

      if( debug_start.gt.0 ) call report_bf(6,'RECOMP')

****  OK For moment; see how it goes
      return
      end

CTITLE DEL_USR_BF

      subroutine del_usr_bf

      implicit none

*     Routine to scan data and remove any bias flags requested by the user
*     usr_delbf command


      include 'track_com.h'

* LOCAL VARIABLES

      integer*4 ns, ep, k   ! Site epoch and channel counter
      integer*4 j  ! Loop over usr_delbf entries

      logical kbit  ! Checks bit status

***** Loop over all the sites and times, and at each time check to see if any
*     bias flag needs to be removed
      if( num_rbf.gt.0 ) write(*,120)
 120  format('CHECKING BIAS FLAGS to be removed')
      do ns = 1, num_site
         do ep = 1, num_epochs

*           Scan channels to see if bias flag is set.
            do k = 1, num_chan_se(ns,ep)

*              See if bias flag on this obs
               if( kbit(data_flag_cse(k,ns,ep),4) ) then
*                  These is a bias flag.  See if it should be
*                  removed
                   do j = 1, num_rbf
*                     Check site and PRN
                      if( ss_rbf(1,j).eq.ns .and.
     .                    ss_rbf(2,j).eq. ctop_cse(k,ns,ep) ) then
*                         Site and PRN match, see if time macthes
                          if( abs(tt_rbf(j)-
     .                        ((ep-1)*usr_interval/86400.d0+ref_start))
     .                        .lt.usr_interval/(2*86400.d0) ) then
*                           TImes match:
                            call sbit(data_flag_cse(k,ns,ep),3,0)
                            call sbit(data_flag_cse(k,ns,ep),4,0)
                            write(*,240) j, ep, 
     .                         site_names(ns), ss_rbf(2,j) 
 240                        format('Deleting bias ',i3,' at Epoch ',
     .                          i5, ' site ',a4,' PRN ',i3.2) 
                           end if
                      endif
                   end do
               endif
            end do
         end do
      end do

*     Thats all
      return
      end

CTITLE GET_AVWL

      subroutine get_avwl( ep, i, kv, av_mwwl, av_exwl)

      implicit none

*     Routine to compute 5-pt average mw-wl and ex-wl

      include 'track_com.h' 
      include '../includes/const_param.h'

      integer*4 ep     ! Epoch 
     .,         i      ! site number
     .,         kv     ! satellite number

      real*8 av_mwwl, av_exwl  ! Average change on MW-wl and EX-wl

      real*8 bmwwl, bexwl  ! Before mw-wl and exwl
     .,      amwwl, aexwl  ! After mw-wl and exwl
     .,      mwwl, exwl    ! Mw and ex widelanes
     .,      data(4)       ! Phase and range data
     .,      elev          ! Elevation angle
     .,      dcberr(2)     ! DCB range contributions

      integer*4 bnum, anum  ! Before and after number of values
     .,         df          ! Data flag
     .,         k           ! Epoch conunter

      logical OK            ! True if data OK

***** Compute the average before
      bmwwl = 0
      bexwl = 0
      bnum  = 0
      k = ep
      do while ( k.gt.1 )
          k = k - 1
          call get_ow(k,i, kv, data,  elev, OK, df)
          call get_dcb(i, kv, dcberr) 
          if( OK ) then
              exwl =  data(1) - (fR1/fR2)*data(2)
              mwwl =  data(1) - data(2) -
     .                dfsf*((data(3)+dcberr(1))*fR1/vel_light+
     .                      (data(4)+dcberr(2))*fR2/vel_light)
              bmwwl = bmwwl + mwwl
              bexwl = bexwl + exwl
              bnum  = bnum + 1
          end if
          if( bnum.gt.5 ) k = 0
      end do
      if( bnum.gt.0 ) then
          bmwwl = bmwwl/bnum
          bexwl = bexwl/bnum
      endif

*     Now do the after average values
      amwwl = 0
      aexwl = 0
      anum  = 0

      do k = ep, min(ep+5,num_epochs)
          call get_ow(k,i, kv, data,  elev, OK, df)
          call get_dcb(i, kv, dcberr) 

          if( OK ) then
              exwl =  data(1) - (fR1/fR2)*data(2)
              mwwl =  data(1) - data(2) -
     .             dfsf*((data(3)+dcberr(1))*fR1/vel_light+
     .                   (data(4)+dcberr(2))*fR2/vel_light)
              amwwl = amwwl + mwwl
              aexwl = aexwl + exwl
              anum  = anum + 1
          end if
      end do
      if( anum.gt.0 ) then
          amwwl = amwwl/anum
          aexwl = aexwl/anum
      endif

***** Compute change
      av_mwwl = amwwl - bmwwl
      av_exwl = aexwl - bexwl

***** Thats all
      return 
      end

CTITLE REMAP_BF

      subroutine remap_bf(bf_name, remap_iter, made_updates)

      implicit none

*     Routine that will test all the 'O' settings (bf_enf(5,na) bit 4) 
*     and remap them to ambqiguities that have been fixed already.
*     Invoked when no new ambiguities can be resolved.

      include 'track_com.h' 
      include '../includes/const_param.h'

* PASSED 
      character*(*) bf_name  ! FLOAT Solution name
      integer*4 remap_iter   ! Count used to test number of data neeeded.
      logical made_updates   ! Set true when updates are made.  If none
                             ! made due to short spans of data, system is
                             ! is iterated.

* LOCAL
      integer*4 best_comb(3) ! Three other entries that give best WL VAR
     .,   best_num           ! Best number of values
     .,   num_wls            ! Number of widelanes
     .,   i, j, ri, rj       ! Ambiguity numbers for 4-dd in WLS
     .,   ol_start, ol_end   ! Start and end of overlap of data
     .,   ns, lv             ! Site and PRN number 
     .,   min_num            ! minimum number of WLS

      real*8 best_var        ! Best WL var computed from RMS values
     .,   rms_wls(2)         ! RMS of MW and EX WLS
     .,   wl_var             ! Variance computed for specific combination
                             ! EW and MW are weighted (wl_fact, lg_fact)
      logical kbit           ! Tests if bit is set
      logical none_left      ! Set true when there are no more O ambiguities
                             ! (causes exit from iterations)
 

****  Loop over all ambiguities and test ones that could not be 
*     resolved due to other combinaion not being resolved.
      none_left = .true.
      min_num = max(1,num_epochs/(2**(remap_iter+1))) 
      do i = 1, num_ambs
*        See if we need to change
         ns = bf_ents(1,i)
         lv = bf_ents(2,i)
         if( kbit(bf_ents(5,i),4) ) then
*            This ambiquity could not bre resolved due to other
*            Find a better combination
             none_left = .false.  ! Found an 'O' ambiguity
             best_var = 1.d20
             best_comb(:) = 0
             do j = 1, num_ambs
*               Make sure correct station and fixed
                if( bf_ents(1,j).eq.ns .and.
     .              kbit(bf_ents(5,j),2) ) then  ! Bit 2: Bias fixed

*                  See if we have over lap
                   ol_start = max(bf_ents(3,i),bf_ents(3,j))
                   ol_end   = min(bf_ents(4,i),bf_ents(4,j))
                   if( ol_end.ge.ol_start ) then
*                     OK, this bias can be tested. Fine other 
*                     combinations to use
                      do ri = 1, num_ambs
                         if( bf_ents(1,ri).ne.ns .and.
     .                      bf_ents(2,ri).eq.lv .and.
     .                      kbit(bf_ents(5,ri),2) ) then
*                           This one is OK for second station, same satellite
*                           List here is: i -- ambiquity to re-ordered
*                                         j -- second satellite at site (1st entry)
*                                        ri -- Different site same SV as i (2nd entry)
*                                        rj -- Different site, same SV as j (3rd entry)
                            do rj = 1, num_ambs
                               if( bf_ents(1,rj).ne.ns .and.           ! Different site
     .                            bf_ents(1,rj).eq.bf_ents(1,ri) .and. ! 3rd is same site as 2nd
     .                            bf_ents(2,rj).eq.bf_ents(2,j) .and.  ! 3rd entry and SV as 1st
     .                            kbit(bf_ents(5,rj),2) ) then
*                                 See if we have over lap
                                  ol_start = max(bf_ents(3,i),
     .                                bf_ents(3,j),bf_ents(3,ri),
     .                                bf_ents(3,rj))
                                  ol_end   = min(bf_ents(4,i),
     .                                bf_ents(4,j),bf_ents(4,ri),
     .                                bf_ents(4,rj))
*                                 If we have over lap, compute the WL RMS
                                  if( ol_end.ge.ol_start ) then
                                     call comp_wlrms( i,j, ri, rj,  
     .                                  ol_start, ol_end,  num_wls, 
     .                                  rms_wls)
*                                    See how var looks
                                     if( num_wls.gt.0 ) then
                                         wl_var = (wl_fact*
     .                                        rms_wls(1)**2+
     .                                        lg_fact*
     .                                        rms_wls(2)**2)/num_wls
                                         if( wl_var.lt.best_var ) then
                                             best_var = wl_var
                                             best_comb(1) = j
                                             best_comb(2) = ri
                                             best_comb(3) = rj
                                             best_num = num_wls
                                         endif
                                     end if  
                                  end if
                                end if
                            end do
                         end if
                      end do
                   end if
*               See if we have passed our station range (if so exit)
                else if( bf_ents(1,j).gt.ns ) then
                   exit
                end if
             end do
*            Now see what the best combination was.  Check that span
*            is not too small
             if( best_num.gt.0 .and. min_num.gt.min_good )
     .            num_tot_resolved = num_tot_resolved + 1 
             if( best_comb(1).ne.0 .and. best_num.ge.min_num ) then
                 write(*,420) i, site_names(bf_ents(1,i)), bf_ents(2,i),
     .               site_names(bf_ents(1,best_comb(3))),
     .               bf_ents(2,best_comb(3)), best_comb, 
     .               best_num, sqrt(best_var)
 420             format('BF_UPDATE: Updating BF ',i4,1x,a4,
     .                  ' PRN ',I3.2,' with ',a4,' PRN ',i3.2,
     .                  ' BFS ',3I5,' Best Num ',i6,' RMS ',F8.3,' cyc')
*                Update entries
                 wls_ref(1,i) = best_comb(1)                 
                 wls_ref(2,i) = best_comb(2)                 
                 wls_ref(3,i) = best_comb(3) 
                 made_updates = .true.
                 num_tot_resolved = num_tot_resolved + 1 ! Set to kept iteratipn
             else
                 write(*,440) i, site_names(bf_ents(1,i)), bf_ents(2,i),
     .               best_num, min_num, remap_iter
 440             format('BF_UPDATE: No Update for ',i4,1x,a4,
     .                  ' PRN ',I3.2,' Num WLS ',I7,' Min ',i7,
     .                  ' Iter ',i3)
                 if( min_num.ge.min_good )
     .               num_tot_resolved = num_tot_resolved + 1  ! Increase number
                         ! so that we continue iterating float solution
             endif

          end if   ! If this one need updating
      end do       ! Looping over all ambiguities

****  If there were no left; set num_tot_resolved to 0.  We are done
      if( none_left .or. min_num.lt.min_good ) then
         num_tot_resolved = 0
         made_updates = .true.
         if( .not.none_left )
     .   write(*,520) min_num, min_good
 520     format('Ending Remapping biases with ',i5,' min number ',
     .          ' with min_good ',I5)
      endif


****  Thats all 
      return
      end 

                                   
CTITLE COMP_WLRMS

      subroutine comp_wlrms( tar_b2, tar_b1, ref_b2, ref_b1, 
     .          ol_start, ol_end, num_wls, rms_wls)

      implicit none 

*     Routine the compute the RMS scatter of MW and EX widelanes
*     for the bf_combinations of tar_b2, tar_b1 with reference
*     biases ref_b2, ref_b1
                           
      include '../includes/const_param.h'
      include 'track_com.h'

* PASSED VARIABLES
* tar_b2 -- Target site bias flag being resolved
* tar_b1 -- Target site bias flag used as reference one for double
*           differences
* ref_b2 -- Reference site bf_ent being used (to get site)
* ref_b1 -- Final bias to which things are done relative
* num_wls -- Number of value (if no overlap then zero).
* ol_start -- First epoch of overlap
* ol_end   -- Last epoch of overlap

      integer*4 tar_b1, tar_b2, ref_b2, ref_b1 
     .,    num_wls, ol_start, ol_end 

      real*8 rms_wls(2) ! RMS scatter of MW and EX

* LOCAL VARIABLES
* i   -- Loop counter
* ns  -- Site number of target site
* ks  -- Site number of reference site
* lv  -- Satellite PRN trying to be resolved
* kv  -- Reference PRN for double differences
* ep  -- Epoch counter
* nwl -- Number of wide lane estimates

      integer*4 i, j, ns, ks, lv, kv, ep, nwl

* D11(4), D12(4), D21(4), D22(4) -- The four
*     one-ways sets of values that make up the double differences
* data_dd(4) -- Double differenced data (L1,L2,P1,P2)
* mw_wl, ex_wl -- MW and Extra-wide lanes estimates
* sum_wl, var_wl, wgh_wl  -- Sum, sum**2 and weights of MW-WL
* sum_ex, var_ex, wgh_wl  -- Sum, sum**2 and wieghts of Extra-wide lane
* weight  -- Weight given to widelane estimates.  Based on inverse of
*            sin(elev)
* mean_wl, rms_wl -- Mean and RMS of MW-WL
* mean_ex, rms_ex -- Mean and RMS of Extra-Wide lane

      real*8 D11(4), D12(4), D21(4), D22(4), weight,
     .       sum_wl, var_wl, wgh_wl, sum_ex, var_ex, wgh_ex,
     .       mean_wl, mean_ex
 
      real*8 mw_wl(4), ex_wl(4), mw_wldd, ex_wldd  ! MW EX WLs by site,
             ! satellite and then double difference version

* elev_11, elev_12, elev_21, elev_22 -- Elevation angles (deg) of the
*     four one-way measurements

      real*8 elev_11, elev_12, elev_21, elev_22

* new_met, save_met -- Metric for judging if wide-lane double difference
*     should be saved.  Current is min(diff,0.1)/sigma
c       real*8 new_met, save_met


* OK  -- Logical to indicate that data is OK

      logical OK

* df  -- Data flag returned by get_ow
      integer*4 df

***** OK, loop over all of the data for this combination
*     computing the mean and rms of the MW-WL and ION-WL
      nwl = 0
      sum_wl = 0
      sum_ex = 0
      var_wl = 0
      var_ex = 0
      wgh_wl = 0
      wgh_ex = 0

      ns  = bf_ents(1,tar_b2)
      ks  = bf_ents(1,ref_b2)
      kv  = bf_ents(2,tar_b2)
      lv  = bf_ents(2,tar_b1)

      do ep = ol_start, ol_end

*        Get each of the oneway values that we need.
         call get_ow(ep,ns,lv, D11, elev_11, OK, df)
         if( OK ) call comp_wl(ep, ns,lv, D11,mw_wl(1), ex_wl(1))

         if( OK ) call get_ow(ep,  ns,kv, D12, elev_12, OK, df)
         if( OK ) call comp_wl(ep, ns,kv, D11,mw_wl(2), ex_wl(2))

         if( OK ) call get_ow(ep,  ks,lv, D21, elev_21, OK, df)
         if( OK ) call comp_wl(ep, ks,lv, D21, mw_wl(3), ex_wl(3))

         if( OK ) call get_ow(ep,  ks,kv, D22, elev_22, OK, df)
         if( OK ) call comp_wl(ep, ks,kv, D22, mw_wl(4), ex_wl(4))

*        If all is still OK, compute the MW-WL and ION
         if( OK ) then

*            For the double difference of the data
!            do i = 1, 4
!               data_dd(i) = (data_11(i)-data_12(i))-
!    .                       (data_21(i)-data_22(i)) +
!    .                       (dcb_11(i)-dcb_12(i))-
!    .                       (dcb_21(i)-dcb_22(i))
!            end do

*            Compute the weight based on the elevation angles
             weight =  1.d0/(1.d0/sin(elev_11*pi/180)**2 +
     .                    1.d0/sin(elev_12*pi/180)**2 +
     .                    1.d0/sin(elev_21*pi/180)**2 +
     .                    1.d0/sin(elev_22*pi/180)**2 )

*            Now compute the ML-WL and ion_wl (extra-wide lane)
             mw_wldd = (mw_wl(1)-mw_wl(2))-(mw_wl(3)-mw_wl(4))
             ex_wldd = (ex_wl(1)-ex_wl(2))-(ex_wl(3)-ex_wl(4))


*            Save the WL's in an array for further processing
             nwl = nwl + 1

*            Sum for mean and RMS
             sum_wl = sum_wl + mw_wldd*weight
             var_wl = var_wl + mw_wldd**2*weight
             wgh_wl = wgh_wl + weight
             sum_ex = sum_ex + ex_wldd*weight
             var_ex = var_ex + ex_wldd**2*weight 
             wgh_ex = wgh_ex + weight
         end if
      end do

****  OK, if we have values compute the RMS values
      num_wls = nwl 
      if( nwl.gt.1 ) then
          mean_wl = sum_wl/wgh_wl
          mean_ex = sum_ex/wgh_ex
          rms_wls(1)  = sqrt(abs((var_wl-mean_wl**2*wgh_wl)/wgh_wl))
          rms_wls(2)  = sqrt(abs((var_ex-mean_ex**2*wgh_wl)/wgh_wl))
      else
*        If only one value, then set large sigma
         rms_wls(1)  = 99.d0
         rms_wls(2)  = 99.d0
      end if
****  Thats all
      return
      end

CTITLE GET_DCB 

      subroutine get_dcb(ns, pn, dcberr) 

      implicit none

*     Routine to return DCB values for site kn, satellite pn.

      include 'track_com.h'

      integer*4 ns, pn  ! Site and PRN number
      real*8 dcberr(2)  ! P1/P2 dcb contibutions (m)

      dcberr(1) = 0   ! For debug output 
      dcberr(2) = 0
      if( rcv_type(ns)(1:1).eq.'C' ) then
          dcberr(1) = dcbs(pn)
      elseif( rcv_type(ns)(1:1).eq.'P' ) then
          dcberr(1) = dcbs(pn)
          dcberr(2) = dcbs(pn)
      end if

****  Thats all
      return
      end

     


