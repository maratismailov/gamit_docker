CTITLE TRIM_ONEWAYS
 
      subroutine trim_oneways( ctol_cse, data_flag_cse, out )

      implicit none
 
*     Routine to trim the one-way data by flagging all data with too
*     short a gap or not enough data between bias flags and at
*     end of the stream of data for a particular satellite.
* MOD TAH 051227: Introduced output option so that routine can be called 
*     more frequently during cleaning.
 
* INCLUDES
 
      include '../includes/kalman_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
 
*   data_flag_cse(num_chan, num_cfiles, num_ep) - Data flag for each
*                   - measurement (same definition as in Gobs File)
*   ctol_cse(num_chan, num_cfiles, num_ep)  - Conversion from
*                   - channel number to satellite number
 
 
      integer*4 data_flag_cse(num_chan, num_cfiles, num_ep),
     .    ctol_cse(num_chan, num_cfiles, num_ep)

*   out  -- Set to Y to output summary
      character*(*) out
 
* LOCAL VARIABLES

*   i,j,k       - Loop counters
*   dnum_good   - Number of good data since last bias flag
*   last_good_ep    - EPoch number of the last good observation
*   ch          - Channel number for current satellite at
*               - current epoch.  (-1 if not observed)
*   bias_deleted_mask	- Mask which does not check if point
*                 hase been flagged with bit 6 (no enough data between
*                 bias flags.  (This bit is set here and during clean_dd)
 
      integer*4 i,j, bias_deleted_mask
 
*   edits       - Set true in trim_ss whenever there is edited
*                 data.
 
      logical edits
 
 
***** Set the standard phase masks (these mask ignore low-el
*     data.
 
      call set_phs_mask( phs_mask, phs_bias_mask )

*     setthe bias_deleted mask to ignore bit 6 and bit 5
      bias_deleted_mask = phs_mask
      call sbit(bias_deleted_mask,6,0)
      call sbit(bias_deleted_mask,5,0)
 
      write(*,100)
 100  format(/' TRIM_ONEWAYS starting -- cleaning up tails etc.')
 
****  Now loop over all stations and satellites
 
	
      do i = 1, num_cfiles
          do j = 1, num_sat

              edits = .true.
              do while ( edits )

*                 Now trim the one-way sequence.  Keep editing until
*                 all data pass the test.
                  call trim_ss(i,j, ctol_cse, data_flag_cse, 
     .                         bias_deleted_mask, edits)
              end do
          end do
 
*         Give stats for stations
          if( i.eq.1 .and. out(1:1).eq.'Y' ) then
*             write the header
              write(uns,300) 
 300          format(/,' DATA AMOUNTS (Good: # good data;',
     .                 ' Gap: # deleted in gaps; BF: # bias flags ',
     .                 ' < 2*max separation)',/,
     .                 ' SITE ',4(' PRN Good  Gap  BF '))
          end if
          if( out(1:1).eq.'Y' )
     .    write(uns,310) cf_codes(i), (prn_list(j),num_good(j,i), 
     .                  num_deleted(j,i),
     .                  num_bias_flag(j,i), j = 1, num_sat)
 310      format(' ',a4,9(4(:,' PN',i2.2,2i5,1x,i3),:,/,5x))
      end do

****  Thats all
      return
      end
 
 
CTITLE FLAG_DATA
 
      subroutine flag_data( ctol_cse, data_flag_cse, start, stop,
     .                    ns, nl, bit )

      implicit none
 
*     Subroutine the flag data by setting bit 'bit' in the data_flag
*     array.
 
* INCLUDES
 
      include '../includes/kalman_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
 
*   data_flag_cse(num_chan, num_cfiles, num_ep) - Data flag for each
*                   - measurement (same definition as in Gobs File)
*   ctol_cse(num_chan, num_cfiles, num_ep)  - Conversion from
*                   - channel number to satellite number
*   start, stop     - Start and stop epoch numbers
*   ns, nl          - site and satellite list number
*   bit             - bit to be set in data_flag
 
 
      integer*4 data_flag_cse(num_chan, num_cfiles, num_ep),
     .    ctol_cse(num_chan, num_cfiles, num_ep),
     .    start, stop, ns, nl, bit
 
* LOCAL VARIABLES
 
*   i,j,k       - Loop counters
*   ltoc        - function to return the channel number
*               - for a particular observation.  If the
*               - satellite is not being oberved then -1
*               - is returned
*   ch          - Channel number for current satellite at
*               - current epoch.  (-1 if not observed)
*   dn          - Number of data deleted (just for debug)
 
      integer*4 i, ltoc, ch, dn

*   data_OK     - Function returns true if data OK
      logical data_OK
 
***** Check to see if valid start epoch, if not probably no
*     data so just return
      if( start.le.0 ) RETURN
 
*     Loop over the epoch range
      dn = 0
      do i = start, stop
 
*         Get the channel number of satellite to be flagd
          ch = ltoc(ctol_cse(1,ns,i), nl, actual_max_chan )
         
          if( ch.gt.0 .and. 
     .        data_OK(data_flag_cse(ch,ns,i),0,phs_mask) ) then
              call sbit(data_flag_cse(ch,ns,i), bit, 1 )
              dn = dn + 1
          end if
      end do
      if( 1.eq.2 )
     .write(*,900) cf_codes(ns), prn_list(nl), start, stop, dn
 900  format(' Flagged data at site ',a4,' PRN ',i2.2,' from epoch ',
     .        i4,' to ',i4,'. ',i4,' data deleted')
 
 
****  Thats all
      return
      end
 
CTITLE TRIM_SS 
 
      subroutine trim_ss(i,j, ctol_cse, data_flag_cse, 
     .                   bias_deleted_mask, edits )

      implicit none
 
*     Routine to trim the one-way data by flagging all data with too
*     short a gap or not enough data between bias flags and at
*     end of the stream of data for a particular satellite.
 
* INCLUDES
 
      include '../includes/kalman_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
 
*   data_flag_cse(num_chan, num_cfiles, num_ep) - Data flag for each
*                   - measurement (same definition as in Gobs File)
*   ctol_cse(num_chan, num_cfiles, num_ep)  - Conversion from
*                   - channel number to satellite number
*   i,j    - Site and satellite number being trimmed.
 
 
      integer*4 data_flag_cse(num_chan, num_cfiles, num_ep),
     .    ctol_cse(num_chan, num_cfiles, num_ep), i, j

*   edits  - Set true if there have been edites to the data.

      logical edits
 
* LOCAL VARIABLES

*   i,j,k       - Loop counters
*   ltoc        - function to return the channel number
*               - for a particular observation.  If the
*               - satellite is not being oberved then -1
*               - is returned
*   dnum_good   - Number of good data since last bias flag
*   prev_bias_ep    - Epoch number of previous bias flag
*   last_good_ep    - EPoch number of the last good observation
*   ch          - Channel number for current satellite at
*               - current epoch.  (-1 if not observed)
*   bias_deleted_mask	- Mask which does not check if point
*                 hase been flagged with bit 6 (no enough data between
*                 bias flags.  (This bit is set here and during clean_dd)
 
      integer*4 k, ltoc, dnum_good, prev_bias_ep, last_good_ep, ch,
     .          bias_deleted_mask
 
*   dtl_bias        - Time difference since last bias flag
*   dtr_end         - Ratio of time difference since last bias
*                   - and total duration of data.
 
      real*8 dtl_bias, dtr_end
 
*   data_ok     - Logical function which returns true if
*               - if none of the bits in the phase mask are
*               - set.
*   first_bias  - Logical to indicate that this is first bias
*               - flag.
*   kbit        - Checks if bit set
 
 
      logical data_ok, first_bias, kbit

*   mess        - Message for warning file
      character*80 mess
 
 
****  Now loop over all stations and satellites
 
      first_bias = .true.
      edits = .false.
      num_good(j,i) = 0
      num_deleted(j,i) = 0
      dnum_good = 0
      last_good_ep = 1

*     Now loop over epochs
      do k = 1, num_ep
          ch = ltoc( ctol_cse(1,i,k), j, actual_max_chan)
*                          !  See if data OK.
          if( ch.gt.0 ) then
              if( data_OK(data_flag_cse(ch,i,k),0,
     .                    phs_mask) )  then
 
*                see if bias flag
                 dnum_good = dnum_good + 1
                 if((kbit(data_flag_cse(ch,i,k),31) .or.
     .                kbit(data_flag_cse(ch,i,k),32)) .and.
     .               .not.first_bias ) then
 
*                   See how long since fast bias flag
                    dtl_bias = (k-prev_bias_ep)*sampling
 
*                   see if passes.  
* MOD TAH 950828: Increased tolerance by 150% to get get rid of
*                 close bias flags
                    if( (dtl_bias.lt. min_dtl_bias*2.0 .or.
     .                  dnum_good.lt. min_good_bias*2.0) .and.
     .                  dnum_good.gt.0  ) then
 
*                       Flag this data as not to be used.
                        edits = .true.
                        call flag_data(ctol_cse,
     .                            data_flag_cse, prev_bias_ep,
     .                             k-1, i, j , 6)

                        num_deleted(j,i) = num_deleted(j,i) + 
     .                                     dnum_good 
*                       By setting dnum_good here to zero, this
*                       will stop data in counted in the total
*                       number of good data.
                        dnum_good = 0
                    end if
 
*                   save the epoch number of this bias
                    prev_bias_ep = k
*                   Reset the number of good data since last bias
*                   after we update the number of data
                    num_good(j,i) = num_good(j,i) + dnum_good
                    dnum_good = 1
                 end if
 
****             If this is the first bias, then initialize
                 if( first_bias ) then
                    first_bias = .false.
                    prev_bias_ep = k
                 end if
                 last_good_ep = k
              end if
          end if
      end do
 
*     Check for the end of data tolerances.  Get total number
*     good data
      num_good(j,i) = num_good(j,i) + dnum_good 
      if( num_good(j,i).gt.0 ) then
          dtr_end = float(last_good_ep-prev_bias_ep)/
     .                    num_good(j,i)
      else
          dtr_end = 1
      end if
 
*     See if passes.
      if( (dtr_end.lt. min_dtr_end .or.
     .    dnum_good.lt.min_good_end).and.
     .    dnum_good.gt.0  )  then
          edits = .true. 
          call flag_data( ctol_cse, data_flag_cse,
     .                    prev_bias_ep, last_good_ep,
     .                    i, j , 5)
          num_deleted(j,i) = num_deleted(j,i)+dnum_good+1
          num_good(j,i) = num_good(j,i) - dnum_good
      end if

****  Now count up the number of bias flags in the data.  Set value
*     -1 first to allow for the bias flag on the first data point
*     (should not be counted, all data have this bias flag)
      num_bias_flag(j,i) =  0
      num_good(j,i) = 0
      num_deleted(j,i) = 0
      last_good_ep = -1000
      do k = 1, num_ep
          ch = ltoc( ctol_cse(1,i,k), j, actual_max_chan)
*                          !  See if data OK.
          if( ch.gt.0 ) then
              if( data_OK(data_flag_cse(ch,i,k),0,
     .                    phs_mask) )  then

*                See if bias flag Only count is separation
*                is less tham twice that allowed for bias flag
*                to be removed.
                 if((kbit(data_flag_cse(ch,i,k),31) .or.
     .                kbit(data_flag_cse(ch,i,k),32)) .and.
     .                (k-last_good_ep)*sampling.le.
     .                      2* dchi2_max_sep) then
                     num_bias_flag(j,i) = num_bias_flag(j,i) + 1
                 end if
                 num_good(j,i) = num_good(j,i) + 1
                 last_good_ep = k
              else

*                 see if we deleted due to too little data i.e.,
*                 would be good except for bit 6
                  if( data_OK(data_flag_cse(ch,i,k),0, 
     .                        bias_deleted_mask) ) then
                      num_deleted(j,i) = num_deleted(j,i)+1
                  end if 
              end if
               
          end if
      end do

* MOD TAH 050217: Final check on data to see that min_ow_data is 
*     OK
      if( num_good(j,i).le. min_ow_data .and. 
     .    num_good(j,i).gt. 0                 ) then
*        Not enough good one-way data so delete whole sequence
         write(*,310) cf_codes(i), prn_list(j), num_good(j,i)
 310     format('REMOVING ',a4,' PRN ',I2.2,': Only ',i4,
     .          ' good one-way data')
         write(mess,310) cf_codes(i), prn_list(j), num_good(j,i)
         call report_stat('status','autcln','trim_oneways',' ',mess,0)
         do k = 1, num_ep
             ch = ltoc( ctol_cse(1,i,k), j, actual_max_chan)
*                          !  See if data OK.
             if( ch.gt.0 ) then
                call sbit(data_flag_cse(ch,i,k),19,1)
                num_deleted(j,i) = num_deleted(j,i)+1
             end if
         end do
         num_good(j,i) = 0
      end if

 
****  Thats all
      return
      end
 
 
