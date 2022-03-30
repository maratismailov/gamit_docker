CTITLE TRIM_SHORTSEG
 
      subroutine trim_shortseg( ctol_cse, data_flag_cse, out, luo )

      implicit none
 
*     Routine to remove small blocks of data with gaps in between
*     the data.  Intial version any single missing data point
*     is considered to have started a gap.  Bit 6 will be set
*     for these data.
 
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
      integer*4 luo   ! Output unit when out == Y (Added 200612).

*   out  -- Set to Y to output summary
      character*(*) out
 
* LOCAL VARIABLES

*   i,j       - Loop counters over site and satellite.

      integer*4 i,j

 
***** Set the standard phase masks (these mask ignore low-el
*     data.
 
      call set_phs_mask( phs_mask, phs_bias_mask )

      write(*,100)
 100  format(/' Trim_shortseg starting')
 
****  Now loop over all stations and satellites
 
	
      do i = 1, num_cfiles
          do j = 1, num_sat
*            Now trim the one-way sequence.  Keep editing until
*            all data pass the test.
             call trim_gap(i,j, ctol_cse, data_flag_cse, 
     .                    phs_bias_mask)
          end do
 
*         Give stats for stations
          if( i.eq.1 .and. out(1:1).eq.'Y' ) then
*             write the header
              write(luo,300) ('PRN   Gap  ',j=1,num_sat)
 300          format(' DATA AMOUNTS: Gap: # deleted in gaps',/,
     .               ' SITE ',50a)
          end if
          if( out(1:1).eq.'Y' )
     .    write(luo,310) cf_codes(i), (prn_list(j), 
     .                  num_deleted(j,i), j = 1, num_sat)
 310      format(' ',a4,50(:,' PN',i2.2,i5,1x))
      end do

****  Thats all
      return
      end
  
CTITLE TRIM_GAP
 
      subroutine trim_gap(i,j, ctol_cse, data_flag_cse, 
     .                   gap_mask)

      implicit none
 
*     Routine to trim the one-way data by flagging all data with too
*     short a gap or not enough data between gaps.
 
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

*   k           - epoch counter.
*   ltoc        - function to return the channel number
*               - for a particular observation.  If the
*               - satellite is not being oberved then -1
*               - is returned
*   first_good  - Epohc number of first good data point.
*   dnum_good   - Number of good data since last bias flag
*   ch          - Channel number for current satellite at
*               - current epoch.  (-1 if not observed)
*   gap_mask	- Mask which does not check if point
*                 hase been flagged with bit 6 (no enough data between
*                 bias flags.  (This bit is set here and during clean_dd)
 
      integer*4 k, ltoc, dnum_good, first_good, ch,
     .          gap_mask

*   data_ok     - Logical function which returns true if
*               - if none of the bits in the phase mask are
*               - set.
*   break       - Set true when missing or bad data
 
      logical data_ok, break

****  Now loop over epochs for this station/satellite
 
      num_deleted(j,i) = 0
      dnum_good = 0

*     Now loop over epochs
      do k = 1, num_ep
          ch = ltoc( ctol_cse(1,i,k), j, actual_max_chan)
*                          !  See if data OK.
          if( ch.gt.0 ) then
              if( data_OK(data_flag_cse(ch,i,k),0,
     .                    gap_mask) )  then
 
*                See if first good data after gaps
                 if( dnum_good.eq.0 ) first_good = k
                 dnum_good = dnum_good + 1
                 break = .false.
              else ! Break in data
                 break = .true.  ! bad data
              endif
           else
              break = .true.     ! No data 
           endif

****       OK; If we have break and some good data, see if
*          habe enough
           if( break .and. dnum_good.lt.min_good_bias .and.
     .                     dnum_good.gt.0  ) then
*              Not enough; so flag data
               call flag_data(ctol_cse,data_flag_cse, first_good,
     .              k-1, i, j , 6)

               num_deleted(j,i) = num_deleted(j,i) + 
     .                            dnum_good 
               dnum_good = 0
           elseif ( break ) then
               dnum_good = 0     ! Break but enough good data so OK.
           endif

      end do
 

****  Thats all
      return
      end
 
 
