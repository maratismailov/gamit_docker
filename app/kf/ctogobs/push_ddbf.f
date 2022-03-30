CTITLE PUSH_DDBF
 
      subroutine push_ddbf( eps, ns, lv, data_flag_cse, bf_type_cse,
     .                      ctol_cse)

      implicit none
 
*     Routine to find next good point after epoch eps and put a
*     bias flag on it.  (In most cases there will already be a bias
*     flag there, but occasionally it has been removed by a previous
*     patch.)
 
* INCLUDES
 
      include '../includes/kalman_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
 
*   eps  - start epoch to start searching
*   ns, lv - site and satellite list number.
*   data_flag_cse(num_chan, num_cfiles, num_ep) - Data flag for each
*                   - measurement (same definition as in Gobs File)
*   bf_type_cse(num_chan, num_cfiles, num_ep) - Bias flag type, records
*                   - why a bias flag was set.  (Never reset even when
*                     the bias flag is removed)

*   ctol_cse(num_chan, num_cfiles, num_ep)  - Conversion from
*                   - channel number to satellite number
 
 
 
      integer*4 eps, ns, lv,
     .    data_flag_cse(num_chan, num_cfiles, num_ep),
     .    bf_type_cse(num_chan, num_cfiles, num_ep),
     .    ctol_cse(num_chan, num_cfiles, num_ep)
 
* LOCAL VARIABLES
 
*   i,j,k   - Counters of stations, satellites and epochs
*   ch      - Channel number at current epoch corresponding to
*           - list number lv.
*   ltoc    - function to convert sv list number to channel
*           - number.  Returns -1 if not observed.
*   save_bf_type - Bias flag type from the pushed bias flag
 
 
      integer*4 k, ch, ltoc, save_bf_type
 
*   data_OK     - Tests to see if data is OK using mask
*   kbit        - Checkes a bit is set
*   finished    - Indicates that we have pushed the bias flag
 
 
      logical data_OK, kbit, finished
 
***** Loop until we find the next good point
      finished = .false.
      k = eps
      ch = ltoc(ctol_cse(1,ns,k), lv, actual_max_chan)
      if( ch.gt.0 ) then
         save_bf_type = bf_type_cse(ch,ns,k)
      else
         save_bf_type = 0
      end if

      do while( .not.finished )
          k = k + 1
          if( k.gt.num_ep ) finished = .true.
          if( .not.finished ) then
              ch = ltoc(ctol_cse(1,ns,k), lv, actual_max_chan)
              if( ch.gt.0 ) then
*                 See if this data point is good
                  if( data_OK(data_flag_cse(ch,ns,k),0,phs_mask) ) then
 
*                     See if bias flag already set.
                      if( .not. kbit(data_flag_cse(ch,ns,k),31) .and.
     .                    .not. kbit(data_flag_cse(ch,ns,k),32) ) then
 
*                         Not set, so tell user and push bias flag
                          call sbit(data_flag_cse(ch,ns,k),31,1)
                          bf_type_cse(ch,ns,k) = save_bf_type
                          if( 1.eq.2 )
     .                    write(*,200) k, cf_codes(ns),
     .                        prn_list(lv)
 200                      format(' Bias pushed to epoch ',i5,
     .                           ' at site ', a,' PRN ',i2.2)
                      end if
                      finished = .true.
                  end if
              end if
          end if
      end do
 
***** Thats all
      return
      end
 
 
