CTITLE REPORT_BF
 
      subroutine report_bf( L1_cyc_cse, L2_cyc_cse, ctol_cse, 
     .        data_flag_cse, bf_type_cse )

      implicit none
 
*     Routine add up each of the types of bias flags that have been
*     set and un-set and report on the results.  These summaries can
*     be used to determine how many of the bias flags set were really
*     needed.
 
* INCLUDES
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include 'ctogobs_com.h'
 
* LOCAL PARAMETERS
 
*   max_bf_types    - Maximum number of bias flag types
 
      integer*4 max_bf_types
 
      parameter ( max_bf_types = 7 )
 
* PASSED VARIABLES
 
*   data_flag_cse(num_chan, num_cfiles, num_ep) - Data flag for each
*                   - measurement (same definition as in Gobs File)
*   bf_type_cse(num_chan, num_cfiles, num_ep) - Bias flag type, records
*                   - why a bias flag was set.  (Never reset even when
*                     the bias flag is removed)
*   ctol_cse(num_chan, num_cfiles, num_ep)  - Satellite list number in
*                     in each channel
 
 
      integer*4 data_flag_cse(num_chan, num_cfiles, num_ep),
     .    bf_type_cse(num_chan, num_cfiles, num_ep),
     .    ctol_cse(num_chan, num_cfiles, num_ep)
 
*   L1_cyc_cse(num_chan, num_cfiles, num_ep) - Number of cycles to be
*         added to L1
*   L2_cyc_cse(num_chan, num_cfiles, num_ep) - Number of cycles to be
*         added to L2
      real*8 L1_cyc_cse(num_chan, num_cfiles, num_ep),
     .    L2_cyc_cse(num_chan, num_cfiles, num_ep)
 
* LOCAL VARIABLES
 
*   i,j,k,l     - Loop counters
*   lv          - List number of the satellite
*   num_bf_set(max_bf_types)    - Number of each type of bias flag
*               - set as data processed
*   num_bf_rem(max_bf_types)    - Number of each type of bias flag
*               - remaining after data cleaned
*   num_bf_edt(max_bf_types)    - Number of each type of bias flag
*               - remaining in edited data
*   num_bf_slp(max_bf_types)    - Number of each type of bias flag
*               - which were removed and for which there was change
*               - in the number of cycles.
 
      integer*4 i,j,k,l, lv, num_bf_set(max_bf_types),
     .    num_bf_rem(max_bf_types), num_bf_edt(max_bf_types),
     .    num_bf_slp(max_bf_types)
 
*   data_OK     - Function which returns true is no bit in the mask
*               - are set
*   good_bf     - Function that returns true if bias flag set on a good
*               - data point
*   kbit        - Tests if bit s set.
 
      logical data_OK, good_bf, kbit
 
*   last_cyc(2,max_gprn)    - last number of cycle slips on each each
*               - list entries for the satellites.
      real*8 last_cyc(2,max_gprn)
 
*   bf_type_codes(max_bf_types) - codes for the types of bias flags
*                   - (ie, reason bias flag was added).
 
      character*4 bf_type_codes(max_bf_types)
 
      data bf_type_codes / ' ORG',' JMP',' ION',' GAP',' DDS',
     .                     ' WLS',' DDC' /
 
***** Write out the heading
      write(uns, 100) ((bf_type_codes(i),i=1,max_bf_types), j=1,5)
 100  format(/' BIAS FLAG REPORT: Types',/,
     .        A4,'-Original',t15,a4,'-Big Jump',t30,a4,'-Ion Jump',t45,
     .        a4,'-Data Gap',t60,a4,'-DD scan',t75,
     .        a4,'-Wide Lane',t90,a4,'-DD cleaning',/,
     .        'SITE',t12,'# Flagged',t33,'|',t40,'# Remaining',
     .        t62,'|',t68,'# Edited',t91,'|',
     .        t96,'# with jump',/,
     .        4x,7a4,'|',7a4,'|',7a4,'|',7a4)
 
*     Loop over the sites
      do i = 1, num_cfiles
          do j = 1, max_bf_types
              num_bf_set(j) = 0
              num_bf_rem(j) = 0
              num_bf_edt(j) = 0
              num_bf_slp(j) = 0
          end do
          do j = 1, 32
              edit_counts(j,i) = 0
          end do

          do j = 1, max_gprn
              last_cyc(1,j) = 0.0d0
              last_cyc(2,j) = 0.0d0
          end do
 
*         Now loop over epochs of data
          do k = 1, num_ep
 
*             Loop over channels
              do j = 1, num_chan
*                 See if we have data in this channel
                  if( .not.kbit(data_flag_cse(j,i,k),30) ) then
                      lv = ctol_cse(j,i,k)
 
*                     See if any of the bias flags were ever set
                      do l = 1,max_bf_types
                          if( kbit(bf_type_cse(j,i,k),l) ) then
                              num_bf_set(l) = num_bf_set(l) + 1
                          end if
*                         See if bias flag still set
                          if( good_bf(data_flag_cse(j,i,k),0,
     .                                phs_mask) .and.
     .                         kbit(bf_type_cse(j,i,k),l) ) then
                              num_bf_rem(l) = num_bf_rem(l) + 1
                          end if
 
*                         See if edited out
                          if( .not.data_OK(data_flag_cse(j,i,k),0,
     .                                    phs_mask) .and.
     .                         kbit(bf_type_cse(j,i,k),l) )then
                              num_bf_edt(l) = num_bf_edt(l) + 1
                          end if
 
*                         If there is no bias flag, see if the number of
*                         cycles changed.
                          if( data_OK(data_flag_cse(j,i,k),0,
     .                        phs_bias_mask) .and.
     .                         kbit(bf_type_cse(j,i,k),l) .and.
     .                        (L1_cyc_cse(j,i,k)-last_cyc(1,lv).ne.0 
     .                          .or.
     .                        L2_cyc_cse(j,i,k)-last_cyc(2,lv).ne.0)
     .                                                           ) then
                              num_bf_slp(l) = num_bf_slp(l) + 1
                          end if
*                             ! Looping over types
                      end do

*                     Now sum up the number of bits sets in data_flag
                      do l = 1, 27
                         if( kbit(data_flag_cse(j,i,k),l) ) then
                             edit_counts(l,i) = edit_counts(l,i) + 1
                         end if
                      end do
*                     See if data is good (ignoring the bias flag)
                      if( data_OK(data_flag_cse(j,i,k),0, 
     .                    phs_mask) ) then
                          edit_counts(30,i) = edit_counts(30,i) + 1
                      end if
                      
 
*                     Save the number of cycles for this satellites
                      if( data_OK(data_flag_cse(j,i,k),0,
     .                    phs_mask) ) then
                          last_cyc(1,lv) = L1_cyc_cse(j,i,k)
                          last_cyc(2,lv) = L2_cyc_cse(j,i,k)
                      end if
*                             ! Data in channell
                  end if
*                             ! Looping over channels
              end do
*                             ! Looping over epochs
          end do
 
*         Writethe results out
          write(uns,300) cf_codes(i),(num_bf_set(j),j=1,max_bf_types),
     .        (num_bf_rem(j),j=1,max_bf_types),
     .        (num_bf_edt(j),j=1,max_bf_types),
     .        (num_bf_slp(j),j=1,max_bf_types)
 300  format(A4,7I4,'|',7i4,'|',7i4,'|',7i4)
      end do

*     Now write out the counts for the editing of the data
      write(uns,400) (edit_names(i),i=1,7), (edit_names(i),i=14,15),
     .               (edit_names(i),i=19,27)
 400  format(/,' EDITING REPORT AND SITE PARAMS',/,
     .         ' SITE  MnCLN  MnOUT  SNR ',1x,18(1x,a4),2x,
     .         ' Good',/,
     .         '       (deg)  (deg)  L1 L2')
      do i = 1, num_cfiles
         if( .not.use_gamit_elc )
     .   write(uns,420) cf_codes(i), site_celev(i)*180.0/pi,
     .                 site_oelev(i)*180.0/pi, site_snr(1,i),
     .                 site_snr(2,i), (edit_counts(j,i), j=1,7),
     .                 (edit_counts(j,i),j=14,15), 
     .                 (edit_counts(j,i),j=19,27), 
     .                 edit_counts(30,i)
 420     format(1x,a4,1x,f6.2,1x,f6.2,1x,i2,1x,i2,1x,18i5,i7)
         if( use_gamit_elc )
     .   write(uns,420) cf_codes(i), min_cfile_elev*180.0/pi,
     .                 min_cfile_elev*180.0/pi, site_snr(1,i),
     .                 site_snr(2,i), (edit_counts(j,i), j=1,7),
     .                 (edit_counts(j,i),j=14,15), 
     .                 (edit_counts(j,i),j=19,27), 
     .                 edit_counts(30,i)
      end do
 
***** Thats all
      return
      end
 
