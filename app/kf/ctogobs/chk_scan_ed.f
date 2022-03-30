CTITLE CHK_SCAN_EDIT
 
      subroutine chk_scan_edit( cpass, scan_edits,
     .            L1_cyc_cse, L2_cyc_cse, ctol_cse, data_flag_cse,
     .            bf_type_cse )

      implicit none

 
*     This routine will check the number of bias flags set during scanning
*     and if these exceed the maximum number allowed then the site/svs
*     combination will be edited.

* INCLUDES
 
      include '../includes/kalman_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
 
*   data_flag_cse(num_chan, num_cfiles, num_ep) - Data flag for each
*                   - measurement (same definition as in Gobs File)
*   bf_type_cse(num_chan, num_cfiles, num_ep) - Bias flag type, records
*                   - why a bias flag was set.  (Never reset even when
*                     the bias flag is removed)

*   ctol_cse(num_chan, num_cfiles, num_ep)  - Conversion from
*                   - channel number to satellite number

*   cpass    - Iteration number
 
      integer*4 data_flag_cse(num_chan, num_cfiles, num_ep),
     .    bf_type_cse(num_chan, num_cfiles, num_ep),
     .    ctol_cse(num_chan, num_cfiles, num_ep), cpass
 
*   L1_cyc_cse(num_chan, num_cfiles, num_ep)    - Number of cycles
*                   - needed for each L1 phase measurement .  May be
*                   - fracttional for half cycle units.
*   L2_cyc_cse(num_chan, num_cfiles, num_ep)    - number of cycles
*                   - needed for each L2 phase measurement
 
 
      real*8 L1_cyc_cse(num_chan, num_cfiles, num_ep),
     .    L2_cyc_cse(num_chan, num_cfiles, num_ep)

*   scan_edits  - Set true if we decide we need to edit based on the
*                 number of bias flags added.

      logical scan_edits
 
* LOCAL VARIABLES
 
*   i,j,k,l - Counters of stations, satellites and epochs
 
      integer*4 i,j,k,l
 

****  Start by reporting the numbers of scanning bias flags
*     added.


      call rep_num_dd_flag( 6, cpass )

****  Now scan to see if we have too many bias flags
      scan_edits = .false.

*     Check to see if the max has been set.  If it is still
*     zero then just return

      if( max_scan_edit.le.0 ) RETURN

      do i = 1, num_cfiles
         do j = 1, num_sat
            if( num_dd_flags(i,j).gt. max_scan_edit) then

*               Report that we are editing this site/sv combination
                scan_edits = .true.
                write(*,120) cf_codes(i), prn_list(j), 
     .               num_dd_flags(i,j)
 120            format(' EDITING Site ',a4,' PRN_',i2.2,' due to ',i4,
     .                 ' DDSCAN bias flags')

*               Now loop over all epochs and edit the data
                do k = 1, num_ep
                   do l = 1, num_chan
                      if( ctol_cse(l,i,k).eq. j ) then
                          call sbit(data_flag_cse(l,i,k),19,1)
                      end if
                   end do
                end do
            end if
         end do
      end do

****  Now if we edited the data, re-set all the cycle adjusments
*     back to zero, clear the num_dd_flags counter and try again
      if( scan_edits ) then
          call init_scan_edits( num_cfiles, num_sat, num_dd_flags,
     .                          max_cfiles )

          do k = 1, num_ep        
             do i = 1, num_cfiles
                do j = 1, num_chan
* MOD TAH 000309: Commented out reseting the cycles 
C                  L1_cyc_cse(j,i,k) = 0.d0
C                  L2_cyc_cse(j,i,k) = 0.d0
*                  Clear the bit for number of DDScan biases 
*                  detected.                  
                   call sbit(bf_type_cse(j,i,k), 5,0)
                end do
             end do
         end do
      end if

***** Thats all
      return
      end

CTITLE REP_NUM_DD_FLAG

      subroutine rep_num_dd_flag( un, cpass )

      implicit none

*      Routine to report the numbers of scanning bias flags added.

* INCLUDES
 
      include '../includes/kalman_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
*  un   - unit number for output
*  cpass - Pass number for iterations in removing data

      integer*4 un, cpass

* LOCAL VARIABLES

*  i,j  - Loop counter

      integer*4 i,j

***** Write out header
      write(un, 120 ) cpass, (prn_list(j), j=1,num_sat)
 120  format(/,' DDScan bias flags added report for pass ',i2,/,
     .         ' SITE PN',i2.2,100(2x,i2.2))
      do i = 1, num_cfiles
         write(un,140) cf_codes(i), (num_dd_flags(i,j), j=1,num_sat)
 140     format(1x,a4,1x,100i4)   
      end do

****  Thats all
      return
      end

      
