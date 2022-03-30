CTITLE INIT_SCAN_EDITS 

      subroutine init_scan_edits( num_cfiles, num_sat, num_dd_flags,
     .                            max_cfiles )

      implicit none

*     Routine to clear the counters of the number of dd bias flags
*     added.

* PASSED Variables
* num_cfiles   - number of cfiles
* num_sat      - Number of satelittes expected
* max_cfiles   - Maximum number of cfiles allowed.
* num_dd_flags(max_cfiles, num_sat)  - number of dd_scan bias 
*                flags added

      integer*4 num_cfiles, num_sat, max_cfiles, 
     .          num_dd_flags(max_cfiles, num_sat)

* LOCAL VARIABLES

      integer*4 i,j


***** Just loop over complete range
      do i = 1, num_cfiles
         do j = 1, num_sat
            num_dd_flags(i,j) = 0
         end do
      end do

***** Thats all
      return
      end



