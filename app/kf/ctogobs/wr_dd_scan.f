CTITLE WRITE_DD_SCAN
 
      subroutine write_dd_scan( unit, ep, ch1, ch2, s1, s2, dL1_slip,
     .            dL2_slip, chiqual, reliable )

      implicit none
 
*     Routine to write the double difference scan results in a
*     format readable by cview.
 
* PASSED VARIABLES
 
*   unit        - Unit number for output file
*   ep          - Epoch number being checked.
*   ch1, ch2    - Channel 1 and 2 of the double difference
*   s1, s2      - Site 1 and 2 of the double difference
*               - (nominally the bias flag is at ch1, s1)
 
      integer*4 unit, ep, ch1, ch2, s1, s2
 
*   reliable    - True if we can remove this bias flag.
 
      logical reliable
 
*   dL1_slip, dL2_slip  - Estimated number of cycle slips
*               - at L1 and L2
*   chiqual     - The bias flag removal chi**2 ratio for this
*               - double difference.
 
      real*8 dL1_slip, dL2_slip, chiqual
 
* LOCAL VARIABLES
 
*   cnt         - Current number of values written to file
*               - (Saved local an increment for each observation)
*   ierr        - IOSTAT erorr.  If we have a problem writing the
*               - file then we stop writing it.  Again ierr is
*               - saved locally.
 
      integer*4 cnt, ierr
 
      data cnt / 0 /, ierr / 0 /
 
****  Increment the count of the number values written
      if( ierr.ne.0 ) RETURN
 
      cnt = cnt + 1
      if( reliable ) then
          write(unit, 100, iostat=ierr) cnt, ep, ch1, ch2, s1, s2,
     .        dL1_slip, dL2_slip,'0000',chiqual
100       format(I4,i5,'       9.99',3x,4I5,2F10.1,7x,a,1x,F8.2)
      else
          write(unit, 100, iostat=ierr) cnt, ep, ch1, ch2, s1, s2,
     .        dL1_slip, dL2_slip,'1000',chiqual
      end if
      call report_error('IOSTAT',ierr,'writ','DD_REPORT File',0,
     .                  'write_dd_scan')
 
****  Thats all
      return
      end
 
 
 
