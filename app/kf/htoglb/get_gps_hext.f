CTITLE GET_GPS_HEXT
 
      subroutine get_gps_hext( line, hf_ext )
 
      implicit none

*     Routine to read the keys line of the hfile entry and
*     get string starting column 75.
*     ** WARNING *** If format changed then this routine will need
*     to be updated.
*
* Example of line:
* keys: DEFLT FULL  DBLE  LC    NOION NOATM FREE  STN   ORB   ZEN   NOCLK  GLR
 
* PASSED VARIABLES:
 
*   line        - Line read from hfile (should have keys: at start)
*   hf_ext      - Extent read from line.  (Returns blank if string not
*               - found)
 
      character*(*) line, hf_ext

* LOCAL VARIABLES:
*   ierr        - IOSTAT error

      integer*4 ierr
 
*     Clear the hf_ext entry
      hf_ext = ' '
 
*     If beginning of line is OK read the extent
      if( line(2:6).eq.'keys:' ) then
          read(line,100, iostat=ierr) hf_ext
 100      format(74x,a4)
          call caseunfold(hf_ext)
          if( ierr.ne.0 ) hf_ext = ' '
      end if
 
****  Thats all
      return
      end
 
 
 
