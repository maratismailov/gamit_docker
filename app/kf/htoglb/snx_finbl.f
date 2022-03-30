CTITLE SNX_FINBL
 
      subroutine snx_finbl(unit)
 
      implicit none

*     Routine to finish reading a sinex block
 
* PASSED VARIABLES
 
*   unit        - Unit number
 
      integer*4 unit
 
* LOCAL VARIABLES
 
*   ierr    - IOSTAT error
 
      integer*4 ierr
 
*   end_found   - Set true when end of block
*           - found
 
      logical end_found
 
*   line    - Line read from file
 
 
      character*80 line
 
****  Loop until end of block is found
      ierr = 0
      end_found = .false.
      do while ( ierr.eq.0 .and. .not.end_found)
          read(unit,'(a)', iostat=ierr) line
          if( line(1:1).eq.'-' ) end_found = .true.
      end do
 
****  Thats all
      return
      end
 
 
