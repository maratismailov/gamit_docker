CTITLE GET_CAPROG
 
      subroutine get_caprog( runstring, caprog_name)

      implicit none
 
*     Routine to find out if ctogobs or autcln is being
*     run.
 
* PASSED VARIABLES
 
*   runstring       - Name of program from runstring
*   caprog_name     - Generated name (either AUTCLN or CTOGOBS)
 
      character*(*) runstring, caprog_name
 
* LOCAL VARIABLES
 
*   slash_pos       - Position of last slash in name
*   i               - Loopcounter
*   trimlen         - Length of string
*   len_run         - Length of runsttring
 
 
      integer*4 slash_pos, i, trimlen, len_run
 
***** Get the length
      len_run = trimlen(runstring)
      slash_pos = index(runstring,'/')
      if( slash_pos.gt.0 ) then
          i = len_run
          do while (runstring(i:i).ne.'/')
              i = i - 1
          end do
          slash_pos = i - 1
      end if
      if( index(runstring(slash_pos+1:),'ctog') .gt.0 ) then
          caprog_name = 'CTOGOBS'
      else
          caprog_name = 'AUTCLN '
      end if
 
****  Thats all
      return
      end
 
