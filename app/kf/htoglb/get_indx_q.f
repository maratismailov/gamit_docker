CTITLE GET_INDX_Q
 
      subroutine get_indx_q( hfile, indx_q )
 
      implicit none

*     Routine to determine the index of the first . in the
*     hfile name after the directory structure. (Thus allowing
*     dots in the directory names)
 
*   indx_q      - Character position of the first . in the
*               - name after the directory.
 
      integer*4 indx_q
 
*   hfile       - Name of h-file
 
      character*(*) hfile
 
* LOCAL VARIABLES
 
*   last_slash  - Position of last / in name.
*   trimlen     - Length of non-blank portion of string.
 
      integer*4 last_slash, trimlen
 
****  Find out where the directoty names end
 
      last_slash = trimlen(hfile)
      do while ( hfile(last_slash:last_slash).ne.'/' .and.
     .            last_slash.gt.2 )
          last_slash = last_slash - 1
      end do
 
*     Now find the position of the . relative to this point
      indx_q = index( hfile(last_slash:),'.')
 
*     See if we found .
      if( indx_q.gt.0 ) then
          indx_q = indx_q + last_slash - 1
      else
*         MOD TAH 980803: See of '_X' in name
          indx_q = index( hfile(last_slash:),'_X')
          if( indx_q.gt.0 ) then
             indx_q = indx_q + last_slash - 1
          else
             indx_q = trimlen(hfile)
          end if
      end if
 
***** If last_slash is less than 2 (i.e., the character we will use
*     is one before this.  Print warning asn return 2
      if( indx_q.lt.2 ) then
          write(*,100) hfile(1:max(1,trimlen(hfile)))
 100  format(' ***WARNING *** Problem with isolating extent in ',
     .            a, '.  First character will be used')
          indx_q = 2
      end if
 
***** Thats all
      return
      end
 
 
 
 
