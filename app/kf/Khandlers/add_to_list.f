CTITLE ADD_NAME_TO_LIST
 
      subroutine add_name_to_list( next_name, lnum, gnames, gnum )
 

      implicit none
 
*     Routine to add the name just read from the to the global
*     list of names.  If next_name is blank then nothing is done.
 
*   gnum        - Total number of names in the list
*   i,j         - Loop counters
*   lnum        - Number of this name in list (will be added
*               - if not already there.)
*   trimlen     - HP function for length of string
 
      integer*4 gnum, i, lnum, trimlen
 
*   gnames(1)   - Global list of names.
*   next_name   - Name to be positioned in the list
 
      character*(*) gnames(1), next_name
 
*   name_found  - Indicates that the name is found.
 
      logical name_found
 
***** Check to see if next_name is blank
 
      lnum = 0
      name_found = .false.
 
*                                             ! No name
      if( trimlen(next_name).eq.0 ) RETURN
 
***** See if we can find name
 
      do i = 1, gnum
          if( next_name.eq.gnames(i) ) then
              lnum = i
              name_found = .true.
          end if
      end do
 
***** If we did not find name add to list
      if( .not.name_found ) then
          gnum = gnum + 1
          gnames(gnum) = next_name
          lnum = gnum
      end if
 
***** Thats all
      return
      end
 
