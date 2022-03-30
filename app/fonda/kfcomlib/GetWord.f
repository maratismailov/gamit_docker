CTITLE 'GETWORD'
 
      subroutine GetWord( line, word, indx )
 
 
*     Routine to get the first set of non--blank characters from
*     the 'line'.  The word is returned with no leading blanks.
 
*         start - Start position in line
*         fin   - stop position in line
*         i     - Loop counter
*         indx  - Pointer to position in line to get word from
*               - Indx is updated to give position of space after
*               - word.
*         len_line  - Total used length of line
*         trimlen   - Get length of string
 
      integer*4 start, fin, i, indx, len_line, trimlen
 
*         found     - Indicates either start or end found
 
      logical found
 
*             line  - Line to be read
*         word      - Word to be pulled from line
 
      character*(*) line, word
 
****  Get the length of line
      len_line = trimlen( line )
 
***** Start at front and find first non--blank character
      i = max(indx,1)
 
      found = .false.
      do while ( i.le.len_line .and. .not.found )
          if( line(i:i).eq.' ' ) then
              i = i + 1
          else
              found = .true.
          end if
      end do
 
*     Save start
      start = i
      found = .false.
      do while ( i.le.len_line .and. .not.found )
          if( line(i:i).ne.' ' ) then
              i = i + 1
          else
              found = .true.
          end if
      end do
 
*     Save stop
      fin = i
 
****  Save word
      word = line(start:fin)
      indx = fin
 
***** Thats all
      return
      end
 
