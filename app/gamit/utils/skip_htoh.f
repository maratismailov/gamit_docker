CTITLE SKIP_HTOH

      subroutine skip_htoh ( unit, line, ierr )

*     Routine to skip over the htoh comment lines inserted
*     when station names are changed.  The comment lines
*     are echoed to the user.
*
* MOD TAH 941223: Read the next line after 'biases' found so
*     that we are at the keys: line.  (After htoh run this
*     ensures that hfile type is returned correctly
* MOD TAH 941223: Made check on keys more specific (to be in
*     column 2) incase keys is used in the htoh comments.

* PASSED VARIABLES

*   unit            - Unit number of file being read
*   ierr            - IOSTAT error on read

      integer*4 unit, ierr

*   line            - Line read from file (checked for word
*                   - biases

      character*(*) line

* LOCAL VARAIBLES

*   indx            - Position of 'biases' line.
*   trimlen         - Length of string

      integer*4 indx, trimlen

*   found           - Set true when 'biases' is found

      logical found

****  echo line already read.
      write(*,'(a)') line(1:max(1,trimlen(line)))

***** Loop in file until 'biases' found or EOF reached.

      found = .false.
      ierr  = 0
      do while (ierr.eq.0 .and. .not.found )
          read(unit,'(a)', iostat=ierr) line
          indx = index(line,'biases')

*         If biases line found, read the next line:
          if( indx.ne.0 ) then
              read(unit,'(a)', iostat=ierr) line
              indx = 0
          end if

*         If the biases line not found try for the keys
*         line for GLOBALT files
          indx = index(line,'keys')

*                                     ! Write comment to user
          if( indx.ne.2 ) then
              write(*,'(a)') line(1:max(1,trimlen(line)))
*                                     ! Set found
          else
              found = .true.
          end if
      end do

***** Thats all
      return
      end

