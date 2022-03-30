CTitle fullfilename
 
      subroutine fullfilename(name, type, size, recl, fullname)

      implicit none

* MOD TAH 190531: Increase format for size to avoid possible overflow.
 
*     This rotuine will construct the full file name of a file given the
*     type and record length.
 
* VARIABLES
 
*   type        - File type
*   recl        - Record length in i*4 words
*   size        - size of file (NOT used in UNIX version).  Added to name just
*                 sowe can see size
 
      integer*4 type, recl, size
 
*   name        - Input name
*   fullname    - the full name with appended information
 
      character*(*) name, fullname
 
* LOCAL VARIABLES
 
*   trimlen     - Returns length of string
*   lenname     - Length of the name
 
      integer*4 trimlen, lenname
 
****  Get the current name length
      lenname = trimlen(name)
*                                 ! Start appending information
      if( lenname.gt.0 ) then
          fullname = name(:lenname) // ':::'
          lenname = lenname + 3
 
*         Add type to file
          write(fullname(lenname+1:),'(I1,'':'')') type
          lenname = lenname + 2
 
*         Now write our the record length
          write(fullname(lenname+1:),'(i4.4,'':'')') recl
          lenname = lenname + 5

*         Add size even though it wont be used
* MOD TAH 190531: Icreased format to avoid possible overflow.
          write(fullname(lenname+1:),'(I9.9,'':'')') size
 
*                                 ! Set name to blank
      else
          fullname = ' '
      end if
 
****  Thats all
      return
      end
 
 
