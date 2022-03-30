 
CTITLE FMPSHORTNAME
 
 
      integer*4 function fmpshortname( dcb, ierr, name )

      implicit none
 
*     This routine will enquire as to the name of the file connected
*     to the dcb buffer, and return the name of the file with the
*     current working directory prepended.
 
*         dcb(16)       - DCB buffer (see FmpOpen for meaning of
*                       - of the words in this routine)
*   ierr                - IOSTAT (with negative sign of the return)
 
      integer*4 dcb(16), ierr
 
*   name                - name of the file connected to the dcb
*                       - buffer
 
      character*(*) name
 
* LOCAL VARIABLES
 
*   trimlen             - HP1000 emulation routine for length of a
*                       - string
*   lendir              - length of the directory return
 
      integer*4 trimlen, lendir,getcwd
 
*   localname           - local name of the file (as used in the open)
*   dirname             - Name of the current working directory
 
      character*132 localname, dirname
 
****  First get the working directory
c      call getcwd(dirname)
      ierr = getcwd(dirname)
      lendir = trimlen( dirname)
 
      inquire( dcb(1), name = localname, iostat=ierr )
 
****  Now construct full name if ierr OK
      if( ierr.eq.0 ) then
          name = dirname(:lendir) // '/' // localname
      else
          name = ' '
          ierr = -ierr
      end if
      
      fmpshortname = -ierr
      
****  Thats all
      return
      end
