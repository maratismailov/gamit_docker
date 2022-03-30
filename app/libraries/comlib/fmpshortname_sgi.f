CTITLE FMPSHORTNAME

      integer*4 function fmpshortname( dcb, ierr, name )
      implicit none

c     This routine will enquire as to the name of the file connected
c     to the dcb buffer, and return the name of the file with the
c     current working directory prepended.

c         dcb(16)- DCB buffer (see FmpOpen for meaning of
c                - of the words in this routine)
c   ierr         - IOSTAT (with negative sign of the return)

      integer*4 dcb(16), ierr

c   name         - name of the file connected to the dcb
c                - buffer

      character*(*) name

c LOCAL VARIABLES

c   trimlen      - HP1000 emulation routine for length of a string
c   lendir       - length of the directory return

      integer*4 trimlen, lendir

c   localname    - local name of the file (as used in the open)
c   dirname      - Name of the current working directory

      character*132 localname, dirname, getcwd

c   getcwd       - system function.  stores current path name in its
c                  argument and returns a status code (0 => ok)

c MOD JFC 970320 - porting to Solaris 2.5: adapt getcwd definition
C MOD SSS 970718 - Ported to IRIX

c***  First get the working directory
      lendir = len(dirname)
      dirname = getcwd(lendir)
      lendir = trimlen(dirname)
      inquire(dcb(1), name = localname, iostat=ierr)

c***  Now construct full name if ierr OK
      if(ierr.eq.0) then
         name = dirname(:lendir) // '/' // localname
      else
          name = ' '
          ierr = -ierr
      end if

      FmpShortName = ierr

      end
