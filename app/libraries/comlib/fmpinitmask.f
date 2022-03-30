      subroutine FmpInitMask(dcb, ierr, mask, curpath)

      implicit none
 
*     Routine to emulate the routine of the same name from the HP1000
*     system.  Here MASK containes a wild card file mask, and this
*     routine will system an ll command with the mask and generate an
*     output file.  FmpNextMask will then read this generated file.
 
*   dcb(*)  - Old DCB buffer, Here just contains the UNIT number
*   ierr    - negative of IOSTAT error
 
      integer*4 dcb(*), ierr, null_terminate, nlen
 
*   mask    - The mask to be searched on
*   curpath - Current path name (although we do not use this yet)
 
 
      character*(*) mask, curpath
 
* LOCAL VARIABLES
 
*   trimlen     - Length of string
*   lenmask     - Length of mask
*   system      - Runs csh command
 
      integer*4 trimlen, lenmask, system  
* PT & HM 960920: make "system" an external variable for g77 compilation 
      external system

*   ll_command  - Command containing the ll command to be shelled.
 
 
      character*128 ll_command
                         

****  Construct the shell command we will need
      lenmask = trimlen(mask)
      ll_command = ' '
      if( lenmask.gt.0 ) then
          write(ll_command,100, iostat=ierr) mask(1:lenmask)
	  nlen = null_terminate(ll_command)
 100      format('ls -lg ',a,' > KDL@DIR@FILE')
      else
          ierr = -1
          return
      end if
 
****  Now send the shell command
      ierr = system( ll_command )
****  Now open the output file created
      if( ierr.eq.0 ) then
          dcb(1) = 690
	  close(dcb(1))
          open(dcb(1),file='KDL@DIR@FILE', iostat=ierr, status='old')
      end if
      ierr = -ierr
      return
      end
 
 
