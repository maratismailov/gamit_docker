CTITLE GEN_HELP_NAME
 
      subroutine gen_help_name( local, full )

      implicit none
 
*     This routine will generate the full name of the
*     help file from its local name using the HELP_DIR
*     environment that should be set in the login file.
 
* PASSED VARIABLES
 
*   local       - short name for file
*   full        - Full name with directory prepended
 
      character*(*) local, full
 
* LOCAL VARIABLES
 
*   help_dir        - Name of help directory.
 
      character*128 help_dir
 
*   trimlen         - LEngth of string
 
      integer*4 trimlen
 
*     Get the help directory
 
      help_dir = ' '
      call getenv('HELP_DIR', help_dir)
      if( trimlen(help_dir).gt.0 ) then
          full = help_dir(1:max(1,trimlen(help_dir))) // local 
      else
          full = local
      end if
 
****  Thats all
      return
      end
 
 
 
