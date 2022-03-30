CTITLE FMPRENAME
 
      integer*4 function fmprename( old, new, option )

      implicit none
 
*     This routine will rename a file using the C subroutine RENAME
*     The option here indicates with an existing name should be
*     overwritten or not.  If option is 'K' (for kill) then an
*     existing file will be replaced.  Any other character will
*     protect existing files.
*
*     The function returns either the C error number or -2 if the file
*     already exists.
 
*     Restriction: File names should be less than 256 characters, and
*         have to be on the same file system.
*
 
* PASSED VARIABLES
 
*   old     - Old file name (must already exist)
*   new     - New file name (must not exist unless option is
*           - is set to K
*   option  - Only one option - K - kill any existing file
*           - with new file name
 
 
      character*(*) old, new, option
 
* LOCAL VARIABLES
 
*   null_terminate  - Function to add null character to the
*           - end of string.
*   rename      - Fortran  routine which calls rename.
*   ierr        - IOSTAT or C error return from RENAME
 
 
      integer*4 null_terminate, ierr
 
*   new_exists  - Return from the INQUIRE statement to see if
*           - file exists.
 
      logical new_exists
 
*   old_null        - Null terminated version of old name
*   new_null        - Null terminated version of new name.
 
 
      character*256 old_null, new_null
 
****  See if the new file exists if the K option has not been given
 
      if( option.ne.'k' .and. option.ne.'K' ) then
          inquire( file=new, exist=new_exists )
*                                 ! Get out, set error to -2
          if( new_exists ) then
              fmprename = -2
              RETURN
          end if
      end if
 
****  Now null terminating the file names.  Copy names first.
 
      old_null = old
      new_null = new
 
      ierr = null_terminate ( old_null )
*                                         ! Check for error
      if( ierr.ne.0 ) then
          fmprename = ierr
          RETURN
      end if
      ierr = null_terminate ( new_null )
*                                         ! Check for error
      if( ierr.ne.0 ) then
          fmprename = ierr
          RETURN
      end if
 
***** Now call the fortran interface routine rename 
c      ierr = rename( old_null, new_null )
      call rename( old_null, new_null,ierr )
 
      fmprename = ierr
 
***** Thats all
      RETURN
      end
