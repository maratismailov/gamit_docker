CTITLE OPEN_SOLVK
 
      subroutine open_solvk(CO_name,ierr)

      implicit none
 
c     Routine to open SOLVK common
C MOD TAH 870128 Modified to use FmpOpen routine rather than OPEN call
 
 
      include '../includes/kalman_param.h'
      include '../includes/solvk_cntl.h'
 
 
      character*(*) CO_name
 
*   FmpShortName    - Returns name of a file connected to a DCB
*   ierr            - FMP error flag
 
      integer*4 FmpShortName, ierr
 
***** Open file
 
      call FmpOpen(isoldc, ierr, CO_name, 'RW', 1)
 
      if( ierr.ne. -6 ) then
          call report_error('FmpOpen', ierr,'open/creat',CO_name,1,
     .                      'OPEN_SOLVK')
      end if
 
***** Now use FmpShortName to get the name of the file.  In particular
*     this will return the cartridge on which this is located
 
      if( ierr.ge.0 ) then
          ierr = FmpShortName( isoldc, ierr, CO_name )
          call report_error('FmpShortName',ierr,'nam',CO_name,0,
     .                      'OPEN_SOLVK')
      end if
 
***** Set value_read false
 
*                                ! Force the first record of file to be
      control_read = .false.
*                             ! read, so that we will get the record
*                             ! lengths and start record numbers in
*                             ! RW_Solvk_Block
      end
 
