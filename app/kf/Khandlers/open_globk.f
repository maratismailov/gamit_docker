CTITLE OPEN_GLOBK
 
      subroutine open_globk(CO_name,ierr)

      implicit none
 
c     Routine to open GLOBK common
C MOD TAH 870128 Modified to use FmpOpen routine rather than OPEN call
 
 
      include '../includes/kalman_param.h'
      include '../includes/globk_cntl.h'
 
 
      character*(*) CO_name
 
*   FmpShortName    - Returns name of a file connected to a DCB
*   ierr            - FMP error flag
*   trimlen         - HP function for length of string
 
      integer*4 ierr, trimlen
 
****  Check to see if we have a file name.  If name is blank then use
*     the default
 
*                                         ! Use default
      if( trimlen(CO_name).eq.0 ) then
          CO_name = glb_com_default
      end if
 
***** Open file
 
      call FmpOpen(glb_com_dcb, ierr, CO_name, 'RW', 1)
 
      call report_error('FmpOpen', ierr,'open/creat',CO_name,1,
     .                  'OPEN_GLOBK')
 
***** Set value_read false
 
*                             ! Force the first record of file to be
      glb_con_read = .false.
*                             ! read, so that we will get the record
*                             ! lengths and start record numbers in
*                             ! RW_Globk_Block
      end
 
