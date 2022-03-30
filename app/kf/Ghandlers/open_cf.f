CTITLE OPEN_CF
 
      subroutine open_cf( unit, name, ierr )

      implicit none 
 
*     This routine will open the c-file as fortran unit number unit.
*     Any error is reported, but it is the users responsibility to
*     process the error.
 
* INCLUDE FILES
*     None
* PASSED VARIABLES
 
*   unit        - Unit number to be attached to c-file
*   ierr        - IOSTAT error on open
 
      integer*4 unit, ierr
 
*   name        - Name of cfile
 
 
      character*(*) name
 
* LOCAL VARIABLES
*     None
 
***** Open the c-file (as old) and report any error.
 
      open(unit, file=name, iostat=ierr, status='old',
     .    form='unformatted', access='sequential' )
 
      call report_error('IOSTAT', ierr,'open', name, 0,
     .                'OPEN_CF')
 
***** Thats all
      return
      end
 
