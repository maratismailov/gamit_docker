CTITLE OPENFILE
 
      subroutine openfile( unit, name, nstatus, naccess, nformat,ierr )

      implicit none
 
*     This routine will open a file as fortran unit number unit.
*     The data form could be ASCII or BINARY
*     Gang Chen, July 3, 1996
*     Any error is reported, but not processed.
 
* INCLUDE FILES
*     None
* PASSED VARIABLES
 
*   unit        - Unit number to be attached to the file
*   ierr        - IOSTAT error on open
 
      integer*4 unit, ierr
 
*   name        - Name of file

      character*(*) name
      
*   nstatus    - status of file: new, unknown, old
*   nformat    - format option of file: ASCII/BINARY (formatted/unformatted)
*   naccess    - access option of file: sequential
 
      character*(*) nstatus, nformat, naccess
* LOCAL VARIABLES
*   nform_real -- Real format type to use.

      character*16 nform_real 
 
***** Open the new file
      if(nformat.eq.'BINARY'.or.nformat.eq.'unformatted') then
           nform_real ='unformatted'
      else
           nform_real = 'formatted'
      endif
      
      if (naccess.ne.'sequential') then
        open(unit, file=name, iostat=ierr, status=nstatus)
c     .    form=nform_real, access=naccess )
      else
        open(unit, file=name, iostat=ierr, status=nstatus,
     .    form=nform_real)
      endif 

***** Report any error.

      If (ierr.ne.0) then
        call report_error('IOSTAT', ierr,nstatus, name, 1,
     .                'openfile')
      endif
        
 
***** Thats all
      return
      end
 
