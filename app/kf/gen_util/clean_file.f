ctitle
 
      subroutine clean_file(file)

      implicit none 
c
c     This routine will purge file 'file'
c
c Variables
c ---------
c file -- the file to be purged
c
      character*(*) file
c
*   FmpPurge    - FMP function to purge file (returns error)
*   ierr        - FMP error message
 
      integer*4 FmpPurge, ierr
c
      ierr = FmpPurge(file)
 
      return
      end
 
c......................................................................
