 
      integer*4 function fmpendmask( dcb )

      implicit none
 
*     Routine to close and delete the file containing the
*     directory information
 
*   dcb(*)  - Dcb buffer only contains the unit number
 
      integer*4 dcb(*)
 
*   ierr    - IOSTAT error in close
 
      integer*4 ierr
 
*     CLose the file and delete
 
      close(dcb(1), status='delete', iostat=ierr)
      ierr = -ierr
      
      fmpendmask = ierr
      
      return
      end
