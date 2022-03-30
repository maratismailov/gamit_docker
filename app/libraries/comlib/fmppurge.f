CTITLE FMPPURGE
 
 
      integer*4 function fmppurge (name)

      implicit none
 
*     Rotuine to emulate HP purge function.  The function returns
*     the error number.
 
* variables
*             name  - Name of file to purge
 
      character*(*) name
 
* Local Variable
*   ierr            - IOSTAT error flag.
*   num    -- Unit number that file is open to.
 
      integer*4 ierr, num
 
****  First we open the file and then we close with delete
      inquire(file=name,iostat=ierr,number=num)
      call report_error('IOSTAT',ierr,'inquir',name,0,'fmppurge') 
      if( num.gt.0 ) close(num)

      open(600, file=name, status = 'unknown', iostat=ierr)
      call report_error('IOSTAT',ierr,'del/open',name,0,'fmppurge') 
      if( ierr.eq.0 ) then
          close(600, status='delete', iostat=ierr)
          call report_error('IOSTAT',ierr,'del/close',name,0,'fmppurge') 
      end if
 
      fmppurge = -ierr
 
****  Thats all
      return
      end
 
 
