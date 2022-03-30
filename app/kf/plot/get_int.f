CTITLE    ................................................................
 
      subroutine get_int(buffer, int, num)
c
c     Routine to get an array of 'num' integers from the buffer
c
c Variables
c ---------
c buffer  -- the user buffer
c int  -- the array in which the intergers will be stored.
c num  -- the number of integers to get
c
      character*(*) buffer
c
      integer*4 int(1), num, ierr
      integer*4 indx
      character*4 cd
c
c
c.... Get the values from the string
c     read(buffer(9:),*, iostat = ierr, err=1000, end=1000 )
c    .  (int(i),i=1,num)
*     Replace with multiread call
      indx = 9 
      call multiread(buffer, indx, 'I4', ierr, int, cd, num)
c
      return
c
c.... Error return
c**   These can no longer be reached--commented by rwk 970920
c 1000 continue
c      call report_error('IOSTAT',ierr,'decod',buffer,0,'GET_INT')
c
c     return
      end
 
