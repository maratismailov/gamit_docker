
      integer*4 function  null_terminate(string)

      implicit none
c
c     This function will put a null at the end of fortran string
*     variable string.  
*     Return is -1 if string is not long enough for the terminator.

      character*(*) string
 
      character*1 blank/' '/
 
      integer*4 leng, i
 
c
      null_terminate = 0
      leng = len(string)
      do i = leng,1,-1
        if(string(i:i).ne.blank) then
          if(leng.gt.i) then
            string(i+1:i+1) = char(0)
          else
            null_terminate = -1
          endif
c                             break out of loop because we are done.
          go to 100
          endif
      enddo
c
      string(1:1) = char(0)
 100  return
      end
 
