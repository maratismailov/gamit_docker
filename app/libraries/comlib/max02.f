      integer*2 function min02(ival1,ival2)

c  function to get around a problem with the g77 compiler which won't allow
c  you to take the minimum value of a I*2. This routine will convert the 
c  variables ival1,ival2 (I*2) into  I*4, find the minimum value and then
c  return it as an I*2
c
c  P. Tregoning
c 24th September 1996

      implicit none

      integer*2 ival1,ival2
      integer*4 temp,temp1,temp2

      temp1 = ival1
      temp2 = ival2
      temp = min0(temp1,temp2)
      min02 = temp

      return
      end


