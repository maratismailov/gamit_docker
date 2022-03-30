      integer*2 function iabs2(ivalin)

c  function to get around a problem with the g77 compiler which won't allow
c  you to take the absolute value of a I*2. This routine will convert the 
c  variable ival (a I*2) into a a I*4, take the absolute value and then
c  return it as an I*2
c
c  P. Tregoning
c 24th September 1996

      implicit none

      integer*2 ivalin
      integer*4 temp

      temp = ivalin
      temp = iabs(temp)
      iabs2 = temp

      return
      end


