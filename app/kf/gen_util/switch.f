  
  
      subroutine switch_ch(a,b,temp)

      implicit none 
  
c
c     Switches character values
c
      character*(*) a, b, temp
c
      temp = a
      a    = b
      b    = temp
c
      end
  
c.......................................................................
  
      subroutine switch_8(a,b,temp)
 
      implicit none 
 
c
c     Switches real*8 values
c
      real*8 a, b, temp
c
      temp = a
      a    = b
      b    = temp
c
      end
  
c.......................................................................
  
      subroutine switch_4(a,b,temp)

      implicit none 
  
c     Switches real*4 values
c
      real*4 a, b, temp
c
      temp = a
      a    = b
      b    = temp
c
      end
  
c.......................................................................
  
      subroutine switch_I4(a,b,temp)
 
      implicit none 
 
c
c     Switches integer*4 values
c
      integer*4 a, b, temp
c
      temp = a
      a    = b
      b    = temp
c
      end
  
