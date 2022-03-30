CTITLE S2MOV
      subroutine s2mov(x,y)
 
 
*     Routine to save the current position of the cursor.  This
*     position will be used to save the label position information.
 
      include 'plot_param.h'
      include 'plot_com.h'
 
*   x,y     - The position to be saved
 
      real*4 x,y
 
*     Save the position
 
      pos(1) = x
      pos(2) = y
 
      return
      end
 
