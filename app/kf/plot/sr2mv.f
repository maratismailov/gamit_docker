CTITLE SR2MV
      subroutine sr2mv(dx,dy)
 
 
*     Routine to save the current relative shift of the cursor.  This
*     position will be used to save the label position information.
 
      include 'plot_param.h'
      include 'plot_com.h'
 
*   dx,dy     - The position shift to be saved
 
      real*4 dx,dy
 
*     Save the position
 
      pos(1) = pos(1) + dx
      pos(2) = pos(2) + dy
 
      return
      end
 
