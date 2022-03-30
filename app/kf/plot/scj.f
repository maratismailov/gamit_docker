CTITLE SCJ
 
      subroutine scj(or_x,or_y, j_x,j_y)
 
 
*     Routine to save the current orienataion and justifcation
*     of the labeling.  These values will be used in save_label.
 
      include 'plot_param.h'
      include 'plot_com.h'
 
*   j_x,j_y     - Justification to be saved
*   or_x, or_y  - Oriennataion to be saved
 
      real*4 j_x,j_y, or_x, or_y
 
*     Save values
 
      cori(1) = or_x
      cori(2) = or_y
 
      just(1) = j_x
      just(2) = j_y
 
      return
      end
 
 
