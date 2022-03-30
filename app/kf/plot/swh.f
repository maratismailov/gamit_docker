CTITLE SWH
 
      subroutine swh( w,h)
 
 
*     Routine to save the width and height of labels on an axis
*     These values are used later when the title of the axis has
*     to be put onto the plot
 
      include 'plot_param.h'
      include 'plot_com.h'
 
*   w,h - The width and height of the label in characters
 
      integer*4 w,h
 
***** Save the values
 
      wh(1) = w
      wh(2) = h
 
***** Thats all
      return
      end
 
