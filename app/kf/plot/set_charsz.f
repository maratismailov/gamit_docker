CTITLE    ..................................................................
 
      subroutine set_charsz(width_mm, height_mm)
c
c     Routine to set the character size in Graphics/1000
c
c Include files
c -------------
*                          ! the parameter file
      include 'plot_param.h'
c
*                          ! the plot_com common block
      include 'plot_com.h'
c
c Variables
c ---------
c width_mm -- the width of the character in mm
c height_mm -- the height of the characters in mm
c
      real*4 width_mm, height_mm
 
c
c Local variables
c ---------------
c mmtow -- conversion from mm to width world coordinates
c mmtoh -- conversion from mm to height world coordinates
c width -- the width of the character  (for jcsiz)
c height -- the height  of the character  (for jcsiz)
c gap -- the gap between characters
c
      real*4 mmtow, mmtoh
 
c
      real*4 width, height, gap
 
c
c.... Compute the width of the characters in world coordinates
c     First get the mm/world coordinates
      mmtow = x_size_mm/( (scale_size(2)-scale_size(1))/
     .         (view_size(2)-view_size(1)) )
      mmtow = 1.0
      width = width_mm/mmtow
c
c.... Now height
      mmtoh = y_size_mm/( (scale_size(4)-scale_size(3))/
     .         (view_size(4)-view_size(3))/aspect_ratio )
      mmtow = 1.0
      height = height_mm/mmtoh
c
      gap = 0.
c
c.... Now set the character size
      call jcsiz(width, height, gap)
c
c.... Thats all
      return
      end
 
