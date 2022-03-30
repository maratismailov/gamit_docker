CTITLE COMLF
 
      subroutine comlf( tlfc, tlf )
 
 
*     This routine will compute the spacing in x and y to get a
*     "line feed" on the plot.  The routine uses the line feed in
*     characters and converts this to world coordinates.  The character
*     size used depends on the orientation of the labels being
*     written.  One of the TLFC values is assumed to zero i.e,
*     this routine will only move parallel to the axes.
 
      include 'plot_param.h'
      include 'plot_com.h'
 
*   cxy(2)  - Size of characters in the x and y directions
*           - rather than in width and height.
*   tlf(2)  - The line feed in world coordinates  (OUTPUT)
*   tlfc(2) - The line feed in characaters. (INPUT)
 
 
      real*4 cxy(2), tlf(2), tlfc(2)
      real*4 S   ! Scale for font
      S = 1.8
 
*     Get correct character sizes
 
*                                        ! label || to x axis
*     if( anint(cori(2)).eq.0 ) then
*         cxy(1) = charsz_x*1.5
*         cxy(2) = charsz_y*1.5
*                                 ! label || to y axis
*     else
*         cxy(1) = charsz_y*1.5
*         cxy(2) = charsz_x*1.5
*     end if
*
**** MOD TAH 131018: Use actual char size ! label || to x axis
      if( anint(cori(2)).eq.0 ) then
          cxy(1) = font_size(1)/3.5*S   ! 72 dpi=3.5 mm/pixel
          cxy(2) = font_size(2)/4.8*S
*                                 ! label || to y axis
      else
          cxy(1) = font_size(2)/4.8*S
          cxy(2) = font_size(1)/3.5*S
      end if

*     Now compute the line feed values
*                                      ! line feed in y direction
      if( anint(tlfc(1)).eq.0 ) then
          tlf(1) = 0.0
          call conv_mmtow(0.0,cxy(2)*tlfc(2), y_size_mm, tlf(2))
      else
          call conv_mmtow(cxy(1)*tlfc(1),0.0, x_size_mm, tlf(1))
          tlf(2) = 0.0
      end if
 
***** Thats all
      return
      end
 
