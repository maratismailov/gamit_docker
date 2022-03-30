 
      subroutine status
c
c     This routine will list the status of the parameters current set
c     PLTSL.
c
c
c Include files
c -------------
*                         ! the parameter file
      include  'plot_param.h'
c
*                         ! the control common block
      include  'plot_com.h'

      integer*4 trimlen
c
c
c.... START status output
      call report(' STATUS SUMMARY OF PLOT PARAMETERS')
 
      write(termlu,100) point_type, line_type, errb_type, 
     .      font_name(1:trimlen(font_name))
 100  format(" Point type ",i2,t20,"Line type ",2i3,t40,
     .       "Error bar type ",i2,t60,"Font name ",a10)
 
      write(termlu,150) errb_scale, pen_type
 150  format(" Error bar scale ",2(f8.3,1x),t40,"Pen number ",i2)
 
      write(termlu,200) view_size, x_size_mm, y_size_mm
 200  format(" VIEW        ",4f6.3,t40,"Paper size    ",2(f7.1," mm"))
c
      write(termlu,300) aspect_ratio, charsz_x, charsz_y
 300  format(" Aspect ratio",f6.3,t40,"Character size ",2(f7.1," mm"))
c
      write(termlu,400) x_field, y_field, p_field
 400  format(" X field columns ",3i3,t40,"Y field columns ",3i3,/,
     .       " P field columns ",2i3)
c
      call report_window(poly_window)
 
      call report_scales(scale_size)
c
      call sum_file
c
c.... Thats all
      return
      end
