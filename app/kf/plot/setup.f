      subroutine setup

      implicit none 
c
 
*     Initializes the GKS system. Using the OPNGKS routine in ncar
*     graphics.
c
c Include files
c -------------
*                            ! the parameter file
      include 'plot_param.h'
c
*                            ! the control common block
      include 'plot_com.h'

      integer*4 x0, y0,  wdt0, hgt0  ! Initial window position and size

c**** Set some values and open GKS
 
      aspect_ratio = 1.0
      x_size_mm = 180.0
      y_size_mm = 180.0
      x0 = 0 
      call jbegn( meta_file, x0, y0,  wdt0, hgt0  )
 
*     Turn off clipping.  We will turn on later during data plotting
      call gsclip(0)
 
      call plot_final_setup
 
      return
      end
 
