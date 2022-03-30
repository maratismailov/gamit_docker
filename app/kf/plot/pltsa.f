CTITLE
 
      subroutine pltsa
c
c     This segment does all of the high quality labeling work for
c     PLOT.  The axis routines and put_label save the label information
c     and then this segment is called to output the labels.
c
c Include files
c -------------
*                          ! the parameter file
      include 'plot_param.h'
c
*                          ! the common block
      include 'plot_com.h'
c
c
      call output_labels
 
*                   ! flush the plot buffers
      call jmcur
*                   ! return to main segment
      return
c
      end
 
