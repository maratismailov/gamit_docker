CTITLE OUTPUT_LABELS
 
      subroutine output_labels
 
 
*     Routine to write out all of the labels which have been saved
*     from axis routines and labeling routines.  This subroutine should
*     be in a segment by itself
 
      include 'plot_param.h'
      include 'plot_com.h'
 
*   end     - End position of current label
*   i       - Loop counter
*   start   - start position of current label minus one
 
      integer*4 end, i, start
 
*   or_x, or_y  - Orienation of baseline
*   pl_x, pl_y  - Orientation of perpendicular to baseline
 
      real*4 or_x, or_y, pl_x, pl_y
 
***** Loop over all of the labels which have been saved
 
      do i = 1, num_labels
 
*         Get the orientation of labels
          or_x =  label_orient(1,i)
          or_y =  label_orient(2,i)
          pl_x = -label_orient(2,i)
          pl_y =  label_orient(1,i)
 
          call jcori(or_x,or_y,0.0, pl_x,pl_y, 0.0)
 
*         Get the justification of the label
          call jjust(label_just(1,i), label_just(2,i))
 
*         Set the cursor position
          call j2mov(label_pos(1,i), label_pos(2,i))
 
*         Now extract the label to be written
          if( i.eq.1 ) then
              start = 0
          else
              start = label_char(i-1)
          end if
 
          end = label_char(i)
 
*                                              ! Really an error somewhere
          if( start+1.gt.end ) start = end - 1
*                                              ! this will stop 601 error.
 
*****     Now write the label
          call write_label( label_all(start+1:end) )
 
      end do
 
****  Now reset the label data
      num_labels = 0
      label_all  = ' '
 
****  Thats all
      return
      end
 
