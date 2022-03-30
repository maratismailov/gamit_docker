CTITLE DECODE_GLFOR_RUN
 
      subroutine decode_glfor_run

      implicit none  
 
*     Routine to decode the glfor runstring.  Only thing passed
*     is the name of the common.
 
      include '../includes/kalman_param.h'
      include '../includes/globk_cntl.h'
 
*   len_run     - Length of runstring
*   rcpar       - Function to return runstring
 
      integer*4 len_run, rcpar
 
****  Get the runstring
 
      len_run = rcpar(1, glb_com_file )
      if( len_run.eq.0 ) then 
          write(*,100) 
 100      format('***DISASTER*** No common file passed')
          stop 'GLFOR: No common file name'
      end if
 
***** Thats all
      return
      end
 
