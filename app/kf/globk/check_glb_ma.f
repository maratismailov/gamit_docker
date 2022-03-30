CTITLE CHECK_GLB_MAX
 
      subroutine check_glb_max( type, num, max )
 
      implicit none 

 
*     Routine to check if we have exceeded the maximum values allowed
*     for various arrays in the Global Kalman filter.
*
*     The program stops, if num is greater than max.
 
*   iout    - User's terminal
*   loglu   - HP function for user's LU
*   id      - Dummy argument for LogLu
*   num     - Current numbers of TYPE
*   max     - Maximum number of TYPE
 
      integer*4 iout, loglu, id, num, max
 
*   type    - Type of quantity we have exceeded.
 
      character*(*) type
 
***** See if we have to many
*                                 ! Too many, report and stop
      if( num.gt.max ) then
 
          iout = loglu(id)
          write(iout,100) type, max, num, type
  100     format(/' **** DISASTER **** Number of ',a,' exceeded.',/,
     .            ' Maxiumim number allowed is ',i8,' and ', i8,
     .            ' requested',/,
     .            ' Either reduce number of ',a,' or modify',
     .            ' Kalman_param.ftni'/)
          stop ' RUN TERMINATED: array sizes exceeed'
      end if
 
***** Thats all
      return
      end
 
