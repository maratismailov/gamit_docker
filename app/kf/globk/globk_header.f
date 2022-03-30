CTITLE GLOBK_HEADER
 
      subroutine globk_header ( icrt, version )

      implicit none 
 
c
c     routine to write the header information to the user terminal
c
*   icrt    - User's LU number
*   version - number of the program (passed from Kalman_param.h
*   hsver   - Function that returns version with H or S added 
*             depending on system type.
 
      integer*4 icrt, trimlen
      character*(*) version
      character*12 hsver, full_ver

      full_ver = hsver(version) 
 
      write(icrt,100) full_ver(1:trimlen(full_ver))
 100  format(///" +++++++++++++++++++++++++++++++++++++++++++++",/,
     .          " + GLOBAL KALMAN FILTER             Ver ",a," +",/
     .          " +++++++++++++++++++++++++++++++++++++++++++++",//)
c
      return
      end
 
 
