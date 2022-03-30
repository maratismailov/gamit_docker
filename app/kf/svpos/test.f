
      implicit none
 
*     Program to compute the rise and set times of the GPS
*     satelites based on the ephemeris in the rinex navigiation
*     files
 
*     The runstring of the program is
*     % svsp3 [sp3 orbit file] [data file]
*     where nav_file is nave file name
*         data_file is the name of the Rinex data file.
 
*     include '../includes/const_param.h'
      include '../includes/sp3_def.h' 
*     include 'svsp3.h'
 
* Main pogram variables

      integer*4 i, ep
 
      print *,' EP ',ep

      end

