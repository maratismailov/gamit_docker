CTITLE INCREMENT_KALVER
 
      subroutine increment_kalver

      implicit none
 
*     This routine will increment the KalObs file version number.
*     It first checks to see if the current version is greater than
*     the database version.  If it is not then the version is intialized.
*     This would be first done in READIN normally, but for older KalObs
*     file, the first time they are updated.
*
      include '../includes/kalman_param.h'
      include '../includes/obs_header.h'
 
*     Check the current value of the version, and see if it needs
*     initialization.
      if( Kal_ver.lt. version*100 ) then
          Kal_ver = version*100
      else
          Kal_ver = Kal_ver + 1
      end if
 
*     Thats all
      return
      end
 
 
