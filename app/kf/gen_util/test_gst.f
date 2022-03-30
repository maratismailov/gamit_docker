      program test_gst

      implicit none 

*     Routine to test the two version of GST I have and see
*     if they differ.
*
*     SUBROUTINE GST_JD(jd, GST)
*
*     SUBROUTINE GSTUT(IM,ID,IY,FRACT)
*
      include '../includes/const_param.h'

      integer*4 im, id, iy, date(5)
      real*8 jd, gst1, gst2, sectag, fract


*     Do a values
      date(1) = 89
      date(2) =  2
      date(3) = 17
      date(4) = 15
      date(5) = 20
      sectag  = 10.d0
      call ymdhms_to_jd( date, sectag, jd )

      fract = date(4)/24.d0 + date(5)/1440.d0 + sectag/86400.d0

      call gstut(date(2), date(3), date(1), fract, gst1)
      call gst_jd( jd, gst2 )
      gst1 = mod(gst1,2*pi)
      gst2 = mod(gst2,2*pi)

      write(*,* ) ' JD, gst_ut, gst_jd ',jd, gst1, gst2
      write(*,*) ' Difference (secs) ', (gst1-gst2)/(2*pi)*86400.d0

      end


      
