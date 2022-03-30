CTITLE 'GSTUT'
 
      SUBROUTINE GSTUT(IM,ID,IY,FRACT,GST)

      implicit none 
C
C     T.HERRING                                 4 MARCH 1981
C
C     SUBROUTINE TO COMPUTE GST (rads) GIVEN UT AND DATE
C
C# LAST COMPC'ED  810311:14:23                  #
C
C
*         im, id, iy              - M D Y
*   ,   date(5)                   - M D Y H M
 
      integer*4 im, id, iy, date(5)
 
 
      real*8 fract, gst
 
 
      real*8 fjd, pi, t, ratio, gstd
 
C     DOUBLE PRECISION FJLDY,FJD,PI,T,RATIO,GSTD
C
      DATA PI / 3.141592654D0 /
C
C**** COMPUTE JULIAN DAT, USING ROUTINE FROM &SUTIL
      date(1) = iy
      date(2) = im
      date(3) = id
      date(4) = 0
      date(5) = 0
      call YMDHMS_to_JD(date,0.0D0,fjd)
      T   = FJD - 2415020.0D0
      RATIO = 1.002737909
C
C**** COMPUTE GST
      GSTD = .173993591D1 + T*(.17202791266D-1 + T*.50640897D-14)
      GST = GSTD + (RATIO*FRACT- IDINT(GSTD/2./PI)) *2.*PI
 
      end
 
