      real*8 function timdif (jd2,sod2,jd1,sod1)
c
c     return difference in seconds between epoch 2 minus epoch 1
c     handles going over days and positive and negative diffs
c
C     Julian Day numbers
      integer*4   jd1,jd2
C     second of day
      real*8      sod1,sod2

      real*8      spd

C     seconds per day
      data        spd/86400.d0/

      timdif = spd*(jd2-jd1) + sod2-sod1

      return
      end
