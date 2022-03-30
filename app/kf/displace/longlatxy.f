************************************************************************

      subroutine transxy(long0, lat0, long, lat, localx, localy, 
     .     number)

C	This subroutine finds the x and y coordinates in kilometers
C	of the points in long and lat relative to the origin
C	long0 and lat0.  Number is the number of points in 
C	long and lat.

      integer
     +         number, n

      real
     +      long0, lat0, long(*), lat(*), localx(*), localy(*),
     +      PI, radius

      data radius, PI  / 6378.153, 3.141592654 /

      do 300 n = 1, number
        localx(n) = (long(n)-long0) * PI/180.0 * radius *
     +                 cos(PI/360.0 * (lat(n)+lat0))
        localy(n) = (lat(n)-lat0) * PI/180.0 * radius
300   continue

      return
      end

************************************************************************

      subroutine transcoord(long0, lat0, long, lat, x, y, number)

C	This subroutine finds the long and lat of the points 
C	in x and y (in kilometers) relative to the origin at
C	long0 and lat0.  Number is the number of points in 
C	x and y.

      integer
     +         n, number

      real 
     +      long0, lat0, long(*), lat(*), x(*), y(*),
     +      PI, radius

      data radius, PI  / 6378.153, 3.141592654 /

      do 300 n = 1, number
        lat(n) = (y(n) * 180 / (PI * radius)) + lat0
        long(n) = ((180 * x(n)) / (PI * radius * cos(PI/360 * 
     +            (lat(n) + lat0)))) + long0
300   continue

      return
      end

************************************************************************

      subroutine translatexy(fx, fy, x, y, tx, ty, number)

C	This subroutine translates the x and y coordinates, expressed
C	in kilometers relative to long0 and lat0 to tx and ty
C	coordinates expressed in kilometers relative to a new origin
C	at long and lat.

      integer
     +         n, number
      real
     +       fx, fy, x(*), y(*), tx(*), ty(*)

      do 310 n = 1, number
        tx(n) = x(n) - fx
        ty(n) = y(n) - fy
310   continue

      return
      end

************************************************************************

      subroutine rotatexy(strike, tx, ty, rx, ry, number)

      integer
     +         n, number

      real
     +       dist, alpha, beta, strike, 
     +       tx(*), ty(*), rx(*), ry(*),
     +       PI

      data PI  / 3.141592654 /

      do 350 n = 1, number
        if ((tx(n).eq.0).and.(ty(n).eq.0)) then 
          rx(n) = 0
          ry(n) = 0

        else
          dist= sqrt(tx(n)**2 + ty(n)**2)
          alpha= atan2(ty(n),tx(n)) * 180.0/PI
          beta= 180.0 - alpha - strike
          rx(n)= dist * sin(beta * PI/180.0)
          ry(n)= dist * cos(beta * PI/180.0)

        endif
350   continue

      return
      end

************************************************************************

      subroutine unrotate(strike, deltan, deltae, number)

      integer
     +         n, number

      real
     +      dist, strike, 
     +      tempn, tempe,
     +      deltan(*), deltae(*),
     +       PI

      data PI  / 3.141592654 /

      do 550 n = 1, number
        dist = sqrt(deltan(n)**2 + deltae(n)**2)
        tempn = dist * sin ( (180.0 - strike) * PI/180.0 
     +             - atan2( deltae(n),deltan(n) ) )
        tempe =  dist * cos ((180.0 - strike) * 
     +               PI/180.0-atan2(deltae(n),deltan(n)))
        deltan(n) = tempn
        deltae(n) = tempe
550   continue
     
      return
      end

************************************************************************
