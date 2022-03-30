      subroutine wdms (latlon,rad,string)

c     write a rad into string with format Xddd:mm:ss.sssss

c     input
c     latlon = 1 for latitude in a system in which N is reckoned positive
c             -1 for latitude in a system in which N is reckoned negative (WEIRD!)
c     latlon = 2 for longitude in a system in which E is reckoned positive
c             -2 for longitude in a system in which E is reckoned negative
c     latlon = 3 for radians converted simply to their ddd:mm:ss.sssss equivelant
c             -3 for radians converted with a sign change to their ddd:mm:ss.sssss equivelant

c     rad     latitude, longitude, or any angular value in radians

c     passed
      real*8 rad
      integer latlon
      character*16 string

c     local
      real*8 radian,absdeg,sec,deg
      integer ideg,imin
      character*1 flag

      radian = 45.0d0/datan(1.0d0)

      deg = rad * radian

      if (latlon .eq. -1) then
         if (deg .lt. 0.0d0) then
            flag = 'N'
         else
            flag = 'S'
         endif
      else if (latlon .eq. +1) then
         if (deg .le. 0.0d0) then
            flag = 'S'
         else
            flag = 'N'
         endif
      else if (latlon .eq. -2) then
c        positive reckoned West
         if (deg .le. -180.0d0) then
            flag = 'W'
            deg = deg + 360.0d0
         else if (deg .le. 0.0d0) then
            flag = 'E'
         else if (deg .le. 180.0d0) then
            flag = 'W'
         else
            flag = 'E'
            deg = deg - 360.0d0
         endif
      else if (latlon .eq. +2) then
c        positive reckoned East
         if (deg .le. -180.0d0) then
            flag = 'E'
            deg = deg + 360.0d0
         else if (deg .le. 0.0d0) then
            flag = 'W'
         else if (deg .le. 180.0d0) then
            flag = 'E'
         else
            flag = 'W'
            deg = deg - 360.0d0
         endif
      else if (latlon .eq. +3) then
c        direct conversion from radians to ddd:mm:ss.ssss
         if (deg.lt.0) then
            flag = '-'
         else
            flag = ' '
         endif
      else if (latlon .eq. -3) then
c        direct change sign of radians when converting to ddd:mm:ss.ssss
         if (deg.lt.0) then
            flag = ' '
         else
            flag = '-'
         endif
      endif
      absdeg = dabs(deg)
      ideg = absdeg
      imin = (absdeg - dble(ideg))*60
      sec  = (absdeg - dble(ideg) - dble(imin)/60.0d0 )*3600.0d0

      if (abs(latlon) .eq. 1) then
         write (string,10) flag,ideg,imin,sec
 10      format (1x,a1,i2.2,':',i2.2,':',f8.5)
      else if (abs(latlon) .eq. 2) then
         write (string,20) flag,ideg,imin,sec
 20      format (a1,i3.3,':',i2.2,':',f8.5)
      else if (abs(latlon) .eq. 3) then
         write (string,30) flag,ideg,imin,sec
 30      format (a1,i3.3,':',i2.2,':',f8.5)
      else
         write (string,40) rad*radian
 40      format (g16.10)
      endif

c     make 0.1 seconds print as 00.1
      if (string(9:9) .eq. ' ') string(9:9) = '0'

      return
      end


