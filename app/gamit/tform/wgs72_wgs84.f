C
C  Converts WGS72 to WGS84 given geodetic coordinates in dms,
C     Source:  DMA TR 8350.2
C
      implicit real*8 (a-h,q-z)
c
      delf = -0.8120450E-07
      a = 6378145.0d0
      dela = -8.0d0
      delr = -3.80d0
C
      radcon = datan(1.d0)/45.d0
      read(5,*) lat,latm,seclat,lon,lonm,seclon,hgt
      asec = radcon/3600.d0
      alat = (seclat/60.d0 + dble(latm))/60.d0 + dble(lat)
      alon = (seclon/60.d0 + dble(lonm))/60.d0 + dble(lon)
      alat = alat*radcon
      alon = alon*radcon
c
      dellat = (4.5d0*cos(alat))/(a*sin(asec)) +
     .         delf*sin(2.d0*alat)/sin(asec)
      dellon = 0.814
      delh = 4.5d0*sin(alat) + a*delf*sin(alat)*sin(alat)
     .       - dela + delr
c
      print*, lat,latm,seclat+dellat,lon,lonm,seclon+dellon,hgt+delh
      end
