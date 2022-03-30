      subroutine conetrans(slong,slat,numsta,sx,sy,long0,lat0)
      
      integer numsta, indx
      real*4 slong(numsta),slat(numsta),sx(numsta),sy(numsta),
     .     maxlat, minlat,maxlon,minlon,long0,lat0
      
      
C     loop through stations, find centroid of the network
      minlon = slong(1)
      maxlon = slong(1)
      minlat = slat(1)
      maxlat = slat(1)

C     write(*,*) 'Translating longitudes and latitudes to x and y'
C     write(*,*) 'coordinates using Murrays''s Polyconic projection'
C     write(*,*) 
C     write(*,*) 'Longitudes and latitudes:'
C     write(*,'(f12.6, f12.6)') slong(1),slat(1)

      do 1001 indx = 2, numsta
         if (slong(indx).gt.maxlon) then
            maxlon = slong(indx)
         elseif (slong(indx).lt.minlon) then
            minlon = slong(indx)
         endif
         if (slat(indx).gt.maxlat) then
            maxlat = slat(indx)
         elseif (slat(indx).lt.minlat) then
            minlat = slat(indx)
         endif
C        write(*,'(f12.6,f12.6)') slong(indx),slat(indx)
 1001 continue

      write(*,*) 'minlon ', minlon, '  maxlon ', maxlon
      write(*,*) 'minlat ', minlat, '  maxlat ', maxlat
         
C     set this equal to long0,lat0
      long0 = (minlon + maxlon)/2
      lat0 = (minlat + maxlat)/2

      write(*,*) 
C     write(*,*) 'X and Y coordinates:'
C     loop through stations, send each to poly4
      do 1002 indx = 1,numsta
         call poly44(long0,lat0,slong(indx),slat(indx),
     .        sx(indx),sy(indx))
 1002 continue

      return
      end

************************************************************************
      subroutine poly44(long0, lat0, longp, latp, px, py)
c     subroutine poly44(p1,p2,il,px,py)
      
c     polyconic projection of point lat=p2, diff long=il from arbitrary
c     central meridian. lat of arbitrary origin is p1. px=dist from cm
c     along lat p2. py=dist from p1 to p2. px,py in meters.
c     p1,p2,and il in seconds.
      
      real*4 long0, lat0, longp, latp, px, py, 
     *     il,la,ip,ipr, p1, p2, theta, sinp2, cosp2, a, cot, 
     *     pr, arcone, esq, a0, a2, a4, a6, a8
      
C     arcone = pi / 180.0 / 60.0 / 60.0
C     esq = ??
      
      data arcone,esq,la,a0,a2,a4,a6,a8/4.8481368e-6,6.7686580e-3,
     *     6378206.4,6367399.7,32433.888,34.4187,.0454,6.0e-5/

      p2 = latp * 3600.0
      p1 = lat0 * 3600.0
      il = (longp - long0) * 3600.0

      ip=p2-p1
      sinp2=sin(p2*arcone)
      cosp2=cos(p2*arcone)
      theta=il*sinp2
      a=sqrt(1.0-(esq*(2.*sinp2)))/(la*arcone)
      cot=cosp2/sinp2
      px=(cot*sin(theta*arcone))/(a*arcone)
      ipr=ip*arcone
      pr=((p2+p1)/2.)*arcone
      py= a0*ipr-(a2*cos(2.*pr)*sin(ipr))+(a4*cos(4.*pr)*sin(2.*ipr))-
     *     (a6*cos(6.*pr)*sin(3.*ipr))+a8*cos(8.*pr)*sin(4.*ipr)

      px = px / 1000.0
      py = py / 1000.0
C     write(*,'(f12.6,f12.6)') px,py

      return
      end
