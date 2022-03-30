      subroutine read_ref_net(ifil1)
c
c     read reference network file (assuming free format for the time being)
c                              Danan Dong 07/24/95
c
c     coordinate system : (comode)
c         1. geocentric spherical coordinate (lat,lon,radius)
c         2. geodetic (ellipsoidal) coordinate (lat,lon,height)
c         3. geocentric Cartesian coordinate (x,y,z)
c         4. topocentric coordinate (x(east),y(north),z(height))
c         Usually, we name the geocentric coordinate as 'global' frame.
c         Meanwhile, the topocentric coordinate is named as 'local' frame.
c
c     unit:
c         x, y, z : m
c        vx,vy,vz : mm/year
c        time     : year
c
      implicit real*8(a-h,o-z)
      include 'maked.fti'

      character*8 name
      character*256 line
      dimension ax1(10)
      integer   ifil1,i,k,k1
      integer   nblen
c     
c     default format: free format (geocentric coordinate)
c      if (netfmt(1:1).eq.'*')
c     . netfmt(1:40) = '(a1,f7.3,1x,f7.3,6f8.1,f7.3,3f8.1,2x,a8)'
c
      n_ref = 0
      do 20 i = 1,10000
c        shift format line by line
         read (ifil1,'(a)',err=20,end=50) line 
         k1 = nblen(line)
c        skip comment line
         if (line(1:1).ne.' '.or.k1.le.8) goto 20
         read (line,*,err=20)
     .        name,(ax1(k),k=1,7)
         n_ref = n_ref+1
         ref_nam(n_ref) = name
         ref_x(n_ref) = ax1(1)
         ref_y(n_ref) = ax1(2)
         ref_z(n_ref) = ax1(3)
         ref_vx(n_ref) = ax1(4)
         ref_vy(n_ref) = ax1(5)
         ref_vz(n_ref) = ax1(6)
         ref_tim(n_ref) = ax1(7)
c
 20   continue
c
c10   format (a1,f7.3,1x,f7.3,6f7.1,f8.3,3f7.1,1x,a8)
c10   format (a1,f7.3,1x,f7.3,6f8.1,f7.3,3f8.1,2x,a8)
c10   format (a1,f8.3,1x,f8.3,1x,6f7.1,f7.3,2x,3f7.1,1x,a8)
 30   format (1x,a8,1x,a12,1x,a1,2(i2,1x),f8.5,1x,
     .        a1,i3,1x,i2,1x,f8.5,f13.4,3f8.4,f9.3,3f8.4)
 40   format (1x,2f11.5,2f9.3,2x,a8)
c
 50   continue
      return
      end
