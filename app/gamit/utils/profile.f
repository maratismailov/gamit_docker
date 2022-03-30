      program profile

*=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
*  Generate output for plotting velocity profile.  Enter longitude,  
*  latitude and azimuth of the profile, and the width of the site swath. 
*
* Author:       R.A. Bennett,  Thu Sep 18 09:49:36 EDT 1997
* Modified by:  S.C. McClusky, Thu Sep 18 13:45:42 EDT 1997
* Modified by:  P. Vernant, Aug 2005, now rho is really used in the computation 
* Modified by:  S.C. McClusky, March 31 2011, now reads input .vel file in free format!
*=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

      implicit none

* Set dimensions
*   nsit        = maximum number of sites

      integer nsit
      parameter(nsit=5000) 

* Station parameters:
*   slat        = site latitude
*   slon        = site longitude
*   x1          = unrotated cartesian site coordinate -- east
*   x2          = unrotated cartesian site coordinate -- north
*   xr1         = rotated cartesian site coordinate -- orthog to Nuvel1
*   xr2         = rotated cartesian site coordinate -- Nuvel1 direction
*   name        = site name
*   nvel        = unrotated north component of site velocity 
*   evel        = unrotated east component of site velocity
*   r1vel       = component of site velocity orthog to Nuvel1 
*   r2vel       = component of site velocity in Nuvel1 direction
*   nsig        = unrotated uncertainty in north comp of vel
*   esig        = unrotated uncertainty in east comp of vel
*   r1sig       = rotated uncertainty orthog to Nuvel1 
*   r2sig       = rotated uncertainty in Nuvel1 direction
*   rho         = correlation coefficient
*   isit        = site index

      character*8 name(nsit),cd
      integer isit, iopt, iargc, ierr, jerr, terr, indx, trimlen
      real*8 theta, uvel, usig, pleng, pwidth, x, y, val
      real*8 rlon1, rlat1, rlon2, rlat2, rlon3, rlat3, rlon4, rlat4
      real*8 rlon5, rlat5, rlon6, rlat6, rlon7, rlat7, rlon8, rlat8
      real*8 slat(nsit), slon(nsit), rho(nsit)
      real*8 nvel(nsit),evel(nsit),hvel(nsit),r1vel(nsit),r2vel(nsit)
      real*8 nsig(nsit),esig(nsit),hsig(nsit),r1sig(nsit),r2sig(nsit)
      real*8 nadj(nsit), eadj(nsit), hadj(nsit)
      real*8 x1(nsit), x2(nsit), xr1(nsit), xr2(nsit)

* Profile coordinate system parameters
*   olat        = lat of origin
*   olon        = lon of origin
*   az          = Nuvel1 azimuth at origin of profile
*   v           = Nuvel1 velocity of origin of profile

      real*8 olat, olon, az, v

* Swath parameters
      real*8 width

* Misc parameters

      integer i_in1, iout
      real*8 pi, d2r
      character*20 arg, outname
      character*80 filename
      character*120 line
      
* Set fixed parameters 

      pi   = 4.0d0*atan(1.0d0)
      d2r  = pi/180.0d0
      olon = 244.0d0*d2r
      outname = 'profile.out'

* Assign file handles

      i_in1 = 12
      iout  = 13

* Get command line argument
                          
      if (iargc().lt.5)  then 
         write(*,*)' '
         write(*,*)' '
         write(*,*)'Program profile:'
         write(*,*)' '
         write(*,*)' Generate data for plotting velocity profiles'
         write(*,*)' '
         write(*,*)'Example:'
         write(*,*)'profile <lon> <lat> <vel> <width> <az>'
         write(*,*)' '
         write(*,*)'where <lon> <lat> = lon & lat of profile origin'
         write(*,*)'      <vel>       = velocity input file'
         write(*,*)'      <width>     = 1/2 width of profile swath (km)'
         write(*,*)'      <az>        = opt azimuth of prof normal'
         write(*,*)'                    0.0 for Nuvel-1'
         write(*,*)' '
         stop
      endif

* Read command line args and open files

      call rcpar(1,arg)
      read(arg,*) olon 
      if(olon.lt.0) olon = 360.0d0 + olon

      call rcpar(2,arg)
      read(arg,*)olat

      call rcpar(3,filename)
      open(i_in1, file = filename)

      call rcpar(4,arg)
      read(arg,*)width

      call rcpar(5,arg)
      read(arg,*)az

      iopt = 1

* Open output file
      open(iout, file = outname)
      write(iout,5)olon, olat,  width, az, iopt, filename
5     format ('profile  olon     olat   width   azmith  iopt filename',/
     .,'profile: ',f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,i2,4x,a80,/,
     .'---------------------------------------------------------------')
      write(iout,*)'  norm_dist   norm_vel   norm_sig   para_vel', 
     .'   para_sig   site'  

* Get Nuvel1 azimuth
      if(az.eq.0.0d0)then
         call nuvel1(olon, olat, v, az)
      endif

      theta = (360.d0-az)*d2r

      write(*,*)' '
      write(*,*)' Theta   =', theta
      write(*,*)' Azimuth =', az
      write(*,*)' '

*     Set site counter to 0
      isit = 0

*     Read input file
      if(iopt.eq.1)then
*     Loop over the file.
        do while ( ierr.eq.0 )
          read(i_in1,'(a)', iostat=ierr, end = 999 ) line
 
*         Process record if no error and the first character is blank
          if( ierr.eq.0  .and. line(1:1).eq.' ' .and.
     .        trimlen(line).gt.0                    ) then

              isit = isit + 1 

*             Initialize record index and error flags
              indx = 0
              terr = 0
              jerr = 0
*             Now each of the values from the line
              call read_line(line, indx, 'R8', jerr, slon(isit), cd)
              terr = terr + jerr
              call read_line(line, indx, 'R8', jerr, slat(isit), cd)
              terr = terr + jerr
              call read_line(line, indx, 'R8', jerr, evel(isit), cd)
              terr = terr + jerr
              call read_line(line, indx, 'R8', jerr, nvel(isit), cd)
              terr = terr + jerr
              call read_line(line, indx, 'R8', jerr, eadj(isit), cd)
              terr = terr + jerr
              call read_line(line, indx, 'R8', jerr, nadj(isit), cd)
              terr = terr + jerr
              call read_line(line, indx, 'R8', jerr, esig(isit), cd)
              terr = terr + jerr
              call read_line(line, indx, 'R8', jerr, nsig(isit), cd)
              terr = terr + jerr
              call read_line(line, indx, 'R8', jerr, rho(isit),  cd)
              terr = terr + jerr
              call read_line(line, indx, 'R8', jerr, hvel(isit), cd)
              terr = terr + jerr
              call read_line(line, indx, 'R8', jerr, hadj(isit), cd)
              terr = terr + jerr
              call read_line(line, indx, 'R8', jerr, hsig(isit), cd)
              terr = terr + jerr
*              Get the site name (tUadj is a dummy here: will not be
*             changed by call.
              call read_line(line, indx, 'CH', jerr, val,  name(isit))
              terr = terr + jerr

*        print waht we read!
100           format(2x,f7.3,2x,f7.3,1x,6f8.2,f7.3,2x,3f8.2,1x,a8)
              write(*, 100) slon(isit), slat(isit), evel(isit),
     .        nvel(isit), eadj(isit), nadj(isit), esig(isit),
     .        nsig(isit), rho(isit), hvel(isit), hadj(isit), 
     .        hsig(isit), name(isit)  

*             Now see if any errors. If there are error assume that
*             this is comment line, and ignore.
              if ( terr.ne.0 ) then
                 isit = isit - 1
                 print*,'invalid velocity record line - ignored:',line
              else

*        Transform coords to local cartesian and rotate into Nuvel-1 system
*        xr1 orthogonal to Nuvel-1 or specified direction
                 call geod2xy(olon,olat,slon(isit),slat(isit),
     .                x1(isit),x2(isit))
                 call vrot(x1(isit),x2(isit),xr1(isit),xr2(isit),theta)
*        Rotate Velocities
                 call vrot(evel(isit),nvel(isit),r1vel(isit),
     .                r2vel(isit),theta)
*        Rotate uncertainties
                 call erot(esig(isit),nsig(isit),rho(isit),r1sig(isit),
     .           r2sig(isit),theta)

*        Output Nuvel-1 parallel and perpendicular velocities
                 if(abs(xr2(isit)).lt.width)then
                    write(iout,150) xr1(isit),r1vel(isit),r1sig(isit), 
     .              r2vel(isit),r2sig(isit),name(isit)
150                 format(f11.4,4f11.5,3x,a8)

                 endif  

              endif

          endif

        enddo

      endif

      if(iopt.eq.2) then
 20      read(i_in1, '(a90)', end = 999) line
         if(line(1:1).ne.' ') goto 20
         isit = isit + 1
         read(line, 200) evel(isit), esig(isit), nvel(isit), nsig(isit),
     .        uvel, usig, rho(isit), slon(isit), slat(isit), name(isit)  
 200     format(6(1x,f7.2),1x,f7.4,2(1x,f9.4),1x,a8)
*        Transform coords to local cartesian and rotate into Nuvel-1 system
*        xr1 orthogonal to Nuvel-1 or specified direction
         call geod2xy(olon,olat,slon(isit),slat(isit),x1(isit),x2(isit))
         call vrot(x1(isit),x2(isit),xr1(isit),xr2(isit),theta)
*        Rotate Velocities
         call vrot(evel(isit),nvel(isit),r1vel(isit),r2vel(isit),theta)
*        Rotate uncertainties
         call erot(esig(isit),nsig(isit),rho(isit),r1sig(isit),
     .        r2sig(isit),theta)

*        Output Nuvel-1 parallel and perpendicular velocities
         if(abs(xr2(isit)).lt.width)then
            write(iout,*) xr1(isit), r1vel(isit), r1sig(isit), 
     .           r2vel(isit), r2sig(isit),'  ',name(isit)

         endif
         goto 20
      endif

 999  continue

*     Build file for plotting profile
      pleng = 200.d0
      pwidth = width 
      x = pleng*dsin((az+90.d0)*d2r)
      y = pleng*dcos((az+90.d0)*d2r)  
      call xy2geod(olon,olat,rlon1,rlat1,x,y)
      call xy2geod(olon,olat,rlon2,rlat2,-x,-y)
      x = pwidth*dsin(az*d2r)
      y = pwidth*dcos(az*d2r)
      call xy2geod(olon,olat,rlon3,rlat3,x,y)
      call xy2geod(olon,olat,rlon4,rlat4,-x,-y)
      x = pleng*dsin((az+90.d0)*d2r)
      y = pleng*dcos((az+90.d0)*d2r)
      call xy2geod(rlon3,rlat3,rlon5,rlat5,x,y)
      call xy2geod(rlon3,rlat3,rlon6,rlat6,-x,-y)
      call xy2geod(rlon4,rlat4,rlon7,rlat7,x,y)
      call xy2geod(rlon4,rlat4,rlon8,rlat8,-x,-y)

      open(25,file='profile_map')
      write(25,*)'>> profile', olon, olat,' ', filename, width, az, iopt
      write(25,*)">"
c      write(25,*)rlon2 , rlat2
      write(25,*)rlon3 , rlat3
      write(25,*)olon  , olat
c      write(25,*)rlon1 , rlat1
      write(25,*)rlon4 , rlat4
      write(25,*)">"
      write(25,*)rlon6 , rlat6
      write(25,*)rlon3  ,rlat3
      write(25,*)rlon5 , rlat5
      write(25,*)">"
      write(25,*)rlon8 , rlat8
      write(25,*)rlon4  ,rlat4
      write(25,*)rlon7 , rlat7

*     The End
      stop
      end

      subroutine geod2xy(o_lon,o_lat,lon,lat,x,y)
*     Subroutine to convert geodetic lat lon coords to x y coords  

      real*8 o_lat,o_lon,lat,lon,x,y,d2r
      real*8 scale

*     Degrees to radians
      d2r = acos(-1.d0)/180.d0

*     Correct For Latitudinal Dependence of Longitude Difference
      scale = d2r*cos(d2r*o_lat)
     .        *(6378.139d0*(1.0d0-sin(d2r*o_lat)
     .        *sin(d2r*o_lat)/298.247d0))

      x = (lon - o_lon)*scale
      y = (lat - o_lat)*111.32d0

      return
      end

      subroutine nuvel1(lon,lat, v, az)

* Initialize parameters and specify planet
      real*8 v,az,w,a,b,latp,longp,latx,longx,f,r,pi,ang1,ang2
      real*8 lon, lat, d2r
      parameter(a=6378.14d0,b=6356.75d0)

* NUVEL-1 numbers
      w=7.829e-07
      latp=48.709d0
      longp=-78.167d0

      pi=4.0d0*atan(1.0d0)
      d2r=pi/180.0d0
      latx=lat*d2r
      longx=lon*d2r
      latp=latp*d2r
      longp=longp*d2r
      w=w*d2r 
      f=(a-b)/a
      r=a*(1-f*sin(latx)*sin(latx))
      ang1=acos(sin(latx)*sin(latp)+
     .cos(latx)*cos(latp)*cos(longx-longp))
      ang2=asin(cos(latp)*sin(longp-longx)/sin(ang1))
      v=w*r*sin(ang1)*1000000.0d0
      az=270.d0+ang2/d2r

      return
      end


      subroutine vrot(v1,v2,vr1,vr2,angle)

      real*8 v1, v2, vr1, vr2, angle
      
      vr1 =  v1*cos(angle) + v2*sin(angle)
      vr2 = -v1*sin(angle) + v2*cos(angle)

      return
      end
      
      subroutine erot(e1,e2,rho,er1,er2,angle)

      real*8 e1, e2, rho, er1, er2, angle
      real*8 r11, r12, r21, r22

      r11 =  cos(angle)
      r12 =  sin(angle)
      r21 = -r12
      r22 =  r11
      er1 = sqrt(r11*r11*e1*e1+2*rho*r11*r12*e1*e2+r12*r12*e2*e2)
      er2 = sqrt(r21*r21*e1*e1+2*rho*r22*r21*e1*e2+r22*r22*e2*e2)
      
      return
      end

      
      subroutine xy2geod(o_lon,o_lat,lon,lat,x,y)
*     Subroutine to convert catesian x y coords to geodetic lat lon 

      real*8 o_lat,o_lon		
      real*8 x, y, lat, lon
      real*8 d2r, xscale, yscale, rad

*     Degrees to radians
      d2r = 4*atan(1.0d0)/180.0d0

*     Earth radius
      rad =  6.3674e+03

*     Correct For Latitudinal Dependence of Longitude Difference 
      xscale = d2r * rad * cos(d2r*o_lat)
      yscale = d2r * rad

      lon = (x/xscale) + o_lon
      lat = (y/yscale) + o_lat

      return
      end

