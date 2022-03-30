      subroutine centroid (dlats,dlons,dlat0,dlon0,n)

*     Kurt Feigl April, 1991
*
*
*     Find the centroid in Mercator coordinates.
*
*     input:  arrays of coordinates in degrees: dlats,dlons
*             number of points, n
*     output: centroid in degrees dlat0,dlon
*
*
*

c     lat,lon of data points in degrees
      real*8 dlats(*),dlons(*)
c     lat,lon of centroid
      real*8 dlat0,dlon0
c     number of points 
      integer n

c     local variables
      real*8 rat,rat0,dlatmin,dlatmax,ron,ron0,dlonmin,dlonmax,
     .     ratc,ronc,degrad,b0,rscale,d0

      integer i

c     Mercator X,Y of data points in rubber units and sums
      real*8 xm,ym,xsum,ysum,xc,yc
      
c     prevent disaster
      if (n .lt. 0) then
         write (*,*) 'CENTROID: n .le. 0. ',n
         return
      endif

      degrad = datan(1.0d0)/45.0d0

c     initialize extrema
      dlatmax = -90.
      dlatmin =  90.
      dlonmax = -360.
      dlonmin =  360.
      xsum = 0.
      ysum = 0.

c     find extreme values
      do i = 1,n
         if (dlats(i) .gt. dlatmax) dlatmax = dlats(i)
         if (dlats(i) .lt. dlatmin) dlatmin = dlats(i)
         if (dlons(i) .gt. dlonmax) dlonmax = dlons(i)
         if (dlons(i) .lt. dlonmin) dlonmin = dlons(i)
      enddo

c     convert to Transverse Mercator coordinates and sum
c     origin is at (rat0,ron0)
      rat0 = degrad*(dlatmin+dlatmax)/2.0d0 
      ron0 = degrad*(dlonmin+dlonmax)/2.0d0 
      rscale = 6378.2064d0
      do i=1,n
         rat = dlats(i)*degrad
         ron = dlons(i)*degrad
         b0 = cos(rat) * sin(ron-ron0)
         xm = 0.50d0 * rscale * log ((1.0d0+b0)/(1.0d0-b0))
         ym = rscale*(datan(dtan(rat)/dcos(ron-ron0))-rat0)
         xsum = xsum + xm
         ysum = ysum + ym
      enddo

c     calculate centroid
      xc = xsum/n
      yc = ysum/n

c     inverse formula for Transverse Mercator Projection
      d0 = yc/rscale + rat0
      ratc = dasin(sin(d0)/cosh(xc/rscale))
      ronc = ron0 + datan(sinh(xc/rscale)/dcos(d0))

c     convert back to degrees
      dlat0 = ratc / degrad
      dlon0 = ronc / degrad

      return
      end


 
 
 
 
