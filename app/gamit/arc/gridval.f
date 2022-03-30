      Subroutine gridval (xgridmin,ygridmin,cellsize
     , ,nlonrow,nlatcol,grid,lat_deg,lon_deg,val)
C     Get the correct 4 values from a UCL SRP grid for a particular 
C     lat/long of radiation and interpolate them to give an output value

C     In:
C         xgridmin      xvalue, lower left corner of entire grid
c         ygridmin      yvalue, llc of entire grid
c         cellsize      size of grid cells, same for x and y
C         nlonrow       number of rows in grid (varying longitude)
C         nlatcol       number of columns in the grid (varying latitude)
c         grid          array holding grid values (currently sized as
c                       (100,120) in arc.h)
C         lat_deg, lon_deg latitude and longitude of the sat-sun vector
C                       in the satellite body fixed system in degrees
C
C     Out:
c         val           value interpolated from the grid

      
      implicit none

      real*8 xgridmin,ygridmin,cellsize,grid(100,120)
     .    ,  lat_deg,lon_deg
     .    ,  val
     .    ,  x0,y0
     .    ,  u1,u2,u3,u4
     .    ,  xlatfactor,ylonfactor 
     .    ,  cellsizey
C               grid row, number of rows
      integer*4 i,i1,nlonrow,ioerr
C               grid column, number of columns
     .          ,j,nlatcol

      val=0.d0

C         Find the coords of the lower left corner of the grid cell$
C         containing the point and get the values at the corners of
C         the grid cells for x, y, and z
C         The UCL grid files are ESRI grid file format. This means that the
C         values for the first column have x coordinates of xgridmin +
C         cellsize/2
c         x0,y0 :  the coords of the lower left point of the points
C                  the point to interpolate is within.
c         u1    : value at lower left point     u4 ----- u3
c         u2    : value at lower right point     |       |
c         u3    : value at upper right point     |       |
c         u4    : value at upper left point   y0 u1 ---- u2
C                                                x0
C         Find x0 (dint truncates to whole number)
Cd          print*,'latdeg',lat_deg          
          if (lat_deg.le.-90.d0 ) then
               lat_deg=lat_deg+180.d0
          elseif (lat_deg.gt.91.d0) then        
               lat_deg=lat_deg-180.d0
          endif
C         Check
          if ((lat_deg.le.-90.d0).or.(lat_deg.gt.91.d0)) then
             call report_stat('FATAL','ARC','gridval',lat_deg,
     .            'Strange value of lat_deg',ioerr)  
          endif        
Cd          print*,'latdeg',lat_deg          
          xlatfactor=(lat_deg-(xgridmin+(cellsize/2.d0)))/cellsize
          x0=xgridmin+cellsize/2.d0+cellsize*dint(xlatfactor)
C         Find y0
C         Find x (lat) values
          j=dint(xlatfactor)+1
C         Find y0 (problem - grid is not actually overlapping for
C         longitude
          cellsizey=cellsize
C         Adjust so within the 360 degrees up from
          if (lon_deg .lt.(ygridmin+cellsize/2.d0)) then
                  lon_deg=lon_deg+360.d0
          elseif (lon_deg .gt.( ygridmin+cellsize/2.d0+360)) then  
                  lon_deg=lon_deg-360.d0
          endif
C         Check
          if ((lon_deg .lt.(ygridmin+cellsize/2.d0)) .or.
     .       (lon_deg .gt.(ygridmin+cellsize/2.d0+360))) then
             call report_stat('FATAL','ARC','gridval',lon_deg,
     .            'Strange value of lon_deg',ioerr)  
          endif   
C         If the value is out of the grid, adjust the cellsize
          ylonfactor=(lon_deg-(ygridmin+(cellsizey/2.d0)))/cellsizey
          if (dint(ylonfactor).ge.nlonrow-1) then 
                 cellsizey=360.d0-cellsize*float(nlonrow-1)
          endif
          y0=ygridmin+cellsize/2.d0+cellsize*dint(ylonfactor)
C         find y (lon) values (smallest (-ve) value at bottom of array)
          i=nlonrow-dint(ylonfactor)
          i1=i-1
          if (i.le.1) then
                  i1=nlonrow
          endif
C         look up the values at the coordinates
          u1=grid(i,j)
          u2=grid(i,j+1)
          u3=grid(i1,j+1)
          u4=grid(i1,j)
Cd          print*,'lon_deg',lon_deg
Cd          print*,'lat_deg',lat_deg
Cd          print*,'maxlon',ygridmin+cellsize/2.d0+360
Cd          print*, 'u1-4',u1,u2,u3,u4
Cd          print*,'x0',x0,'y0',y0
Cd          print*,'cellsizey',cellsizey
Cd          print*,'xlatfactor',xlatfactor,dint(xlatfactor)
Cd          print*,'ylonfactor',ylonfactor,dint(ylonfactor)
Cd          print*,'i',i,'j',j,'i1',i1
C         Call the bilin8 (copied from ../grdtab/bilin4 ) but real*8) bilinear 
C         interpolation subroutine
Cd          print*,lat_deg,lon_deg,x0,y0,cellsize,cellsizey
Cd     .     ,u1,u2,u3,u4,val
         call bilin8(lat_deg,lon_deg,x0,y0,cellsize,cellsizey
     .     ,u1,u2,u3,u4,val)
         print*,'val',val

      end

c common/uclgrd/
c ,busX3xmin,busX3ymin,busX3cellsize,busX3grid(100,120)
c ,busY3xmin,busY3ymin,busY3cellsize,busY3grid(100,120)
c ,busZ3xmin,busZ3ymin,busZ3cellsize,busZ3grid(100,120)
c ,busX4xmin,busX4ymin,busX4cellsize,busX4grid(100,120)
c ,busY4xmin,busY4ymin,busY4cellsize,busY4grid(100,120)
c ,busZ4xmin,busZ4ymin,busZ4cellsize,busZ4grid(100,120)
c ,pnlX3xmin,pnlX3ymin,pnlX3cellsize,pnlX3grid(100,120)
c ,pnlY3xmin,pnlY3ymin,pnlY3cellsize,pnlY3grid(100,120)
c ,pnlZ3xmin,pnlZ3ymin,pnlZ3cellsize,pnlZ3grid(100,120)
c ,pnlX4xmin,pnlX4ymin,pnlX4cellsize,pnlX4grid(100,120)
c ,pnlY4xmin,pnlY4ymin,pnlY4cellsize,pnlY4grid(100,120)
c ,pnlZ4xmin,pnlZ4ymin,pnlZ4cellsize,pnlZ4grid(100,120)
