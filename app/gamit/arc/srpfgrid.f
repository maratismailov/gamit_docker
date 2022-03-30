       Subroutine srpfgrid( ssvec,srpxacc,srpyacc,srpzacc)
C      Calculate accels on satellite in a satellite XYZ frame
C      using the UCL grid models, at a nominal distance to the sun, and 
C      no allowance for eclipsing. The six grid files for each satellite type 
C      (x,y,z satellite bus, x,y,z solar panels) will be read into memory
C      earlier by arc, to avoid opening the files at each inegration step.

C      Called by ertorb.f
C      EJP 13 April 2012

C      IN:  
C        ssvec    satellite->sun unit vector in the satelite frame
c        These in ../includes/arc.h:
c        antbody  SV body type (works for all GNSS)  (arc.h)
C        twopi    two*pi 
C        sbmass   satellite mass /kg?

           


C       OUT:
C         srpxacc, srpyacc, srpzacc   /km s**-2

C       Parameters:
           

       implicit none
       include '../includes/dimpar.h'
       include '../includes/arc.h'
                                    
       real*8 ssvec(3),srpxacc,srpyacc,srpzacc
     .    , lat_deg,lon_deg,lat_rad,lon_rad
     .    , lat_deg_pnl,lon_deg_pnl,lat_rad_pnl,lon_rad_pnl
     .    , xaccbus,yaccbus,zaccbus
     .    , xaccpnl,yaccpnl,zaccpnl

C       Apply the correct accelerations for satellite type

       if( antbody(1:9).eq.'BLOCK IIA'.or.
     .      antbody(1:9).eq.'BLOCK IIR ') then
C         From satellite sun BFS frame unit vector, work out the 
C         latitude and longitude of the vector in the satellite XYZ frame  
          lat_rad  = atan2(ssvec(3),sqrt(ssvec(1)**2.d0+
     .           ssvec(2)**2.d0))
          lon_rad = atan2(ssvec(2),ssvec(1))
C         convert to degrees to be compatible with grid file
          lat_deg = lat_rad*360.d0/twopi
          lon_deg = lon_rad*360.d0/twopi

C         Pick the correct grids for the satellite and find the value
C          If assume cellsize and origin same for x,y,z grids could
C          ondense this
          if( antbody(1:9).eq.'BLOCK IIA' ) then 
C          Block IIA (3) bus      
            call gridval(busX3xmin,busX3ymin,busX3cellsize
     .        ,bus3nrows,bus3ncols,busX3grid
     .        ,lat_deg,lon_deg,xaccbus)    
            call gridval(busY3xmin,busY3ymin,busY3cellsize
     .        ,bus3nrows,bus3ncols,busY3grid
     .        ,lat_deg,lon_deg,yaccbus) 
            call gridval(busZ3xmin,busZ3ymin,busZ3cellsize
     .        ,bus3nrows,bus3ncols,busZ3grid
     .        ,lat_deg,lon_deg,zaccbus)
C          Block IIA (3) panels
            call gridval(pnlX3xmin,pnlX3ymin,pnlX3cellsize
     .        ,pnl3nrows,pnl3ncols,pnlX3grid
     .        ,lat_deg_pnl,lon_deg_pnl,xaccpnl)    
            call gridval(pnlY3xmin,pnlY3ymin,pnlY3cellsize
     .        ,pnl3nrows,pnl3ncols,pnlY3grid
     .        ,lat_deg_pnl,lon_deg_pnl,yaccpnl) 
            call gridval(pnlZ3xmin,pnlZ3ymin,pnlZ3cellsize
     .        ,pnl3nrows,pnl3ncols,pnlZ3grid
     .        ,lat_deg_pnl,lon_deg_pnl,zaccpnl)
          else 
	     call report_stat('WARNING','ARC','srpfgrid',
     .         ' ','Grid files available only for GPS Block IIA'
     .            ,0 )
C        Don't have grid files yet so set accels to zero
                srpxacc=0.d0
                srpyacc=0.d0
                srpzacc=0.d0
          endif  

C         Combine the x accels for the bus and solar panels, repeat for y,
C         z, adjust for mass of satellite and return the values to ertorb
C         Grid files have values in units of metres per sec squared (check
C         this again with Stuart), so convert to km per sec squared.
C         Grids were calculated with a nominal mass of ?? kg and a mean solar irradience of ???
C          srpxacc=1.d-3*(xaccbus+xaccpnl)*nominal.d0/sbmass
C          srpyacc=1.d-3*(yaccbus+yaccpnl)*nominal.d0/sbmass
C          srpzacc=1.d-3*(zaccbus+zaccpnl)*nominal.d0/sbmass       
                    srpxacc=1.d-3*(xaccbus+xaccpnl)
                    srpyacc=1.d-3*(yaccbus+yaccpnl)
                    srpzacc=1.d-3*(zaccbus+zaccpnl)    
        else 
                srpxacc=0.d0
                srpyacc=0.d0
                srpzacc=0.d0
        endif
        return
        end
