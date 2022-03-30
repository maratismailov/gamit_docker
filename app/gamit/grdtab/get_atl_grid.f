      Subroutine get_atl_grid( slat,slon,sdoy )

c     Interpolate a binary grid file to get values for atmospheric tidal loading 
c     for a particular station.  Adapted from P. Tregoning utils/interp_atmtide.f.
c     Written by R. King  10 August 2006, modified by R. King 31 January 2013 to 
c     read a time-dependent grid file written at MIT from Tonie Van Dam's 12 monthly 
c     grid files using atl_asc2bin.       

c     The ANU grid is arranged in 73 latitude rings at 2.5 deg spacing 
c     from 90N to 90S and 720 longitude values within each ring from 
c     0E to 360-dlatE.  There are 12 real*4 coefficients in each record.

c     The LU/MIT grid is arranged in 181 latitude rings at 1 deg spacing
c     from 90N to 90S and 361 longitude values within each ring from
c     0 to 360 E.  There are 144 I*4 coefficients in each record, 12 for
c     each month of the year.

c     This routine assumes that atl.grid has been opened and the header
c     records read in rd_atl_grid.f.  The ANU file is version 1.0, the
c     LU/MIT file version 2.0.

      implicit none      
                      
      include '../includes/dimpar.h'
      include '../includes/grdtab.h'   
      include '../includes/model.h'

c Input calling arguments
c   slat     r*4  latitude of site (decimal deg)
c   slon     r*4  longitude of site (decimal deg)
c   sdoy     i*4  day-of-year
      integer*4 sdoy
      real*4 slat,slon
 
c Input from grdtab.h  
c   atlgver   r*4      version number of grid file
c   atlgswap  logical  true if the grid file needs to be byte-swapped
c   atlgmod   c*8      model name 
c   atlglat   i*4      number of latitude values on grid
c   atlglon   i*4      number of longitude values on grid 
c   natlg     i*4      number of tidal components available by interpolating
                 
c  Output in common /ufcom/ in model.h
c  atides(maxatl,6) r*4  tidal components for site

c   1st subscript (col) gives the tidal constituent (2, for S1 S2)
c   2nd subscript (row) gives the component: cos U  sin U  cos N sin N  cos E sin E (mm)
c
c Local
                                          
      integer*4 maxgval
      parameter(maxgval=144)
        
      integer*2 ival(maxgval)
      integer*4 irecbox(2,2),nhead,ngval,i,j
      real*4   rval(maxgval,4),tval(maxgval)
     .       , alat1,alon1,dlat,dlon,dx,dx1,dy,dy1

c  convert longitude range from -180/180 to 0/360
      if(slon.lt.0)slon = 360.0 + slon

c     latitude and longitude intervals (as per the file header )
c PT040106: divide by atlglat-1 because the files include 90S and 90N
      dlat = 180./float(atlglat-1)
      dlon = 360./float(atlglon)
c     lon, lat, of lower left of box
      alon1 = float(int(slon/dlon))*dlon
      if( slat.ge.0.) then
         alat1=float(int(slat/dlat))*dlat
      else
         alat1 = (float(int((90.0+slat)/dlat)))*dlat - 90.0
      endif
c  special case when slat is exactly zero
      if(slat.eq.0.0)alat1 = -dlat
c  special case when slat is exactly -90
      if(slat.eq.-90.0)alat1 = slat
      
c Get record numbers to be read from the grid file
c        (for linear interpolation, need 4)
            
      if( atlgver.eq.1.0 ) then
         nhead = 3                       
         ngval = 12
      elseif( atlgver.ge.2.0.and.atlgver.lt.3.0 ) then
         nhead = 1 
         ngval = 144                                   
      else
        call report_stat('FATAL','GRDTAB','get_atl_grid',' '
     .           ,'Unknown version for ATL grid file',' ')
      endif
c     normal case (longitudes within one latitude ring)
c       lower left (lon1,lat1)
      irecbox(1,1)= int((90.0-alat1)/dlat)*atlglon + int(alon1/dlon)
     .    + nhead
c       lower right (lon2,lat1)
      irecbox(2,1) = irecbox(1,1) + 1
c       if lower right longitude becomes 360 deg, need to use the value
c       for zero; therefore, remove the records of a whole row of longitude
      if(alon1.eq.(360.0-dlon))then
        irecbox(2,1) = irecbox(2,1) - atlglon
      endif
c       upper left (lon1,lat2)
c PT050511: with three header lines, add "3" to the computation 
c        irecbox(1,2) = (int((90.0-alat1-dlat)/dlat))*atlglon
      irecbox(1,2) = (int((90.0-alat1-dlat)/dlat))*atlglon
     .                               +int(alon1/dlon) + nhead
c       upper right (lon2,lat2)
      irecbox(2,2) = irecbox(1,2) + 1
      if(alon1.eq.(360.0-dlon))irecbox(2,2) = irecbox(2,2) - atlglon
          
cd       print*,'box corners ',irecbox(1,1),irecbox(2,1)
cd     .                       ,irecbox(1,2),irecbox(2,2)

c Get the coefficients for the four corners.  ANU is real*4 with only 12
c values per record.  LU/MIT is integer*4 with 144 values per record.

cd        print *,'luatlg ngval atlgver',luatlg,ngval,atlgver
cd        read(luatlg,rec=2,err=200) ival(1)
cd        print *,'rec 2 1st value ',ival(1) 
cd        read(luatlg,rec=2,err=200) (ival(i),i=1,ngval)
cd        print *,'rec 2 ival ',(ival(i),i=1,ngval)
            
      if( atlgver.eq.1.0 ) then
            
        read(luatlg,rec=irecbox(1,1),err=200) (rval(i,1),i=1,ngval)
        if( atlgswap ) call swap_bytes(4,rval(1,1),ngval) 
        read(luatlg,rec=irecbox(2,1),err=200) (rval(i,2),i=1,ngval)
        if( atlgswap ) call swap_bytes(4,rval(1,2),ngval)
        read(luatlg,rec=irecbox(2,2),err=200) (rval(i,3),i=1,ngval)    
        if( atlgswap ) call swap_bytes(4,rval(1,3),ngval)
        read(luatlg,rec=irecbox(1,2),err=200) (rval(i,4),i=1,ngval) 
        if( atlgswap ) call swap_bytes(4,rval(1,4),ngval)
 
                                                     
      elseif( atlgver.eq.2.0) then
                                       
cd        print *,'again ngval ',ngval
        read(luatlg,rec=irecbox(1,1),err=200) (ival(i),i=1,ngval)
        if( atlgswap ) call swap_bytes(2,ival,ngval) 
        do i=1,ngval
          rval(i,1) = float(ival(i))*1.e-4
        enddo
        read(luatlg,rec=irecbox(2,1),err=200) (ival(i),i=1,ngval)
        if( atlgswap ) call swap_bytes(2,ival,ngval) 
        do i=1,ngval
          rval(i,2) = float(ival(i))*1.e-4
        enddo
        read(luatlg,rec=irecbox(2,2),err=200) (ival(i),i=1,ngval)      
        if( atlgswap ) call swap_bytes(2,ival,ngval) 
        do i=1,ngval
          rval(i,3) = float(ival(i))*1.e-4
        enddo
        read(luatlg,rec=irecbox(1,2),err=200) (ival(i),i=1,ngval) 
        if( atlgswap ) call swap_bytes(2,ival,ngval) 
        do i=1,ngval
          rval(i,4) = float(ival(i))*1.e-4
        enddo
      endif 
cd      print *,'ngval, values at corners ',ngval 
cd      do j=1,4
cd        write(*,'(12f9.4)') (rval(i,j),i=1,ngval)
cd      enddo 

                                                                  
c Perform a bilinear spatial interpolation for all coefficients
      dy  = abs(slat-alat1)/dlat
      dy1 = (dlat - dy*dlat)/dlat
      dx  = (slon - alon1)/dlon
      dx1 = (dlon - dx*dlon)/dlon
      do i=1,ngval
        tval(i)= dx1*dy1*rval(i,1) + dx*dy1*rval(i,2)
     .          + dx1*dy*rval(i,4) + dx*dy *rval(i,3)
     .         
      enddo

             
c Transfer the interpolated values to the /ufcom/ array

      if( atlgver.eq.1.0 ) then
                                                             
c       Coefficient order: cS1n sS1n cS1e sS1e cS1u sS1u   cS2n sS2n cS2e sS2e cS2u sS2u
c       For ANU files, need to swap the order from N E U to U N E 
c         S1 cos/sin  U N E
        atides(1,1) = tval(5)
        atides(1,2) = tval(6)
        atides(1,3) = tval(1)
        atides(1,4) = tval(2)
        atides(1,5) = tval(3)
        atides(1,6) = tval(4)
c         S2 cos/sin U N E
        atides(2,1) = tval(11)
        atides(2,2) = tval(12)
        atides(2,3) = tval(7)
        atides(2,4) = tval(8)
        atides(2,5) = tval(9)
        atides(2,6) = tval(10)  
                
      elseif( atlgver.eq.2.0 ) then

cd        print *,'tval(1-12) ',(tval(i),i=1,12)
c       LU coefficient order is cS1u sS1u cS2u sS2u   cSs1n sS1n cS2n sS2n   cS1e sS1e cS2e sS2
c           x 12 months  -- assigned in interp_months
        call interp_months( ngval,tval,sdoy,atides )
      
      endif
      
      return

200   continue
c  we will have branched here only if there was an error reading the grid file
           call report_stat('FATAL','GRDTAB','get_atl_grid',' '
     .                    ,'Error reading grid file',' ')

      end
