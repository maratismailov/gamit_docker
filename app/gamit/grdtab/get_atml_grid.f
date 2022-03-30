      Subroutine get_atml_grid( syear,sdoy,sdur,slat,slon )

c     Interpolate a binary grid file to get values for non-tidal atmospheric
c     loading for a particular station.  Adapted from P. Tregoning utils/interp_atm.f.
c     Written by R. King  9 August 2006       
        

c     The grid is arranged in 360 latitude rings at unknown spacing 
c     from 90N to 90S and 720 longitude values within each ring from 
c     0E to 360-dlatE.

      implicit none   
     
      include '../includes/dimpar.h'
      include '../includes/grdtab.h'
      include '../includes/model.h'
         
c Input calling arugments
c   syear  year of session
c   sdoy   day-of-year of session   
c   sdur   duration of session (days)
c   slat    latitude reguested
c   slon    longitude requested  
      integer*4 syear,sdoy
      real*4 sdur,slat,slon
                            
c Input from grdtab.h                                        
c   luprnt     i*4 print file  
c   luatmlg    i*4 ATML grid file
c   atmlver    r*4      version number of grid file
c   atmlgswap  logical  true if the grid file needs to be byte-atmlgswap
c   atmlmod    c*8      model name
c   atmlgrecl  i*4      record length of grid file (binary, direct-access)
c   atmlglat   i*4      number of latitude values on grid
c   atmlglon   i*4      number of longitude values on grid 
c   natmlg     i*4      number of components available from file (always 3)
c   atmlg_int  R*4      epoch interval of grid file (in decimal days)

c Output in common /ufcom/ in model.h, to be written on the u-file 
c   atml_time(maxatml)  times in decimal day-of-year of output values
c   atml_val(maxatml,3) N E U (changed to U N E for the u-file) in mm for site 
c   natml number of components used  (always 3)
c   ntatml number of time values in atml_time amd atml_val arrays (usually 4/day x 3 days = 12)

  
c     Record #s and values read from the grid file (linear interpolaton)
      integer*4 irecbox(2,2)
      real*4   n(4),e(4),u(4) 

c  flag to tell grdtab.f whether the read went beyond the end of the year
      character*1 eoy

c      Local variables
      integer*4 nprior,ncurrent,nfinish,i,j,ndays
      real*4 alat1,alon1,dlat,dlon,displ(3),time,recsperday     
      integer*2 yr_start2,yr_end2
      integer*4 yr_start,yr_end
      real*4 doy_start,doy_end,start_decyr,end_decyr
      character*256 message

c  Check whether date within range
      read(luatmlg,rec=3) yr_start2,doy_start,yr_end2,doy_end
      if(atmlgswap)then
          call swap_bytes(2,yr_start2,1)
          call swap_bytes(4,doy_start,1)
          call swap_bytes(2,yr_end2,1)
          call swap_bytes(4,doy_end,1)
      endif         
      if(yr_start2.lt.50)then
        yr_start = yr_start2 + 2000 
      else
        yr_start = yr_start2 + 1900 
      endif
      if(yr_end2.lt.50)then
        yr_end = yr_end2 + 2000 
      else
        yr_end = yr_end2 + 1900 
      endif                                     
      if( mod(yr_start,4).eq.0 )then
        ndays = 366
      else
        ndays = 365
      endif
      start_decyr = yr_start + doy_start / ndays
      end_decyr = yr_end + doy_end /ndays

      if(syear + (sdoy + sdur) / ndays .lt.start_decyr .or.
     .   syear + (sdoy + sdur) / ndays .gt.end_decyr )then
         write(message,'(a,a,i5,a1,i3.3,f9.5,2(a,i5,f10.5))')
     .     ' Requested time span exceeds limits of input grid file.'
     .     ,'Requested start and duration (days): ',syear,':',sdoy,sdur
     .     ,' Grid start:',yr_start,doy_start
     .     ,' Grid end: ',yr_end,doy_end 
         call report_stat('FATAL','GRDTAB','get_atml_grid'
     .    ,'atml.grid',message,0)
       endif


c Compute the number of values to read. Add one to cover the 2400 UT value
      ntatml = int( sdur/atmlg_int + 0.5) + 1
      if(ntatml.gt.maxatml)then
        write(message,'(a,i5,a,i5)')
     .   'Maximum number of atml values (',maxatml,') exceeded:',ntatml
         call report_stat('FATAL','GRDTAB','get_atml_grid'
     .    ,' ',message,0)
      endif
  
c  convert longitude range from -180/180 to 0/360
      if(slon.lt.0)slon = 360.0 + slon

c     latitude and longitude intervals (as per the file header )
c PT040106: divide by atmlglat-1 because the files include 90S and 90N
      dlat = 180./float(atmlglat-1)
      dlon = 360./float(atmlglon)
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

      
c      Get record numbers to be read from the grid file
c       (for linear interpolation, need 4)
c       (these calculations count from beginning of grid, assuming
c        that there is one header record stating the number of 
c        records/line and the grid dimensions)
                                                                        
c       normal case (longitudes within one latitude ring)
c       lower left (lon1,lat1)
c PT040106: NOTE that the ATM grids go from 90N to 90S 
c  lower left record       
c PT050110: with an additional header line, add "2" to the computation rather than "1"
c PT060825: now there are 3 header lines so add "3"
        irecbox(1,1) = int((90.0-alat1)/dlat)*atmlglon+int(alon1/dlon)+3
c       lower right (lon2,lat1)
        irecbox(2,1) = irecbox(1,1) + 1
c  if lower right longitude becomes 360 deg, need to use the value
c  for zero; therefore, remove the records of a whole row of longitude
        if(alon1.eq.(360.0-dlon))then
          irecbox(2,1) = irecbox(2,1) - atmlglon
        endif
                                                                        
c       upper left (lon1,lat2)
c PT050110: with an additional header line, add "2" to the computation rather than "1"
c PT060825: now there are 3 header lines so add "3"
        irecbox(1,2) = (int((90.0-alat1-dlat)/dlat))*atmlglon
     .                               +int(alon1/dlon)+3
c       upper right (lon2,lat2)
        irecbox(2,2) = irecbox(1,2) + 1
        if(alon1.eq.(360.0-dlon))irecbox(2,2) = irecbox(2,2) - atmlglon


c  Bon, maitenant on a tous qu'il faut pour tout faire. Ca ne reste
c  que voyager dedans le grand fichier, sauter jusqu'au le premier
c  jour, calculer heuresement jusqu'au le dernier jour, arreter la
c  bas et prendre un apero. Sante!
c
c  there are 73 x 144 x 4 records per day (=42048 records/day)
c  PT040108: the day counter here must start at zero but start
c            has the value of 1 for doy 001. 
c  PT050104: now that there is a header line in the file, we must
c            add one record to the counter anyway .... 
c  PT060824: make the computation of records per day generic
      recsperday =  atmlglon*atmlglat/atmlg_int   
c  Get the last record of the previous epoch (always an even day)
      nprior = (sdoy-1)*int(recsperday)
c  Now add 1 to set the pointer to the first record of the next epoch            
      ncurrent = nprior + 1
c  we don't know how long the observation session is. Calculate it
      nfinish = (sdoy+sdur-1)*recsperday
                         
      natml = 0  
      time = sdoy 
        do i=1,ntatml
          natml = natml + 1
c  correct the record indicator for prior days and increment it for
c  days/times already read
c  and  read the three displacement values at these four corner points.
          read(luatmlg,rec=irecbox(1,1)+ncurrent,err=200)n(1),e(1),u(1)
          read(luatmlg,rec=irecbox(2,1)+ncurrent,err=200)n(2),e(2),u(2)
          read(luatmlg,rec=irecbox(2,2)+ncurrent,err=200)n(3),e(3),u(3)
          read(luatmlg,rec=irecbox(1,2)+ncurrent,err=200)n(4),e(4),u(4)

*     If file is byte-atmlgswap, then switch the values just read
          if( atmlgswap ) then
            call swap_bytes(4,n,4)
            call swap_bytes(4,e,4)
            call swap_bytes(4,u,4)
          end if

c  now perform a bi-linear interpolation
      call bilin4(slon,slat,alon1,alat1,dlon,dlat,n(1),n(2),n(3),n(4)
     .          ,displ(1))
      call bilin4(slon,slat,alon1,alat1,dlon,dlat,e(1),e(2),e(3),e(4)
     .          ,displ(2))
      call bilin4(slon,slat,alon1,alat1,dlon,dlat,u(1),u(2),u(3),u(4)
     .          ,displ(3))
      
          atml_time(natml) = time 
          time = time + atmlg_int

c PT060825: these are in order N, E, U - compatible with the order in
c           the binary grid files and as extracted by atmtoasc. 
          do j=1,3
            atml_val(natml,j) = displ(j)
          enddo
c          write(*,'(4f10.2)')atml_time(natml),displ

c increment ncurrent by the number of records in a single grid (=atmlglon*atmlglat)
          ncurrent = ncurrent + atmlglon*atmlglat
        enddo

c print what was used to grdtab.out

      write(luprnt,'(a,a8)') '  ATML  grid ',atmlgmod              

c'est tout!
      return

200   continue
c  we will have branched here only if there was an error reading the grid file
c  PT050418: this can happen if we try to read beyond the end of the file
c            (ie day 366 for non-leap years and day 367 for leap years). GRDTAB
c            doesn't know what year it is (!!) so we can't discriminate.
c            The next version of the grid files should contain a line stating
c            the time interval contained in the file!
c
c            For now, end fatally only if the day of year is less than 366
           if((float(ncurrent)/42048.0).ge.366.0)then
             eoy = 'Y'
           else
             call report_stat('FATAL','GRDTAB','get_atml_grid',' '
     .                    ,'Error reading grid file',0)
           endif


      end


