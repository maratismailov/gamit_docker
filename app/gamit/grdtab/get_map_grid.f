      Subroutine get_map_grid( syear,sdoy,sdur,slat,slon,ell_ht )

c     Interpolate a binary grid file to get values for mapping function 
c     coefficients by Joannes Broehm of IGG/Vienna.

c     Adapted from P. Tregoning utils/interp_vmf1grd.f
c     Written by R. King  10 August 2006       

c   MODS
c   
c   PT061213: reduce each node ZHD value to the station ellipsoidal
c             height, then perform the bilinear interpolation. This
c             now invokes new routine correct_zhd2        

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
c   ell_ht  ellipsoidal height of the site requested  
      integer*4 syear,sdoy
      real*4 sdur,slat,slon,ell_ht

c Input from grdtab.h                                        
c   luprnt    i*4 unit number of print file  
c   lumapg    i*4  unit number of mapping function grid file
c   mapver    r*4      version number of grid file
c   mapgswap  logical  true if the grid file needs to be byte-swapped
c   mapmod    c*8      model name
c   mapgrecl  i*4      record length of grid file (binary, direct-access)
c   mapglat   i*4      number of latitude values on grid
c   mapglon   i*4      number of longitude values on grid 
c   nmapg     i*4      number of components available from file (3: Ah Aw ZHD)
c   mapg_int  r*4      epoch interval (decimal days) of grid file

c Output in common /ufcom/ in model.h, to be written on the u-file 
c   map_time(maxmap)  times in decimal day-of-year of output values
c   map_val(maxmap,3) Ah Aw ZHD values of mapping function (and met) values for site 
c   nmap number of values in map_time amd map_val arrays (usually 4/day x 3 days = 12)

                                                             
c    Record #s and values read from the grid file (linear interpolaton)
      integer*4 irecbox(2,2)
      real*4   val1(4),val2(4),val3(4) 

c      Local variables
      integer*4 nprior,ncurrent,nfinish,ioerr,i,j,nrec
     .         ,nrecsperday,iorog(4)
      real*4 alat1,alon1,dlat,dlon,map(3),orog(4)
      real*4 deg2rad,tmp,time
      character*256 message

c PT060531: variables specific for reading the VMF1 grid
c RWK130116; rename the coefficients with 'i' to avoid conflict real*8 variables
c            in model.h
      integer*4 iAh(4),iAw(4)
      real*4 ZHD(4)

c PT061213: variables for correct_zhd2
      real*8 dmjd
      real*4 dzhd,pi
      parameter(pi = 3.141592653)


c  flag to tell grdtab.f whether the read went beyond the end of the year


      deg2rad = 3.14159265/180.0
     

c Compute the number of values to read. Add one to cover the 2400 UT value
      ntmap = int( sdur/mapg_int + 0.5) + 1
      if(ntmap.gt.maxmap)then
        write(message,'(a,i5,a,i5)')
     .   'Maximum number of mapf values (',maxmap,') exceeded:',ntmap
         call report_stat('FATAL','GAMIT','GRDTAB/get_map_grid'
     .    ,' ',message,0)
      endif

c  convert longitude range from -180/180 to 0/260
      if (slon.lt.0) slon = 360.0 + slon

c     latitude and longitude intervals (as per the file header )
c PT040106: divide by mapglat-1 because the files include 90S and 90N
      dlat = 180./float(mapglat-1)
      dlon = 360./float(mapglon)
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
c PT040106: NOTE that the VMF1 grids go from 90N to 90S 
c  lower left record
c PT060531: the VMF1 grid has three header record, so add 3
        irecbox(1,1) = int((90.0-alat1)/dlat)*mapglon+int(alon1/dlon)+3
c       lower right (lon2,lat1)
        irecbox(2,1) = irecbox(1,1) + 1
c  if lower right longitude becomes 360 deg, need to use the value
c  for zero; therefore, remove the records of a whole row of longitude
      if(alon1.eq.(360.0-dlon))then
        irecbox(2,1) = irecbox(2,1) - mapglon
      endif
                                                                        
c     upper left (lon1,lat2)
      irecbox(1,2) = (int((90.0-alat1-dlat)/dlat))*mapglon
     .                               +int(alon1/dlon)+3
c     upper right (lon2,lat2)
      irecbox(2,2) = irecbox(1,2) + 1
      if(alon1.eq.(360.0-dlon))irecbox(2,2) = irecbox(2,2) - mapglon

c        print*,'alat1 and alon1',alat1,alon1,slat,slon
c        print*,'row number of LL = ', (int((90.0-alat1)/dlat))*mapglon
c        print*,'colum number of LL = ',int(alon1/dlon) + 1

c  read the orography to get the heights of the four box corners. The grid of this
c  information starts at the first record of what would be day 368.
      nrec = int(368 * mapglon*mapglat/mapg_int)
 
c PT060830: need to subtract "3" from irecbox because there are no "headers" to
c           consider when reading the orography records.
      read(lumapg,rec=nrec+irecbox(1,1)-3,iostat=ioerr)
     .      tmp,tmp,iorog(1)
      read(lumapg,rec=nrec+irecbox(2,1)-3,iostat=ioerr)
     .      tmp,tmp,iorog(2)
      read(lumapg,rec=nrec+irecbox(2,2)-3,iostat=ioerr)
     .      tmp,tmp,iorog(3)
      read(lumapg,rec=nrec+irecbox(1,2)-3,iostat=ioerr)
     .      tmp,tmp,iorog(4)
c      print*,'Read Orography values are:',iorog
c      print*,'box corners for orography:',nrec,irecbox(1,1)-3
c     .   ,nrec+irecbox(2,1)-3,nrec+irecbox(2,2)-3,irecbox(1,2)-3
       if(mapgswap)then
         call swap_bytes(4,iorog,4)
       endif
       do i=1,4
         orog(i) = iorog(i)*1.0
       enddo

c PT061213: no longer required
c       call bilin4(slon,slat,alon1,alat1,dlon,dlat,orog(1),orog(2)
c     .          ,orog(3),orog(4),orog_hgt)
cc      print*,'Orography values are:',orog
cc      print*,'Interpolated height is:',orog_hgt,' m. Ell_ht =',ell_ht

c PT061213: get the MJD for use in the GPT model for pressure and temperature
      call yds_to_jd ( syear, sdoy, 0, dmjd )
      dmjd = dmjd - 2.4d6

c  Bon, maitenant on a tous qu'il faut pour tout faire. Ca ne reste
c  que voyager dedans le grand fichier, sauter jusqu'au le premier
c  jour, calculer heuresement jusqu'au le dernier jour, arreter la
c  bas et prendre un apero. Sante!
c
c  there are 91 x 144 x 4 records per day (=52416 records/day). This is
c  the same for the IMF and the VMF1 grids.
c  PT040108: the day counter here must start at zero but start
c            has the value of 1 for doy 001.
      nrecsperday = int(mapglon*mapglat/mapg_int)
      nprior = (sdoy-1)*nrecsperday  

c  we don't know how long the observation session is. Calculate it
      nfinish = (sdoy+sdur-1)*nrecsperday

c  adjust the counter for the header record
      ncurrent = nprior + 1
               
c  start the time counter for the values to be saved and written on the u-file
      nmap = 0
      time = sdoy 
      do i=1,ntmap
        nmap = nmap + 1
                                                               
c  correct the record indicator for prior days and increment it for
c  days/times already read
c  and  read the three mapping values at these four corner points.
        read(lumapg,rec=irecbox(1,1)+ncurrent,err=200,iostat=ioerr)
     .           iAh(1),iAw(1),ZHD(1)
        read(lumapg,rec=irecbox(2,1)+ncurrent,err=200,iostat=ioerr)
     .           iAh(2),iAw(2),ZHD(2)
        read(lumapg,rec=irecbox(2,2)+ncurrent,err=200,iostat=ioerr)
     .           iAh(3),iAw(3),ZHD(3)
        read(lumapg,rec=irecbox(1,2)+ncurrent,err=200,iostat=ioerr)
     .           iAh(4),iAw(4),ZHD(4)
                                                                        
*     If file is byte-swapped, then switch the values just read
        if( mapgswap ) then
          call swap_bytes(4,iAh,4)
          call swap_bytes(4,iAw,4)
          call swap_bytes(4,ZHD,4)
        end if

c PT060531: convert the integer values of coefficients into real*4
          do j=1,4
            val1(j) = iAh(j)  / 1.d8
            val2(j) = iAw(j)  / 1.d8
            val3(j) = ZHD(j) * 1.d3
c PT061213: now we reduce each of the ZHD values to the height
c           of the geodetic site
c            print*,'reduce_zhd input',slat,slon,dmjd,orog(j),ell_ht
            call reduce_zhd(slat*pi/180.0,slon*pi/180.0,dmjd
     .         ,orog(j),ell_ht,dzhd)
c            print*,'node zhd and correction',val3(j),dzhd
            val3(j) = val3(j) - dzhd
          enddo

c      print*,'z200' ,n(1),n(2),n(3),n(4)
c      print*,'swpf3',e(1),e(2),e(3),e(4)
c      print*,'temp' ,u(1),u(2),u(3),u(4)
c  now perform a bi-linear interpolation
      call bilin4(slon,slat,alon1,alat1,dlon,dlat,val1(1),val1(2)
     .          ,val1(3),val1(4),map(1))
      call bilin4(slon,slat,alon1,alat1,dlon,dlat,val2(1),val2(2)
     .          ,val2(3),val2(4),map(2))
      call bilin4(slon,slat,alon1,alat1,dlon,dlat,val3(1),val3(2)
     .          ,val3(3),val3(4),map(3))

        map_time(nmap) = time 
        time = time + mapg_int

c      print*,'interpolated values: ',map
c      stop
      
      do j=1,nmapg
        map_val(nmap,j) = map(j)
      enddo


c increment ncurrent by the number of records in a single 6 hr grid 
c  (=10512 for 73 x 144 grid; 13104 for 91 x 144 grid)
        ncurrent = ncurrent + mapglon*mapglat
      enddo
            

c write what was used to grdtab.out

      write(luprnt,'(a,a8)') '  MAP  grid ',mapgmod              


c'est tout!
      return

200   continue
c  we will have branched here only if there was an error reading the grid file
           call report_stat('FATAL','GRDTAB','get_map_grid','map.grid'
     .               ,'Error reading mapping function grid file',ioerr)

      end

















