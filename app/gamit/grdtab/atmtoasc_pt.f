      Program ATMTOASC

c   Program to read a binary atmospheric grid file and extract
c   the displacements for a single geographic location for a
c   particular time interval. This may be a useful utility but
c   I am writing it primarily to see whether I can actually
c   read and interpolate Tonie van Dam's binary grid files
c   correctly
c
c   P. Tregoning  6 January 2004
c   R. King 9 October 2006: mods for versions 2.0 and 2.1

      implicit none

      real*4 slat,slon,alat1,alon1,dlat,dlon
     .      ,n(4),e(4),u(4),displ(3),dx,dx1,dy,dy1
     .      ,tcurrent
      integer*4 ngval,nglon,nglat,recl,ioerr,ts,tend
     .         ,irecbox(2,2),i,ncurrent,nfinish,nprior
      character atmfile*15,arg1*30,arg2*3 
                       

c PT050104: variables for the header lines of the atmdisp files
      integer*2  version,nhead,i2val,i2lon,i2lat
     .         , yr_start,yr_end
      real*4 interval,doy_start,doy_end,fversion                 

      character*5 press_src
      character*3 ref_frame          
      character*2 spare

      logical swapped

c RWK DEBUG
c      real*4 values(3)
c      integer*4 irec
c end DEBUG 



      swapped = .false.

c   the input atmospheric loading file contains displacements
c   (not corrections) of the Earth as a result of the pressure
c   loading. The first line of the file contains three values:
c
c   number of values per line, number of grid nodes in long, num in lat
c
c   The file then contains the displacements in N, E, U (in mm) for
c   each grid node in rows of longitude, starting at 90N, 00E then
c   running across to 357.5E before going to 87.5N, 00E. The first
c   set of grid entries correspond to decimal day 001.00, followed
c   by 001.25 etc through until the last day of the year. Thus, the
c   last row of the file will contain NEU values for grid point
c   90S, 357.5E at decimal day 365.75 (or 366.75 if it is a leap year).
c
c   Clear as mud? 

c  each line has 3 x real*4 values. Therefore, the record length is 12
      recl = 12

c   decode command line
      call getarg(1,atmfile)
      if(atmfile(1:1).eq.' ')then
         print*,'Runstring: atmtoasc atmfile lat long doy_st doy_end'
         stop 'No command line arguments found. Program stopped'
      else
         open(unit=10,file=atmfile,form='unformatted',access='direct'
     .        ,status='old',recl=12,iostat=ioerr)
         if(ioerr.ne.0) then 
           write(*,'(3a,i4)') 'Error opening ',atmfile,' iostat=',ioerr  
           stop
         endif
      endif

c  get the site coordinates
      call getarg(2,arg1)
      read(arg1,*)slat
      call getarg(3,arg1)
      read(arg1,*)slon

c  get the start and end days
      call getarg(4,arg2)
      read(arg2,'(i3.3)') ts
      call getarg(5,arg2)
      read(arg2,'(i3.3)') tend
             
      write(*,'(a)') 'ATMTOASC input arguments: '
      write(*,'(2a)') '  Binary atml file: ',atmfile
      write(*,'(a,2f9.3)') '  Lat lon : ',slat,slon
      write(*,'(a,2i4)') '  Start, end days ',ts,tend
            
c  Record 1:   
c
c  Version 1.0 or 2.0       version  spare pressure-source  reference-frame
c                             I*2     C*2     C*5              C*3
c  
c  Version 2.1 (same except:          nhead
c                                     I*2

      read(10,rec=1,iostat=ioerr) version 
      print*,'Initial read ',version
      if(version.gt.100.or.version.lt.0)then
*         If the grid size does not seem be correct, swap the bytes
*         and see if OK.   
        write(*,'(a,i4,a)') '**version =',version,' swapping**'
        call swap_bytes(2,version,1)
*         Now see if OK size.
        if(version.gt.100.or.version.lt.0) then 
           write(*,'(a,i4,a)') 'Invalid version (',version
     .            ,') of ATML grid file '
*             Still a problem so kill.   
           stop
        else
*         File seems to be swapped.  Set the swapped status and continue.
*         (We now need to swap every record read)  
          write(*,'(a)') 'Byte-swapping input grid file '
          swapped = .true.    
        end if
      endif  
      if( version.eq.10 ) then  
        read(10,rec=1,iostat=ioerr) version,spare,press_src,ref_frame  
        if( swapped ) call swap_bytes(2,version,1)
        nhead= 2    
        write(*,'(/,a,1x,i4,1x,a5,1x,a3)') 
     .    'Read Record 1 Version 1.0  version pressure-source frame '
     .            ,version,press_src,ref_frame 
        write(*,'(a)') '--assume nhead = 2'
      elseif( version.eq.20 ) then   
        read(10,rec=1,iostat=ioerr) version,spare,press_src,ref_frame    
       if( swapped ) call swap_bytes(2,version,1)
        nhead = 3      
        write(*,'(/,a,1x,i4,1x,a5,1x,a3)') 
     .    'Read Record 1 Version 2.0  version pressure-source frame '
     .            ,version,press_src,ref_frame 
        write(*,'(a)') '--assume nhead = 3'
      elseif( version.eq.21 ) then 
        read(10,rec=1,iostat=ioerr) version,nhead,press_src,ref_frame 
        if( swapped ) then          
           call swap_bytes(2,version,1)
           call swap_bytes(2,nhead,1)
        endif 
        write(*,'(a,2i3,1x,a5,1x,a3)') 
     . 'Read Record 1 Version 2.1  version nhead pressure-source frame '
     .            ,version,nhead,press_src,ref_frame  
      else
         write(*,'(a,i2)') 'Unknown version ',version
         stop
      endif
       
c  Record 2:  array limits  

c     Version 1.0 :  ngval   nglon   nglat
c                     i4       i4     i4

c     Version 2.0 :  interval  ngval   nglon   nglat
c                      r4        i2      i2      i2

      if( version.eq.10 ) then 
         read(10,rec=2,iostat=ioerr) ngval,nglon,nglat
         if(swapped)then
            call swap_bytes(4,ngval,1)
            call swap_bytes(4,nglon,1)
            call swap_bytes(4,nglat,1)
        endif
        write(*,'(a,3i5,/,a)') 
     .       'Read Record 2 Version 1.0 ngval nglon nglat '
     .            ,ngval,nglon,nglat
     .           ,'--assume interval = 0.25 d'
        interval = 0.25 

      elseif( version.ge.20 ) then  
        read(10,rec=2,iostat=ioerr) interval,i2val,i2lon,i2lat  
        if( ioerr.ne.0 ) then  
           write(*,'(a)') 'Error reading record 2 of grid file '
           stop    
        endif
        if(swapped)then
            call swap_bytes(4,interval,1)
            call swap_bytes(2,i2val,1)
            call swap_bytes(2,i2lon,1)
            call swap_bytes(2,i2lat,1)   
        endif    
        ngval = i2val
        nglon = i2lon
        nglat = i2lat
        fversion = float(version)/10.d0                 
        write(*,'(/,a,f3.1,a,1x,f4.2,1x,3i5)')
     .   'Read Version ',fversion,' Record 2 interval ngval nglon nglat'
     .                     ,interval,ngval,nglon,nglat        
      endif   

     
c  Record 3 (version 2.0ff only): start/stop times

      if( version.ge.20 ) then
        read(10,rec=3,iostat=ioerr) yr_start,doy_start,yr_end,doy_end
        if( ioerr.ne.0 ) then
           write(*,'(a,i4)') 
     .             'Error reading record 3 of grid file, iostat=',ioerr
           stop   
        endif
        if( swapped ) then
           call swap_bytes(2,yr_start,1)
           call swap_bytes(4,doy_start,1)
           call swap_bytes(2,yr_end,1)
           call swap_bytes(4,doy_end,1) 
        endif                                
        write(*,'(a,1x,i2,1x,f6.2,2x,i2,1x,f6.2)') 
     .        'Record 3: yr/doy start yr/doy end'
     .       , yr_start,doy_start,yr_end,doy_end
        else
         print*,'Not reading header line 3'
      endif
                      
c**** TEMP DEBUG
c   last 365 record 
c      irec = 73*144*4*365 + nhead
c      read(10,rec=irec,iostat=ioerr) values 
c      if( ioerr.ne.0 ) write(*,'(a,i5)') '365 iostat ',ioerr
c      if( swapped ) call swap_bytes(4,values,3)
c      write(*,'(a,i10,3f6.2)') 'Last  365 rec ',irec,values
c      irec = irec + 1
c      read(10,rec=irec,iostat=ioerr) values
c      if( ioerr.ne.0 ) write(*,'(a,i5)') 'iostat ',ioerr
c      if( swapped ) call swap_bytes(4,values,3)
c      write(*,'(a,i10,3f6.2)') 'First day 366 rec ',irec,values    
c      irec = 73*144*4*366 + nhead
c      read(10,rec=irec,iostat=ioerr) values     
c      if( ioerr.ne.0 ) write(*,'(a,i5)') 'last iostat ',ioerr
c      if( swapped ) call swap_bytes(4,values,3)
c      write(*,'(a,i10,3f6.2)') 'Last day 366 rec ',irec,values    
c      stop                
c*****  end DEBUG


c  Now find the four grid points that are closest
c  to the input coordinates.
c       normal case (longitudes within one latitude ring)
c       lower left (lon1,lat1)
c PT040106: NOTE that the ATM grids go from 90N to 90S 
c PT040106: divide by nglat-1 because the files include 90S and 90N
      dlat = 180./float(nglat-1)
      dlon = 360./float(nglon)
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
c  special case when slat is exactly 90
      if(slat.eq.90.0)alat1 = slat - dlat

c debug
c      print*,'LL grid lat long:',alat1,alon1

c  lower left record
        irecbox(1,1) = int((90.0-alat1)/dlat)*nglon+int(alon1/dlon) 
     .                  + nhead          
c        print *,'nglon nhead ',nglon,nhead
c       lower right (lon2,lat1)
        irecbox(2,1) = irecbox(1,1) + 1
c  if lower right longitude becomes 360 deg, need to use the value
c  for zero; therefore, remove the records of a whole row of longitude
        if(alon1.eq.(360.0-dlon))then
          irecbox(2,1) = irecbox(2,1) - nglon
        endif

c       upper left (lon1,lat2)
        irecbox(1,2) = (int((90.0-alat1-dlat)/dlat))*nglon
     .                      +int(alon1/dlon) + nhead
c       upper right (lon2,lat2)
        irecbox(2,2) = irecbox(1,2) + 1
        if(alon1.eq.(360.0-dlon))irecbox(2,2) = irecbox(2,2) - nglon

c  debug
c        print*,'Box corner record numbers LL LR UR UL:'
c     .   ,irecbox(1,1),irecbox(2,1),irecbox(2,2),irecbox(1,2)

c  Bon, maitenant on a tous qu'il faut pour tout faire. Ca ne reste
c  que voyager dedans le grand fichier, sauter jusqu'au le premier
c  jour, calculer heuresement jusqu'au le dernier jour, arreter la
c  bas et prendre un apero. Sante!
c
c  there are 73 x 144 x 4 records per day (=42048 records/day)
c  PT040108: for the first day, ts is 1 but we don't want to
c            skip through the file because the first record is
c            for day 1 not day 0. So subtract 1 from ts.
      nprior = (ts-1)*nglon*nglat/interval
      nfinish = (tend+1)*nglon*nglat/interval
      ncurrent = nprior + 1

c debug 
c      print*,'ts,nprior,nfinish,ncurrent',ts,nprior,nfinish,ncurrent

      do while (ncurrent.lt.nfinish)
                
c  need to read four full grids per day
        do i=1,4
c  correct the record indicator for prior days and increment it for
c  days/times already read
c  and  read the three displacement values at these four corner points.
          read(10,rec=irecbox(1,1)+ncurrent,iostat=ioerr)n(1),e(1),u(1)
          read(10,rec=irecbox(2,1)+ncurrent,iostat=ioerr)n(2),e(2),u(2)
          read(10,rec=irecbox(2,2)+ncurrent,iostat=ioerr)n(3),e(3),u(3)
          read(10,rec=irecbox(1,2)+ncurrent,iostat=ioerr)n(4),e(4),u(4)
          if( ioerr.ne.0 ) then  
            tcurrent = 1.0+float(ncurrent)/42048.0
            write(*,'(a,i8,a,f6.2)') 
     .         'Attempt to read beyond EOF, record = '
     .           ,ncurrent,'   t=',tcurrent
            stop
          else
            if( swapped ) then
              call swap_bytes(4,n,4)
              call swap_bytes(4,e,4)
              call swap_bytes(4,u,4)
            endif  
          endif
c          print*,'LL values:',n(1),e(1),u(1)


c  now perform a bi-linear interpolation
c PT040106: I found a simple, one-line formula that seems to do this!
c                  VAL = DX1*DY1*ARRAY(ILO,JLO) + DX*DY1*ARRAY(IHI,JLO) +
c     &               DX1*DY*ARRAY(ILO,JHI) + DX*DY*ARRAY(IHI,JHI)

          dy  = abs(slat-alat1)/dlat
          dy1 = (dlat - dy*dlat)/dlat
          dx  = (slon - alon1)/dlon
          dx1 = (dlon - dx*dlon)/dlon
c then
          displ(1) = dx1*dy1*n(1)+dx*dy1*n(2)+dx1*dy*n(4)+dx*dy*n(3)
          displ(2) = dx1*dy1*e(1)+dx*dy1*e(2)+dx1*dy*e(4)+dx*dy*e(3)
          displ(3) = dx1*dy1*u(1)+dx*dy1*u(2)+dx1*dy*u(4)+dx*dy*u(3)

c the day counter here assumes that 0 h on day 1 is day 0. We need to
c add 1 to the doy calc so that we can write out that 0h on doy 001 is 1.00
          write(*,'(f6.2,3f11.4)') 1.0+float(ncurrent)/42048.0,displ

c increment ncurrent by the number of records in a single 6 hr grid (=10512)
          ncurrent = ncurrent + nglon*nglat
        enddo
      enddo

c eh voila, c'est finit. Sante!
 
      close(10)
      end


