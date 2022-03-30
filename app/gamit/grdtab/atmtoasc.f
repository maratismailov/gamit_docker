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
     .      ,t
      integer*4 ngval,nglon,nglat,recl,ioerr,ts,tend
     .         ,irecbox(2,2),nrec,ngrid
      character atmfile*128,arg1*30,arg2*3 
                       

c PT050104: variables for the header lines of the atmdisp files
      integer*2  version,nhead,i2val,i2lon,i2lat
     .         , yr_start,yr_end
      real*4 interval,doy_start,doy_end,fversion                 

      character*5 press_src
      character*3 ref_frame          
      character*2 spare

      logical swapped

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
      write(*,'(2a)') '  Binary atml file: ',trim(atmfile)
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
      ngrid = nglon*nglat
     
c  Record 3 (version 2.0ff only): start/stop times 

c     Version 2.0 : yr_start  doy_start  yr_end  doy_end
c                      i2       r4        i2      r4

      if( version.ge.20 ) then
        read(10,rec=3,iostat=ioerr) yr_start,doy_start,yr_end,doy_end
        if( ioerr.ne.0 ) then
           write(*,'(a,i4)') 
     .             'Error reading record 3 of grid file, iostat=',ioerr
           stop   
        endif
        if( swapped ) then
cd           print *,'swapped '
           call swap_bytes(2,yr_start,1)
           call swap_bytes(4,doy_start,1)
           call swap_bytes(2,yr_end,1)
           call swap_bytes(4,doy_end,1) 
        endif                                
        write(*,'(/,a,1x,i2,1x,f6.2,2x,i2,1x,f6.2)') 
     .        'Record 3: yr/doy start yr/doy end'
     .       , yr_start,doy_start,yr_end,doy_end
      else  
        write(*,'(a,1x,i2,1x,f6.2,2x,i2,1x,f6.2)') 
     .        'Version 1.0, no Record 3; assume start/stop 1.0 365.75'
        doy_start = 1.0
        doy_end = 365.
      endif   
         
cd      do nrec=4,50                       
cd       read(10,rec=nrec) n(1),e(1),u(1) 
cd       print *,'record ',nrec,n(1),e(1),u(1)  
cd      enddo
cd      stop


c  Now find the four grid points that are closest  to the input coordinates.
 
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
c       special case when slat is exactly zero
      if(slat.eq.0.0)alat1 = -dlat
c       special case when slat is exactly -90
      if(slat.eq.-90.0)alat1 = slat
c        special case when slat is exactly 90
      if(slat.eq.90.0)alat1 = slat - dlat

c debug
      print*,'LL grid lat long:',alat1,alon1
  
c  Calculate the record-number offsets from the 1st point of the epoch
c     (i.e., upper left = 0) 
c  lower left record
        irecbox(1,1) = int((90.0-alat1)/dlat)*nglon+int(alon1/dlon) 
        print *,'alat1 dlat alon1 dlon nhead LL '
     .         , alat1,dlat,alon1,dlon,nhead,irecbox(1,1)
c       lower right (lon2,lat1)
        irecbox(2,1) = irecbox(1,1) + 1
c  if lower right longitude becomes 360 deg, need to use the value
c  for zero; therefore, remove the records of a whole row of longitude
        if(alon1.eq.(360.0-dlon))then
          irecbox(2,1) = irecbox(2,1) - nglon
        endif
c       upper left (lon1,lat2)
        irecbox(1,2) = (int((90.0-alat1-dlat)/dlat))*nglon
     .                      +int(alon1/dlon) 
c       upper right (lon2,lat2)
        irecbox(2,2) = irecbox(1,2) + 1
        if(alon1.eq.(360.0-dlon))irecbox(2,2) = irecbox(2,2) - nglon

c  debug
        print*,'Box corner record numbers LL LR UR UL:'
     .   ,irecbox(1,1),irecbox(2,1),irecbox(2,2),irecbox(1,2)

c Now loop through the requested epochs
                                           
      t = ts
      do while(t.le.(tend+interval))    
         call grid_time('T',doy_start,interval,ngrid,nhead,t,nrec )
c        read values for the four corners at this epoch
         read(10,rec=irecbox(1,1)+nrec,iostat=ioerr) n(1),e(1),u(1)
         read(10,rec=irecbox(2,1)+nrec,iostat=ioerr) n(2),e(2),u(2)
         read(10,rec=irecbox(2,2)+nrec,iostat=ioerr) n(3),e(3),u(3)
         read(10,rec=irecbox(1,2)+nrec,iostat=ioerr) n(4),e(4),u(4)
         if( ioerr.ne.0 ) then  
           print *,'ioerr rec ',ioerr,irecbox(1,1)+nrec
            write(*,'(a,i8,a,f6.2)') 
     .         'Attempt to read beyond EOF, record = '
     .           ,nrec,'   t=',t
           stop
         else
           if( swapped ) then
             call swap_bytes(4,n,4)
             call swap_bytes(4,e,4)
             call swap_bytes(4,u,4)
           endif  
         endif
c        perform a bi-linear interpolation
         dy  = abs(slat-alat1)/dlat
         dy1 = (dlat - dy*dlat)/dlat
         dx  = (slon - alon1)/dlon
         dx1 = (dlon - dx*dlon)/dlon
         displ(1) = dx1*dy1*n(1)+dx*dy1*n(2)+dx1*dy*n(4)+dx*dy*n(3)
         displ(2) = dx1*dy1*e(1)+dx*dy1*e(2)+dx1*dy*e(4)+dx*dy*e(3)
         displ(3) = dx1*dy1*u(1)+dx*dy1*u(2)+dx1*dy*u(4)+dx*dy*u(3)
cd DEBUG
cd       write(*,'(a)') 'irec U '
cd       write(*,'(i5,f11.4)') irecbox(1,1),u(1)
cd       write(*,'(i5,f11.4)') irecbox(2,1),u(2)
cd       write(*,'(i5,f11.4)') irecbox(2,2),u(3)
cd       write(*,'(i5,f11.4)') irecbox(1,2),u(4)
c        write out the time and displacements
         write(*,'(f6.2,3f11.4)') t,displ
         t = t + interval
      enddo

c eh voila, c'est finit. Sante!
 
      close(10)
      end


