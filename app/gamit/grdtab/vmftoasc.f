      Program vmftoasc
c                  
c   List values in a VMF1 grid file by lat lon and date
c   Paul Tregoning September 2005.
c   Last modified by R. King August 2007
c         

c   The Version 1.0  VMF1 grid is a direct-access file containing
c   the hydrostatic (Ah)and wet (Aw) mapping function coefficients 
c   and the hydrostatic zenith delay (ZHD, meters) for a global grid
c   at 6-hr intervals. c   It also contains, at the end the orography 
c   used by the Vienna group to calculate the values.  Since the orography 
c   is at the end (position of day 368), a file that is incomplete for
c   the year will still be full size (231626292 bytes).  Future versions
c   of the file may have more values (adding the lapse rate, e.g.) and
c   place the orography at the beginning of the file.  The grid is
c 
c   The grid has rows of rows of longitude, starting at 90N, 00E then
c   running across to 357.5E before going to 87.5N, 00E (73 rows x
c   144 columns). The first set of grid entries correspond to decimal 
c   day 001.00, followed by 001.25 etc through until the last day of the 
c   year, usually padded through day 367 for interpolation. Thus, the last 
c   row of the grid will contain values for 90S, 357.5E at decimal day 367.75.

c   Version 1.0 has a record length of 12 bytes, filled by integers (later
c   (converted to real*4) for the met values and orography but a combination
c   of integer, character, and real for the header (see below).

c  PT060829: the header lines are as follows:
c Header line 1: 10VMF1JB (version 1.0, VMF1, source: JB = Johannes Boehm: C*2 C*4 C*2)
c Header line 2: 06214491 (06 hour epochs, 2 values/record, nglon=144, nglat=91: C*2 C*1 C*3 C*2)
c Header line 3: start_time,stop_time (R*4 R*4)

      implicit none

c  Variables from the VMF file
      integer*4 Ah(4),Aw(4),orog(4)
      integer*4 vmf(2)
      real*4 ZHD(4),tmp
      real*8    zhdval,height

c  Variables from the header recordss
      character  mapf*4,source*2
      integer*4 nglon,nglat,nrec
      integer*2 version,interval,ngval   
      real*4 start_time,stop_time

c  Other variables
      real*4 slat,slon,alat1,alon1,dlat,dlon
     .      ,dx,dx1,dy,dy1,intervald,start_doy,stop_doy,t
      integer*4 recl,ngrid,ierr,ts,tend,nhead
     .         ,irecbox(2,2),start_yr,stop_yr,ioerr  
c**   debug
cd      integer*4 i,irec
      character mapgrid*32,arg1*30,arg2*3
      logical swapped

c  Function

      integer*4 nydays

c  Initialize record size and # header records for Version 1.0
      recl = 12  
      nhead = 3

c  Decode the command line

      call getarg(1,mapgrid)
      if(mapgrid(1:1).eq.' ')then
         write(*,'(a)') 
     .       'Runstring: vmftoasc vmffile lat long doy_st doy_end'
         stop 'No command line arguments found. Program stopped'
      else
         open(unit=10,file=mapgrid,form='unformatted',access='direct'
     .        ,status='old',recl=recl,iostat=ierr)
         if(ierr.ne.0)then
            write(*,'(2a)') 'Error opening file ',mapgrid
            stop 'Stopped with errors'
         endif
       endif

c  Get the site coordinates
      call getarg(2,arg1)
      read(arg1,*)slat
      call getarg(3,arg1)
      read(arg1,*)slon

c  Get the start and end days
      call getarg(4,arg2)
      read(arg2,'(i3.3)') ts
      call getarg(5,arg2)
      read(arg2,'(i3.3)') tend

cd      print*,'Input parameters:',mapgrid,slat,slon,ts,tend
                         
c  Read the headers     

c    Use detection of a valid version number to decide whether integer
c    and real numbers need to be byte swapped to convert from BIG ENDIAN
c    to little endian representation.

c    Record 1:   
c
c    Version 1.0     version  map id   source   spare
c                      I*2     C*4     C*2       C*4
c
c             e.g  10VMF1JB (version 1.0, VMF1, source: JB = Johannes Boehm)

      read(10,rec=1,iostat=ioerr) version,mapf,source 
      if(version.lt.10.or.version.gt.100)then
cd        print*,'calling swap_bytes'
        call swap_bytes(2,version,1)
        if(version.ge.10.or.version.lt.100) then
           call report_stat('STATUS','GRDTAB','vmftosac'
     .            ,mapgrid, 'Byte-swapping input grid file',0)
           swapped = .true.
        else
          call report_stat('FATAL','GRDTAB','vmftoasc',
     .        mapgrid,'Invalid values read from binary header',0)
        endif
       else
         swapped = .false.
       endif      


c    Record 2:
c
c    Version 1.0   interval  #values  #lon   #lat
c                    I*2       I*2     I*4    I*4
c            e.g      6         3      144    91  
c                 (6h intervals, 3 values, long, lat at 2.5-deg spacing)
      read(10,rec=2,iostat=ioerr) interval,ngval,nglon,nglat
      if(swapped)then
        call swap_bytes(2,interval,1)
        call swap_bytes(2,ngval,1)
        call swap_bytes(4,nglon,1)
        call swap_bytes(4,nglat,1)
      endif

c     Record 3: 

c     Version 1.0    start       stop
c                     R*4         R*4
c             e.g.   2006.0      2006.2157
c              ( 2007 day 1.0 thru day 78.75)
      read(10,rec=3) start_time,stop_time
      if(swapped)then
        call swap_bytes(4,start_time,1)
        call swap_bytes(4,stop_time,1)
      endif
      print*,'Map.grid  first line: ',version," ",mapf," ",source
      print*,'Map.grid second line:',interval,ngval,nglon,nglat
      print*,'Map.grid  third line:',start_time,stop_time    
      start_yr = int(start_time)
      start_doy=amod(start_time,1.0)*float(nydays(start_yr))+1.0    
      call round6h(start_doy) 
      stop_yr = int(stop_time)
      stop_doy=amod(stop_time,1.0)*float(nydays(stop_yr))+1.0   
      call round6h(stop_doy)             
      write(*,'(22x,a,2(i6,f7.2),a)') '( = ',start_yr,start_doy
     .     , stop_yr,stop_doy,' )'

      
c  Find the spatial grid points closest to the input coordinates

c       normal case (longitudes within one latitude ring)
c       lower left (lon1,lat1)
c       divide by nglat-1 because the files include 90S and 90N
      dlat = 180./float(nglat-1)
      dlon = 360./float(nglon)
c     lon, lat, of lower left of box
      alon1 = float(int(slon/dlon))*dlon
      if( slat.ge.0.) then
         alat1=float(int(slat/dlat))*dlat
      else
         alat1 = (float(int((90.0+slat)/dlat)))*dlat - 90.0
      endif
c        special case when slat is exactly zero
      if(slat.eq.0.0)alat1 = -dlat
c        special case when slat is exactly -90
      if(slat.eq.-90.0)alat1 = slat
c       special case when slat is exactly 90
      if(slat.eq.90.0)alat1 = slat - dlat   

c* debug
cd      print*,'LL grid lat long:',alat1,alon1
     
c  Calculate the record-number offsets from the 1st point of the epoch
c     (i.e., upper left = 0) 
c         lower left record
      irecbox(1,1) = int((90.0-alat1)/dlat)*nglon+int(alon1/dlon) 
cd        print *,'alat1 dlat alon1 dlon nhead LL '
cd     .         , alat1,dlat,alon1,dlon,nhead,irecbox(1,1)
c          lower right (lon2,lat1)
      irecbox(2,1) = irecbox(1,1) + 1
c        if lower right longitude becomes 360 deg, need to use the value
c        for zero; therefore, remove the records of a whole row of longitude
      if(alon1.eq.(360.0-dlon))then
        irecbox(2,1) = irecbox(2,1) - nglon
      endif
c        upper left (lon1,lat2)
      irecbox(1,2) = (int((90.0-alat1-dlat)/dlat))*nglon
     .                    +int(alon1/dlon) 
c         upper right (lon2,lat2)
      irecbox(2,2) = irecbox(1,2) + 1
      if(alon1.eq.(360.0-dlon))irecbox(2,2) = irecbox(2,2) - nglon
          
c  debug
cd        print*,'Box corner record numbers LL LR UR UL:'
cd     .   ,irecbox(1,1),irecbox(2,1),irecbox(2,2),irecbox(1,2)  
       

c  Read the orography to get the heights of the four box corners. 
c     The intent was to write the values at the first record of day
c     368 ( 367*91*144*4 + 3[head] +1 = 19 236 676 ); in fact in 
c     Version 1.0 they are written 4 records prior to the beginning 
c     of day 369 (368*91*144*4 + 3 +1 - 4 = 19 289 088), so put in 
c     special code to handle this.
cd      print *,'DEBUG 19 289 085 - 19 289 094 '
cd      do i=1,10 
cd       irec = 19289085+i 
cd       read(10,rec=irec) tmp,tmp,orog(1)
cd       print *,irec,tmp,tmp,orog(1)
cd      enddo
      ngrid = nglat*nglon                           
      intervald = interval/24.       
      call grid_time('T',start_doy,intervald,ngrid,nhead,369.,nrec ) 
      if( version.eq.10 ) nrec = nrec -4 
      read(10,rec=nrec+irecbox(1,1),iostat=ioerr)
     .      tmp,tmp,orog(1)
      read(10,rec=nrec+irecbox(2,1),iostat=ioerr)
     .      tmp,tmp,orog(2)
      read(10,rec=nrec+irecbox(2,2),iostat=ioerr)
     .      tmp,tmp,orog(3)
      read(10,rec=nrec+irecbox(1,2),iostat=ioerr)
     .      tmp,tmp,orog(4)
      if(swapped)  call swap_bytes(4,orog,4)
cd      print*,'nrec, box corners for orography:',nrec,nrec+irecbox(1,1)
cd     .   ,nrec+irecbox(2,1),nrec+irecbox(2,2),nrec+irecbox(1,2)
cd      print*,'Orography values are:',orog
      dy  = abs(slat-alat1)/dlat
      dy1 = (dlat - dy*dlat)/dlat
      dx  = (slon - alon1)/dlon
      dx1 = (dlon - dx*dlon)/dlon 
c then
      height = dx1*dy1*orog(1)+dx*dy1*orog(2)+dx1*dy*orog(4)
     .                 +dx*dy*orog(3) 
      write(*,'(a,f10.3,a)')'Interpolated height is: ',height
     .  ,' m. Interpolated ZHD value refers to this ellipsoidal height.'
      

c  Bon, maitenant on a tous qu'il faut pour tout faire. Ca ne reste
c  que voyager dedans le grand fichier, sauter jusqu'au le premier
c  jour, calculer heuresement jusqu'au le dernier jour, arreter la
c  bas et prendre un apero. Sante!

c  That is to say, Get the values for the position and time requeted


c Now loop through the requested epochs
 
      t = ts  
      do while(t.le.(tend+intervald))   
        call grid_time('T',start_doy,intervald,ngrid,nhead,t,nrec )   
cd        print *,'start_doy intervald ngrid nhead t nrec '
cd     .        , start_doy,intervald,ngrid,nhead,nrec
c       VMF1 Ver 1.0 has Ah Aw I*4, ZHD R*4
        read(10,rec=irecbox(1,1)+nrec,iostat=ierr) Ah(1),Aw(1),ZHD(1)
        read(10,rec=irecbox(2,1)+nrec,iostat=ierr) Ah(2),Aw(2),ZHD(2)
        read(10,rec=irecbox(2,2)+nrec,iostat=ierr) Ah(3),Aw(3),ZHD(3)
        read(10,rec=irecbox(1,2)+nrec,iostat=ierr) Ah(4),Aw(4),ZHD(4)   
cd        print *,'Rec # LL LR UR UL ',irecbox(1,1)+nrec,irecbox(2,1)+nrec
cd     .               ,irecbox(2,2)+nrec,irecbox(1,2)+nrec
cd        print *,'Ah Values LL LR UR UL :',Ah  
        if( ioerr.ne.0 ) then  
           print *,'ioerr rec ',ioerr,irecbox(1,1)+nrec 
           write(*,'(a,i8,a,f6.2)') 
     .       'Attempt to read beyond EOF, record = '
     .         ,nrec,'   t=',t
          stop
        else
          if( swapped ) then
            call swap_bytes(4,Ah,4)
            call swap_bytes(4,Aw,4)
            call swap_bytes(4,ZHD,4)
          end if
        endif
c      perform a bi-linear interpolation
c      PT040106: I found a simple, one-line formula that seems to do this!
c               VAL = DX1*DY1*ARRAY(ILO,JLO) + DX*DY1*ARRAY(IHI,JLO) +
c               DX1*DY*ARRAY(ILO,JHI) + DX*DY*ARRAY(IHI,JHI)
        vmf(1) = dx1*dy1*Ah(1)+dx*dy1*Ah(2)+dx1*dy*Ah(4) + dx*dy*Ah(3) 
        vmf(2) = dx1*dy1*Aw(1)+dx*dy1*Aw(2)+dx1*dy*Aw(4) + dx*dy*Aw(3) 
        zhdval = dx1*dy1*ZHD(1)+dx*dy1*ZHD(2)+dx1*dy*ZHD(4)+dx*dy*ZHD(3)

c PT060926:  correct the zhd value for any height difference
c            In fact, we can't do this because we don't know
c            the ellipsoidal height of the input coords!
c       call correct_zhd(slat,height,ell_ht,zhdval)
          
        write(*,'(f6.2,2i10,f11.4)') t,vmf,zhdval


        t = t + intervald

      enddo

c eh voila, c'est finit. Sante!
 
      close(10)
      end




