      program imftoasc
c
c  PT050908: this is program atmtoasc but modified to read the
c            imf grid files instead. The only differences are:
c
c 1. The IMF grid file doesn't have a header line. It starts
c    with ngval nlon nlat (ie line 2 of the atmgrid file)
c 2. Therefore, the counter for the line in the file needs to 
c    be one less than for atmtoasc (ie back to 1 from 2)

c   program to read a binary atmospheric grid file and extract
c   the displacements for a single geographic location for a
c   particular time interval. This may be a useful utility but
c   I am writing it primarily to see whether I can actually
c   read and interpolate Tonie van Dam's binary grid files
c   correctly
c
c   P. Tregoning
c   6 January 2004 
c
C  

      implicit none

      real*4 slat,slon,alat1,alon1,dlat,dlon
     .      ,n(4),e(4),u(4),displ(3),dx,dx1,dy,dy1
      integer*4 ngval,nglon,nglat,recl,ierr,ts,tend
     .         ,irecbox(2,2),i,ncurrent,nfinish,nprior
     .         , ioerr

      character atmfile*12,arg1*30,arg2*3
                       

c PT050104: variables for the header line of the atmdisp files
      integer*2  version
      character*2 spare
      character*5 press_src
      character*3 ref_frame

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
         print*,'Runstring: imftoasc imffile lat long doy_st doy_end'
         stop 'No command line arguments found. Program stopped'
      else
         open(unit=10,file=atmfile,form='unformatted',access='direct'
     .        ,status='old',recl=12,iostat=ierr)
         if(ierr.ne.0)then
            print*,'Error opening file',atmfile
            stop 'Stopped with errors'
         endif
       endif

c  get the site coordinates
      call getarg(2,arg1)
      read(arg1,*)slat
      call getarg(3,arg1)
      read(arg1,*)slon

c  get the start and end days
      call getarg(4,arg2)
      read(arg2,'(i3.3)')ts
      call getarg(5,arg2)
      read(arg2,'(i3.3)')tend

c      print*,'Input parameters:',atmfile,slat,slon,ts,tend
                         
c  PT050104: first line is now a header line for atmtoasc 
c  PT050908: but not for imftoasc at this stage!
c        read(10,rec=1,iostat=ioerr) version,spare,press_src,ref_frame 
c        print*,version,spare,press_src,ref_frame
c       read the next record to get array limits
        read(10,rec=1,iostat=ioerr) ngval,nglon,nglat
      print*,'imfgrid first line',ngval,nglon,nglat

        if(ngval.ne.3)then
*           If the grid size does not seem be correct, swap the bytes
*           and see if OK.
          call swap_bytes(4, ngval, 1)
*           Now see if OK size.
          if(ngval.ne.3) then
*               Still a problem so kill.
             call report_stat('FATAL','UTILS','imftoasc',
     .          atmfile,'Invalid number of variables on imfgrid file',0)
          else
*               File seems to be swapped.  Set the swapped status and continue.
*               (We now need to swap every record read)
                call report_stat('STATUS','UTILS','imftosac'
     .             ,atmfile, 'Byte-swapping input grid file',0)
                swapped = .true.
                call swap_bytes(4,ngval,1)
                call swap_bytes(4,nglon,1)
               call swap_bytes(4,nglat,1)
          end if
        end if

c  Now find the four grid points that are closest
c  to the input coordinates.
c       normal case (longitudes within one latitude ring)
c       lower left (lon1,lat1)
c PT040106: NOTE that the ATM grids go from 90N to 90S whereas
c           the tide loading models go from 90S to 90N
c     latitude and longitude intervals (as per the file header )
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
c PT050110: with an additional header line, add "2" to the
c           computation rather than "1"
c PT050908: for the IMF grids, there currently aren't two header lines, so make this 1
        irecbox(1,1) = int((90.0-alat1)/dlat)*nglon+int(alon1/dlon) + 1
c       lower right (lon2,lat1)
        irecbox(2,1) = irecbox(1,1) + 1
c  if lower right longitude becomes 360 deg, need to use the value
c  for zero; therefore, remove the records of a whole row of longitude
        if(alon1.eq.(360.0-dlon))then
          irecbox(2,1) = irecbox(2,1) - nglon
        endif

c       upper left (lon1,lat2)
c PT050110: with an additional header line, add "2" to the
c           computation rather than "1"
c PT050908: for the IMF grids, there currently aren't two header lines, so make this 1
        irecbox(1,2) = (int((90.0-alat1-dlat)/dlat))*nglon
     .                               +int(alon1/dlon) + 1
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
      nprior = (ts-1)*42048
      nfinish = (tend+1)*42048
      ncurrent = nprior + 1

c debug 
c      print*,'ts,nprior,nfinish,ncurrent',ts,nprior,nfinish,ncurrent

      do while (ncurrent.lt.nfinish)

c  need to read four full grids per day
        do i=1,4
c  correct the record indicator for prior days and increment it for
c  days/times already read
c  and  read the three displacement values at these four corner points.
          read(10,rec=irecbox(1,1)+ncurrent,iostat=ierr)n(1),e(1),u(1) 
c      print*,'LL values:',n(1),e(1),u(1)
          read(10,rec=irecbox(2,1)+ncurrent,iostat=ierr)n(2),e(2),u(2)
          read(10,rec=irecbox(2,2)+ncurrent,iostat=ierr)n(3),e(3),u(3)
          read(10,rec=irecbox(1,2)+ncurrent,iostat=ierr)n(4),e(4),u(4)

*     If file is byte-swapped, then switch the values just read
          if( swapped ) then
            call swap_bytes(4,n,4)
            call swap_bytes(4,e,4)
            call swap_bytes(4,u,4)
          end if

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
          ncurrent = ncurrent + 10512
        enddo
      enddo

c eh voila, c'est finit. Sante!
 
      close(10)
      end


