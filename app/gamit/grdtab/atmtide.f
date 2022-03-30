      program atmtide


c  program to extract the cosine and sine coefficients for the S1 and S2 tides
c  for a given lat/long. This will be called from script sh_plotatmtide   
c
c  The code has been cut/pasted largely from grdtab.f and interp_atmtide.f 
c
c  P. Tregoning
c  21 September 2005
c
c PT050928: modified to read a phase shift from line 2 of the atmtide gridfile

      implicit none
             
      integer ngval,nglon,nglat,ioerr,lu6,ngval0,recl
      integer*2 version
      real*4 slat,slon,fversion
      real*8 xyz(3),geod_pos(3),rot_mat(3,3)    
      logical swapped
      character arg*20, press*5,frame*2,source*3,tide*5,message*256
     .         ,atmgrid*12


c     Record #s and values read from the grid file (linear interpolaton)
      integer*4 irecbox(2,2)
      real*4   cS1(3,4),cS2(3,4),sS1(3,4),sS2(3,4),val(12) 

c      Local variables
      integer*4 i
      real*4 alat1,alon1,dlat,dlon,phase

      lu6 = 10
      atmgrid = "atmtide.grid"
c  PT050510: the program tvd2grid converts files of cosine/sine amplitudes of loading (in mm)
c            for N, E, U. 
c            There are now 3 components x 2 amplitudes x 2 tides x R*4 values = recl 48
      ngval0 = 12
      swapped = .false.
      recl = 48

c  decode command line arguments
      call getarg(1,arg)
      if(arg(1:1).eq.' ')then
         print*,'Runstring: atmtide < X >  < Y > < Z > '
         stop 'No command line arguments found. Program stopped'
      else
c       get the XYZ site coordinates 
         read(arg,*)xyz(1)
         call getarg(2,arg)
         read(arg,*)xyz(2)
         call getarg(3,arg)
         read(arg,*)xyz(3)

c  convert XYZ to lat/long/height
         call xyz_to_geod(rot_mat, xyz, geod_pos )
         slat = 90.0 - geod_pos(1) * 180.0/3.141529653
         slon = geod_pos(2) * 180.0/3.141529653
      endif

c  open the atm tide grid file
       open(unit=lu6,file="atmtide.grid",form='unformatted'
     .        ,access='direct'
     .        ,status='old',recl=recl,iostat=ioerr)
       if(ioerr.ne.0)then
          print*,'Error opening file atmtide.grid. Does it exist?'
          stop 'Stopped with errors'
       endif

c  PT050104: first line is now a header line
        read(lu6,rec=1,iostat=ioerr) version,press,frame,source,tide
        if(version.gt.100.or.version.lt.0)then
*           If the grid size does not seem be correct, swap the bytes
*           and see if OK.
          call swap_bytes(2, version, 1)
*           Now see if OK size.
          if(version.gt.100.or.version.lt.0)then
*               Still a problem so kill.
             call report_stat('FATAL','GRDTAB','atmtide',
     .               atmgrid,
     .              'Invalid version number of atmgrid file',0)
          else
*               File seems to be swapped.  Set the swapped status and continue.
*               (We now need to swap every record read)
                call report_stat('STATUS','GRDTAB','atmtide'
     .             ,atmgrid, 'Byte-swapping input grid file',0)
                swapped = .true.      
          end if
        end if

c PT050104: write out header information to the output atmtid.DDD file
        fversion = version
        fversion = fversion/10.  
        write(message,'(a,f5.1,1x,a5,1x,a2,1x,a2,1x,a5)') 
     .   'Grid file version tides field, frame '
     .   ,fversion,press,frame,source,tide
        call report_stat('STATUS','GRDTAB','atmtide',' ',message,0)

c PT050505: the second record is blank - left for future expansion if necessary
c PT050928: first value on second line is now a phase shift for the whole file
c           (to account for the 12 hr difference between Leonid's and Tonie's files
        read(lu6,rec=2,iostat=ioerr) phase

c       read the 3rd record to get array limits
        read(lu6,rec=3,iostat=ioerr) ngval,nglon,nglat
        if(swapped)then
            call swap_bytes(4,ngval,1)
            call swap_bytes(4,nglon,1)
            call swap_bytes(4,nglat,1)
        end if      

        call report_stat('STATUS','GRDTAB','atmtide',' '
     .                ,'Interpolating atm tide grid file',' ')




c Now, all the same logic from interp_atmtide.f can be used to extract the
c cosine and sine coefficients - including the bilinear interpolation

c  convert longitude range from -180/180 to 0/360
      if(slon.lt.0)slon = 360.0 + slon

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

c      Get record numbers to be read from the grid file
c       (for linear interpolation, need 4)
c       (these calculations count from beginning of grid, assuming
c        that there is one header record stating the number of 
c        records/line and the grid dimensions)
                                                                        
c       normal case (longitudes within one latitude ring)
c       lower left (lon1,lat1)
c PT040106: NOTE that the ATM grids go from 90N to 90S 
c  lower left record       
c PT050110: with an additional header line, add "2" to the
c           computation rather than "1"
c PT050505: with a spare header line added into the atm grid file, add "3" rather than "2"
        irecbox(1,1) = int((90.0-alat1)/dlat)*nglon + int(alon1/dlon)+3

c       lower right (lon2,lat1)
      irecbox(2,1) = irecbox(1,1) + 1
c  if lower right longitude becomes 360 deg, need to use the value
c  for zero; therefore, remove the records of a whole row of longitude
      if(alon1.eq.(360.0-dlon))then
        irecbox(2,1) = irecbox(2,1) - nglon
      endif
                                                                        
c     upper left (lon1,lat2)
c PT050511: with three header lines, add "3" to the computation 
c        irecbox(1,2) = (int((90.0-alat1-dlat)/dlat))*nglon
        irecbox(1,2) = (int((90.0-alat1-dlat)/dlat))*nglon
     .                               +int(alon1/dlon)+3
c       upper right (lon2,lat2)
      irecbox(2,2) = irecbox(1,2) + 1
      if(alon1.eq.(360.0-dlon))irecbox(2,2) = irecbox(2,2) - nglon


c  There is only one epoch of the atm tide grid. So simply extract from the grid the
c  indices as calculated above. Note that there are 12 values per line, being
c  cS1 sS1 for north, then east, then up, followed by the same for S2.


c  need to read only one full grids per day
      read(lu6,rec=irecbox(1,1),err=200)
     .     (cS1(i,1),sS1(i,1),i=1,3),(cS2(i,1),sS2(i,1),i=1,3)   

c      print*,'cosine/sine for S1 and S2 N,E,U for lleft corner'
c     . ,(cS1(i,1),sS1(i,1),i=1,3),(cS2(i,1),sS2(i,1),i=1,3)   


      read(lu6,rec=irecbox(2,1),err=200)
     .     (cS1(i,2),sS1(i,2),i=1,3),(cS2(i,2),sS2(i,2),i=1,3)   
      read(lu6,rec=irecbox(2,2),err=200)
     .     (cS1(i,3),sS1(i,3),i=1,3),(cS2(i,3),sS2(i,3),i=1,3)   
      read(lu6,rec=irecbox(1,2),err=200)
     .     (cS1(i,4),sS1(i,4),i=1,3),(cS2(i,4),sS2(i,4),i=1,3)   

c      print*,'box corners ',irecbox(1,1),irecbox(2,1)
c     .                       ,irecbox(1,2),irecbox(2,2)
                                                                  
*     If file is byte-swapped, then switch the values just read
      if( swapped ) then
        call swap_bytes(4,cS1,12)
        call swap_bytes(4,cS2,12)
        call swap_bytes(4,sS1,12)
        call swap_bytes(4,sS2,12)
      end if

c  now perform a bi-linear interpolation
c
c S1 north Cosine and Sine
      call bilin4(slon,slat,alon1,alat1,dlon,dlat
     .          ,cS1(1,1),cS1(1,2),cS1(1,3),cS1(1,4)
     .          ,val(1))
      call bilin4(slon,slat,alon1,alat1,dlon,dlat
     .          ,sS1(1,1),sS1(1,2),sS1(1,3),sS1(1,4)
     .          ,val(2))  


c S1 east Cosine and Sine
      call bilin4(slon,slat,alon1,alat1,dlon,dlat
     .          ,cS1(2,1),cS1(2,2),cS1(2,3),cS1(2,4)
     .          ,val(3))
      call bilin4(slon,slat,alon1,alat1,dlon,dlat
     .          ,sS1(2,1),sS1(2,2),sS1(2,3),sS1(2,4)
     .          ,val(4))

c S1 Up Cosine and Sine
      call bilin4(slon,slat,alon1,alat1,dlon,dlat
     .          ,cS1(3,1),cS1(3,2),cS1(3,3),cS1(3,4)
     .          ,val(5))
      call bilin4(slon,slat,alon1,alat1,dlon,dlat
     .          ,sS1(3,1),sS1(3,2),sS1(3,3),sS1(3,4)
     .          ,val(6))

c S2 north Cosine and Sine
      call bilin4(slon,slat,alon1,alat1,dlon,dlat
     .          ,cS2(1,1),cS2(1,2),cS2(1,3),cS2(1,4)
     .          ,val(7))
      call bilin4(slon,slat,alon1,alat1,dlon,dlat
     .          ,sS2(1,1),sS2(1,2),sS2(1,3),sS2(1,4)
     .          ,val(8))

c S2 east Cosine and Sine
      call bilin4(slon,slat,alon1,alat1,dlon,dlat
     .          ,cS2(2,1),cS2(2,2),cS2(2,3),cS2(2,4)
     .          ,val(9))
      call bilin4(slon,slat,alon1,alat1,dlon,dlat
     .          ,sS2(2,1),sS2(2,2),sS2(2,3),sS2(2,4)
     .          ,val(10))

c S2 Up Cosine and Sine
      call bilin4(slon,slat,alon1,alat1,dlon,dlat
     .          ,cS2(3,1),cS2(3,2),cS2(3,3),cS2(3,4)
     .          ,val(11))
      call bilin4(slon,slat,alon1,alat1,dlon,dlat
     .          ,sS2(3,1),sS2(3,2),sS2(3,3),sS2(3,4)
     .          ,val(12))

c  write out the values to the ascii file
c  PT050928: .... including the phase shift value
      write(*,100)slat,slon,val,phase
100   format("ATMTIDE ",2f12.7,1x,12f12.8,f8.4)

      stop 'Ended normally'


200   continue
c  we will have branched here only if there was an error reading the grid file
           call report_stat('FATAL','GRDTAB','atmtide',' '
     .                    ,'Error reading grid file',' ')


      end


      subroutine bilin4(x,y,x0,y0,xstep,ystep,u1,u2,u3,u4,val)

c   computes a bilinear interpolation of the value of a single
c   point within a rectangle. I found this formula via google
c   and it seems to work!
c
c   INPUT:
c         x,y  :  the coordinates of the point to be interpolated
c        x0,y0 :  the coords of the lower left corner
c     xstep,ystep : the grid step sizes
c         u1   : value at lower left corner
c         u2   : value at lower right corner
c         u3   : value at upper right corner
c         u4   : value at upper left corner
c
c    OUTPUT:
c        val   : the interpolated value at the requested point
c
c   P. Tregoning
c   7 January 2004

      implicit none

c  argument variables
      real*4 x,y,x0,y0,xstep,ystep,u1,u2,u3,u4,val

c  local variables
      real*4 dy,dx,dy1,dx1

          dy  = abs(y-y0)/ystep
          dy1 = (ystep - dy*ystep)/ystep
          dx  = (x - x0)/xstep
          dx1 = (xstep - dx*xstep)/xstep
c then
          val = dx1*dy1*u1+dx*dy1*u2+dx1*dy*u4+dx*dy*u3

      return
      end






