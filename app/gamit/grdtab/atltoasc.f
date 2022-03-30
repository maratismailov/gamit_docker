      Program ATLTOASC

c   Program to read a binary atmospheric tidal loading grid file and 
c   extract the displacements for a single geographic location. This 
c   utility may be used to check the writing of the binary atl.grid
c   file (program ATL_ASCTOBIN) or to check a particular loading
c   displacement.  
c
c   R. King  23 January 2013, modeled on ATMTOASC.
          
c Version 1.0 (ANU file) has constant coefficients, 3 header records, and
c    a record length of 48 bytes      
c    Record 1     
c      nversion,press_src,ref_frame,source,tide
c        I*2        C*5       C*3      C*3   C*5 (these not used?)      
c     Record 2 is blank - left for future expansion 
c     Record 3  -  array limits 
c        ngtid   nglon   nglat
c        I*4      I*4     I*4
c     Data records
c       rval(i),i=12
c         R*4  

c Version 2.0 (LU file ) has monthly coefficients, 1 header record, and
c    a record length of 288
c   Header record
c     version  nhead  model frame interval  ngval  nglat nglon 
c       I*2     I*2    C*8   C*3    R*4      I*4    I*4   I*4
c
c   Data records
c      144 I*2 values representing 12 coefficients for 12 months at
c      a particular lat and lon value.  The order of the records
c      (lon,lat) is 0,90 0,89, ..0,-90,...1,90,....360,-90  

      implicit none
               
      real*4 day,slat,slon,alat1,alon1,dlat,dlon
     .      ,dx,dx1,dy,dy1,tyr,atides(2,6),pi,fs1,fs2,time
      integer*4 grecl,ioerr,irecbox(2,2),nrec,ntid,nepc,iday,ihr,i,j
      character atlfile*128,arg*30 

c        header values
      integer*2 version,nhead
      integer*4 nglat,nglon,ngval
      character*8 atidemod
      character*3 ref_frame,tide 
      character*5 press_src
             
c        data values        
      integer*4 maxgval
      parameter(maxgval=144)
      integer*2 ival(maxgval,4)    
      real*4 rval(maxgval,4)

      real*4 tval(maxgval),disp(6)
      logical swapped,debug
                                             
      data pi/3.14159265/

      swapped = .false.

c  each line (except the header) has 12x12 integer*2 values. 
      grecl = 288

c   Decode the command line  

c     input grid file name
      call getarg(1,atlfile)
      if(atlfile(1:1).eq.' ')then
         print*,'Runstring: atltoasc atlfile lat long doy [debug] '
         stop 'No command line arguments found. Program stopped'
      endif
c     site coordinates
      call getarg(2,arg)
      read(arg,*) slat
      call getarg(3,arg)
      read(arg,*) slon
c     epoch doy 
      call getarg(4,arg)
      read(arg,*) day                 
c     see if debug
      call getarg(5,arg)
      if( arg(1:5).eq.'debug') then
        debug = .true.
      else
        debug = .false.
      endif
      write(*,'(a)') 'ATLTOASC input arguments: '
      write(*,'(2a)') '  Binary atl file: ',trim(atlfile)
      write(*,'(a,2f9.3)') '  Lat lon : ',slat,slon
      write(*,'(a,f5.0)') '  Day-of-year : ',day
      write(*,'(a,l1)') '  Debug = ',debug

   
c  Provisionally open the grid file to determine the version

      open(unit=10,file=atlfile,form='unformatted',access='direct'
     .        ,status='old',recl=grecl,iostat=ioerr)
      if(ioerr.ne.0) then 
        write(*,'(3a,i4)') 'Error opening ',atlfile,' iostat=',ioerr  
        stop
      endif
      read(10,rec=1,iostat=ioerr) version 
      if(debug) print*,'DEBUG Initial read ',version
      if(version.gt.100.or.version.lt.0)then
*         If the grid size does not seem be correct, swap the bytes
*         and see if OK.   
        write(*,'(a,i4,a)') '**version =',version,' swapping**'
        call swap_bytes(2,version,1)      
        if( debug ) print *,'Swapped version ',version
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
            
c  Now we know the version and record length, read the headers for real

      if( version.eq.10 ) then   
        close(10)                   
        grecl = 48
        open(unit=10,file=atlfile,form='unformatted',access='direct'
     .        ,status='old',recl=grecl,iostat=ioerr)
c        nversion,press_src,ref_frame,source,tide
c          I*2        C*5       C*3      C*3   C*5 (these not used?)      
        read(10,rec=1,iostat=ioerr) version,press_src,ref_frame,tide
        if(ioerr.ne.0 ) then
         write(*,'(a,i4)') 'Error reading Record 1 of Version 1.0 file '
     .      ,ioerr 
         stop
        endif 
        if( swapped ) call swap_bytes(2,version,1)  
        write(*,'(/,2a,1x,i4,1x,a5,1x,a3,1x,a3,1x,a5)') 
     .       'Read Record 1 Version 1.0 '
     .      ,' version press_src ref_frame tide '
     .      ,version,press_src,ref_frame,tide 
        write(*,'(/,a)') 'Record 2 unused'
        read(10,rec=3,iostat=ioerr) ngval,nglon,nglat  
        if( ioerr.ne.0 ) then
          write(*,'(a,i4)') 'Error reading record 3 of Version 1.0 file'
     .      ,ioerr
          stop
        endif
        if( debug ) print *,'unswapped ngval nglon nglat '
     .        ,ngval,nglon,nglat
        if( swapped ) then
           call swap_bytes(4,ngval,1)
           call swap_bytes(4,nglon,1)
           call swap_bytes(4,nglat,1)                         
           if( debug ) print *,'swapped ngval,nglon,nglat '
     .                        , ngval,nglon,nglat
        endif                        
        write(*,'(/,a,3i4)') 'Read Record 3: ngval,nglon,nglat '
     .                        ,ngval,nglon,nglat
 
      elseif( version.eq.20 ) then

c         version  nhead  model frame  ngval  nglat nglon 
c           I*2     I*2    C*8   C*3    I*4    I*4   I*4
        read(10,rec=1,iostat=ioerr) version,nhead,atidemod,ref_frame
     .      ,ngval,nglat,nglon                       
        if(ioerr.ne.0 ) then
         write(*,'(a,i4)') 'Error reading Record 1 of Version 2.0 file '
     .      ,ioerr 
         stop
        endif 
        if( debug ) print *,'DEBUG swapped ',swapped    
        if(swapped)then
           call swap_bytes(2,nhead,1)
           call swap_bytes(4,ngval,1)
           call swap_bytes(4,nglon,1)
           call swap_bytes(4,nglat,1)
        endif
        write(*,'(/,2a,1x,i4,1x,i2,1x,a8,1x,a3,3i4)')
     .       'Read Record 1 Version 1.0 '
     .      ,' version nhead atidemod ngval nglat nglon'
     .      ,version,nhead,atidemod,ref_frame,ngval,nglat,nglon 
      else
         write(*,'(a,i2)') 'Unknown version ',version
         stop
      endif
      if( ngval.gt.maxgval ) then
        write(*,'(a,i5,a,i5,a)') 
     .     'ngval=',ngval,' > maxgval=',maxgval,' stop'
        stop
      endif
                     
c*** FOR DEBUG simply read the first 200 datarecords
cd      do i=2,200
cd        read(10,rec=i,iostat=ioerr) (ival(j,1),j=1,ngval)
cd        print *,i,(ival(j,1),j=1,ngval)
cd      enddo
cd      stop

c  Now find the four grid points that are closest to the input coordinates.
 
c The grid files go from 90N to 90S, and , and within each lat group, from 0 to 360E.
c (opposite from Tonie's ascii files but consistent with the ATML, OTL, and VMF1
c grid files).
c PT040106: divide by nglat-1 because the files include 90S and 90N
      dlat = 180./float(nglat-1)
      dlon = 360./float(nglon-1)
      if( debug ) print *,'DEBUG nglat dlat nglon dlon '
     .    ,nglat,dlat,nglon,dlon
      if( debug ) print *,'Grid intervals (lon,lat) are ',dlon,dlat
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

c  Calculate the record-number offsets 
c     lower left record
      irecbox(1,1) = int((90.0-alat1)/dlat)*nglon +
     .      int((alon1+dlon)/dlon) + nhead 
      if( debug ) 
     .   print *,'alat1 dlat alon1 dlon nhead LL '
     .         , alat1,dlat,alon1,dlon,nhead,irecbox(1,1)
       
c     lower right (lon2,lat1)
      irecbox(2,1) = irecbox(1,1) + 1
c     if lower right longitude becomes 360 deg, need to use the value
c     for zero; therefore, remove the records of a whole row of longitude
      if(alon1.eq.(360.0-dlon))then
        irecbox(2,1) = irecbox(2,1) - nglon
      endif
c     upper left (lon1,lat2)
      irecbox(1,2) = (int((90.0-alat1-dlat)/dlat))*nglon
     .                      +int(alon1/dlon) + 1 + nhead
c     upper right (lon2,lat2)
      irecbox(2,2) = irecbox(1,2) + 1
      if(alon1.eq.(360.0-dlon))irecbox(2,2) = irecbox(2,2) - nglon
      if(debug)  print*,'Box corner record numbers LL LR UR UL:'
     .   ,irecbox(1,1),irecbox(2,1),irecbox(2,2),irecbox(1,2)


C  Read the records at the  four corners

      if( version.eq.10 ) then

        read(10,rec=irecbox(1,1),iostat=ioerr) (rval(i,1),i=1,ngval)
        if( swapped ) call swap_bytes(4,rval(1,1),ngval)
        read(10,rec=irecbox(2,1),iostat=ioerr) (rval(i,2),i=1,ngval) 
        if( swapped ) call swap_bytes(4,rval(1,2),ngval)
        read(10,rec=irecbox(2,2),iostat=ioerr) (rval(i,3),i=1,ngval) 
        if( swapped ) call swap_bytes(4,rval(1,3),ngval)
        read(10,rec=irecbox(1,2),iostat=ioerr) (rval(i,4),i=1,ngval)     
        if( swapped ) call swap_bytes(4,rval(1,4),ngval)
        if( debug ) 
     .   print *,'S1cos U values at 4 corners : ',(rval(1,j),j=1,4) 
        
      elseif (version.eq.20 ) then


        read(10,rec=irecbox(1,1),iostat=ioerr) (ival(i,1),i=1,ngval)
cd      print *,'DEBUG irec ival ',irecbox(1,1),(ival(i,1),i=1,ngval)
cd      print *,'stop'
cd      stop
        read(10,rec=irecbox(2,1),iostat=ioerr) (ival(i,2),i=1,ngval)
        read(10,rec=irecbox(2,2),iostat=ioerr) (ival(i,3),i=1,ngval)
        read(10,rec=irecbox(1,2),iostat=ioerr) (ival(i,4),i=1,ngval)
        if( ioerr.ne.0 ) then  
          print *,'ioerr rec ',ioerr,irecbox(1,1)
          write(*,'(a,i8)') 
     .      'Attempt to read beyond EOF, record = ',nrec
          stop
        else
          if( swapped ) then 
            call swap_bytes(2,ival(1,1),ngval)    
            call swap_bytes(2,ival(1,2),ngval)  
            call swap_bytes(2,ival(1,3),ngval)
            call swap_bytes(2,ival(1,4),ngval)
          endif
        endif   
        do i=1,ngval 
          do j=1,4                  
c           convert to mm 
            rval(i,j) = float(ival(i,j))*1.e-4
          enddo
        enddo 
      endif

      if( debug .and. version.eq.20 ) then    
        print *,'Monthly S1cos U values at corners '
        do j=1,4
          write(*,'(11f9.4)') (rval(i,j),i=1,ngval,12)
        enddo
      endif 


c Perform a bi-linear interpolation 

      dy  = abs(slat-alat1)/dlat
      dy1 = (dlat - dy*dlat)/dlat
      dx  = (slon - alon1)/dlon
      dx1 = (dlon - dx*dlon)/dlon
      if( debug ) then 
        print *,'dy dy1 dx dx1 ',dy,dy1,dx,dx1 
cd       print *,'rval(1) ',(rval(i,1),i=1,ngval)
      endif
      do i=1,ngval
        tval(i)=  dx1*dy1*rval(i,1) + dx*dy1*rval(i,2)
     .          + dx1*dy *rval(i,4) + dx*dy *rval(i,3) 
      enddo
cd    if( debug ) then 
cd      write(*,'(a)') 'irec  corners for 1st coefficient'
cd      write(*,'(i5,f7.4)') irecbox(1,1),rval(1,1)
cd      write(*,'(i5,f7.4)') irecbox(2,1),rval(1,2)
cd      write(*,'(i5,f7.4)') irecbox(2,2),rval(1,3)
cd      write(*,'(i5,f7.4)') irecbox(1,2),rval(1,4)                   
cd      if( version.eq.10 ) then
cd        write(*,'(a,f7.4)') 'Interpolated 1st coeff ',tval(1)
cd      elseif (version.eq.20) then
cd         write(*,'(a,f7.4)') 'Interpolated Jan 1st coeff ',tval(1)
cd         write(*,'(a,f7.4)') 'Interpolated Dec 1st coeff ',tval(133)
cd      endif  
cd    endif

c** extra debug: read the value for record 1085, which I think should
cd     be lon 0, lat 87 
cd      read(10,rec=1085,iostat=ioerr) ival(1,1)
cd      print *,'0,87 irec 1085 ',ival(1,1)
cd      read(10,rec=726,iostat=ioerr) ival(1,1)
cd      print *,'2,88 irec 726 ',ival(1,1)
 
c  Assign the coefficients in the order expected by GAMIT grdtab

      if( version.eq.10 ) then        
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
        write(*,'(/,a,2f10.2)') 'Coefficents for lat/lon ',slat,slon
        
      elseif( version.eq.20 ) then   
c       if time-dependent, interpolate from monthly values
        iday = nint(day)                
        call interp_months( ngval,tval,iday,atides )
        tyr = day/365.25                                     
        write(*,'(/,a,2f10.2,a,i4,a,f6.3)') 'Coefficents for lat/lon '
     .     ,slat,slon,' on DOY ',iday,' FracYr ',tyr   
      endif
      write(*,'(a,6f8.3)') ' S1 UNE (mm); ',(atides(1,i),i=1,6)
      write(*,'(a,6f8.3)') ' S2 UNE (mm); ',(atides(2,i),i=1,6)     
      write(*,'(/,a)') 'Up displacement (mm) by hours' 
      write(*,'(a)') '   0      4      8     12      16      20 '
      i = 0      
      fs1=2.*pi
      fs2=4.*pi
      do ihr = 0,24,4
        i = i+1
        time = float(ihr)/24.
        disp(i) = atides(1,1)*cos(fs1*time) + atides(1,2)*sin(fs1*time)
     .     +      atides(2,1)*cos(fs2*time) + atides(2,2)*sin(fs2*time)
cd       print *,'S1 ihr freq t phi ',ihr,fs1,time,fs1*time
cd       print *,'S2 ihr freq t phi ',ihr,fs2,time,fs2*time
      enddo
      write(*,'(6f7.3)') (disp(i),i=1,6)

      end
                            


