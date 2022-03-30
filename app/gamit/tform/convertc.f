      Program CONVERTC

 
c     Generalized version of gamit/utils/convertc to convert coordinates 
c     from one type to another; replaces gamit/utils/convertc and largely
c     supersedes gamit/tform/tform, except for differencing files, conversions
c     to cylindrical coordinates, and geodetic datums other than WGS84.

c     R. King 27 February 2008

C     CONVERTC can read input either three input coordinates (XYZ or lat/lon/ht)
c     or a file of one of following types:
c
c       GAMIT l-file  (spherical, deg/min/sec, radius)
c       GLOBK apr file (Cartesian position and velocity)
c       GLOBK vel file (geodetic lon/lat/ht and ENU velocity)
c       GLOBK glist file (geodetic lon/lat)

c     The output file types supported are:
c
c       GAMIT old-style L-file (spherical lat/lon/radius)
c       GLOBK apr file (new GAMIT L-file, XYZ and velocities)  
c       GLOBK vel file (lon/lat ENU velocities)
c       KML file for Google Earth (geodetic lat/lon/ht)

      implicit none

      include '../../libraries/includes/const_param.h'   

      integer*4 iarg,iclarg,ioerr,luin,luout,nblen,i
         
      real*8 inval(6),x(3),xdot(3),tstrad,latd,lond,ht,date
     .     , geod_pos(3),rot(3,3),rottr(3,3),neudot(3),loc_coord(3)
     .     , lfile_epoch
                           
      character*3 intype,outtype
      character*8 site 
      character*32 arg,infile,outfile   
      character*256 line

      logical infile_flag,vel_flag,values_ok,finished,comment

c  Output program header

c**      call tversn  

c  Assign unit numbers  (all in common in tform.h)
              
      luin = 1
      luout= 2
                    
c  Initialize values that might not be set
      do i=1,3
        xdot(i) = 0.d0
      enddo                                            
      date = 0.d0
   
c  Get the input values
             
      iarg = iclarg(1,arg) 
      if( iarg.le.0 ) then
        call convertc_help 
        stop
      endif
c     read the first argument to see whether coordinates or a file name
        read(arg,*,iostat=ioerr) inval(1)          
        if(ioerr.eq.0) then     
          infile_flag = .false.
c         argument is a numerical value     
          iarg = iclarg(2,arg) 
          read(arg,*,iostat=ioerr) inval(2)   
          if( ioerr.ne.0) then
            write(*,'(a)') 'Error reading second coordinate'
            stop
          endif  
          iarg = iclarg(3,arg)
          read(arg,*,iostat=ioerr) inval(3)
          if( ioerr.ne.0) then
            write(*,'(a)') 'Error reading third coordinate'
            stop
          endif   
          iarg = iclarg(4,arg) 
c         see if the 4th token is a velocity or output type
          read(arg,*,iostat=ioerr) inval(4)
          if( ioerr.eq.0 ) then
c           argument is numeric, must be velocity
            vel_flag = .true. 
            iarg = iclarg(5,arg) 
            read(arg,*,iostat=ioerr) inval(5)
            iarg = iclarg(6,arg)
            read(arg,*,iostat=ioerr) inval(6)    
            if( ioerr.ne.0) then
             write(*,'(a)') 'Error reading velocities from command line'
             stop  
            endif
            iarg = iclarg(7,arg)   
            if( iarg.gt.0 ) then 
              outtype = arg(1:3)
            else
              write(*,'(a)') 'Missing output type in command-line'
              stop
            endif 
          else
c           argument is non-numeric, must be output type
            vel_flag = .false.
            outtype = arg 
          endif   

        else                   
c         argument is name of input file   
          infile_flag = .true.
          infile = arg        
          iarg = iclarg(2,arg)
          outfile = arg       
          iarg = iclarg(3,arg)
          outtype = arg(1:3)
        endif               
        

c  Open the files
      if( infile_flag ) then
        open(unit=luin,file=infile,status='old',iostat=ioerr)
        if( ioerr.ne.0 ) then
          write(*,'(a,i4)') 'Error opening input file ',ioerr
          stop
        endif  
        open(unit=luout,file=outfile,status='unknown',iostat=ioerr)
        if( ioerr.ne.0 ) then
          write(*,'(a,i4)') 'Error opening output file ',ioerr
          stop
        endif   
      endif
                                  

c  Determine the input file (or command-line coordinate) type
c  and position the file for reading values

      if( .not.infile_flag ) then              
        tstrad = dsqrt(inval(1)**2+inval(2)**2+inval(3)**2)
        if( tstrad.gt.6.d6 ) then
           intype = 'XYZ'
        else
           intype = 'GEO'
        endif
      else
        call get_filetype(luin,intype,lfile_epoch)
        date = lfile_epoch
      endif  
      if( infile_flag ) write(*,'(a,a3)') 'Input file type is ',intype

c  Write the header of the output file according to type   
       
      if( infile_flag ) then
        if( outtype.eq.'LFI' ) then   
          write(luout,'(3a)') 'L-file from ',infile(1:nblen(infile))
     .                   ,' by CONVERTC'
        elseif( outtype.eq.'APR' ) then   
          write(luout,'(3a)') '* Apr file from ',infile(1:nblen(infile))
     .                   , ' by CONVERTC'         
          if( date.eq.0.d0 ) then 
            write(*,'(a)') 
     .       'WARNING: No date available for apr-file coordinates'
          endif
        elseif( outtype.eq.'VEL' ) then 
          write(luout,'(3a)') 
     .       '*SUMMARY VELOCITY ESTIMATES converted from '
     ,         ,infile(1:nblen(infile)),' by CONVERTC'
          write(luout,'(3a)') '* Long.     Lat.        E & N Rate '
     .     ,'     E & N Adj.      E & N +-   RHO        H Rate  '
     .     ,' H adj.    +-  SITE'
          write(luout,'(2a)') '* (deg)    (deg)          (mm/yr)'
     .    ,'       (mm/yr)       (mm/yr)                 (mm/yr)'
        elseif( outtype(1:2).eq.'GE' ) then   
           write(luout,'(3a)') 'Geodetic coordinates from '
     .       ,infile(1:nblen(infile)),' by CONVERTC'
        elseif (outtype.eq.'KML') then   
          call write_kmlf(luout,1,site,x)      
        endif
      endif

       
 
c  Now loop over the sites, reading the input file and writing the output file   
                                
      finished = .false.
      do while ( .not.finished )    
                                            
c       this false so far only in old-style l-file case
        values_ok = .true.    
      
        if( intype.eq.'XYZ' ) then   
          if( infile_flag ) then
            write(*,'(a)') 'XYZ not allowed as input file (use APR)'
            stop   
          else
            x(1) = inval(1)
            x(2) = inval(2) 
            x(3) = inval(3)
            if( vel_flag ) then
              xdot(1) = inval(4)
              xdot(2) = inval(5)
              xdot(3) = inval(6)
            endif
          endif   
        
        elseif( intype.eq.'GEO' ) then
          if( infile_flag ) then  
            read(luin,'(a)',iostat=ioerr) line
            if( ioerr.eq.-1 ) then 
              finished = .true.
            elseif( ioerr.ne.0 ) then
              write(*,'(a,i4)') 'Error reading GEO or GEP line ',ioerr
              stop
            else
              read(line,'(2f11.5,f7.1,2x,a4)',iostat=ioerr) 
     .           (inval(i),i=1,3),site(1:4)
              if(ioerr.eq.-1) then
                finished = .true.
              elseif(ioerr.ne.0) then
               write(*,'(a,i4)') 'Error decoding GEO or GEP line ',ioerr
               write(*,'(a)') line 
               stop 
              endif
            endif  
          endif
          geod_pos(1) = (90.d0 - inval(2))*pi/180.d0 
          geod_pos(2) = inval(1)*pi/180.d0 
          geod_pos(3) = inval(3)     
          call GEOD_to_XYZ(geod_pos,x)  
          if( vel_flag ) then 
            neudot(1) = inval(5)
            neudot(2) = inval(4)
            neudot(3) = inval(6)
            call XYZ_to_GEOD(rot,x,loc_coord)
            call transp(rot,rottr,3,3)
            call mmply(rottr,neudot,xdot,3,3,1) 
          endif  
                                   
        elseif( intype.eq.'APR' ) then 
          comment =.true.
          do while( comment )
            read(luin,'(a)',iostat=ioerr) line  
            if( ioerr.eq.-1 ) then
              finished = .true.   
              comment = .false.
            elseif( ioerr.eq.0 ) then
              if( line(1:1).eq.' ' ) then
                comment = .false.
                read(line,'(1x,a8)',iostat=ioerr) site
                read(line(10:nblen(line)),*,iostat=ioerr) x,xdot,date   
cd                print *,'ioerr read site x xdot date '
cd     .                 , ioerr,site,x,xdot,date
                if( ioerr.ne.0 )  then  
                  write(*,'(a,i3)') 'Error reading x xdot date ',ioerr
                  write(*,'(a)') 'LINE: ',line
                  comment = .true.
                endif
              endif
            else
              write(*,'(a,i4)') 'Error reading apr-file line ',ioerr
              stop
            endif                                   
          enddo

        elseif( intype.eq.'VEL' ) then 
           call read_velf( luin,site,x,xdot,finished)
          
        elseif( intype.eq.'GLI' ) then
          call read_glistf( luin,site,x,finished)      

        elseif( intype.eq.'LFI' ) then
          call read_lfile( luin,site,x,date,values_ok,finished)  
          if( date.eq.0.d0 ) date = lfile_epoch

        else
          write(*,'(a,a3)') 
     .       'Input file coordinate type not recognized ',intype
          stop
        endif
        
        if( .not.finished .and. values_ok ) then

          if( outtype.eq.'XYZ' ) then
            if( .not.infile_flag ) then
              if( vel_flag ) then 
                write(*,'(3f15.5,3f10.5)') x,xdot
              else
                write(*,'(3f15.5)') x    
              endif                             
              finished = .true.   
            else
              write(*,'(a)') 
     .            'XYZ not allowed as output file type (use APR)'
              stop
            endif

          elseif( outtype.eq.'APR' ) then 
            write(luout,'(1x,a8,3f15.5,3f10.5,f9.3)') site,x,xdot,date
  
          elseif( outtype.eq.'VEL' ) then
            call write_velf(luout,site,x,xdot)

          elseif( outtype.eq.'LFI' ) then
            call write_lfile(luout,site,x) 
        
          elseif( outtype(1:2).eq.'GE' ) then
            call XYZ_to_GEOD(rot,x,geod_pos)  
            latd = 90.d0 - geod_pos(1)*180.d0/pi 
            lond = geod_pos(2)*180.d0/pi   
            if( outtype.eq.'GEO'.and.lond.gt.180.d0 ) then
              lond = lond -360.d0  
            elseif( outtype.eq.'GEP'.and.lond.lt.0.d0) then
              lond = lond + 360.d0
            endif 
            ht = geod_pos(3)                   
cd            print *,'converted ',lond,latd,ht  
cd            print *,'infile_lfag vel_flag ',infile_flag,vel_flag
            if( .not.infile_flag ) then  
c             RWK 110701: Note that both pos and vel are now lon lat up, not NEU
              if( vel_flag ) then
                write(*,'(2f11.5,f7.1,3f10.5)') 
     .            lond,latd,ht,neudot(2),neudot(1),neudot(3)  
              else   
                write(*,'(2f11.5,f7.1)') lond,latd,ht     
              endif
              finished = .true.
            else
              write(luout,'(2f11.5,f7.1,2x,a4)') lond,latd,ht,site(1:4)
            endif

          elseif( outtype.eq.'KML' ) then
            call write_kmlf(luout,2,site,x) 
        
          else
            write(*,'(a,a3)') 
     .         'Output file coordinate type not recognized ',outtype
            stop
          endif
        endif
      enddo
c Write the kml footer
      if( outtype.eq.'KML' ) then
            call write_kmlf(luout,3,site,x)   
      endif   

      if( infile_flag ) then
        write(*,'(4a)') 'Wrote output file ',outfile(1:nblen(outfile))
     .   ,' of type ',outtype 
      endif
      stop
      end

