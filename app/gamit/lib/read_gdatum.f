      Subroutine READ_GDATUM( iud,datum,semi,finv,shft,shftdot
     .                      , rot,rotdot,scale,scaledot,epoch )

c     Read the datum parameters form gdetic.dat

c     R. King 27 Sept 2002
         
c     Input
c       iud     unit number (integer)
c       datum   5-character datum name; if blank, assume WGS84

c     Output (all real*8)
c       semi        semi-major axis  (m)
c       finv       inverse flattening
c       shft(3)    translation, x y z (m) 
c       shftdot(3) translation rate  xdot ydot zdot ( m/yr)
c       rot(3)     rotation (radians) 
c       rotdot(3)  rotation rate (radians/yr)
c       scale      scale
c       scaledot   scale rate (/yr) 
c       epoch      epoch for rates (decimal years)
                                                                      
      character*3 gformat
      character*5 datum,adatum,fdatum 
      character*10 lowerc
      character*256 line,prog_name,message

      integer*4 iud,ioerr,rcpar,len,i

      real*8 semi,finv,shft(3),shftdot(3),rot(3),rotdot(3)
     .     , scale,scaledot,epoch,pi 

      logical eof,found                 

c       Get the constants for conversion

      pi = 4.0d0*datan(1.0d0)    
        
c       Get the name of the calling program for report_stat calls

      len = rcpar(0,prog_name)

c        If datum blank (geocentric), set to WGS84
                  
      adatum = datum
      if( adatum.eq.'    ' ) adatum = 'WGS84'       

c       Open the file and see if old or new-style format
            
      open(unit=iud,file=lowerc('gdetic.dat'),status='old',iostat=ioerr)
      if(ioerr.ne.0 ) call report_stat('FATAL',prog_name,'read_gdatum'
     .                   ,'gdetic.dat','Error opening datum file',ioerr)
      read(iud,'(a)',iostat=ioerr) line
      if( ioerr.ne.0) call report_stat('FATAL',prog_name,'read_gdatum'
     .              ,' ','Error reading first line of gdetic.dat',ioerr)
c     assume old-style has 'Geodetic' in first 8 columns of header, but
c     allow for reconsideration if no rotation parameters found
      gformat = 'new'
      if( line(1:4).eq.'Geod'.or.line(2:5).eq.'Geod' ) gformat = 'old'  
      
c       Read an old-style file (2 header lines, then datum in column 1, no rotation)
       
      if( gformat.eq.'old' ) then 
c       second header line
        read(iud,'(a)',iostat=ioerr)
        if( ioerr.ne.0) call report_stat('FATAL',prog_name,'read_gdatum'
     .   ,' ','Error reading second line of old-style gdetic.dat',ioerr)
        eof = .false.  
        found = .false.
        do while( .not.found .and. .not.eof ) 
          read(iud,'(a)',iostat=ioerr) line  
          if( ioerr.eq.0 ) then
            read(line(1:5),'(a5)',iostat=ioerr) fdatum  
            if( ioerr.ne.0) call report_stat('FATAL',prog_name
     .        ,'read_gdatum',' ','Error reading datum name (old)',ioerr)
            if( fdatum.eq.adatum ) then 
              found = .true.
              read(line,'(9x,f11.3,1x,f13.9,3f9.3)',iostat=ioerr) 
     .              semi,finv,shft  
              if( ioerr.ne.0) call report_stat('FATAL',prog_name
     .      ,'read_gdatum',' ','Error reading datum values (old)',ioerr) 
              scale = 0.d0
              scaledot = 0.d0
              do i=1,3        
                shftdot(i) = 0.d0
                rot(i) = 0.d0    
                rotdot(i) = 0.d0 
              enddo
            endif
          elseif( ioerr.eq.-1 ) then
            eof = .true. 
          else
            call report_stat('FATAL',prog_name,'read_gdatum',' '
     .              ,'Error reading line in old-style gdetic.dat',ioerr)
          endif
        enddo 

c         Read a new-style file (assumes at least the first line is a comment)

      else 
                                          
        eof = .false.
        found = .false.
        do while(.not.found .and. .not.eof)
          read(iud,'(a)',iostat=ioerr)  line    
          if( ioerr.eq.0 .and. line(1:1).eq.' ') then 
            read(line(2:6),'(a5)',iostat=ioerr) fdatum
            if( ioerr.ne.0) call report_stat('FATAL',prog_name
     .     ,'read_gdatum','gdetic.dat','Error reading datum name',ioerr)
            if( fdatum.eq.adatum ) then 
              found = .true.
              read(line(7:256),*,iostat=ioerr) semi,finv
     .               ,shft,shftdot,rot,rotdot,scale,scaledot,epoch  
              if( ioerr.ne.0) call report_stat('FATAL',prog_name
     .             ,'read_gdatum',' '
     .             ,'Error reading datum values--old-style file?',ioerr) 
              scale = scale*1.d-9
              scaledot = scaledot*1.d-9
              do i=1,3        
c               units of table are mas
                rot(i) = rot(i)*pi/180.d0/3600.d0/1.d3
                rotdot(i) = rotdot(i)*pi/180.d0/3600.d0/1.d3
              enddo
            endif
          elseif( ioerr.eq.-1 ) then
             eof = .true.      
          elseif( ioerr.ne.0 ) then
            call report_stat('FATAL',prog_name,'read_gdatum',' '
     .              ,'Error reading line in gdetic.dat',ioerr)  
c         else comment, so continue
          endif
        enddo 
          
      endif   
      
      if( .not.found ) then
        write(message,'(a,a5,a)') 'Datum (',datum
     .     ,') not found in gdetic.dat'
        call report_stat('FATAL',prog_name,'read_gdaum',' ',message,0)
      endif

      return
      end


     






                 

