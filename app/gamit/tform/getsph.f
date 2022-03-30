      Subroutine GETSPH( idms,numsit,x,sitnam,ifile )

c      Read in a set or file of spherical coordinates

c      IDMS = 1 decimal degrees
c      IDMS = 2 L-file format (deg/min/sec)


      implicit none

      include '../includes/tform.h'
        
      logical eof 

      integer*4 idms,numsit,ifile,latd,latm,lond,lonm,badline,ioerr

      real*8 x(3,nsdim),slat,slon,alat,along,rad   
                                      
      character*1 latflag,lonflag 
      character*8 siteid  
      character*12 sname
      character*16 sitnam(nsdim),blank    
      character*60,lfmt1
      character*63 lfmt2
      character*256 line,message 

      data blank/'                '/  

      data lfmt1/     
     .'(a4,1x,a12,a1,i2,1x,i2,1x,f8.5,1x,a1,i3,1x,i2,1x,f8.5,f13.4)'/  
      data lfmt2/
     .'(1x,a4,1x,a12,a1,i2,1x,i2,1x,f8.5,1x,a1,i3,1x,i2,1x,f8.5,f13.4)'/
          
c     'sitnam' stores a 4-character id followed by a 12-character descriptor.
c     For apr-type files, these are output as an 8-character name SITE_GPS;
c     for l-files, the output is 'SITE description '
                 
                           
c     Allowed formats
c     ---------------
c     if screen input (one site): free-format decimal degrees, no site name

c     if file input with decimal degrees:  8-character site-name beginning in column 2, 
c                                           free-format coordinates
c     if file input with deg/min/sec:  L-file format, 4-char id in column 1, name in columns 5-16
                


      if( ifile.eq.0 ) then   
c       read from the terminal
        sitnam(1)= blank
        numsit=1
        if( idms.eq.1 ) then
          write(iscrn,'(/,a)') 
     .        'Enter Lat Lon(E) in decimal degrees, Radius (in meters)' 
          read(iterm,*,iostat=ioerr) alat,along,rad    
          if(ioerr.ne.0) call report_stat('FATAL','TFORM','getsph',' '
     .          ,'Error reading input coordinates',ioerr)
        elseif( idms.eq.2 ) then
          write(iscrn,'(/,a)') 
     .           'Enter coordinates in L-file format (complete line)'
          read(iterm,lfmt1,iostat=ioerr) siteid(1:4),sname
     .              ,latflag,latd,latm,slat,lonflag,lond,lonm,slon,rad     
          if(ioerr.ne.0 ) call report_stat('FATAL','TFORM','getsph'
     .         ,' ','Error reading coordinates in L-file format',ioerr)  
          call getsname(siteid,sname,sitnam(1),-1) 
          call dmsdeg(latflag,latd,latm,slat,alat)   
          call dmsdeg(lonflag,lond,lonm,slon,along) 
        endif 
        call sphxyz(alat,along,rad,x(1,1))   


      else     
c       read from file 
        numsit = 0  
        eof = .false. 
        badline = 0
        do while (.not.eof)
          read(ifile,'(a)',iostat=ioerr) line        
          if( ioerr.eq.0 ) then  
            if(idms.eq.1 ) then     
              if(line(1:1).eq.' '.and.line(1:2).ne.'  ') then
c             decimal degrees follows GLOBK conventions, non-blank first column is comment 
c             blank second column implies a blank line, so skip 
              numsit = numsit + 1  
              if( numsit.gt.nsdim ) call report_stat('FATAL','TFORM'
     .                ,'getsph',' '
     .                ,'Number of sites in file exceeds dimensions',0)   
              read(line(1:9),'(1x,a8)',iostat=ioerr) siteid 
              call getsname(siteid,sname,sitnam(numsit),-2)
              if(ioerr.ne.0 ) call report_stat('FATAL','TFORM','getsph'
     .           ,' ','Error reading site name from input file',ioerr)   
              read(line(10:256),*,iostat=ioerr) alat,along,rad 
              if(ioerr.ne.0) 
     .            call report_stat('FATAL','TFORM','getsph',' '
     .            ,'Error reading decimal degrees and radius input file'
     .            , ioerr)   
              if( iprnt.gt.0 ) write(iprnt,'(a8,2f16.10,f14.4)')
     .                             sitnam(numsit)(1:8),alat,along,rad 
              call sphxyz(alat,along,rad,x(1,numsit))  
              endif
            elseif( idms.eq.2 ) then
c             L-file format      
              if( line(1:4).ne.'   ' ) then 
                read(line,lfmt1,iostat=ioerr)  siteid(1:4),sname
     .                ,latflag,latd,latm,slat,lonflag,lond,lonm,slon,rad 
                if(ioerr.eq.0 ) then 
                  numsit = numsit + 1    
                  call getsname(siteid,sname,sitnam(numsit),-1)
                  if(iprnt.gt.0) write(iprnt,lfmt2) siteid,sname
     .              ,latflag,latd,latm,slat,lonflag,lond,lonm,slon,rad
                  call dmsdeg(latflag,latd,latm,slat,alat)
                  call dmsdeg(lonflag,lond,lonm,slon,along) 
                  call sphxyz(alat,along,rad,x(1,numsit))  
                else   
c                 if read-error, assume that this is a comment line and keep going,
c                 but count the lines and report at the end  
                  badline = badline + 1
                endif
              endif    
            else
               call report_stat('FATAL','TFORM','getsph',' '
     .                      ,'Unrecognized idms value',0)   
            endif
          elseif( ioerr.eq.-1) then 
              eof = .true.   
          else
            call report_stat('FATAL','TFORM','getsph',' '
     .                      ,'Error reading data line ',0)
          endif  
        enddo 
        if( badline.gt.0 ) then                                                        
          write(message,'(i4,a,a,i4)') badline
     .      ,' non-standard lines read from L-file, ok if comments;'
     .      ,' numsit = ',numsit
          call report_stat('WARNING','TFORM','getsph',' '
     .      ,message,0)   
        endif
      endif
           
      return
      end

      Subroutine SPHXYZ(alat,along,rad,x)

      implicit none 
      
      real*8 x(3),pi,alat,along,rad
C
      PI= 4.D0*DATAN(1.D0)
C
      ALAT = ALAT*PI/180.D0
      ALONG = ALONG*PI/180.D0
      X(1)= RAD*DCOS(ALAT)*DCOS(ALONG)
      X(2)= RAD*DCOS(ALAT)*DSIN(ALONG)
      X(3)= RAD*DSIN(ALAT)
C     write(*,120) alat, along, rad, x
C120  format('Coords ',2F13.10,F14.4,3F15.4)
C
      RETURN
      END

