      SUBROUTINE GETGEO( IDMS,NUMSIT,X,SITNAM,IFILE)

C       Read in a set or file of geodetic coordinates.

      implicit none

      include '../includes/tform.h'

      logical fcheck,eof

      integer*4 idms,numsit,ifile,ioerr,latd,latm,lond,lonm,badline,i

      real*8 x(3,nsdim),semi,finv,shft(3),shftdot(3),rot(3),rotdot(3)
     .     , scale,scaledot,refepoch,epoch,ht,slat,slon,alat,along
     .     , geodrad,pi,convmas
                       
      character* 1 latflag,lonflag
      character* 5 datum   
      character* 8 siteid
      character*12 sname
      character*16 sitnam(nsdim) 
      character*60,lfmt1
      character*63 lfmt2
      character*256 message,line
                        
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
                



 
c       Constant

       pi = 4.d0*datan(1.d0)   
c      mas to radians
       convmas = 4.84813481d-9
                     
 
c       Open the datum file to display its contents, then close and call 
c         library routine read_gdatum to get the values.  

      if (fcheck('gdetic.dat')) then
        open(unit=idatum,file='gdetic.dat',iostat=ioerr,status= 'OLD')
        if (ioerr .ne. 0)  then
          call report_stat('FATAL','TFORM','getgeo',' '
     .                             ,'Error opening gdetic.dat',ioerr) 
        else  
           eof = .false.
           write(iscrn,'(/,a)') 'Datums available from gdetic.dat:'
           do while (.not.eof )           
             read(idatum,'(1x,a)',iostat=ioerr) line
             if( ioerr.eq.0 ) then 
               write(iscrn,'(a)') line
             elseif( ioerr.eq.-1 ) then
               eof = .true.
             else
               call report_stat('FATAL','TFORM','getgeo',' '
     .                          ,'Error reading gdetic.dat',ioerr)
             endif  
          enddo
        endif    
        close( idatum )
    5   write(iscrn,'(a)') 
     .    'Enter 5-character code for datum (case insensitive):' 
        read(iterm,'(a5)',iostat=ioerr) datum  
        call uppers(datum)
        if( ioerr.ne.0 ) go to 5
        call read_gdatum( idatum,datum,semi,finv,shft,shftdot
     .                  , rot,rotdot,scale,scaledot,refepoch )
        write(iscrn,'(a,a5,a,f11.3,a,f13.9)') 
     .       'Selected ',datum,' a=',semi,' 1/f=',finv  
        if( refepoch.eq.0.d0 ) then
          write(iscrn,'(/,a,3f9.3)') '  dX=',shft
        else     
          if( refepoch.lt.1900.d0 .or. refepoch.gt.2100.d0 ) then
             write(message,'(a,f20.1)') 'Invalid epoch from gdetic.dat:'
     .               ,refepoch
             call report_stat('FATAL','TFORM','getgeo',' ',message,0) 
          endif    
    6     write(iscrn,'(a)') 
     ,       'Enter epoch for transformation (decimal yr, e.g. 1997.0):'
          read(iterm,*,iostat=ioerr) epoch
          if( ioerr.ne.0 ) then
             write(iscrn,'(a,f20.0)') 'Error reading epoch: ',epoch
             goto 6                                   
          elseif( refepoch.lt.1900.d0 .or. refepoch.gt.2100.d0 ) then
              write(iscrn,'(a,f20.0)') 'Invalid epoch: ',epoch
              goto 6
          else      
           write(iscrn,'(/,2(a,3f9.3),a)') '  dX=',shft,'  dXdot='
     .            ,shftdot,' (m m/yr)'
           write(iscrn,'(2(a,3f9.3),a)') '  R =',(rot(i)/convmas,i=1,3)
     .             ,'  Rdot =',(rotdot(i)/convmas,i=1,3),' (mas mas/yr)'
           write(iscrn,'(2(a,d9.2),a)') '  Scale=',scale*1.d9
     .         ,'  Sdot=',scaledot*1.d9,' (ppb ppb/yr)'
           write(iscrn,'(a,f7.2)')    '  Ref Epoch=',refepoch   
           shft(1) = shft(1) + shftdot(1)*(epoch-refepoch)
           shft(2) = shft(2) + shftdot(2)*(epoch-refepoch)
           shft(3) = shft(3) + shftdot(3)*(epoch-refepoch)
           rot(1)  = rot(1)  + rotdot(1) *(epoch-refepoch) 
           rot(2)  = rot(2)  + rotdot(2) *(epoch-refepoch) 
           rot(3)  = rot(3)  + rotdot(3) *(epoch-refepoch) 
           scale = scale + scaledot *(epoch-refepoch)
           write(iscrn,'(/,a,f7.2)') 'Transformation at epoch ',epoch
           write(iscrn,'(a,3f12.3,a)') 'dX=',shft,' m'
           write(iscrn,'(a,1p3d12.4,a)') 'dR=',rot,' radians'
           write(iscrn,'(a,d12.4)') 'dS=',scale 
          endif                                
        endif 

      else    
        call report_stat('WARNING','TFORM','getgeo',' '
     .  ,'Missing gdetic.dat, use NAD83 with time-dependent corrections'
     .       ,0)
        semi=6378137.d0
        finv=298.257222101d0   
        scale = 0.d0 
        scaledot = 0.d0
        refepoch = 0.d0
        do i=1,3      
          shft(i) = 0.d0
          shftdot(i) = 0.d0
          rot(i) = 0.d0 
          rotdot(i) = 0.d0 
        enddo
      endif


c      Enter the coordinates          
                
      if( ifile.eq.0 ) then
        sitnam(1) = ' '
        numsit = 1
        if( idms.eq.1 ) then  
          write(iscrn,'(/,1x,a)') 'Enter Lat Long(E) Ht in deg, m' 
          read(iterm,*,iostat=ioerr) alat,along,ht   
          if( ioerr.ne.0) call report_stat('FATAL','TFORM','getgeo',' '
     .         ,'Error reading input coordinates',0)
        elseif( idms.eq.2) then   
          write(iscrn,'(1x,a)') 
     .      'Enter coordinates in L-file format (complete line)'
          read(iterm,lfmt1,iostat=ioerr) siteid(1:4),sname
     .        ,latflag,latd,latm,slat,lonflag,lond,lonm,slon,ht   
          siteid(5:8) = '_GPS'
          if( ioerr.ne.0) call report_stat('FATAL','TFORM','getgeo',' '
     .         ,'Error reading input coordinates',0)   
          call getsname(siteid,sname,sitnam(1),-1)  
c         get hemisphere     
          call dmsdeg(latflag,latd,latm,slat,alat)
          call dmsdeg(lonflag,lond,lonm,slon,along)
        else
           call report_stat('FATAL','TFORM','getgeo',' '
     .         ,'Unrecognized coordinate type',0)
        endif
        if( iprnt.gt.0 ) write(iprnt,'(/,a,/,1x,2f14.8,f9.3)') 
     .     'Input Geodetic Coordinates (lat lon ht(m) : '
     .      ,alat,along,ht 
        alat = alat*pi/180.d0
        along = along*pi/180.d0  
        call geoxyz( semi,finv,alat,along,ht,geodrad
     .              , X(1,1),X(2,1),X(3,1),1 )     
        if( refepoch.ne.0.d0.and.iprnt.gt.0 ) 
     .      write(iprnt,'(a,3f14.4)') 
     .      'Cartesian coordinates before applying datum shift: '
     .        ,(x(i,1),i=1,3)
           x(1,1) = -shft(1) + (1.d0-scale)*x(1,1) 
     .                 - rot(3)*x(2,1) + rot(2)*x(3,1) 
           x(2,1) = -shft(2) + (1.d0-scale)*x(2,1)
     .                 + rot(3)*x(1,1) - rot(1)*x(3,1)
           x(3,1) = -shft(3) + (1.d0-scale)*x(3,1)
     .                 - rot(2)*x(1,1) + rot(1)*x(2,1)

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
c              decimal degrees follows GLOBK conventions, non-blank first column 
c              is comment  blank second column implies a blank line, so skip 
               numsit = numsit + 1  
               if( numsit.gt.nsdim ) call report_stat('FATAL','TFORM'
     .                  ,'getgeo',' '
     .                  ,'Number of sites in file exceeds dimensions',0)
               read(line(1:9),'(1x,a8)',iostat=ioerr) siteid 
               call getsname(siteid,sname,sitnam(numsit),-2)
               if(ioerr.ne.0 ) call report_stat('FATAL','TFORM','getgeo'
     .            ,' ','Error reading site name from input file',ioerr)   
               read(line(10:256),*,iostat=ioerr) alat,along,ht 
               if(ioerr.ne.0) 
     .             call report_stat('FATAL','TFORM','getgeo',' '
     .            ,'Error reading decimal degrees and radius input file'
     .             , ioerr)   
               if( iprnt.gt.0 ) write(iprnt,'(a8,2f16.10,f14.4)')
     .                             sitnam(numsit)(1:8),alat,along,ht 
              endif
            elseif( idms.eq.2 ) then
c             L-file format      
              if( line(1:4).ne.'   ' ) then        
                read(line,lfmt1,iostat=ioerr)  siteid(1:4),sname
     .                ,latflag,latd,latm,slat,lonflag,lond,lonm,slon,ht  
                if(ioerr.eq.0 ) then 
                  numsit = numsit + 1    
                  call getsname(siteid,sname,sitnam(numsit),-1)    
                  if(iprnt.gt.0) write(iprnt,lfmt2) siteid,sname
     .              ,latflag,latd,latm,slat,lonflag,lond,lonm,slon,ht
                  call dmsdeg(latflag,latd,latm,slat,alat)
                  call dmsdeg(lonflag,lond,lonm,slon,along) 
                else   
c                 if read-error, assume that this is a comment line and keep going,
c                 but count the lines and report at the end
                  badline = badline + 1
                endif
              endif    
            else
               call report_stat('FATAL','TFORM','getgeo',' '
     .                      ,'Unrecognized idms value',0)   
            endif
          elseif( ioerr.eq.-1) then 
              eof = .true.   
          else
            call report_stat('FATAL','TFORM','getgeo',' '
     .                      ,'Error reading data line ',0)
          endif 
          if( numsit.le.nsdim.and..not.eof ) then
            alat = alat*pi/180.d0
            along = along*pi/180.d0     
            call geoxyz( semi,finv,alat,along,ht,geodrad
     .                 , x(1,numsit),X(2,numsit),X(3,numsit),1 )
            if( refepoch.ne.0.d0.and.iprnt.gt.0 ) 
     .        write(iprnt,'(a,3f14.4)') 
     .         'Cartesian coordinates before applying datum shift: '
     .        ,(x(i,numsit),i=1,3)
            x(1,numsit) = -shft(1) + (1.d0-scale)*x(1,numsit) 
     .                 - rot(3)*x(2,numsit) + rot(2)*x(3,numsit) 
            x(2,numsit) = -shft(2) + (1.d0-scale)*x(2,numsit)
     .                 + rot(3)*x(1,numsit) - rot(1)*x(3,numsit)
            x(3,numsit) = -shft(3) + (1.d0-scale)*x(3,numsit)
     .                 - rot(2)*x(1,numsit) + rot(1)*x(2,numsit) 
          elseif( numsit.gt.nsdim) then
            call report_stat('FATAL','TFORM','getgeo',' '
     .        ,'Number of sites in input file exceeds dimensions',0)
          endif
        enddo 
        if( badline.gt.0 ) then                                                        
          write(message,'(i4,a,a,i4)') badline
     .      ,' non-standard lines read from L-file, ok if comments;'
     .      ,' numsit = ',numsit
          call report_stat('WARNING','TFORM','getgeo',' '
     .      ,message,0)   
        endif
      endif
            
      return
      end
      
