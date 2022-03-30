c Convert a GMT velocity file to GLOBK format

      implicit none

      integer ioerr,iugmt,iuvel,item,indx
      integer*2 idum2(8)

      real*4 lat,long,evel,nvel,hvel,evadj,nvadj,hvadj,evsig,nvsig,hvsig
     .     , rho,sig_scale,unit_scale
                     
      character*2 velunit
      character*4 site4,site_ext
      character*8 site 
      character*16 string
      character*30 gmtfile,velfile
      character*256 line
              
      logical endfile
                          
c     initialize endfile
      endfile = .false.

c     initialize units and strings
      iugmt = 1
      iuvel = 2
      gmtfile = ' '
      velfile = ' '
      line = ' '
     
c     set the horizontal adjusts and height values 
      evadj = 0.
      nvadj = 0.
      hvel = 0.
      hvsig = 0.1
      hvadj = 0.

c     get and open the input files
 
      write(*,*) 'Enter name of the input GMT velocity file:'
      read(*,*) gmtfile
      open(unit=1,status='old',file=gmtfile,iostat=ioerr) 
      if( ioerr.ne.0 ) then
        write(*,'(a,a30)')  'Error opening input GMT file ',gmtfile
        stop
      endif
      write(*,'(a,a30)') 'Input GMT file: ',gmtfile

      write(*,*) 'Enter name of the output GLOBK vel file:'
      read(*,*) velfile
      write(*,*) 'Opening output file ',velfile
      open(unit=iuvel,status='unknown',file=velfile)
      if( ioerr.ne.0 ) then
        write(*,'(a,a30)')  'Error opening output vel file ',velfile
        stop
      endif
c     write the header lines
      write(iuvel,'(a)') ' SUMMARY VELOCITY ESTIMATES FROM GMT2VEL'
      write(iuvel,'(3a)') '  Long.     Lat.        E & N Rate '
     .   ,'     E & N Adj.      E & N +-   RHO        H Rate  '
     .   ,' H adj.    +-  SITE'
      write(iuvel,'(2a)') '  (deg)    (deg)          (mm/yr)'
     .  ,'       (mm/yr)       (mm/yr)                 (mm/yr)'

          
c     input the sigma scale factor and site name extent

      write(*,*) 'Enter the scale factor for sigmas '
      read(*,*) sig_scale
      write(*,*) 'Velocity sigmas being scaled by ',sig_scale  
      site_ext = ' ' 
      write(*,*) 'Enter the site-name extent (e.g. _GPS)'
      read(*,*) site_ext   
      velunit = '  ' 
      write(*,*) 'Enter units of velocities (mm, cm, or m) '
      read(*,*) velunit
      if( velunit.eq.'mm' ) then
         unit_scale = 1.   
         write(*,*) 'Input units are mm -- no conversion'
      elseif ( velunit.eq.'cm' ) then
         unit_scale = 10.  
         write(*,*) 'Input units are cm -- convert to mm'
      elseif ( velunit.eq.'m ' ) then
         unit_scale = 1000.   
         write(*,*) 'Input units are m -- convert to mm'
      else
         write(*,'(2a)') 'Invalid units: ',velunit
         stop 2
      endif


c     read and write each line of the files

      if( ioerr.ne.0 ) stop 1   
      do while (.not.endfile ) 
        read(iugmt,'(a)',iostat=ioerr) line 
        if( ioerr.eq.-1 ) then
          endfile = .true.
          print *,'Normal stop at EOF '
          stop
        elseif( ioerr.ne.0 ) then
          stop 3
        endif   
        indx = 1
        do item = 1,8  
           call read_line(line,indx,'CH',ioerr,idum2,string) 
           if( item.eq.1 ) read(string,*,iostat=ioerr) long
           if( item.eq.2 ) read(string,*,iostat=ioerr) lat
           if( item.eq.3 ) read(string,*,iostat=ioerr) evel
           if( item.eq.4 ) read(string,*,iostat=ioerr) nvel
           if( item.eq.5 ) read(string,*,iostat=ioerr) evsig
           if( item.eq.6 ) read(string,*,iostat=ioerr) nvsig
           if( item.eq.7 ) read(string,*,iostat=ioerr) rho
           if( item.eq.8 ) read(string,'(a4)',iostat=ioerr) site4 
           if( ioerr.ne.0 ) then
             write(*,'(a,i1,a)') 'Error reading item ',item,' from line'    
             stop 4
           endif
        enddo   
        site(1:4) = site4
        site(5:8) = site_ext
        evel = evel * unit_scale
        nvel = nvel * unit_scale
        evsig = evsig * unit_scale * sig_scale
        nvsig = nvsig * unit_scale * sig_scale
        write(iuvel,'(2f9.3,1x,4f8.2,2f8.2,f7.3,2x,f8.2,2f8.2,1x,a8)'
     .             ,iostat=ioerr)
     .              long,lat,evel,nvel,evadj,nvadj,evsig,nvsig,rho
     .             ,hvel,hvadj,hvsig,site
        if(ioerr.ne.0 ) then
           write(*,*) 'Error writing output line'
           stop 5 
        endif
      enddo

      stop
      end

