      Program MAKE_STNFOLIST

c     Using site name lists from a grep of sites.defaults ftprnx lines and 
c     and ls of ../rinex, create a single-column list of unique site names
c     for input to mstinf for its -l command of sites to use

c     R. King 091230
            
c     The program is command-line driven with the form:

c         make_stnfolist [sdlistfile] [rxlistfile] [outfile] 

c       where

c          sdlistfile   list from sites.defaults, single-column 4-char codes
c          rxlistfile   list from 'ls ../rinex', of form
c            ../rinex/albh2210.93o
c            ../rinex/albh2220.93o
c        ...
c           ../rinex/10rd2440.93o.gz
c           ../rinex/10rd2450.93o.gz
c        ...
c           ../rinex/alcc1970.93o.Z
c           ../rinex/alcc1980.93o.Z

    
c       One of sdlistfile and rxlistfile can be a blank enclosed in quotes.  
c       The file can also exist but be empty.

      implicit none
                
      character*4 site,sites(1000)
      character*80 line,sdfilef,rxfilef,outfilef

      integer*4 ioerr,nsite,nsiteold,nsitesd,iclarg,iarg,i
                                                        
      logical eof,dup


c       Read the run-string  and open the files
                            
      iarg = iclarg(1,sdfilef)
      if( iarg.gt.0 ) then
        open(unit=1,file=sdfilef,status='old',iostat=ioerr) 
        if( ioerr.ne.0 ) then
          write(*,'(a,a,i5)') 'Error opening sdfile in MAKE_STNFOLIST'
     .        ,sdfilef,ioerr
          stop
        endif  
      else                     
        write(*,'(/,2a)') 'Runstring: '
     .        ,' make_stnfolist [sdlistfile] [rxlistfile] [outfile]'
        stop
      endif 
      iarg = iclarg(2,rxfilef)
      if( iarg.gt.0 ) then 
        open(unit=2,file=rxfilef,status='old',iostat=ioerr) 
        if( ioerr.ne.0 ) then
          write(*,'(a,a,i5)') 'Error opening rxfile in MAKE_STNFOLIST '
     .        ,rxfilef,ioerr
          stop
        endif  
      else
        write(*,'(a)') 'Incomplete runstring for MAKE_STNFOLIST'
        stop
      endif     
      iarg = iclarg(3,outfilef)
      if( iarg.gt.0 ) then 
        open(unit=3,file=outfilef,status='unknown',iostat=ioerr) 
        if( ioerr.ne.0 ) then
          write(*,'(a,a,i5)') 'Error opening outfile in MAKE_STNFOLIST '
     .        ,outfilef,ioerr
          stop
        endif  
      else
        write(*,'(a)') 'Incomplete runstring for MAKE_STNFOLIST'
        stop
      endif     
          
c       Initialize the number of site entries 

      nsiteold = 0    
      nsite = 0  
    
c       Read the sites.defaults list and add any non-duplicate site names to the list
             
      eof = .false.
      do while (.not.eof)
        read(1,'(a4)',iostat=ioerr) site
        if(ioerr.eq.-1.or.site(1:1).eq.' ') then
          eof = .true.  
        elseif( ioerr.ne.0 ) then
          write(*,'(a,i5)') 'Error reading sdlistfile in MAKE_STNFO '
     .       ,ioerr
           stop
        else    
          if( nsiteold.eq.0 ) then
            sites(1) = site
            nsiteold = 1
          else
            dup = .false.
            do i=1,nsiteold
              if( site.eq.sites(i) ) dup = .true.
            enddo
            if( .not.dup ) then
              nsite = nsiteold + 1
              sites(nsite) = site
              nsiteold = nsite
            endif
          endif
        endif
      enddo  
      nsitesd = nsite

      
c       Read the RINEX list  and add any non-duplicate site names to the list
             
      eof = .false.
      do while (.not.eof)
        read(2,'(a)',iostat=ioerr) line
        if(ioerr.eq.-1.or.line(1:1).eq.' ') then
          eof = .true.
        elseif( ioerr.ne.0 ) then
          write(*,'(a,i5)') 
     .       'Error reading rxlistfile line in MAKE_STNFO ',ioerr
           stop
        else        
c         decode the line: assumes form ../rinex/ssssddds.yyo, with possible extents
          read(line(10:13),'(a4)',iostat=ioerr) site
          if( ioerr.ne.0 ) then
            write(*,'(a,i5)') 
     .         'Error decoding RINEX file name in MAKE_STNFO',ioerr
            stop
          else
            if( nsiteold.eq.0 ) then
              sites(1) = site
              nsiteold = 1
            else
              dup = .false.
              do i=1,nsiteold
                if( site.eq.sites(i) ) dup = .true.
              enddo
              if( .not.dup ) then
                nsite = nsiteold + 1
                sites(nsite) = site
                nsiteold = nsite
              endif
            endif
          endif
        endif
      enddo   

c      Possibly add an alphabetical sort here


c      Write out the new file
                     
      write(3,'(a)') '** Site list file for mstinf '
      write(3,'(a)') '**   Entries from sites.defaults:'
      do i=1,nsitesd
       write(3,'(1x,a4)') sites(i)
      enddo
      write(3,'(a)') '**   Entries from rinex directory:'
      do i=nsitesd + 1, nsite
        write(3,'(1x,a4)') sites(i)
      enddo                               

      stop
      end




            
          
            

    
