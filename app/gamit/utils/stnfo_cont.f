      Program stnfo_cont

c     Reads a station.info file made from a sparse collection of RINEX files
c     and enforces continuity by making the stop times for entry n-1 the
c     same as the start time ofr entry n for each station, with the final
c     entry 9999 999 0 0 0 .    R. King 170920

c      Command-line:  stnfo_cont [filename] 
c
c        where [filename] is the input station.info file. The output
c        file will be [filename.NEW].

c     The program assumes that the first 61 columns are SITE, Station Name, 
c     Session Start, and Session Stop. The line length (entries plus comments)
c     can be up to 256 characters.  Only the site name, start time, and stop
c     time are decoded, the rest of the line copied intact.
                    
      implicit none

      integer*4 start(5),stop(5),startlast(5),stoplast(5)
     .        , ioerr,iarg,luin,luout
                   
      logical eoh,eof

c  Functions
       integer*4 iclarg,nblen,i
       character*4 site,sitelast  
       character*32 infile,outfile
       character*256 line,linelast

       luin = 1
       luout = 2    
               
c Get the run-string arguments and open the files

      iarg = iclarg(1,infile)
      if( iarg.gt.0 ) then
        open(unit=luin,file=infile,status='old',iostat=ioerr) 
        if( ioerr.ne.0 ) then
          write(*,'(a,a,i5)') 'Error opening input file ',infile,ioerr
          stop
        endif  
        write(*,'(a,a)') ' Opened input file ',infile
      endif     
      outfile = infile(1:nblen(infile))//'.NEW'
      open(unit=luout,file=outfile,status='unknown',iostat=ioerr)
      if( ioerr.ne.0 ) then
        write(*,'(a,a,i5)') 'Error opening output ',outfile,ioerr   
        stop
      endif  
      write(*,'(a,a)') ' Opened output file ',outfile

c Set the start time-of day always to 0 0 0 
c (stop time will be either 23 59 30 or 0 0 0 : see below)
   
      start(1) = 0
      start(2) = 0 
      start(3) = 0

                                                     
c Write the header and find the starting entry

      eoh = .false.
      do while (.not.eoh) 
        read(luin,'(a)',iostat=ioerr) linelast       
        if( ioerr.ne.0 ) then
           write(*,'(a,i4)') 'Error reading header lines, ioerr=',ioerr
           stop
        endif
        if(linelast(1:1).ne.' ' ) then
          write(luout,'(a)') linelast(1:nblen(linelast)) 
        else
          eoh = .true.
          sitelast = linelast(2:5)   
          read(linelast(26:61),'(2i4,3i3,4x,2i4,3i3)',iostat=ioerr) 
     .          (startlast(i),i=1,5),(stoplast(i),i=1,5)
           if(ioerr.ne.0 ) then
             write(*,'(a,i4)') 'Error decoding 1st site line, ioerr='
     .                   ,ioerr
             write(*,'(a)') linelast
             write(*,'(a4,1x,a)') sitelast,linelast(26:61)
             stop                   
           endif
        endif
      enddo
     
c Read the entries 
             
      eof = .false.
      do while( .not.eof ) 
        read(luin,'(a)',iostat=ioerr) line
        if(ioerr.eq.-1) then
          eof = .true.      
c         write out the last line of the file before exiting
          stoplast(1) = 9999  
          stoplast(2) = 999
          stoplast(3) = 0
          stoplast(4) = 0
          stoplast(5) = 0  
          write(linelast(45:61),'(2i4,3i3)') (stoplast(i),i=1,5)
          write(luout,'(a)') linelast(1:nblen(linelast))
        elseif( ioerr.ne.0 ) then
          write(*,'(a,i4)') 'Error reading site line, ioerr=',ioerr
          write(*,'(a)') line
          stop
        else
          if(line(1:1).eq.' ') then
            site = line(2:5)
            read(line(26:61),'(2i4,3i3,2x,2i4,3i3)',iostat=ioerr) 
     .          (start(i),i=1,5),(stop(i),i=1,5)
            if(ioerr.ne.0 ) then
              write(*,'(a,i4)') 'Error decoding site line, ioerr=',ioerr
              write(*,'(a)') line
              write(*,'(a4,1x,a)') site,line(26:61)
              stop                   
            else     
              if(site.eq.sitelast) then
                stoplast(1) = start(1)
                stoplast(2) = start(2) - 1  
                stoplast(3) = 23
                stoplast(4) = 59
                stoplast(5) = 30
              else
                stoplast(1) = 9999  
                stoplast(2) = 999
                stoplast(3) = 0
                stoplast(4) = 0
                stoplast(5) = 0   
              endif         
              write(linelast(45:61),'(2i4,3i3)') (stoplast(i),i=1,5)
              write(luout,'(a)') linelast(1:nblen(linelast))
              linelast = line
              sitelast = site 
            endif
          else      
c           skip comments
            continue     
          endif
        endif
      enddo  

      write(*,'(a)') 'Normal end of stnfo_cont'
      stop
      end


                           

