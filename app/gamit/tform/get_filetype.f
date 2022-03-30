      Subroutine get_filetype(luin,intype,lfile_epoch)

c     Read the input file for convertc and determine whether it is
c     a GAMIT l-file or a GLOBK apr, vel (or org), or glist file.

c     R. King 28 February 2008

      implicit none

      integer*4 luin,ioerr
      
      real*8 rad,xyzpos(3),lfile_epoch
      character*3 intype     
      character*256 line

      logical finished

      intype = '   '

c  See if an org or vel file
                            
      finished = .false.
      do while(.not.finished)
        read(luin,'(a)',iostat=ioerr) line   
        if( ioerr.eq.0 ) then 
          if( line(2:27).eq.'SUMMARY VELOCITY ESTIMATES' .or.
     .          line(3:8) .eq.'SYSTEM') then
            intype = 'VEL'      
c           skip the two header lines
            read(luin,'(a)') line
            read(luin,'(a)') line
            finished = .true. 
            return
          endif   
        elseif( ioerr.eq.-1 ) then
          finished = .true.
          rewind(luin) 
        else
          write(*,'(a,i4)') 'Error reading input file ',ioerr
        endif
      enddo

c  See if a glist file   
        
      finished = .false.
      do while(.not.finished)
        read(luin,'(a)',iostat=ioerr) line  
        if( ioerr.eq.0 ) then 
          if( line(19:41).eq.'position and occurences') then 
            intype = 'GLI'
c           skip the two header lines
            read(luin,'(a)') line 
            read(luin,'(a)') line   
            finished = .true. 
            return
          endif   
        elseif( ioerr.eq.-1 ) then
          finished = .true.
          rewind(luin)                       
        else
          write(*,'(a,i4)') 'Error reading input file ',ioerr
        endif
      enddo

c  See if a GAMIT old-style  L-file
                        
      finished = .false.
c     look for a date on line 1
      read(luin,'(a)',iostat=ioerr) line   
      if( line(1:5).eq.'Epoch' ) then
        read(line(7:15),'(f9.4)',iostat=ioerr) lfile_epoch
        if( ioerr.ne.0 ) then
          write(*,'(a,i4)') 
     .   'GET_FILETYPE: Error decoding epoch on old-style l-file ',ioerr
          write(*,'(a)') 'Set epoch = 0.0'
          lfile_epoch = 0.0
        endif  
        intype = 'LFI' 
        finished = .true.  
        return 
      else 
c       if 'Epoch not on first line, could still be an l-file, look for N or S
c       or blank in column 18, W or E or blank in col 34, and a valid radius
        do while(.not.finished)
          read(luin,'(a)',iostat=ioerr) line   
          if( ioerr.eq.0 ) then   
            if( ( line(18:18).eq.'N'.or.line(18:18).eq.'S'
     .        .or.line(18:18).eq.' ') .and.
     .        ( line(34:34).eq.'E'.or.line(34:34).eq.'W'
     .         .or.line(34:34).eq.' ') ) then
              read(line(51:62),'(f12.4)',iostat=ioerr) rad
              if(ioerr.eq.0.and.rad.gt.6.3d6.and.rad.lt.6.4d6) then
                intype = 'LFI'
                rewind(luin)  
                finished = .true. 
                return
              endif   
            endif
          elseif( ioerr.eq.-1 ) then
            finished = .true.
          else
            write(*,'(a,i4)') 'Error reading input file ',ioerr
          endif
        enddo
        rewind(luin) 
      endif

c  See if a CONVERTC geodetic coordinates file
         
      read(luin,'(a)',iostat=ioerr) line 
      if( ioerr.eq.0 ) then  
        if( line(1:8).eq. 'Geodetic' ) then  
          intype = 'GEO'
          return
        endif    
      else
       write(*,'(a,i4)') 'Error reading 1st line of input file ',ioerr
      endif


c  See if a GLOBK apr file

      finished = .false.
      do while(.not.finished)
        read(luin,'(a)',iostat=ioerr) line 
        if( ioerr.eq.0 ) then  
          if( line(1:1).eq. ' ' ) then  
            read(line(10:62),*,iostat=ioerr) xyzpos
            if( ioerr.ne.0 ) then
              write(*,'(a,i4)') 
     .            'GET_FILETYPE: Error decoding apr values ',ioerr
              stop
            else
              rad = dsqrt(xyzpos(1)**2+xyzpos(2)**2+xyzpos(3)**2)
              if(rad.gt.6.3d6.and.rad.lt.6.4d6) then
                intype = 'APRL'
                finished = .true. 
                rewind(luin)
                return
              endif    
            endif  
          endif
        elseif( ioerr.eq.-1 ) then
          finished = .true.
          rewind(luin)
        else
          write(*,'(a,i4)') 'Error reading input file ',ioerr
        endif
      enddo      

      if( intype.eq.'   ') then
         write(*,'(a)') 'Cannot identify input file type'
         stop
      endif

      return
      end
