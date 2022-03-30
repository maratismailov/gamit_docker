      Subroutine read_antex_head( luin,antex_vers,pcvtype,refant
     .                          , sinex_code,ncomments,comments )
      
      implicit none


      include '../includes/dimpar.h'
      
      integer*4 luin,ncomments,ioerr,len,rcpar
                   
      real*4 antex_vers

      character*1 pcvtype,svsys 
      character*10 sinex_code
      character*20 refant,label
      character*60 comments(maxtxt)
      character*256 line,prog_name,message

      logical eoh
  
      logical first_call  ! Used to stop repeated errors in case of
                          ! a bad comment line in the header.

      data first_call /.true./
      save first_call
 
      eoh = .false. 

c**   Get the calling module name for report_stat 

      len = rcpar(0,prog_name)
                     

c**   Read the header lines of the ANTEX file
              
c    --version line 
      read(luin,'(a)',iostat=ioerr) line
      if( ioerr.ne.0 ) call report_stat('FATAL',prog_name
     .  ,'lib/read_antex_head',' '
     .  ,'Error reading first record of ANTEX antmod.dat file',ioerr)
c      print *,'LINE',line
c      read(line,'(f8.1,12x,a1,39x,a20)',iostat=ioerr) 
      read(line,'(f8.1,12x,a1,39x,a20)',iostat=ioerr) 
     .     antex_vers,svsys,label  
      if( ioerr.ne.0 .or. label(1:5).ne.'ANTEX')    
     .   call report_stat('FATAL',prog_name,'lib/read_antex_head',' '
     .     ,'Error reading first line of ANTEX file,',ioerr)
      if( antex_vers.gt.1.4 ) 
     .   call report_stat('FATAL',prog_name,'lib/read_antex_head',' '
     .     ,'ANTEX version > 1.4',0)   
c     --type of table line
      read(luin,'(a1,19x,a20,20x,a20)',iostat=ioerr)pcvtype,refant,label
c      print *,'Line 2: pcvtype,refant,label ',pcvtype,refant,label
      if( ioerr.ne.0 .or. label(1:3).ne.'PCV' )  
     .   call report_stat('FATAL',prog_name,'lib/read_antex_head',' '
     .     ,'Error reading second line of ANTEX file,',ioerr)   
c    --SINEX code or comments (both optional) or END OF HEADER
      ncomments=0
      do while (.not.eoh ) 
        read(luin,'(a256)',iostat=ioerr)  line 
        if( line(61:70).eq.'SINEX CODE') then
          read(line,'(a10)',iostat=ioerr) sinex_code
          if( ioerr.ne.0 )  call report_stat('FATAL',prog_name
     .       ,'lib/read_antex_head',' '
     .     ,'Error reading SINEX code from ANTEX header,',ioerr)   
c           print *,'Line3 sinex_code ',sinex_code
        elseif( line(61:67).eq.'COMMENT' ) then
          ncomments = ncomments + 1
          if( ncomments.gt.maxtxt ) then
            write(message,'(a,i3,a)') 
     .         'Number of comments in ANTEX header (',ncomments
     .              ,') > MAXTXT in dimpar.h'
             call report_stat('FATAL',prog_name,'lib/read_antex_head'
     .              ,' ',message,ioerr)
          endif
          comments(ncomments) = line(1:60)
        elseif( line(61:73).eq.'END OF HEADER') then
          eoh = .true.
        else 
* MOD TAH 200209: Tell user more about what was just read (bad comment
*         is example of error
          if( first_call ) then 
             write(*,100) trim(line)
             len = index(line,'COMMENT')
 100         format('Bad Line in Antex header ',a)
             len = index(line,'COMMENT')
             if ( len.gt.0 ) then
                write(*,110) len
 110            format('COMMENT in line, column ',I3,'. Should be 61. ',
     .                 'Continuing to read header') 
             else 
                call report_stat('FATAL',prog_name,
     .              'lib/read_antex_head',' ' ,
     .               'EOF or unexpected line in ANTEX header',ioerr)
             endif
             first_call = .false.
           end if
        endif
      enddo 

      return
      end     


