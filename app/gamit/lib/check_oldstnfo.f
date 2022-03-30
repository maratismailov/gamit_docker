      Subroutine check_oldstnfo( lun, old_stinf ) 

c     Check station.info to see if it's old-style in order to decide
c     write routines to call in gamit and kf.  This routine is temporary
c     for the transition.  --rwk 020308
c     Modified to issue fatal if old-style (no longer supported). rwk 080508

      implicit none

      logical old_stinf   

      integer* 4 lun,ioerr,len,rcpar,i 
               
      character*80  prog_name
      character*256 line,message
          
c     get the main program name calling this routine
      len = rcpar(0,prog_name)
         
      rewind( lun, iostat=ioerr )
      if( ioerr.ne.0 ) call report_stat('FATAL'
     .           ,prog_name,'lib/check_oldstnfo',' ' 
     .           ,'Error rewinding station.info',ioerr) 
      do i=1,3
          read( lun,'(a)',iostat=ioerr) line
           if( ioerr.ne.0 ) call report_stat('FATAL'
     .           ,prog_name,'lib/check_oldstnfo',' ' 
     .           ,'Cannot read line from station.info',ioerr) 
      enddo
      if( line(1:4).eq.'TRCK' .or. line(2:5).eq.'TRCK' ) then 
          old_stinf = .true. 
c         old-style -- stop the program  
          write(message,'(2a)') 'Old-style station.info incompatible'
     .         ,' with GAMIT 10.34, convert with conv_stnfo'
          call report_stat('FATAL',prog_name,'lib/rstnfo',' ',message,0)   
      endif    
      rewind lun 

      return
      end
