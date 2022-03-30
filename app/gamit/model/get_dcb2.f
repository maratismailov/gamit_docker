      Subroutine get_dcb2( iudcb,jd,gnss,isvn,dcb )

c     Read a Version 2 dcb.dat file and return the DCB value for the
c     date requested.  Differs from the original (Version 1) get_dcb.f
c     in that this routine returns a scalar value for the PRN (and date)
c     requested  rather than the full array of values from the monthly
c     CODE solution.  The new format avoids the problem of mid-month
c     changes in PRN assignment

      implicit none

c     Input
c     -----   
c      iudcb  unit number for dcb.dat
c      jd     PEP Julian Day of observation 
c      gnss   GNSS code for requested PRN
c      isvn   SVN requested
                                         
      character*1 gnss
      integer*4 iudcb,jd,isvn

c     Output
c     ------                  
c       dcb   correction to C1 or P2' in ns
      real*8 dcb 
                        
c     Local variables
c     ---------------
                 
      character*1 gnss_tmp
      character*80 line   
      character*256 message
           
      integer*4 iyr,iday,time(5),start(5),stop(5)
     .        , isvn_tmp,iprn_tmp,ioerr,i

      real*8 rms 

      logical eoh,eof,found,warning
      data warning/.false./
            
c     Function
c     ---------
      integer*8 itimdif


c Put the requested time into an array for comparing with file entries
      call dayjul( jd,iyr,iday)
      time(1) = iyr
      time(2) = iday
      time(3) = 0
      time(4) = 0
      time(5) = 0                            
c     times read from file have no seconds
      start(5) = 0
      stop(5) = 0 

c Check the version number and then read but ignore the rest of the header
      read( iudcb,'(a)',iostat=ioerr ) line 
      if( line(11:21).ne.'Version 2.0' ) then 
         call report_stat('FATAL','MODEL','get_dcb2'
     .     ,'dcb.dat','Only Version 2.0 supported',0)
      endif
      eof =.false.
      do while( .not.eoh )
        read(iudcb,'(a)',iostat=ioerr ) line
        if( ioerr.ne.0 ) call report_stat('FATAL','MODEL','get_dcb'
     .     ,' ','Error reading header of V2 dcb.dat file',ioerr)
        if( line(1:4).eq.'*SYS' ) then
          eoh = .true.
        endif
      enddo                                                                              
                                 
      
c  Loop over the file file to find the entry that matches the time tag
      eof = .false.    
      found = .false.
      do while (.not.found .and. .not.eof )
        read(iudcb,'(a)',iostat=ioerr) line
        if( ioerr.eq.-1 ) then
          eof = .true.
        elseif( ioerr.ne.0 ) then
           call report_stat('FATAL','MODEL','get_dcb'
     .         ,' ','Error reading V2 dcb.dat file',ioerr)
        elseif (line(1:1).ne.' ') then
c         comment 
          continue 
        else         
          read(line,'(1x,a1,i6,i4,2(i6,i4,2i3),2f10.3)')
     .          gnss_tmp,isvn_tmp,iprn_tmp
     .        , (start(i),i=1,4),(stop(i),i=1,4),dcb,rms
          if ( gnss_tmp.eq.gnss .and. isvn_tmp.eq.isvn ) then
            if( itimdif(start,time).le.0 .and. 
     .          itimdif(stop,time).ge.0 ) then           
              found = .true. 
cd                print *,'found',found    
cd                print *,'isvn isvn_tmp ',isvn,isvn_tmp
cd                print *,'time,start,stop dcb rms '
cd     .                 , time,start,stop,dcb,rms
            endif
          endif
        endif
      enddo 
                                                                        
c  For GNSS other than GPS, we have no P1-C1 DCBS
      if( gnss.ne.'G' ) then
         dcb = 0.d0
      
c  For GPS stop if not found and warn if we are extrapolating the monthly values  
      else 
        if( .not.found ) then
          write(message,'(a,i3,a)') 'DCB value for SVN ',isvn
     .           ,' not found on dcb.dat'
          call report_stat('FATAL','MODEL','get_dcb2',' ',message,0)
        elseif( itimdif(time,start).gt.86400*31.and..not.warning) then
          write(message,'(a,i3,a)') 'Requested time for SVN ',isvn
     .        ,' > 31 days after the last entry: need to update dcb.dat'
          call report_stat('WARNING','MODEL','get_dcb2','dcb.dat'
     .          , message,0)
          warning = .true.
        endif
      endif
      
      return
      end

