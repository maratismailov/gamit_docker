      Program DCBTAB

c     Update the GAMIT file dcb.dat with new monthly values from AIUB

c     Command-line arguents:
c           dcbtab [old-file] [monthly file] [new-file]
c     e.g.  dcbtab dcb.dat.gps  p1c1.gpsdcb dcb.dat.gps.new

c     R. King 1 October 2004
c     Last modified, R. King 3 August 2015

      implicit none

      character*20 infile,monfile,outfile
      character*80 line 
      character*256 message

      integer*4 ioerr,idoy,year,month,doy,iprn,iarg,iclarg
      integer*4 iday,imon,julian,julday,nblen                

      real*4 dcb,rms

      logical eof
      
c  Read the command-line to get the input file name

      iarg = iclarg(1,infile)
      if( iarg.le.0 ) then  
         write(*,'(a)') 'Missing arguments for dcbtab '
         write(*,'(a)') '  dcbtab [infile] [COD-file] '
         stop
      endif 
      monfile = ' '   
      iarg = iclarg(2,monfile)
      if( iarg.le.0 ) then
         write(*,'(a )') 'Translating a v1 file to v2'
      endif       
        
c  Open the original file

      open(unit=1,file=infile,status='old',iostat=ioerr)
      if( ioerr.ne.0 )  call report_stat('FATAL','DCBTAB','Main',infile,
     .                    'Cannot open input file ',ioerr)

c  Open the new monthly file

      open(unit=2,file=monfile,status='old',iostat=ioerr)
      if( ioerr.ne.0 )  call report_stat('FATAL','DCBTAB','Main'
     .    ,monfile,'Cannot open the monthly AIUB file',ioerr)
                                                                 
c   Open the new updated file
      
      outfile = infile(1:nblen(infile))//'.new'
      print *,'infile outfile ',infile,outfile
      open(unit=3,file=outfile,status='unknown',iostat=ioerr)
      if( ioerr.ne.0 )  call report_stat('FATAL','DCBTAB','Main'
     .      ,outfile,'Failure in opening the new dcb.dat file',ioerr)
                                                                      
c   Copy the original file into the new one
                         
      eof = .false.
      do while (.not.eof) 
        read(1,'(a80)',iostat=ioerr) line   
        if( ioerr.eq.-1.or.line(1:3).eq.'   '  ) then
          eof = .true.
        elseif( ioerr.ne.0 ) then
          call report_stat('FATAL','DCBTAB','Main','dcb.dat',
     .           'Error reading original DCB file ',ioerr)
        else
           write(3,'(a80)') line
        endif    
      enddo

c   Read the new monthly values and add them to the end
c   of the first the head/date line   

c       this line is a mess: there are at least four different formats for the date
c       the one constant seems to be the time-date-stamp for the file, which, unfortunately
c       has the date in letters not numbers.  The three most recent formats seem to be
c       either 
c (V1)   CODE'S 30-DAY GPS P1-C1 DCB SOLUTION, ENDING D274, 2004          12-OCT-04 10:57
c       or 
c (V2)   CODE'S MONTHLY GPS P1-C1 DCB SOLUTION, YEAR-MONTH 04-12          04-JAN-05 10:49  
c       or
c (V3)   CODE'S 30-DAY GPS P1-C1 DCB SOLUTION, ENDING DAY 226, YEAR 2     17-AUG-12 11:31
c       or
c (V4)   CODE'S MONTHLY GPS P1-C1 DCB SOLUTION, YEAR 2012, MONTH 01       07-SEP-12 11:50      
c       or
c (V5)   CODE'S 30-DAY GPS P1-C1 DCB SOLUTION, ENDING DAY 281, 2012       11-OCT-12 11:13

   
c       Use the presence of 'DAY' in columns 46-48 or 'D' in column 46 to detect the 
c       third or first format.  If the dates are unreasonable, stop.

      read(2,'(a80)',iostat=ioerr) line   
      if( ioerr.ne.0 )  call report_stat('FATAL','DCBTAB','Main'
     .      ,monfile,'Error reading line of AIUB DCB file ',ioerr)   
      if( line(46:48).eq.'DAY' .and. line(55:58).ne.'YEAR' ) then
c       format is 'ENDING DAY 281, 2012', Version 5
        read(line,'(49x,i3,2x,i4)',iostat=ioerr) doy,year
        if( ioerr.ne.0 ) 
     .       call report_stat('FATAL','DCBTAB','Main',monfile
     .      ,'Error reading doy, year from V5 AIUB  DCB file',ioerr) 
c       date is end of solution, subtract 15 days to get midpoint
        call monday( doy, imon, iday, year )
        julian = julday( imon, iday, year )
        julian = julian - 15
c       doy = doy - 15  
        call dayjul( julian, year, doy )
        write(*,'(a,a)') 'AIUB monthly file header line : ',line
        write(*,'(a,i4,1x,i3)') 'Read year, doy : ',year,doy
      elseif( line(51:55).eq.'MONTH') then
c       format is 'YEAR yyyy, MONTH mm' Version 4
        read(line,'(44x,i4,8x,i2)',iostat=ioerr) year,month    
        if( ioerr.ne.0 ) 
     .     call report_stat('FATAL','DCBTAB','Main',monfile
     .       ,'Error reading year,month from V4 AIUB  DCB file',ioerr)
        doy = idoy( year,month,1)
c       add 15 days to get the middle of the month
        call monday( doy, imon, iday, year )
        julian = julday( imon, iday, year )
        julian = julian +15
c       doy = doy + 15 
        call dayjul( julian, year, doy )
      elseif( line(55:58).eq.'YEAR' ) then
c       format is 'DAY 226, YEAR 2' Version 3
        read(line,'(49x,i3,7x,i1)',iostat=ioerr) doy,year   
        if( ioerr.ne.0 ) 
     .     call report_stat('FATAL','DCBTAB','Main',monfile
     .    ,'Error reading doy, year from V3 AIUB  DCB file',ioerr) 
c       check for reasonableness
        if( year.lt.0.or.year.gt.9 ) then
          write(message,'(a,i4)') 'Unreasonable year ',year
          call report_stat('FATAL','DCBTAB','Main',' ',message,0) 
        elseif( doy.lt.1.or.doy.gt.366 ) then
          write(message,'(a,i3)') 'Unreasonable doy ',doy
          call report_stat('FATAL','DCBTAB','Main',' ',message,0) 
        else
          year = year + 2010  
c         date is end of solution, subtract 15 days to get midpoint
          call monday( doy, imon, iday, year )
          julian = julday( imon, iday, year )
          julian = julian - 15
c         doy = doy - 15  
          call dayjul( julian, year, doy )
          write(*,'(a,a)') 'AIUB monthly file header line : ',line
          write(*,'(a,i4,1x,i3)') 'Read year, doy : ',year,doy
        endif
      elseif ( line(46:46).eq.'D' ) then  
c       format is 'D274, 2004' Version 1 
        read(line,'(46x,i3,2x,i4)',iostat=ioerr) doy,year   
        if( ioerr.ne.0 )  then 
          write(*,'(a80)') line
          call report_stat('FATAL','DCBTAB','Main',monfile
     .      ,'Error reading doy, year from V1 AIUB  DCB file',ioerr) 
c       check for reasonableness
        elseif( year.lt.2000.or.year.gt.2100 ) then
          write(message,'(a,i4)') 'Unreasonable year ',year
          call report_stat('FATAL','DCBTAB','Main',' ',message,0) 
        elseif( doy.lt.1.or.doy.gt.366 ) then
          write(message,'(a,i3)') 'Unreasonable doy ',doy
          call report_stat('FATAL','DCBTAB','Main',' ',message,0) 
        else
c         date is end of solution, subtract 15 days to get midpoint
          call monday( doy, imon, iday, year )
          julian = julday( imon, iday, year )
          julian = julian - 15
c         doy = doy - 15  
          call dayjul( julian, year, doy )
          write(*,'(a,a)') 'AIUB montly file header line : ',line
          write(*,'(a,i4,1x,i3)') 'Read year, doy : ',year,doy
        endif
      else             
c       Assume format is Version 2
        read(line,'(50x,i2,1x,i2)',iostat=ioerr)  year,month
        if( ioerr.ne.0 )  then 
          write(*,'(a80)') line
           call report_stat('FATAL','DCBTAB','Main',monfile
     .      ,'Error reading year, month from V2 AIUB  DCB file',ioerr) 
c       check for reasonableness
        elseif( year.lt.0.or.year.gt.50 ) then
          write(message,'(a,i2)') 'Unreasonable year ',year
          call report_stat('FATAL','DCBTAB','Main',' ',message,0) 
        elseif( month.lt.1.or.month.gt.12 ) then
          write(message,'(a,i2)') 'Unreasonable month ',month
          call report_stat('FATAL','DCBTAB','Main',' ',message,0) 
        else
c        date is beginning of month  
         write(*,'(a,a)') 'AIUB monthly file header line : ',line
         write(*,'(a,i4,1x,i3)') 'Read month, year: ',month,year
         call fix_y2k(year)   
         doy = idoy( year,month,1)
c        add 15 days to get the middle of the month
         call monday( doy, imon, iday, year )
         julian = julday( imon, iday, year )
         julian = julian + 15
c        doy = doy + 15 
         call dayjul( julian, year, doy )
         if( doy.gt.365 ) then 
c          don't worry about leap years in this file: solution is a 30d average
           doy = doy - 365    
           year = year + 1
         endif 
        endif           
      endif
c     write the time line for the full file   
      print *,'Writing year doy ',year,doy
      write(3,'(a,2i4)',iostat=ioerr) ' Epoch ',year,doy
c     skip 6 lines of header/comments
      read(2,'(/////)') 
c     now read and write the values for each SV
      eof = .false.
      do while (.not.eof) 
         read(2,'(1x,i2,24x,f9.3,3x,f9.3)',iostat=ioerr) iprn,dcb,rms
         if( ioerr.eq.-1 .or. iprn.eq.0 ) then
           eof = .true.
         elseif( ioerr.ne.0 ) then
           call report_stat('FATAL','DCBTAB','Main',monfile,
     .            'Error reading values from monthly DCB file ',ioerr)
         else
           write(3,'(1x,i2,2f9.3)') iprn,dcb,rms
         endif
      enddo
                    
      stop
      end  
                                



      
