      SUBROUTINE ATMRED( IUW,ISCRN,iprnt,JD,T,FJDW1,FJDW2
     1                 , TEMP,PRES,WETVAR,H2OTYP )

c     Interpolate atmospheric values from table
c     R.King  June 1987 from old ATMRED by R. Abbot
c       which incorporated header reads (now in WHDRED)
c     Modified by S.Shimada  May 1990 to have compatibility with
c       SONAP ATMRED program on W-file format   
c     Modified by R. King Dec 2002 to clean up logic and avoid
c       problems at ends of tables

c    Input
c       IUW  :  Unit number of met (W-) file    
c       IUWD :  Unit number of zenith-delay print (WD-) file
c       ISCRN:  Unit number of screen output (for debug)
c       iprnt:  Unit number of print output
c       JD, T:  Epoch of observation (Julian day, seconds of day)
c       FJDW1, FJDW2 :  Start, stop times of table (JD.FRACT)
c       H2OTYP: C*1  Surface measurement of moisture is either relative
c                    humidity ('R') or dew point temperature ('D')

c    Output
c       PRESS:      Input pressure (sealevel or local) (millibars)
c       TEMP:   Site temperature ( degrees Celsius)
c       WETVAR: Surface measurement of moisture, either relative humidity
c               (fraction) or dew point temperature (degrees Celsius)

      implicit none

      character*1  h2otyp,upperc
      character*256 line,message
      integer*4  jd,julday,iscrn,ioerr,iprnt,iyr,idyoyr,month,iday
     .           , ihr,min,iuw

      real * 8  temp,pres,fjdw1,fjdw2,t,tabint,wetvar
     .        , fjd,ftol,p1,p2,t1,t2,w1,w2,fjd1,fjd2

      logical first_call,first_format, debug   
                  
c     save the values to avoid too many rereads
      save fjd1,fjd2,p1,p2,t1,t2,w1,w2,ftol

      data first_call/.true./,first_format/.true./

c     set debug for part of span
      debug = .false.
c**      if( t.gt.7200.d0 ) debug = .true.
     
c     Convert JD, T to fraction of day
      fjd = dble(jd) + t/86400.d0         

c     set tolerance value for 15s to avoid logic problems at edges of table span, 
c     individual entries, or slight changes in delay during the computation
      if( first_call )  ftol = 1.01d0 / 15.d0 / 24.d0

C     Check if time is within table interval
c     allow one minute before or after
      if( fjd.lt.fjdw1 ) then
        if( fjd.lt.(fjdw1-ftol)) then   
          write(message,'(a,f13.4,a,f13.4,a)') 
     .      'Time (',fjd,') before start of weather table (',fjdw1,')'
          call report_stat('FATAL','MODEL','atmred',' ',message,0)
        else
          fjd = fjdw1
        endif
      endif
      if( fjd.gt.fjdw2 ) then
        if( fjd.gt.(fjdw2+ftol)) then   
          write(message,'(a,f13.4,a,f13.4,a)') 
     .      'Time (',fjd,') after end of weather table (',fjdw2,')'
          call report_stat('FATAL','MODEL','atmred',' ',message,0)
        else
          fjd = fjdw2
        endif
      endif

                   
c     To start, read in two table entries
            
      if( first_call ) then
        read(iuw,'(a)',iostat=ioerr) line   
        if( ioerr.eq.-1 ) then
          call report_stat('FATAL','MODEL','atmred',' '
     .           ,'Unexpected EOF on W-file ',ioerr) 
        elseif( ioerr.ne.0 ) then
          call report_stat('FATAL','MODEL','atmred',' '
     .           ,'Error reading line of W-file ',ioerr) 
        endif   
c       preferred format is free, with HH and MM separated; if this fails, try
c       the old HHMM format
        read(line,*,iostat=ioerr) iyr,idyoyr,ihr,min,t1,p1,w1  
        if( ioerr.ne.0) then     
          read( iuw, '(1x,i2,1x,i3,1x,2i2,4x,f5.1,4x,f7.2,4x,f5.1)'
     .        ,iostat=ioerr) iyr, idyoyr, ihr, min, t1, p1, w1    
          if( ioerr.ne.0 )  call report_stat('FATAL','MODEL','atmred'
     .           ,' ','Error reading values from W-file ',ioerr) 
          call report_stat('WARNING','MODEL','atmred',' '
     .            ,'Old-format W-file; change when convenient',0) 
          first_format = .false.
        endif  
        call fix_y2k(iyr) 
        if( debug) print*,'read 1 ',iyr,idyoyr,ihr,min,t1,p1,w1 
        call monday( idyoyr,month,iday,iyr)
        fjd1= julday(month,iday,iyr) + (dble(ihr)+dble(min)/60.d0)/24.d0  
        read(iuw,'(a)',iostat=ioerr) line   
        if( ioerr.eq.-1 ) then
          call report_stat('FATAL','MODEL','atmred',' '
     .           ,'Unexpected EOF on W-file ',ioerr) 
        elseif( ioerr.ne.0 ) then
          call report_stat('FATAL','MODEL','atmred',' '
     .           ,'Error reading line of W-file ',ioerr) 
        endif   
        read(line,*,iostat=ioerr) iyr,idyoyr,ihr,min,t2,p2,w2
        if( ioerr.ne.0) then     
          read( iuw, '(1x,i2,1x,i3,1x,2i2,4x,f5.1,4x,f7.2,4x,f5.1)'
     .        ,iostat=ioerr) iyr, idyoyr, ihr, min, t2, p2, w2    
          if( ioerr.ne.0 )  call report_stat('FATAL','MODEL','atmred'
     .              ,' ','Error reading values from W-file ',ioerr) 
         if( first_format) call report_stat('WARNING','MODEL','atmred'
     .       ,' ','Old-format W-file; change when convenient',0) 
          first_format = .false.
        endif   
        call fix_y2k(iyr)
        if( debug) print *,'read 2 ',iyr,idyoyr,ihr,min,t2,p2,w2  
        call monday( idyoyr,month,iday,iyr)
        fjd2= julday(month,iday,iyr) + (dble(ihr)+dble(min)/60.d0)/24.d0
        first_call = .false.
      endif

c     Check for table values out of order

      if( fjd1.gt.fjd2 ) then
        write(message,'(a,f13.4,a,f13.4,a)') 
     .      'W-file table values out of order, fjd1=',fjd1,' fjd2=',fjd2
        call report_stat('FATAL','MODEL','atmred',' ',message,0)
      endif

c     If time is between the two read-in points, interpolate
              
      if( debug ) write(*,'(a,3f18.9)') 'fjd,fjd1,fjd2 ',fjd,fjd1,fjd2     
      if( fjd.ge.(fjd1-ftol) .and. fjd.le.(fjd2+ftol) ) then
        tabint = fjd2 -fjd1
        if( tabint.le.0 ) call report_stat('FATAL','MODEL'
     .     ,'atmred',' ','Time interval from W-file zero or negative',0)  
        temp = t1 + (t2-t1)/tabint*(fjd-fjd1)
        pres = p1 + (p2-p1)/tabint*(fjd-fjd1)
        wetvar = w1 + (w2-w1)/tabint*(fjd-fjd1)
        if( h2otyp.eq.upperc('R') ) wetvar = wetvar/100.
        if( debug ) print *,'Interpolating T P W',temp,pres,wetvar   
      
c     if time is too late, shift storage and read another value, then interpolate

      elseif( fjd.ge.(fjd2+ftol) ) then  
c       (allow 15s slop to avoid shift or reading beyond table when delay changes slightly)
        fjd1 = fjd2
        t1 = t2
        p1 = p2
        w1 = w2  
        read(iuw,'(a)',iostat=ioerr) line   
        if( ioerr.eq.-1 ) then
          call report_stat('FATAL','MODEL','atmred',' '
     .           ,'Unexpected EOF on W-file ',ioerr) 
        elseif( ioerr.ne.0 ) then
          call report_stat('FATAL','MODEL','atmred',' '
     .           ,'Error reading line of W-file ',ioerr) 
        endif   
        read(line,*,iostat=ioerr) iyr,idyoyr,ihr,min,t2,p2,w2  
        if( debug .and. ioerr.eq.0 ) 
     .     print *,'read 2 ',iyr,idyoyr,ihr,min,t2,p2,w2 
        if( ioerr.ne.0) then     
          read( iuw, '(1x,i2,1x,i3,1x,2i2,4x,f5.1,4x,f7.2,4x,f5.1)'
     .       ,iostat=ioerr) iyr, idyoyr, ihr, min, t2, p2, w2    
          if( debug) print *,'read 2 ',iyr,idyoyr,ihr,min,t2,p2,w2 
          if( ioerr.ne.0 ) call report_stat('FATAL','MODEL','atmred',' '
     .           ,'Error reading values from W-file ',ioerr)  
          if( first_format)  call report_stat('WARNING','MODEL','atmred'
     .     ,' ','Old-format W-file; change when convenient',0) 
          first_format = .false.
        endif 
        call fix_y2k(iyr)
        call monday( idyoyr,month,iday,iyr)
        fjd2= julday(month,iday,iyr) + (dble(ihr)+dble(min)/60.d0)/24.d0    
        if( debug ) write(*,'(a,3f18.9)') 
     .       'shifting fjd,fjd1,fjd2 ',fjd,fjd1,fjd2
  
        tabint = fjd2 -fjd1
        if( tabint.le.0 ) call report_stat('FATAL','MODEL','atmred',' '
     .     ,'Time interval from W-file zero or negative ',0 )  
        temp = t1 + (t2-t1)/tabint*(fjd-fjd1)
        pres = p1 + (p2-p1)/tabint*(fjd-fjd1)
        wetvar = w1 + (w2-w1)/tabint*(fjd-fjd1)
        if( h2otyp.eq.upperc('R') ) wetvar = wetvar/100.
        if( debug ) print *,'Interpolating T P W',temp,pres,wetvar   
      else    
        write(message,'(a,f13.4,a,2f13.4,a)') 
     .  'Requested time (',fjd,') illogical for current tabular times ('
     .        ,fjd1,fjd2,')'
         call report_stat('FATAL','MODEL','atmred',' ',message,0)
      endif
          
      return
      end
