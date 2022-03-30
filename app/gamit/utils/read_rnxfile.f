       subroutine read_metrnx( iurnx,maxmet,times,pres,temp )

c Subroutine to read pressure and temperature values from RINEX met file  
c Adapted from model/read_metrnx.  R. King 26 January 2007
          
c Number, times, and values of met data stored in common /ufcom/ in model.h
c Arrays used for values either from a RINEX met file or those interpolated
c from a global grid and stored on the u-file

      implicit none
       
c  Input
c   iurnx   unit number for RINEX met file
      integer*4 iurnx
c   maxmet  dimensions of time and values arrays
      integer*4 maxmet
c   times    epochs of values values from RINEX file (decimal years)
      real*8 times(maxmet)           
c   pres, temp   values from RINEX file (hPa, degK)
      real*8 pres(maxmet),temp(maxmet)
                                      
      integer*4 iuw,date(5),hr,min,i,j,ii,ioerr,nobs  
      real*4 rxver 
      real*8 val(3),jd,jan1_jd,sec 
      character*2 obstyp(3)
      character*256 line,message
      logical eof,found_type,found_sensor
               
c  Initialization
      line = ' ' 
c      # of values per epoch (P, T, H)_
      nmet = 3                        
c      # of epochs, incremented as the file is read
c      stored with times and values in common /ufcom/ in ../includes/model.h
      ntmet = 0 

            
c  Read the header (order of entries is arbitrary except version number)
                                         
      found_type = .false.
      found_sensor = .false.
      read(iuw,'(a)',iostat=ioerr) line  
      if( line(61:73).ne.'RINEX VERSION' ) then
           call report_stat('FATAL','MODEL','read_metrnx',' '
     .                     ,'Bogus first line of RINEX met file',ioerr) 
      else
c          If the version number is an integer written under version 1 format
c          ( I6 ) it will not be read correctly by the version 2 format ( f9.2)
c          so read the version number free-format
         read( line(1:9),*,iostat=ioerr) rxver  
         if( ioerr.ne.0 ) call report_stat('FATAL','MODEL','read_metrnx'
     .      ,' ' ,'Error reading RINEX met file eversion number ',ioerr) 
      endif
      do while(line(61:73).ne.'END OF HEADER')
        read(iuw,'(a)',iostat=ioerr) line  
        if( ioerr.eq.-1 ) then
           call report_stat('FATAL','MODEL','read_metrnx',' '
     .      ,'End-of-file before END OF HEADER in RINEX met file',ioerr)
        elseif( ioerr.ne.0 ) then
          call report_stat('FATAL','MODEL','read_metrnx',' '
     .           ,'Error reading RINEX met file header',ioerr)
        else 
          if( line(61:69).eq.'# / TYPES' ) then
            read(line(1:24),'(i6,3(4x,a2))',iostat=ioerr) 
     .           nobs,(obstyp(i),i=1,3)  
            if( ioerr.ne.0 )  then
              call report_stat('FATAL','MODEL','read_metrnx',' '
     .             ,'Error reading the observation order',ioerr)   
            else
              found_type = .true.
            endif
          elseif( line(58:70).eq.'PR SENSOR POS' ) then
            read(line(43:56),'(f14.4)',iostat=ioerr) sensor_ht
            if( ioerr.ne.0 ) then
              call report_stat('FATAL','MODEL','read_metrnx',' '
     .           ,'Error reading the pressure sensor position',ioerr) 
            else
              found_sensor = .true.
            endif
          endif
        endif
      enddo
      if( .not.found_type ) then
         call report_stat('FATAL','MODEL','read_metrnx', ' '
     .            ,'Missing # / TYPE line in met rinex file',ioerr)
      endif
      if( .not.found_sensor ) then
        call report_stat('WARNING','MODEL','read_metrnx', ' '
     .     ,'Missing pressure sensor position info  ',0)
      endif                                                     
c      print *,'height of pressure sensor is',sensor_ht 
      

c  now we just read all the data through to the end of the file. Format is:
c  05  5 26  0  0  0  969.3   -8.3   82.7
c  05  5 26  0  0 30  969.3   -8.4   82.7 
             
      ntmet = 0 
      eof =.false.
      do while( .not.eof ) 
        read(iuw,'(5i3,f3.0,3f7.1)',iostat=ioerr) 
     .              yr,mon,day,hr,min,sec,(val(j),j=1,3)
        if( ioerr.eq.-1 ) then
          eof = .true.
        elseif( ioerr.ne.0 ) then
          call report_stat('FATAL','MODEL','read_metrnx', ' '
     .    ,'Error reading values from RINEX met file',ioerr)  
        elseif( ntmet+1.gt.maxmet ) then                           
          write(message,'(a,i4)') 
     .      'Number of met values on RINEX file > maxmet, =',maxmet
          call report_stat('FATAL','MODEL','read_metrnx', ' ',message,0)
        else
          ntmet = ntmet + 1
c         convert the year, month,day, hr, min, sec into decimal year   
          call fix_y2k(yr)
          jdoy = idoy(yr,mon,day)
          sod = ihr*3600. + min*60. + sec
c         time tag is decimal years 
          times(ntmet)= decyrs( iyr,jdoy,sod )  
c         Now, assign the values as indicated by the headers
          if(obstyp(1).eq."PR") pres(ntmet)=val(1)
          if(obstyp(1).eq."TD") temp(ntmet)=val(1)
          if(obstyp(2).eq."PR") pres(ntmet)=val(2)
          if(obstyp(2).eq."TD") temp(ntmet)=val(2)
          if(obstyp(3).eq."PR") pres(ntmet)=val(3)
          if(obstyp(3).eq."TD") temp(ntmet)=val(3)
c          print*,met_time(ntmet),met_val(ntmet,1)
        endif
      enddo
          
c  Augment the values at the beginning and the end to minimize the
c  situation in which we need to extrapolate or substitute nominal
c  values.  Three hours seems reasonably safe and will produce
c  less discontinuity than reverting to the nominal value, which
c  might be the RINEX start, midpoint, or end, or a value from 
c  GPT.   rwk 070110
      
      ntmet = ntmet + 2   
      if (ntmet.gt.maxmet ) then    
        write(message,'(a,i4)') 
     .    'Number of met values on RINEX file > maxmet =',maxmet
          call report_stat('FATAL','MODEL','read_metrnx', ' ',message,0)
      endif
      do i=1,ntmet-2
        ii = ntmet-i
        times(ii) = times(ii-1)  
        pres(ii) = pres(ii-1)
        temp(ii) = temp(ii-1)
      enddo
      times(1) = times(2) - 3./24.        
      times(ntmet) = times(ntmet-1) + 3./24.  
      pres(1) = pres(2)
      temp(1) = temp(2)

      return
      end






