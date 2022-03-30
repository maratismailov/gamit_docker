       subroutine read_metrnx

c Subroutine to read met rinex file and extract all the required information
c
c P. Tregoning c 12 June 2006
c Modified by R. King  23 July 2006
          
c Number, times, and values of met data stored in common /ufcom/ in model.h
c Arrays used for values either from a RINEX met file or those interpolated
c from a global grid and stored on the u-file

      implicit none

      include '../includes/dimpar.h'
      include '../includes/model.h'

      integer*4 date(5),hr,min,i,j,ii,ioerr,nobs  
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
c** rwk 111214: The following is the authorized format for RINEX v 1 or 2, but
c               some providers are using more than 7 columns for the values, so
c               substitute a free-format read to compensate.
c       read(iuw,'(5i3,f3.0,3f7.1)',iostat=ioerr) (date(j),j=1,3)
c    .                                 , hr,min,sec,(val(j),j=1,3)     
        read(iuw,'(a)',iostat=ioerr) line  
        if( ioerr.eq.-1 ) then
          eof = .true.
        elseif( ioerr.ne.0 ) then
          call report_stat('FATAL','MODEL','read_metrnx', ' '
     .    ,'Error reading values line from RINEX met file',ioerr)  
        elseif( ntmet+1.gt.maxmet ) then                           
          write(message,'(a,i4)') 
     .      'Number of met values on RINEX file > maxmet, =',maxmet
          call report_stat('FATAL','MODEL','read_metrnx', ' ',message,0)
        else                        
         read(line,'(5i3,f3.0)',iostat=ioerr) (date(j),j=1,3),hr,min,sec
         if( ioerr.ne.0 ) then
           call report_stat('FATAL','MODEL','read_metrnx', ' '
     .     ,'Error reading time from RINEX met file',ioerr)  
         endif
         read(line(19:60),'(10f7.1)',iostat=ioerr) (val(j),j=1,3) 
c         print *,'ntmet val ',ntmet,val
         if( ioerr.ne.0 ) then
           call report_stat('FATAL','MODEL','read_metrnx', ' '
     .     ,'Error reading values from RINEX met file',ioerr)  
         endif
         ntmet = ntmet + 1
c        convert the year, month,day, hr, min, sec into decimal day of year
         date(4) = 0
         date(5) = 0
         call ymdhms_to_jd(date,0.d0,jd)
*        Now do Jan 1
         date(2) = 1
         date(3) = 1
         call ymdhms_to_jd( date, 0.d0,Jan1_jd)
*        then ...                     
         met_time(ntmet)  = jd - Jan1_jd + 1.d0 + hr*1.d0/24.d0
     .                     + min*1.d0/1440.d0+sec/86400.d0
c         print*,jd,jan1_jd,ntmet,met_time(ntmet)

c  Now, assign the values in the right order: pres, temp, humidity
c    rwk: index order reversed from PT original to be consistent with other arrays in /ufcom/
          if(obstyp(1).eq."PR") met_val(ntmet,1)=val(1)
          if(obstyp(1).eq."TD") met_val(ntmet,2)=val(1)
          if(obstyp(1).eq."HR") met_val(ntmet,3)=val(1)

          if(obstyp(2).eq."PR") met_val(ntmet,1)=val(2)
          if(obstyp(2).eq."TD") met_val(ntmet,2)=val(2)
          if(obstyp(2).eq."HR") met_val(ntmet,3)=val(2)

          if(obstyp(3).eq."PR") met_val(ntmet,1)=val(3)
          if(obstyp(3).eq."TD") met_val(ntmet,2)=val(3)
          if(obstyp(3).eq."HR") met_val(ntmet,3)=val(3)
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
        met_time(ii) = met_time(ii-1)  
        do j=1,3
          met_val(ii,j) = met_val(ii-1,j)
        enddo
      enddo
      met_time(1) = met_time(2) - 3./24.        
      met_time(ntmet) = met_time(ntmet-1) + 3./24.  
      do j = 1,3 
        met_val(1,j)  = met_val(2,j)
        met_val(ntmet,j) =  met_val(ntmet-1,j)
      enddo                         
c      do i=1,ntmet
c        print *,'met_time met_val ',met_time(i),met_val(i,1)
c      enddo

      return
      end






