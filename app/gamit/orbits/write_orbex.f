      Subroutine WRITE_ORBEX(lu,itype)

*     Write SV position and/or attitude values into an ORBEX file.
*     For GAMIT ORBEX, we can list data for every SV at every epoch since
*     the origins of the data are the t-file and y-files, based on models
*     not observations. 

*     R. King September 2019
                               
      implicit none

*      lu     : Logical unit number
*      itype  : Indicates what information to write on the file
*                = 1  SV attitude in quaternions
*                = 2  SV TRF position in meters
*                + 3  Both atitutde and position
                         
      integer*4 lu,itype

      include '../includes/dimpar.h'
      include 'orbex.h'

* Local information for the ORBEX header                       
     
        
      integer*4 start(5),end(5),nq
     .       , iyr,imonth,iday,ihr,imn,isec,ihnsec,iepc,isat,i

      real*4 version
      real*8 sec
      
      character*20 user,rec_types
      character*80 descrip

      logical debug/.true./                                   
                    
* Write the file description parts of the header
                                           
      version = 0.09 
      write(lu,'(a,1x,f5.2)') '%=ORBEX',version
      write(lu,'(a)') '%%'                           
      write(lu,'(a)') '+FILE/DESCRIPTION'
      if(itype.eq.1) then    
        descrip = 'Attitude quaternions for '//org//' products'
      elseif(itype.eq.2) then
        descrip = 'Orbital position (m) for '//org//' products'
      elseif(itype.eq.3) then
        descrip='Attitude quaternions and orbits for '//org//' products'
      endif               
      call getusr(user)
      write(lu,'(1x,a,8x,a)') 'DESCRIPTION',descrip
          write(lu,'(1x,a,9x,a)') 'CREATED BY',user
      call getdat(iyr,imonth,iday)
      call gettim(ihr,imn,isec,ihnsec)                             
      sec = dfloat(isec)
      write(lu,'(1x,a,6x,i4,5i3)') 
     .     'CREATION_DATE',iyr,imonth,iday,ihr,imn,isec
      write(lu,'(1x,a,9x,a)') 'INPUT_DATA',input_data 
      if(input_data(1:1).eq.'x') 
     .       write(lu,'(a,19x,a)') '*',input_data_comment
      write(lu,'(1x,a,12x,a)') 'CONTACT',user                   
      write(lu,'(1x,a,8x,a)') 'TIME_SYSTEM','GPS'              
      if(debug) print *
     .  ,'WRITE_ORBEX obx_jdstart obx_tstart ',obx_jdstart,obx_tstart
      call jdt2cdate(obx_jdstart,obx_tstart,iyr,imonth,iday,ihr,imn,sec)
      if(debug) print *
     . ,'WRITE_ORBEX obx_jdstart obx_tstart iyr imonth iday ihr imn sec'
     .             , obx_jdstart,obx_tstart,iyr,imonth,iday,ihr,imn,sec 
      write(lu,'(1x,a,9x,i4,4i3,f16.12)') 
     .    'START_TIME',iyr,imonth,iday,ihr,imn,sec
      call jdt2cdate( obx_jdstop,obx_tstop,iyr,imonth,iday,ihr,imn,sec)      
      write(lu,'(1x,a,11x,i4,4i3,f16.12)') 
     .    'END_TIME',iyr,imonth,iday,ihr,imn,sec
      write(lu,'(1x,a,5x,f9.3)') 'EPOCH_INTERVAL',obx_inter
      write(lu,'(1x,a,10x,a)') 'COORD_SYS',obx_frame
      write(lu,'(1x,a,9x,a)') 'FRAME_TYPE','ECEF'
      if(itype.eq.1) then
        rec_types = 'ATT'
      elseif(itype.eq.2) then 
        rec_types = 'ORB'
      elseif(itype.eq.3) then
        rec_types = 'ATT ORB'
      endif
      write(lu,'(1x,a,2x,a)') 'LIST_OF_REC_TYPES',rec_types                           
      write(lu,'(a)') '-FILE/DESCRIPTION'

* Write the satellite description part of the header

      write(lu,'(a)') '+SATELLITE/ID_AND_DESCRIPTION'
      do i=1,nosat
        write(lu,'(1x,a3)') prn(i)
      enddo                    
      write(lu,'(a)') '-SATELLITE/ID_AND_DESCRIPTION'

* Write the start line for the data records

      write(lu,'(a)') '+EPHEMERIS/DATA'

* Write the attitude records (optional comments and column headers first)

      if(itype.eq.1.or.itype.eq.3) then
        write(lu,'(2a)') 
     .    '*ATT RECORDS:  TRANSFORMATION FROM TERRESTRIAL FRAME '
     .    ,'COORDINATES (T) to SAT. BODY FRAME ONES (B) SUCH AS '
        write(lu,'(a)') '*            (0,B) = q.(0.T).trans(q)'
        write(lu,'(2a)') 
     .    '*REC ID_              N ____q0_(scalar)_____ _____q1__x____'
     .   ,'______ __q2___y___________ __q3__z_______'
        do iepc = 1,noepoch                                         
          call jdt2cdate( jdobx(iepc),tobx(iepc),iyr,imonth,iday
     .                 , ihr,imn,sec)
          write(lu,'(a,1x,i4,4i3,f16.12,i4)') 
     .       '##' ,iyr,imonth,iday,ihr,imn,sec,nosat
          nq = 4           
          if(debug) then
            print *,'WRITE_ORBEX iepc jdobx tobx '
     .              ,iepc,jdobx(iepc),tobx(iepc)
            print *,'  iyr,imonth,iday,ihr,imn,sec '
     .           ,iyr,imonth,iday,ihr,imn,sec
          endif
          do isat=1,nosat
            write(lu,'(1x,a,1x,a3,14x,i1,4f20.15)') 
     .         'ATT',prn(isat),nq,(quatern(i,isat,iepc),i=1,nq)
          enddo
        enddo      

      else
        call report_stat('FATAL','WRITE_ORBEX','orbits/write_orbex',' '
     .               ,'Only attitude quaternions currently supported',0)
      endif

      write(lu,'(a)') '-EPHEMERIS/DATA'
      write(lu,'(a)') '%END_ORBEX'

      return
      end 


*-----------------------------------------------------------------------------------------

      Subroutine JDT2CDATE(jd,t,yr,mon,day,hr,min,sec)
                      
*     For ORBEX, round to the nearest second

      implicit none

      integer*4 jd,yr,mon,day,hr,min,doy,isec
      real*8 t,tround,sec 
      logical debug/.true./
          
      if(debug) print *,'JDT2CDATE jd t ',jd,t
      call dayjul(jd,yr,doy)
      call monday(doy,mon,day,yr) 
      isec = dfloat(nint(t)) 
 
      tround = dble(isec)
      if(debug) print *,'JDT2CDATE isec tround',isec,tround
      call ds2hms(yr,doy,tround,hr,min,sec)
      if(debug) print *,'JDT2CDATE yr doy mon day hr min sec '
     .                           , yr,doy,mon,day,hr,min,sec 

      return
      end 


