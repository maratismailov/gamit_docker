      Program Y2ORBEX
                     
*     Convert the yaw angles on the y-file to quaternions and write them
*     into an ORBEX file, reading an sp3 file to get the terrestrial positions
*     of the satellites.   King September 2019

      implicit none

      include '../includes/dimpar.h'
      include 'orbex.h'

* Other variables
      integer*4 itype,iepc,isat,ioerr,iarg,start(5),stop(5)
     .        , im,id,i  
      real*8 sun(6),xhat(3),yhat(3),zhat(3),fjd,tdtoff,ytime 
     .     , yaw_s,yaw_e
      character*16 ytfiln 
      character*256 message 

* Function
      integer*4 iclarg,julday
      logical fcheck
                                
      logical debug/.true./
       
* Read the command line

      iarg = iclarg(1,yfiln)
      if ( iarg.le.0 ) then    
        write(*,'(4(/,a))')
     .       'y2orbex <y-file> <t-file> <ORBEX file> '
        stop
      endif
      iarg = iclarg(2,tfiln)
      if( iarg.le.0 ) then   
         call report_stat('FATAL','Y2ORBEX','y2orbex'
     .      ,' ','Missing input t-file name',0)
      endif  
      iarg = iclarg(3,obxfiln)
      if( iarg.le.0 ) then   
         call report_stat('FATAL','Y2ORBEX','y2orbex'
     .      ,' ','Missing output ORBEX file name',0)
       endif
       iarg = iclarg(4,start(1))
       if(debug) print *,'Command line 4 iarg ',iarg 
       if( iarg.ne.0 ) then
         iarg = iclarg(5,start(2))
         iarg = iclarg(6,start(3))
         iarg = iclarg(7,start(4))
         iarg = iclarg(8,start(5))
         iarg = iclarg(9,stop(1))
         if( iarg.eq.0 ) then
           call report_stat('FATAL','Y2ORBEX','y2orbex'
     .      , ' ','If start time given, stop time must be also',0)
         else
           iarg = iclarg(10,stop(2))
           iarg = iclarg(11,stop(3))
           iarg = iclarg(12,stop(4))
           iarg = iclarg(13,stop(5))
         endif 
       else
*        zero out for testing later to set defaults
         do i=1,5
           start(i) = 0
           stop(i) = 0
         enddo
       endif
       if(debug) print *,'Y2ORBEX start stop ',start,stop 
                                                           
* Assign the unit numbers

      luy = 1
      lut = 2
      luobx = 3
      lunut = 20
      luut1 = 21
      lupole = 22
     
* Open the files

*      Y-file                                 
      open(unit=luy,file=yfiln,form='unformatted' 
     .    ,access='sequential',status='old',iostat=ioerr)
      if(ioerr.ne.0) call report_stat('FATAL','Y2ORBEX'
     .   ,'y2orbex',yfiln,'Error opening y-file ',ioerr)

*       T-file
      open(unit=lut,file=tfiln,form='unformatted' 
     .    ,access='sequential',status='old',iostat=ioerr)
      if(ioerr.ne.0) call report_stat('FATAL','Y2ORBEX'
     .   ,'y2orbex',tfiln,'Error opening t-file ',ioerr)

*       ORBEX file
      open(unit=luobx,file=obxfiln,status='unknown',iostat=ioerr)
      if (ioerr .ne. 0) call report_stat('FATAL','Y2ORBEX'
     .    ,'y2orbex',obxfiln,'Error opening ORBEX file:',ioerr)
                       
* EOP files

      if( fcheck('nbody') ) then 
        lunut = 0
        call report_stat('STATUS','Y2ORBEX','y2orbex',' ' 
     .            ,'nbody file exists, using MHB_2000 for nutations',0)
      else
        open( unit=lunut,file='nutabl.',status='old',iostat=ioerr)
        if (ioerr .ne. 0 ) then
            call report_stat('FATAL','Y2ORBEX','y2orbex','nutabl.'
     .                      ,'Error opening nutation table: ',ioerr)
        endif
      endif    
      open(unit=luut1,file='ut1.',status='old',iostat=ioerr)
      if (ioerr .ne. 0 ) then
          call report_stat('FATAL','Y2ORBEX','y2orbex','ut1.'
     .                    ,'Error opening UT1 table: ',ioerr)
      endif
      open(unit=lupole,file='pole.',status='old',iostat=ioerr)
      if (ioerr .ne. 0 ) then
          call report_stat('FATAL','Y2ORBEX','y2orbex','pole.'
     .                    ,'Error opening pole table: ',ioerr)
      endif
                           
* Get the yaw angles from the y-file

      read(luy,iostat=ioerr) nyversn
      if( ioerr.ne.0 ) then
           call report_stat('FATAL','Y2ORBEX','y2orbex',yfiln
     .      ,'Error reading version from y-file header--old-style file?'
     .      ,ioerr)  
      else
        if(nyversn.lt.1051.or.nyversn.gt.1061  ) then
          write(message,'(a,i4,a)') 'Incompatible y-file version ('
     .      ,nyversn,')  File: '
          call report_stat('FATAL','Y2ORBEX','y2orbex',yfiln,message,0)
        endif
      endif
      rewind(luy)
      if( nyversn.eq.1051 ) then
        read(luy) nyversn,ytfiln,yaw_s,yaw_e,nyinter,nyepoch
     .          , nysat,(iprn(isat),ysvn(isat),yblk(isat),isat=1,nysat)
      elseif( nyversn.eq.1061 ) then
        read(luy) nyversn,ytfiln,yaw_s,yaw_e,nyinter,nyepoch,nysat
     .           , (iprn(isat),ysvn(isat),svantbody(isat),isat=1,nysat)
      endif       
      if(debug) print *,'nyversn nyinter nyepoch nysat '
     .                 , nyversn,nyinter,nyepoch,nysat           
      gnss = yfiln(5:5)
      call uppers(gnss)
      do isat=1,nysat
        prn(isat)(1:1) = gnss  
        write(prn(isat)(2:3),'(i2.2)') iprn(isat)
      enddo
      if(debug) write(*,'(a,50(1x,a3))')  'PRNs: ',(prn(i),i=1,nysat)
      do iepc = 1,nyepoch
        read(luy) ytime,(yawang(i,iepc),ievent(i,iepc),i=1,nysat)  
        jdy(iepc) = int(ytime)
        ty(iepc) = (ytime - dfloat(jdy(iepc)))*86400.d0
        if(debug) then 
           if(iepc.le.20) print *,'iepc ytime jdy ty '
     .                    ,iepc,ytime,jdy(iepc),ty(iepc)
        endif
      enddo                         
         
*  Get the Earth-fixed SV coordinates at the y-file epochs

      call get_svcoords
              
*  Set the ORBEX header values and start/stop times 
        
      org = 'MIT'  
      obx_inter = dfloat(nyinter) 
      if(start(1).eq.0) then 
*       obx_jdstart default set in get_svcoords, times are 0h 0m and 23h 59m 30s 
        obx_jdstop = obx_jdstart
        obx_tstart = 0.d0
        obx_tstop = 86370.d0 
      else
        call monday(start(2),im,id,start(1))
        obx_jdstart = julday(im,id,start(1))
        obx_tstart = 3600.d0*start(3) + 60.d0*start(4) + start(5)
        call monday(stop(2),im,id,stop(1))
        obx_jdstop = julday(im,id,stop(1))
        obx_tstop = 3600.d0*stop(3) + 60.d0*stop(4) + stop(5)
      endif           
      if(debug) then 
        print *,'Y2ORBEX start obx_jdstart obx_tstart '
     .    ,start,obx_jdstart,obx_tstart
        print *,'Y2ORBEX stop obx_jdstop obx_tstop '
     .    ,stop,obx_jdstop,obx_tstop
      endif
      obx_frame = ' ' 
      input_data = 'x'
      input_data_comment = 'Values derived from GAMIT yaw models'
      contact = 'tah@mit.edu'
             

*  Convert the yaw angle to quaternions and reduce the arrays to the span requested 
                           
      noepoch = nyepoch        
      nosat = nysat 
      if(debug) print *,'Y2ORBEX nysat nosat noepoch '
     .    ,nysat,nosat,noepoch
      call quaternion
      if(debug) print *,'Y2ORBEX quatern(1,1,1) ',quatern(1,1,1)

*  Write the ORBEX file
            
*     write only attitude values 
      itype = 1    
      if(debug) print *,'Y2ORBEX obx_inter obx_jdstart obx_tstart '
     .                 ,         obx_inter,obx_jdstart,obx_tstart 
      call write_orbex(luobx,itype)

      stop
      end                     


