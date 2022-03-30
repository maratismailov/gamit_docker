c   Extract zenith delays and precipitable water from GAMIT output 

c     R. King 24 June 2004; modified 15 December 2004 to add P. Tregoning code 
c     to read a RINEX met file for P, T, RH and/or a SINEX files w/ TZD values.
c     Also add ability to read/use gradient values from o-file.

c The calling sequence is

c      metutil [ZTD-file] [interval] [scalesig] [z-file] [ RINEX m-file] 
                        
c     [ZTD-file]  File containing the GPS estimates and uncertainties.  
c                 It can be either a GAMIT o-file  or a SINEX zpd-file. 

c     [interval]  Interval at which to evaluate the estimates.  The entry
c                 can be either E, to use the interval of the estimates
c                 provided by the ZTD file, O to use the interval provided
c                 by the observations (z-file), or integer seconds to use
c                 evenly spaced times from the start.  The start time will
c                 be the first epoch of the GAMIT or SINEX file session but
c                 epochs missing observations or surface met data will be
c                 skipped in the output.  The default is to use the epoch of 
c                 the estimates (E).

c     [scalesig]  Factor by which to scale the uncertainties (default 1.0) 


c     [z-file]    from MODEL, containing a priori T, P, RH, DZD, WZD, TZD, DMap, 
c                 WMap, DSD, WSD, TSD at each epoch; name must contain site id 
c                 as characters 2-5    

c     [RINEX m-file]  RINEX met file containing T, P, RH (or ZWD from WVR)               

c     Notes:
c     ------  

c      (1) Either a z-file or a RINEX met file must be present to provide the site 
c          name and day-of-year.  These files also provide the a priori and/or measured 
c          values of pressure for the hydrostatic delay.  If both files are present the 
c          RINEX file values take precedence over the apriori values from MODEL recorded
c          in the z-file.

c      (2) The output file will be named met_[site].[yyddd].

      implicit none        

      include '../includes/dimpar.h'
                                    
      character*2 ayr
      character*3 sampling,adoy       
      character*4 site 
      character*10 zfile
      character*12 rnxfile
      character*14 outfile
      character*16 ztdfile,arg 

      integer*4 iyr2,numatm,numgrad
     .        , iclarg,iarg,ioerr,interval,idate(5)
     .        , iuztd,iuz,iurnx,iuout,doy,j,i
            


c     set maxatm to allow for estimated values and half-hour intervals
c      parameter(maxatm=49)
cRWK 170901: Now set in dimpar.h

    
      real*8 Tm
     .     , tatm(maxatm),estatm(maxatm),adjatm(maxatm),sigatm(maxatm)
     .     , scalesig,estztd,aptzd,adjztd,sigztd,sigpw,zhd,zwd,pw,presr
     .     , t,tstart,tend,dt,temp,sec
     .     , ht,sig_ht,lat,pi,rho,rv,k2p,k3
     .     , tgrad(maxatm)
     .     , estnsgrad(maxatm),signsgrad(maxatm),nsgrad,signs
     .     , estewgrad(maxatm),sigewgrad(maxatm),ewgrad,sigew

      logical finished
                     

c Print help if no arguments

      iarg = iclarg(1,arg)
      if( iarg.eq.0 ) then
        print *,    
     .   'Example:  metutil [ZTD-file] [interval] [scalesig] [z-file]'
     .   ,' [ RINEX m-file]'
        stop
      endif

c Define the constants for converting ZWD to PW
      pi = 3.141526954d0
c       density of liquid water = 1.000 g/cm**3
      rho = 1.000
c      specific gas constant for water vapor = R/Mw = 8.314/18.0152 10**7 dyne-cm/(g K)
c      from Askne et al. [Radio Sci. 22, 379, 1987] 
      Rv = (8.314/18.0152)*1.e7
c      k2 prime from Bevis et al. Eq 4  = 22.1 K/mb = 22.1/(1.013 10**3) K/dyne-cm**2
      k2p = 22.1/1.013*1.e-3     
c      k3 from Bevis et al Table 1 = 3.739*10**-5 K**2/mb = (3.379 10**-5/1.013 10**3) K**2/(dyne cm**2)
      k3 = 3.739e5/1.013e3
        

c Print some stuff to the screen

      print *,' '
      print *,' Starting METUTIL '
      print *,' '  

c Get the run-string arguments and open the files

      call get_input( ztdfile, zfile, rnxfile, iuztd, iuz, iurnx
     .              , sampling, interval, site, doy, scalesig )  
cd      print*,ztdfile, zfile, rnxfile, iuztd, iuz, iurnx
cd     .              , sampling, interval, site, doy, scalesig

c Read the estimated ZD parameters into storage   

      call read_estimates( iuztd,ztdfile,site,doy,maxatm
     .                  , numatm,tatm,estatm,adjatm,sigatm
     .                  , numgrad,tgrad
     .                  , estnsgrad,signsgrad,estewgrad,sigewgrad
     .                  , lat,ht,sig_ht )   
cd       print *,'DEBUG numatm ',numatm
cd       do j=1,numatm
cd         print *,tatm(j),estatm(j),sigatm(j)  
cd       enddo
            
 
c Set the start time and interval    
                       
      tstart = tatm(1)
      tend = tatm(numatm)
      call get_interval( sampling,interval,tstart,tend,numatm,iuz,dt )

c Open the output file        
                  
      iuout = 4    
      outfile = 'met_site.yyddd'  
      call decyrs_to_ydhms(tstart,idate)
      iyr2 = mod(idate(1),100)
      write(ayr,'(i2)') iyr2
      if(ayr(1:1).eq.' ' ) ayr(1:1) = '0' 
      write(adoy,'(i3)') idate(2)
      if(adoy(1:1).eq.' ') adoy(1:1) = '0'
      if(adoy(2:2).eq.' ') adoy(2:2) = '0'  
      write(outfile(5:8),'(a4)') site
      write(outfile(10:11),'(a2)') ayr
      write(outfile(12:14),'(a3)') adoy
      open(unit=iuout,file=outfile,status='unknown',iostat=ioerr)
      if( ioerr.ne.0 ) then
        print*, 'Error opening output',outfile,ioerr
      else
        print *,'Opened output file ',outfile  
      endif
               
c Write the header to the output met file
                                              
      call uppers(site)
      write(iuout,'(a,a4,a,f10.4,a,f7.4,a)') 
     .      '* Estimated atmospheric values for ',site
     .      ,'. Height estimate: ',ht*1.d3,' +/-',sig_ht,' m.'
      write(iuout,'(a)') '* METUTIL Version 3.0  2009-08-27'
      write(iuout,'(a,a16,1x,a10,1x,a12,a,f5.1)') '* Input files: '
     .   ,ztdfile,zfile,rnxfile
     .   ,'  ZTD-file sigmas scaled by ',scalesig
c** old METUTIL format
c      write(iuout,'(2a)') '* Yr  Doy Hr Mn Sec  '
c     .     ,' Total Zen   Wet Zen    Sig Zen       PW       Sig PW (m) '
c** new met_ file format
      write(iuout,'(3a)') '* Yr  Doy Hr Mn Sec  '
     . ,'Total Zen  Wet Zen  Sig Zen    PW   Sig PW (mm)  Press (hPa)'
     . ,'  Temp (K) ZHD (mm)   Grad NS   Sig NS  Grad EW   Sig EW (mm)'
c
             
c  Compute ZWD, PW and write the output file

      t = tstart    
      finished = .false.
      do while( .not.finished )  
        call get_estimates( t,maxatm,numatm,tatm,estatm,adjatm,sigatm  
     .        , numgrad,tgrad,estnsgrad,signsgrad,estewgrad,sigewgrad
     .        , estztd,adjztd,sigztd,ewgrad,sigew,nsgrad,signs )   
        call get_zhd(t,iuz,iurnx,lat,ht,zhd,presr,temp,aptzd)  

c PT060303: if reading from a z-file, use the apriori tzd from
c           there + adjztd to get the correct value of estztd
c           (this avoids problems of SOLVE not writing out the
c           correct ZTDs when time-varying a priori ZTDs are used
c           in MODEL).
        if (iuz.gt.0) then
          estztd = aptzd + adjztd
        endif
        zwd = estztd - zhd    
c       get weighted mean temperature Tm from site temperature using 
c       approximation from Bevis et al. [J. Appl. Metero. ,33, 379, 1994]
        Tm = 70.2 + 0.72*temp
c        then calculate the PW
        pw = zwd * 1.e6/(rho*Rv*(k3/Tm + k2p))  
        sigpw = sigztd * 1.e6/(rho*Rv*(k3/Tm + k2p)) 
        call decyrs_to_ydhms( t, idate )   
        call fix_y2k(idate(1))            
        sec = idate(5)    
c**     temporary write in the old metutil format  (prior to 041229?)
c        write(iuout,'(1x,2i4,2i3,f4.0,5f11.4)')
c     .       (idate(j),j=1,4),sec,estztd/1.d3,zwd/1.d3
c     .       ,sigztd/1.d3,pw/1.d3,sigpw/1.d3      
c**     new format  041229?; gradient columns added 090827 Version 3.0 
        if( zhd.eq.0.d0 ) then
          print *,'ZHD not available, skipping ',(idate(j),j=1,4),sec
        else
          write(iuout,'(1x,2i4,2i3,f4.0,f11.2,1x,4f8.2,7x,f8.2
     .                  ,4x,f8.2,3x,f7.2,2(1x,4f9.2))')
     .       (idate(j),j=1,4),sec,estztd,zwd
     .       ,sigztd,pw,sigpw,presr,temp,zhd
     .       ,nsgrad,signs,ewgrad,sigew
        endif
        t = t + dt     
        if( (t-1.17d-7).gt.tend ) finished = .true.    
c debug        write (*,'(a,f15.8,e10.3,f15.8,l1)') 't dt tend finished '
c     .     ,t,dt,tend,finished 
      enddo

      print *,'Normal end of METUTIL '
      stop
      end      

c---------------------------------------------------------------------------

    
      Subroutine get_input( ztdfile, zfile, rnxfile, iuztd, iuz, iurnx
     .                    , sampling, interval, site, doy, scalesig )    

      implicit none
                     
      character*3 sampling    
      character*4 site     
      character*10 zfile  
      character*12 rnxfile
      character*16 ztdfile,arg
      character*256 line
   
      integer*4 iuztd,iuz,iurnx,iarg,interval,doy,iclarg,ioerr

      real*8 scalesig 
         
      logical found 

c Set the defaults for the file names 

      ztdfile = ' '
      zfile = ' ' 
      rnxfile = ' ' 

c Set the I/O unit number --changed to 0 if file not available
             
      iuztd = 1 
      iuz   = 2
      iurnx = 3
             

c Get the run-string arguments and open the files     

      iarg = iclarg(1,ztdfile)
      open(unit=iuztd,file=ztdfile,status='old',iostat=ioerr) 
      if( ioerr.ne.0 ) then
        print *,'Error opening input o- or SINEX file',ztdfile,ioerr
        stop
      else
        print *,'Opened input file of estimates ',ztdfile                        
      endif
            
      sampling = 'est'             
      interval = 0 
      iarg = iclarg(2,arg)   
      if( iarg.gt.0 ) then 
        if( arg(1:1).eq. 'E') then
          sampling = 'est'  
          print *,'Epochs of estimated values used for output'
        elseif ( arg(1:1).eq.'O' ) then
          sampling = 'obs' 
         print *,'Epochs of observations used for output'
        else
          read(arg,'(i4)',iostat=ioerr) interval
          if( ioerr.ne.0 ) then
            print *,'Unrecognized input interval '
            stop 
          else    
            sampling = 'inp'
            print *,'Interval for output is ',interval,' sec'
          endif
        endif 
      else
        sampling = 'est'
        print *,'Epochs of estimated values used for output' 
      endif                         
  
      scalesig = 1.0
      iarg = iclarg(3,arg)   
      if( iarg.gt.0 ) then
        read(arg,'(f8.0)') scalesig
      endif                         
      print *,'Estimated zenith delay sigmas scaled by ',scalesig
  
      iarg = iclarg(4,zfile)  
      if( iarg.gt.0 ) then 
        open(unit=iuz,file=zfile,status='old',iostat=ioerr)
        if( ioerr.ne.0 ) then
          print*, 'Error opening z-file ',zfile,ioerr
        else
          print *,'Opened z-file ',zfile
     .        , ' --allow 3-hr extrapolation of ZHD'   
        endif
      else
        iuz = 0
      endif  

      iarg = iclarg(5,rnxfile)  
      if( iarg.gt.0 ) then 
        open(unit=iurnx,file=rnxfile,status='old',iostat=ioerr)
        if( ioerr.ne.0 ) then
          print*, 'Error opening RINEX met file',rnxfile,ioerr
        else
          print *,'Opened RINEX met file ',rnxfile
     .        , ' --allow 3-hr extrapolation of ZHD'     
        endif
      else
        iurnx = 0 
      endif  
        

c Get the 4-character site name and day-of-year from the file names or headers

c     if a z-file is available, read the name from the header; otherwise
c     use the name of the RINEX file
      if( iuz.gt.0  ) then
        read(iuz,'(34x,a4)',iostat=ioerr) site  
        if( ioerr.ne.0 ) then
         print *,'Error reading site name from header of z-file ',ioerr
         stop 
        else
          call lowers(site)
c         got site ok, now get day-of-year
c         found = .false.
          do while (.not.found)
            read(iuz,'(a)',iostat=ioerr) line
            if( ioerr.ne.0 ) then
              print *,'EOF before day-of-year read in Z-file'
              stop
            endif
            if( line(1:1).eq.' ' ) then
              read(line,'(6x,i4)',iostat=ioerr) doy
              if( ioerr.ne.0 ) then
                print *,'Error reading day-of-year from Z-file '
                stop
              else
                found = .true.
             endif
            endif
          enddo 
        endif
          
      elseif( iurnx.gt.0 ) then
        read(rnxfile(1:4),'(a4)',iostat=ioerr) site
        if( ioerr.ne.0 ) then
          print *,'Error reading site name from RINEX met file name '
     .       ,ioerr                                                  
           stop
        endif   
        read(rnxfile(5:7),'(i3)',iostat=ioerr) doy
        if( ioerr.ne.0 ) then
          print *,'Error reading day-of-yr from RINEX met file name '
     .      ,ioerr
          stop
        endif

      else
         print *,'Neither z-file nor RINEX met file available, '
     .          ,'cannot get site name '
         stop
      endif     

      return
      end
                 
c---------------------------------------------------------------------------    

      Subroutine read_estimates( iuztd,ztdfile,site,doy,maxatm
     .                         , numatm,tatm,estatm,adjatm,sigatm
     .                         , numgrad,tgrad
     .                         , estnsgrad,signsgrad,estewgrad,sigewgrad
     .                         , lat, ht, sig_ht )

      implicit none

      integer*4 iuztd,maxatm,numatm,deg,iyr,ihr,min,doy,jyr,jdoy,jsod
     .        , idoy,mon,iday,ioerr,i,numgrad

      real*8 tatm(maxatm),estatm(maxatm),adjatm(maxatm),sigatm(maxatm)
     .     , est,adj,sig,sec,semi,finv,clatd,radius,glatd,ht,lat,pi
     .     , rot_mat(3,3),xyz_pos(3),geod_pos(3),sod,decyrs,sig_ht
     .     , tgrad(maxatm)
     .     , estnsgrad(maxatm),signsgrad(maxatm)
     .     , estewgrad(maxatm),sigewgrad(maxatm)

      character*1 nos   
      character*4 site,siteuc,site1
      character*5 filetype
      character*16 ztdfile           
      character*256 line

      logical found,endatm,endgrad,eof

      data pi/3.14159265d0/,semi/6378.137d0 /,finv/298.257222101d0/
             

c See if files of estimates is a GAMIT o-file or a SINEX file
                       
c     if the file name starts with the site-name, assume it's a SINEX file
      if( ztdfile(1:4).eq.site ) then
         filetype = 'sinex'
      else
         filetype = 'ofile'  
cd	 print *,'filetype ofile '
c        check to be sure
         if( ztdfile(1:1).ne.'o' ) then
           print *,'Name of estimates file does not match '
     .           ,'SINEX or o-file convention'
           stop
         endif
      endif

      
c Read the estimates from the file and convert to mm
    
      siteuc = site
      call uppers(siteuc)
      if( filetype.eq.'ofile' ) then 
        call uppers(siteuc)   
c       first look for site coordinates get get latitude and height
        found = .false.     
        do while (.not.found)
          read(iuztd,'(a)',iostat=ioerr) line 
          if( ioerr.eq.-1 ) then
             print *,'Site coordinates not found on o-file ',ioerr
             stop
          elseif( ioerr.ne.0 ) then
             print *,'Error reading coordinates from o-file ',ioerr
             stop
          elseif(line(7:10).eq.siteuc.and.line(12:19).eq.'GEOC LAT')
     .      then
            found =.true. 
            read(line(29:40),'(a1,i2,1x,i2,1x,f5.2)',iostat=ioerr) 
     .           nos,deg,min,sec       
cd              print *,'read nos deg min sec ',nos,deg,min,sec       
            if( ioerr.ne.0 ) then
              print *,'Error decoding lat line of o-file ',ioerr
              stop
            endif
            call dmsdeg( nos,deg,min,sec,clatd)  
c           skip longitude line and read radius
            read(iuztd,'(/,a)',iostat=ioerr ) line
	    if( ioerr.ne.0 ) then
	       print*,'Error reading radius line from o-file'
	    elseif( line(12:17).ne.'RADIUS' ) then
	       print *,'Unexpected string reading radius from o-file'
	    else
              read(line(55:88),'(e11.4,10x,f13.8)',iostat=ioerr)
     .             sig_ht,radius 
cd              print *,'read full sig_ht rad ',sig_ht,radius
              if( sig_ht.eq.0. ) then
c                coordinates not adjusted, use a priori value
                 read(line(28:43),'(f16.10)') radius	
cd		 print *,'read rad only ',radius       
	       endif
            endif
c           units of radius and height are kilometers
c           convert to geodetic
cd             print *,'lat rad semi,finv ',clatd,radius,semi,finv
            call geoc_to_geod( clatd,radius,glatd,ht,semi,finv )
            lat = glatd*pi/180d0     
cd            print *,'radius clatd,glatd,ht,semi,finv '
cd    .             , radius,clatd,glatd,ht,semi,finv
          endif
        enddo                  
c       now get the zenith delay values
        found = .false.     
        do while (.not.found)
         read(iuztd,'(a)',iostat=ioerr) line  
c         print*,line
         if( ioerr.eq.-1 ) then
            print *,'ATM_ZEN not found on o-file ',ioerr
            stop
         elseif( ioerr.ne.0 ) then
            print *,'Error reading ATM_ZEN from o-file ',ioerr
            stop
         elseif(line(1:9).eq.'ATM_ZEN X'.and.line(11:14).eq.siteuc) then
            found =.true.
         endif
        enddo                                
        endatm = .false.       
        numatm = 0                         
        do while (.not.endatm )  
          site1 = ' ' 
          read(line,'(10x,a4,4x,i4,4i3,f8.0,4x,f8.0,f10.4)'
     .       ,iostat=ioerr) site1,iyr,mon,iday,ihr,min,adj,sig,est
          if( ioerr.ne.0 ) then
            print *,'Error decoding line of o-file ',ioerr
            stop     
          endif   
          call uppers(site1)
          if( site1.ne.siteuc ) then
             endatm = .true.   
          else
            numatm = numatm + 1  
            if( numatm.gt.maxatm ) then
              print *,'Number of o-file paramters exceeds maxatm ='
     ,                ,maxatm
              stop
            endif
            jdoy = idoy(iyr,mon,iday) 
            sod = ihr*3600. + min*60.          
c           time tag is decimal years
            tatm(numatm) = decyrs( iyr,jdoy,sod )   
            estatm(numatm) = est*1000.d0
            adjatm(numatm) = adj*1000.d0
            sigatm(numatm) = sig*1000.d0    
            read(iuztd,'(a)',iostat=ioerr ) line
            if( ioerr.eq.-1 ) then
              endatm =.true.
            elseif (ioerr.ne.0 ) then
              print *,'Error reading ATM_ZEN from o-file ',ioerr
              stop 
            endif
          endif
        enddo

c       now get the gradient values, if they were estimated
c       NS Gradient
        found = .false.     
        eof = .false.
        do while ( .not.found .and. .not.eof )
         read(iuztd,'(a)',iostat=ioerr) line  
c         print*,ioerr,line
         if( ioerr.eq.-1 ) then
            print *,'ATM_GRAD not found on o-file ',ioerr
            eof = .true.
         elseif( ioerr.ne.0 ) then
            print *,'Error reading ATM_GRAD from o-file ',ioerr
            stop
         elseif(line(1:9).eq.'NS_GRAD X'.and.line(11:14).eq.siteuc) then
            found =.true.
         endif
        enddo                                
        endgrad = .false.   
        numgrad = 0    
        do while (found .and. .not.endgrad )
          site1 = ' ' 
          read(line,'(7x,3x,a4,1x,2x,2x,i4,4i3,f8.4,4x,f8.4,2x,f8.4)'
     .       ,iostat=ioerr) site1,iyr,mon,iday,ihr,min,adj,sig,est
          if( ioerr.ne.0 ) then
            print *,'Error decoding line of o-file ',ioerr
            stop     
          endif   
          call uppers(site1)
          if( site1.ne.siteuc ) then
             endgrad = .true.   
          else
            numgrad = numgrad + 1  
            if( numgrad.gt.maxatm ) then
              print *,'Number of o-file grad paramters exceeds maxatm ='
     ,                ,maxatm
              stop
            endif
            jdoy = idoy(iyr,mon,iday) 
            sod = ihr*3600. + min*60.          
c           time tag is decimal years
            tgrad(numgrad) = decyrs( iyr,jdoy,sod )   
            estnsgrad(numgrad) = est*1000.d0
c           adjnsgrad(numgrad) = adj*1000.d0
            signsgrad(numgrad) = sig*1000.d0
            read(iuztd,'(a)',iostat=ioerr ) line
            if( ioerr.eq.-1 ) then
              endgrad =.true.
            elseif (ioerr.ne.0 ) then
              print *,'Error reading NS_GRAD from o-file ',ioerr
              stop 
            endif
          endif
        enddo

c       EW Gradient
        found = .false.     
        do while ( .not.found .and. .not.eof )
         read(iuztd,'(a)',iostat=ioerr) line  
c         print*,line
         if( ioerr.eq.-1 ) then
            print *,'ATM_GRAD not found on o-file ',ioerr
c           stop
         elseif( ioerr.ne.0 ) then
            print *,'Error reading ATM_GRAD from o-file ',ioerr
            stop
         elseif(line(1:9).eq.'EW_GRAD X'.and.line(11:14).eq.siteuc) then
            found =.true.
         endif
        enddo                                
        endgrad = .false.       
        numgrad = 0                         
        do while (found .and. .not.endgrad )
          site1 = ' ' 
          read(line,'(7x,3x,a4,1x,2x,2x,i4,4i3,f8.4,4x,f8.4,2x,f8.4)'
     .       ,iostat=ioerr) site1,iyr,mon,iday,ihr,min,adj,sig,est
          if( ioerr.ne.0 ) then
            print *,'Error decoding line of o-file ',ioerr
            stop     
          endif   
          call uppers(site1)
          if( site1.ne.siteuc ) then
             endgrad = .true.   
          else
            numgrad = numgrad + 1  
            if( numgrad.gt.maxatm ) then
              print *,'Number of o-file grad paramters exceeds maxatm ='
     ,                ,maxatm
              stop
            endif
c           jdoy = idoy(iyr,mon,iday) 
c           sod = ihr*3600. + min*60.          
c           time tag is decimal years
c           tgrad(numgrad) = decyrs( iyr,jdoy,sod )   
            estewgrad(numgrad) = est*1000.d0
c           adjewgrad(numgrad) = adj*1000.d0
            sigewgrad(numgrad) = sig*1000.d0
            read(iuztd,'(a)',iostat=ioerr ) line
            if( ioerr.eq.-1 ) then
              endgrad =.true.
            elseif (ioerr.ne.0 ) then
              print *,'Error reading EW_GRAD from o-file ',ioerr
              stop 
            endif
          endif
        enddo

      elseif( filetype.eq.'sinex' ) then   

c       first look for site coordinates get get latitude and height
        found = .false.     
        do while (.not.found)
          read(iuztd,'(a)',iostat=ioerr) line  
          if( ioerr.eq.-1 ) then
             print *,'Site coordinates not found on SINEX file ',ioerr
             stop
          elseif( ioerr.ne.0 ) then
             print *,'Error reading coordinates from SINEX file ',ioerr
             stop
          elseif(line(19:23).eq.'STA_X') then   
            read(iuztd,'(a)',iostat=ioerr) line  
            if( ioerr.ne.0 ) then
              print *,'Error reading site coordinates from SINEX file'
              stop
            else
              read(line(2:5),'(a4)',iostat=ioerr) site1    
              if( ioerr.ne.0 ) then
                print *,'Error decoding site name for coordinates'
                stop
              elseif( site1.ne.siteuc ) then
                print *,'Site ID on SINEX file (',site1
     .             ,') does not match ',siteuc
                stop
              else  
                found = .true.
                read(line(16:54),'(3(1x,f12.3))',iostat=ioerr) 
     .              (xyz_pos(i),i=1,3)
                if( ioerr.ne.0 ) then
                  print *,'Error reading coordinates from SINEX file'
                  stop
                else    
c*debug           print *,'xyz_pos ',xyz_pos
                  call xyz_to_geod(rot_mat,xyz_pos,geod_pos)
                  lat = geod_pos(1)
c                 convert to km for zhd calculation
                  ht  = geod_pos(3)/1.d3
                endif
              endif
            endif
          endif
        enddo
c       end extraction of site coordinates, now get atmospheric parameters
        found = .false.
        do while ( .not.found )
        read(iuztd,'(a)',iostat=ioerr) line 
          if( ioerr.eq.-1 ) then
             print *,'TROP SOLUTION not found on SINEX file',ioerr
             stop
          elseif( ioerr.ne.0 ) then
             print *,'Error reading TPOP SOLUTION from SINEX file',ioerr
             stop
          elseif(line(2:5).eq.'TROP'.and. line(7:14).eq.'SOLUTION') then   
             found =.true.      
             numatm = 0     
             endatm = .false.
             do while (.not.endatm) 
               read(iuztd,'(a)',iostat=ioerr) line  
               if( line(1:5).eq.'-TROP' ) then
                 endatm = .true.
               elseif( ioerr.eq.-1 ) then
                 print *,'EOF before -TROP on SINEX file'
                 stop
               elseif( ioerr.ne.0 ) then
                 print *,'Error reading tropo values from SINEX file ',ioerr
                 stop
               elseif(line(2:5).eq.siteuc ) then
                 read(line(7:30),'(i2,1x,i3,1x,i5,1x,f6.1,1x,f4.1)'
     .             ,iostat=ioerr) jyr,jdoy,jsod,est,sig   
                 call fix_y2k(jyr)   
                 if( jdoy.eq.doy ) then
                   numatm = numatm+1  
                   if( numatm.gt.maxatm ) then
                     print *,'Number of trop values exceeds maxatm '
     .                      ,maxatm
                     stop
                   else 
c                    time tag is decimal years    
                     sod = dfloat(jsod ) 
                     tatm(numatm) = decyrs( jyr,jdoy,sod ) 
c debug              print *,'numatm tatm ',numatm,tatm(numatm)  
                     estatm(numatm) = est
                     sigatm(numatm) = sig
                   endif
                 endif
               endif
             enddo
          endif
        enddo

c     end of branch on file type
      endif

      return
      end

c------------------------------------------------------------------------------

 
      Subroutine get_interval( sampling,interval,tstart,tend,numatm
     .                       , iuz,dt )

c     Set the interval from either the estimates, observations, or input

c     Input
c        sampling    : 'est'  'obs' or 'inp'
c        interval    : input interval in integer seconds, or 0 if 'est' or 'obs'
c        tatm(maxatm): epochs of estimates in decimal years
c        numatm      : number of estimates
c        iuz         : unit number of z-file for reading observation times

c     Output
c        dt          : interval in decimal years

      implicit none
                        
      character*3 sampling 
      character*256 line 

      integer*4 numatm,iuz,iyr,jdoy,ihr,min,interval,ioerr

      real*8 tstart,tend,sod,sec,decyrs,zt1,zt2,dt,dtsec,yrdays 

      logical leapyr,found                   

c Set days in year for screen echo of interval 
       
      
      if( leapyr(idint(tstart)) ) then
        yrdays=366.d0
      else
        yrdays=365.d0 
      endif
                        
               
c The start time is always the beginning of the estimates but
c epochs without observations or met data will be skipped in output

      if( sampling.eq.'est' ) then
        dt = (tend - tstart)/dfloat(numatm-1)   
c        print *,'DEBUG tstart tend dt yrs ',tstart,tend,dt
        dtsec = dnint(dt*86400d0*yrdays)
        write(*,'(a,f6.1,a)') 
     .   ' Sampling interval from ZTD-file is ',dtsec,' s'
 
      elseif( sampling.eq.'obs') then

c       read the first two epochs of the z-file to find the observation interval
        found = .false.
        if( iuz.gt.0 )  then
          do while( .not.found )
            read(iuz,'(a)',iostat=ioerr ) line 
            if( ioerr.ne.0 ) then      
              print *,'Error reading first epoch of z-file',ioerr
              stop
            elseif( line(1:1).eq.' ' ) then
              read(line,'(1x,2i4,2i3,f4.0,i4,1x,2f9.4,2x,f8.1,f7.1,f7.1
     .                   ,8f11.4)',iostat=ioerr) iyr,jdoy,ihr,min,sec
              if( ioerr.ne.0 ) then 
                print *,'Error decoding z-file line (ioerr) ',ioerr
                print *, line
                stop
              else               
                sod = ihr*3600. + min*60. + sec          
c               time tag is decimal years
                zt1 = decyrs( iyr,jdoy,sod )
                found = .true. 
              endif
            endif
          enddo
c         now find the second epoch
          found = .false.
          do while ( .not.found )
            read(iuz,'(1x,2i4,2i3,f4.0,i4,1x,2f9.4,2x,f8.1,f7.1,f7.1
     .              ,8f11.4)',iostat=ioerr) iyr,jdoy,ihr,min,sec
            if( ioerr.ne.0 ) then 
               print *,'Error reading 2d epoch on z-file(ioerr) ',ioerr
               stop
             endif  
             sod = ihr*3600. + min*60. + sec         
             zt2 = decyrs( iyr,jdoy,sod )   
c            slop is 10 seconds
             if( (zt2-zt1).gt.3.17d-7 ) then 
                dt= zt2-zt1    
                found = .true.      
                dtsec = dnint(dt*86400d0*yrdays)
                write(*,'(a,f6.1,a)') 
     .             ' Sampling interval from z-file is ',dtsec,' s'
             endif  
          enddo 
        else 
          print *,'Sampling is by GAMIT observation epoch but no '
     .           ,'z-file input'
          stop
        endif    
        rewind( iuz ) 
    

      elseif( sampling.eq.'inp') then
        dt = dfloat(interval)/86400.d0/yrdays

      else
        print *,'Invalid sampling input: ',sampling
        stop
      endif
     
      return
      end

c----------------------------------------------------------------------------
       
      Subroutine get_estimates(t,maxatm,numatm,tatm,estatm,adjatm,sigatm
     .        , numgrad,tgrad,estnsgrad,signsgrad,estewgrad,sigewgrad
     .        , estztd,adjztd,sigztd,ewgrad,sigew,nsgrad,signs )   

       
c Interpolate the eatimated values to the current epoch

      implicit none

      integer*4 maxatm,numatm,numgrad,index,i

      real*8 t,tatm(maxatm),estatm(maxatm),adjatm(maxatm),sigatm(maxatm)
     .     , tgrad(maxatm),estnsgrad(maxatm),signsgrad(maxatm)
     .     , estewgrad(maxatm),sigewgrad(maxatm)
     .     , estztd,adjztd,sigztd,ewgrad,sigew,nsgrad,signs
     .     , td,ff,t1,t2,var


c Assume that the estimated zenith delays are knots of a linear spline 
c and interpolate to the current epoch     
                       
      if( numatm.le.1 ) then
        write(*,'(a)') 
     .    'Metutil not coded for 0 or 1 zenith delay parameter '
        stop
      endif
      index = 0
      do i=1,numatm-1   
c       add 10s of slope at beginning and end to avoid losing data
        t1 = tatm(i)
        t2 = tatm(i+1)
        if( i.eq.1 ) then
          t1 = tatm(1)-3.17d-7   
        elseif( i.eq.numatm-1) then
          t2 = tatm(i+1)+3.17d-7
        endif
        if( t.ge.t1.and.t.le.t2 ) then 
          index = i
          goto 10
        endif
      enddo  
   10 if( index.eq.0 ) then
        print *,'Observation epoch ',t,' not in o-file table interval ',
     .           tatm(1),tatm(numatm) 
        estztd = 0.d0
        adjztd = 0.d0
        sigztd = 0.d0
      else
        ff = (t-tatm(index))/(tatm(index+1)-tatm(index)) 
        estztd = estatm(index) + ff*(estatm(index+1)-estatm(index))
        adjztd = adjatm(index) + ff*(adjatm(index+1)-adjatm(index))
c       assume the variance is linear between knots (approximation? ) 
        var= sigatm(index)**2 + ff*(sigatm(index+1)**2-sigatm(index)**2)
        sigztd = dsqrt(var)
      endif      

c Same for gradients

      if( numgrad.eq.1 ) then
        nsgrad = estnsgrad(1)
        signs = signsgrad(1)
        ewgrad = estewgrad(1)
        sigew = sigewgrad(1)
        return
      endif

      index = 0
      do i=1,numgrad-1   
c       add 10s of slope at beginning and end to avoid losing data
        t1 = tgrad(i)
        t2 = tgrad(i+1)
        if( i.eq.1 ) then
          t1 = tgrad(1)-3.17d-7   
        elseif( i.eq.numgrad-1) then
          t2 = tgrad(i+1)+3.17d-7
        endif
        if( t.ge.t1.and.t.le.t2 ) then 
          index = i
          goto 20
        endif
      enddo  
   20 if( index.eq.0 ) then
        print *,'Observation epoch ',t,' not in o-file table interval ',
     .           tgrad(1),tgrad(numatm) 
        nsgrad = 0.d0
        signs = 0.d0
        ewgrad= 0.d0
        sigew = 0.d0 
      else
        ff = (t-tgrad(index))/(tgrad(index+1)-tgrad(index)) 
        nsgrad=estnsgrad(index)+ff*(estnsgrad(index+1)-estnsgrad(index))
        ewgrad=estewgrad(index)+ff*(estewgrad(index+1)-estewgrad(index))
c       assume the variance is linear between knots (approximation? ) 
        var=  signsgrad(index)**2 
     .           + ff*(signsgrad(index+1)**2-signsgrad(index)**2)
        signs = dsqrt(var)
        var= sigewgrad(index)**2 
     .           + ff*(sigewgrad(index+1)**2-sigewgrad(index)**2) 
        sigew = dsqrt(var)
      endif      

      return
      end

c------------------------------------------------------------------------------

      Subroutine get_zhd ( t,iuz,iurnx,lat,ht,zhd,pres,temp,ztd )

c     Get the zenith hydrostatic delay from the GAMIT (MODEL) z-file or from
c     pressure values from a RINEX met file.  If zhd returned as 0., values 
c     are either non sensical  or too far from requested epoch.

      implicit none

c        passed arguments  
      integer*4 iuz,iurnx
      real*8   t,lat,ht,zhd,pres,temp,ztd
c        stored values from RINEX or z-file
      integer*4 maxmet,ntmet  
      parameter(maxmet=10000)
      real*8 timr(maxmet),presr(maxmet),tempr(maxmet)
     .     ,  zhdr(maxmet),ztdr(maxmet) 
c        other variables  
c   debug
c      integer*4 i
      real*8  exttim,fn
      logical first_call   

      save timr,presr,tempr,zhdr

c     warn if a large time difference between requested and available (set to 0.5 hr)
      data first_call/.true./
          

c Read the net data into storage
                                
      if( first_call ) then
        if( iurnx.gt.0 ) then
          call read_rnxfile( iurnx,maxmet,ntmet,timr,presr,tempr)
        elseif( iuz.gt.0 ) then
          call read_zfile( iuz,maxmet,ntmet,timr,presr,tempr,zhdr,ztdr )
        else
          print *,'Cannot get hydrostatic delay, iuz=iurnx=0'
          stop
        endif                    
        first_call = .false.                                             
      endif
                                                        
c Initialize to default values

      zhd = 0.d0
      pres = 0.d0 
      temp = 99.d0
      ztd = 0.d0

c Interpolate the arrays to get the values
                           
c      print *,'timr ',(timr(i),i=1,ntmet)
c      print *,'presr ',(presr(i),i=1,ntmet)   
c      print *,'tempr ',(tempr(i),i=1,ntmet)
c      print *,'zhdr ',(zhdr(i),i=1,ntmet)
c      print *,'ztdr ',(ztdr(i),i=1,ntmet)
      if( iurnx.gt.0 ) then
        call lininterp( maxmet,ntmet,timr,presr,t,pres,exttim )
        call lininterp( maxmet,ntmet,timr,tempr,t,temp,exttim )
      elseif( iuz.gt.0 ) then
        call lininterp( maxmet,ntmet,timr,presr,t,pres,exttim )
        call lininterp( maxmet,ntmet,timr,tempr,t,temp,exttim )
        call lininterp( maxmet,ntmet,timr,zhdr,t,zhd,exttim )
        call lininterp( maxmet,ntmet,timr,ztdr,t,ztd,exttim )
      endif       
c      print *,'t ntmet exttim pres zhd ',t,ntmet,exttim,pres,zhd
c     --if time too far off or pressure or zhd  nonsensical, do not use
      if( exttim.lt.0. ) then
        write(*,'(a,f14.8,a,f14.8)') ' Obs time ',t
     .       ,' too early for met file start ',timr(1)
        zhd = 0.d0
      elseif( exttim.gt.0. ) then
        write(*,'(a,f14.8,a,f14.8)') ' Obs time ',t
     .        ,' too late for met file end ',timr(ntmet)   
        zhd = 0.d0
      elseif( iuz.gt.0 ) then 
        if( zhd.lt.1000.d0.or.zhd.gt.2500.d0 ) then
          print *,'Invalid ZHD on z-file ',t,zhd,', skip ' 
          zhd = 0.d0
        else 
          continue
        endif
      elseif( iurnx.gt.0 ) then
        if( pres.gt.1500.d0.or.pres.lt.300.d0 ) then
          print *,'Invalid P on RINEX met file '
     .           ,t,presr,' cannot compute ZHD'    
          zhd = 0.d0
        else     
c         From Bevis et al (1992) quoting Elgered et al (1991)   
          fn = 1.d0 - 0.00266 * dcos (2.d0*lat) - 0.00028 * ht 
          zhd = 2.2779d0 * pres/fn          
        endif
      endif                     
c      print *,'returning zhd ',zhd
      return
      end                                                                 

c-------------------------------------------------------------------------

      Subroutine read_rnxfile( iurnx,maxmet,ntmet,times,pres,temp )

c Subroutine to read pressure and temperature values from RINEX met file  
c Adapted from model/read_metrnx.  R. King 26 January 2007
          
      implicit none
       
c  Input
c   iurnx   unit number for RINEX met file
      integer*4 iurnx
c   maxmet ntmet   dimensions and number of time and values arrays
      integer*4 maxmet,ntmet
c   times    epochs of values values from RINEX file (decimal years)
      real*8 times(maxmet)           
c   pres, temp   values from RINEX file (hPa, degK)
      real*8 pres(maxmet),temp(maxmet)

c  Other variables                                      
      integer*4 yr,mon,day,hr,min,jdoy,i,j,ii,ioerr,nobs  
      real*4 rxver 
      real*8 val(3),sec,sensor_ht,sod
      character*2 obstyp(3)
      character*256 line
      logical eof,found_type,found_sensor

c Functions
      integer*4 idoy 
      real*8 decyrs
               
c  Initialization
      line = ' ' 
c      # of epochs, incremented as the file is read
c      stored with times and values in common /ufcom/ in ../includes/model.h
      ntmet = 0 

            
c  Read the header (order of entries is arbitrary except version number)
                                         
      found_type = .false.
      found_sensor = .false.
      read(iurnx,'(a)',iostat=ioerr) line  
      if( line(61:73).ne.'RINEX VERSION' ) then
         print *,'Bogus first line of RINEX met file',ioerr
         stop
      else
c          If the version number is an integer written under version 1 format
c          ( I6 ) it will not be read correctly by the version 2 format ( f9.2)
c          so read the version number free-format
         read( line(1:9),*,iostat=ioerr) rxver  
         if( ioerr.ne.0 ) then
           print *,'Error reading RINEX file version number ',ioerr
           stop   
         endif
      endif
      do while(line(61:73).ne.'END OF HEADER')
        read(iurnx,'(a)',iostat=ioerr) line  
        if( ioerr.eq.-1 ) then
          print *,'EOF before END OF HEADER in RINEX met file ' 
          stop
        elseif( ioerr.ne.0 ) then
          print *,'Error reading RINEX met file header',ioerr 
          stop
        else 
          if( line(61:69).eq.'# / TYPES' ) then
            read(line(1:24),'(i6,3(4x,a2))',iostat=ioerr) 
     .           nobs,(obstyp(i),i=1,3)  
            if( ioerr.ne.0 )  then
               print *,'Error reading the observation order',ioerr 
               stop
            else
              found_type = .true.
            endif
          elseif( line(58:70).eq.'PR SENSOR POS' ) then
            read(line(43:56),'(f14.4)',iostat=ioerr) sensor_ht
            if( ioerr.ne.0 ) then  
              print *,'Error reading the pressure sensor position',ioerr 
              stop
            else
              found_sensor = .true.
            endif
          endif
        endif
      enddo
      if( .not.found_type ) then
        print *,'Missing # / TYPE line in met rinex file',ioerr
        stop
      endif
c      print *,'height of pressure sensor is',sensor_ht 
      

c  now we just read all the data through to the end of the file. Format is:
c  05  5 26  0  0  0  969.3   -8.3   82.7
c  05  5 26  0  0 30  969.3   -8.4   82.7 
             
      ntmet = 0 
      eof =.false.
      do while( .not.eof ) 
        read(iurnx,'(5i3,f3.0,3f7.1)',iostat=ioerr) 
     .              yr,mon,day,hr,min,sec,(val(j),j=1,3)
c        print *,'Read RNX yr mon day hr min sec val '
c     .         , yr,mon,day,hr,min,sec,(val(j),j=1,3)
        if( ioerr.eq.-1 .or.(mon.eq.0)) then
          eof = .true.
        elseif( ioerr.ne.0 ) then  
          print *,'Error reading values from RINEX met file',ioerr
          stop  
        elseif( ntmet+1.gt.maxmet ) then                           
          print *,'Number of met values on RINEX file > maxmet =',maxmet
          stop   
        else
          ntmet = ntmet + 1
c         convert the year, month,day, hr, min, sec into decimal year   
          call fix_y2k(yr)
          jdoy = idoy(yr,mon,day)
          sod = hr*3600. + min*60. + sec
c         time tag is decimal years 
          times(ntmet)= decyrs( yr,jdoy,sod )  
c         Now, assign the values as indicated by the headers
          if(obstyp(1).eq."PR") pres(ntmet)=val(1)
          if(obstyp(1).eq."TD") temp(ntmet)=val(1) + 273.16d0
          if(obstyp(2).eq."PR") pres(ntmet)=val(2)
          if(obstyp(2).eq."TD") temp(ntmet)=val(2) + 273.16d0
          if(obstyp(3).eq."PR") pres(ntmet)=val(3)
          if(obstyp(3).eq."TD") temp(ntmet)=val(3) + 273.16d0
c          print *,'ntmet time pres ',ntmet,times(ntmet),pres(ntmet)
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
        print *,'Number of met values on RINEX file > maxmet =',maxmet
        stop
      endif
      do i=1,ntmet-2
        ii = ntmet-i
        times(ii) = times(ii-1)  
        pres(ii) = pres(ii-1)
        temp(ii) = temp(ii-1)  
      enddo
      times(1) = times(2) - 3.d0/24.d0/365.25d0        
      times(ntmet) = times(ntmet-1) + 3.d0/24.d0/365.25d0    
      pres(1) = pres(2)
      temp(1) = temp(2)
      pres(ntmet) = pres(ntmet-1)
      temp(ntmet) = temp(ntmet-1)

      return
      end

c---------------------------------------------------------------------------------

      Subroutine read_zfile( iuz,maxmet,ntmet,times,pres,temp,zhd,ztd )

c Subroutine to read pressure and temperature values from RINEX met file  
c Adapted from model/read_metrnx.  R. King 26 January 2007
          
c Number, times, and values of met data stored in common /ufcom/ in model.h
c Arrays used for values either from a RINEX met file or those interpolated
c from a global grid and stored on the u-file

      implicit none
       
c  Input
c   iuz  unit number for z-file file
      integer*4 iuz
c   maxmet  ntmet dimensions and number of time and values arrays
      integer*4 maxmet,ntmet
c   times    epochs of values values from z-file (decimal years)
      real*8 times(maxmet)           
c   pres, temp, zhd, ztd values from z-file (hPa, degK, m)
      real*8 pres(maxmet),temp(maxmet),zhd(maxmet),ztd(maxmet)
                               
c  Other variables       
      integer*4 yr,hr,min,jdoy,prn,i,ii,ioerr  
      real*8 pres1,temp1,zhd1,zwd1,ztd1,azim,elev,e,sec,sod,t,tlast
      character*256 line
      logical eof,found    

c Functions
      real*8 decyrs
               
c  Initialization
      line = ' ' 
c      # of epochs, incremented as the file is read
      ntmet = 0 
            
c  Read the header 

      found = .false.
      do while( .not.found )
        read(iuz,'(a)',iostat=ioerr ) line 
        if( ioerr.ne.0 ) then      
           print *,'Error reading first epoch of z-file',ioerr
           stop
        elseif( line(1:1).eq.' ' ) then     
          found = .true.
          backspace(iuz)
        endif 
      enddo

c  now we just read all the data through to the end of the file. Format is: 
c* Yr  Doy Hr Mn Sec  PRN  Azimuth  Elevation  Pres   Temp    WV Pres  Dry Zen (m) Wet Zen   Total Zen ...
c 2000 220  0  0  0.   4   49.2475  54.3248     927.3  296.1   14.0     2.1122     0.1368     2.2490   ...
c  save only one value per epoch                  

      ntmet = 0 
      eof =.false.     
      tlast = 0.d0
      do while( .not.eof )    
       read(iuz,'(1x,2i4,2i3,f4.0,i4,1x,2f9.4,2x,f8.1,f7.1,f7.1,8f11.4)'
     .   ,iostat=ioerr) yr,jdoy,hr,min,sec,prn,azim,elev,pres1,temp1
     .                   ,e,zhd1,zwd1,ztd1  
c        print *,'READ Z t p zhd ztd ',yr,jdoy,hr,min,sec,pres1,zhd1,ztd1
        if( ioerr.eq.-1 ) then
          eof = .true.    
        elseif (ioerr.ne.0 ) then
          print *,'Error reading data from z-file, ioerr= ',ioerr
          stop   
        else   
          call fix_y2k(yr)  
          sod = hr*3600. + min*60. + sec         
c         time tag is decimal years  
          t =  decyrs( yr,jdoy,sod )
c debug
c         write(*,'(a,2f16.8)' ) 't tlast ',t,tlast
          if( t.gt.tlast ) then 
            ntmet = ntmet + 1   
            if( ntmet.gt.maxmet ) then
               print *,'Number of met values on z-file > maxmet ',maxmet
               stop
            endif
            times(ntmet) = t 
            pres(ntmet)=pres1
            temp(ntmet)=temp1
            zhd(ntmet) = zhd1*1.d3
            ztd(ntmet) = ztd1*1.d3   
            tlast = t
          endif
        endif
      enddo           
c      write(*,'(a,i5,2880f16.8)') 'ntmet times zhd '
c     .     ,ntmet,(times(i),zhd(i),i=1,ntmet)

c  Augment the values at the beginning and the end to minimize the
c  situation in which we need to extrapolate or substitute nominal
c  values.  Three hours seems reasonably safe and will produce
c  less discontinuity than reverting to the nominal value, which
c  might be the RINEX start, midpoint, or end, or a value from 
c  GPT.   rwk 070110
      
      ntmet = ntmet + 2   
      if (ntmet.gt.maxmet ) then    
        print *,'Number of met values on z-file > maxmet =',maxmet 
        stop
      endif     
      do i=1,ntmet-2
        ii = ntmet-i
        times(ii) = times(ii-1)  
        pres(ii) = pres(ii-1)
        temp(ii) = temp(ii-1)
        zhd(ii) = zhd(ii-1)
        ztd(ii) = ztd(ii-1)
      enddo
      times(1) = times(2) - 3.d0/24.d0/365.25d0
      pres(1) = pres(2)
      temp(1) = temp(2) 
      zhd(1) = zhd(2)
      ztd(1) = ztd(2)
      times(ntmet) = times(ntmet-1) + 3.d0/24.d0/365.25d0
      pres(ntmet) = pres(ntmet-1)
      temp(ntmet) = temp(ntmet-1)
      zhd(ntmet) = zhd(ntmet-1)
      ztd(ntmet) = ztd(ntmet-1)                
      
c      write(*,'(a,i5,2880f16.8)') 'ntmet times zhd '
c     .     ,ntmet,(times(i),zhd(i),i=1,ntmet)

      return
      end

c------------------------------------------------------------------------------------------

      Subroutine LININTERP( ndim,ntab,val_time,valtab,obstime
     .                    , valinterp,exttim )

c  Performs a linear interpolation to calculate the value
c  at a particular time. Modified model/linterp.f to interpolate
c  a single variable. 
c
c  INPUT:
c    ndim        : dimensions of the array in the calling program 
c    ntab        : number of tabulated values to be interpolated 
c    val_time    : array of time tags of tabulated values
c    valtab      : array of corresponding tabulated values
c    obstime     : time to which the tabulated values should be
c                  interpolated

c
c  OUTPUT:
c    valinterp   : interpolated value
c     exttim     : if requested interpolation is outside the range
c                  of the table, this variable will have the time
c                  (neg or pos) by which the range is exceeded;
c                  if within range, extflg = 0.
c                  
c
c  P. Tregoning  8 January 2004; R. King August 2018 

      implicit none

      integer*4 ndim,ntab,t1,t2,i

      real*8 val_time(ndim),valtab(ndim),exttim,valinterp
      real*8 obstime,dval,dt,tstep,delta,eps

      logical found,debug/.false./

      if( debug) then 
        print *,'METUTIL/lininterp'
        print*,'ndim,obstime',ndim,obstime
        print*,'val_time',val_time
        print*,'ntab',ntab
        print*,'valtab: ',(valtab(i),i=1,ntab)
      endif
                                    
c     eps is a small number tolerance to allow use of the end values
c     currently set to be 90s, with times in years
      data eps/2.85d-6/
                  
c      if( obstime.gt.220.082 ) debug = .true.
c      if( obstime.gt.220.0826 ) debug = .false.

c   see if the requested time is outside the range of the array   
      exttim = 0.        
      if( debug )
     .   print *,'obstime nrow val_time(n) eps '
     .     ,obstime,ntab,val_time(ntab),eps
      if( obstime.lt.(val_time(1)-eps) )  then
        exttim = obstime - val_time(1)  
        valinterp = valtab(i)
        return
      elseif( obstime.gt.(val_time(ntab)+eps) ) then
        exttim =  obstime - val_time(ntab) 
        if(debug)  print *,'exttim ',exttim 
        valinterp = valtab(i)
        return
      endif

c   run through the time array to find which two times the requested
c   time lies between
      found = .false.
      t1 = 0
      do while (.not.found)
        t1 = t1+1
        if(t1+1.gt.ntab)then
          write(*,'(a,f15.7,a,2f15.7)')
     .     'Requested time',obstime,' outside tabulated values. Min/max'
     .    ,val_time(1),val_time(ntab) 
          stop
        endif
c PT050415: there can be a roundoff problem if obstime is infinitessimally
c           less that the first val_time(t1) .... ! 
c RWK070111: Modified to allow interpolation within eps of last value
        if((dabs(obstime-val_time(t1)).lt.eps.or.
c     .   obstime.ge.val_time(t1)).and.obstime.lt.val_time(t1+1))then
     .    obstime.ge.val_time(t1)).and.obstime.lt.(val_time(t1+1)+eps)) 
     .          then
          found = .true.
        endif    
      enddo

      t2 = t1 + 1

c  now linearly interpolate between these two times
      dt = obstime-val_time(t1)
      tstep = val_time(t2)-val_time(t1)
      dval = valtab(t2)-valtab(t1)
      delta = dval * dt / tstep
      valinterp = valtab(t1)+delta
      if(debug) then
        write(*,'(a,2i3,4f16.8)') 
     .       't1 t2 dt tstep dv delta ',t1,t2,dt,tstep,delta
        write(*,'(a,3f10.2)') 
     .     'val1 val2 val ',valtab(t1),valtab(t2),valinterp
      endif

      return
      end



