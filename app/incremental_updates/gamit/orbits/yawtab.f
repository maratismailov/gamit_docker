      Program YAWTAB

c  Program to read a t-file and svnav.dat file (for maximum yaw rate)
c  and create a binary file of yaw attitude values for all satellites 
c  at all epochs using a modified version of Jan Kouba's  eclips_dec2013.f.
c  The old Bar-Sever code implement by Paul Tregoning has  been removed 
c  and yawtab no longer reads the ascii yawfile created by ARC. MODEL 
c  will read the table created by yawtab to to get the attitude of each 
c  satellite for each epoch.  The binary y-file can be  converted to ASCII 
c  for checking using program ytoasc.  When yawtab runs, it now creates a 
c  short yawtab.out file giving for each satellite at the midpoint epoch: prn, 
c  snv, antbody, beta, and sun-wrt-svnode.
c
c  P. Tregoning  26 November  1997
c  Mods:  R. King 8 January   1998   
c  Mods:  P. Tregoning 27 Feb 1998 
c  Mods:  R. King January     2014
c  Mods:  M. Floyd   14 April 2020

c  Calling sequence:  yawtab <in-tfile> <out-yfile> <tab_interval>
c
c                e.g. yawtab  tvent7.278 yventb.278 30 
c

      implicit none 

      include '../includes/dimpar.h' 
      include '../includes/orbits.h'
      include '../includes/units.h'  
      include '../includes/global.h'
      include '../includes/arc.h'
       
c Unit numbers and file names not in arc.h 
      integer*4 iyawb,iyawp,iout
      character*16 yfile,nutfil,ut1fil,polfil,print_file
            
c Quantities read or derived from the t-file
      integer*4 nsat,itsat(maxsat),nepcht,nintrs
      real*8 start_time,end_time,sdelt

c     jdb tb  jdf tf : start/stop times of t-file    
c     jde te are the IC epoch of the t-file

c Quantities read or derived from svnav.dat
      integer*4 isvn(maxsat),frqchn
     .        , svnstart(5),svnstop(5)
      real*8 yrate(maxsat),sbmassx,ybias(maxsat)
      real*8 antpwrx   ! Antenna transmit power (W) 
                       ! (Added x as with sbmassx above)
         

      character*1 aybias(maxsat)
      character*20 svantbody(maxsat)

c  Satellite/Sun positions   
      real*8 svec(6),sun(6)

c  Interpolation variables and current coordinates for the satellite
      real*8 ytrp(5,2,maxyt2,maxsat),yy(10,maxytp,maxsat)  
      integer*4 ji0,jil,iy1,iy2,jlast,jnow,iendf

c  Saved array of event flags and sun angle at the start of the event
      integer*4 ievent(maxsat)
      real*8 svbcos_start(maxsat)
       
c Constants used internally
c    pi   : 3.1415926
c    dtr  : pi/180.d0  (deg-to-radians)  

c    tdtoff: TDT - TAI in fraction of day for ephemeris interpolation
c    (in arc.h)
      real*8 pi,dtr
          
c Other variables
      integer*4 i,j,kk,ioerr,isat,iepoch,jds,nepoch,xfile_int
     .        , calc_int,jdstart,jdend,prn,yr,month,doy,day,hr
     .        , min,sec,yr2,iarg,iclarg,nversn
     .        ,nepoch_tab,ivel 
c     temporary for debug
     .       ,week

      real*8 tsatic(maxorb,maxsat),ts,etime,yaw_angle(maxsat)
     .     , tstart,tend,satarg,attit(maxsat),tjdb,tjdf
     .     , interval_factor,tjdmid,fjd 
     .     , beta,betar,betadg,ideal_yaw, sc_x_t(3)
     .     , xsv(3),vsvc(3),xsun(3),svbcos,sunlon
c       returned from lib/yaw_attit but not used
     .     , sc_y_t(3),sc_z(3),sc_W(3)
c     temporary for debug
     .     ,sow

c     temporary
      character*1  asys
      character*2  ayr2 
      character*3  adoy
      character*4  buf4
      character*6 module                       
      character*256 message

      logical first_call(maxsat),debug/.false./,bias_warning(maxsat)
     .      ,newcode/.true./,wang/.false./  
      integer*4 iprndb/0/

c Function to return the difference between two time tags
      real*8 timdif    

c Function to check if a file exists
      logical fcheck 

c Function to return the length of a string
      integer*4 nblen
    
c Function to compute a dot produce
      real*8 dot

c  YAWTAB version variables
      character versn*120,uname*16
      integer*4 iyear,imonth,iday,ihr,imn,isec,ihnsec  


c ---------------------------------------------------------------

c Set module for report_stat calls (may change later)
      module = 'YAWTAB'
c     exit if a previous step has failed           
      if( fcheck('GAMIT.fatal') )
     .  call report_stat('FATAL',module,'orbits/yawtab',' '
     .                  ,'GAMIT.fatal exists: YAWTAB not executed',0)
 
      call oversn(versn)  
      write(message,50) versn(1:nblen(versn))
50    format('Program YAWTAB Version ',A)  
      call report_stat('STATUS',module,'orbits/yawtab',' ',message,0)
      call getdat(iyear,imonth,iday)
      call gettim(ihr,imn,isec,ihnsec)
      call getusr(uname)
      write(message,60) iyear,imonth,iday,ihr,imn,isec,uname
60    format(' YAWTAB Run on ',I4,'/',I3,'/',I3,2X,I2,':',I2,':',I2
     .      ,' by ',A16)
      call report_stat('STATUS',module,'orbits/yawtab', ' ',message,0)


c Define variable to correct GPST to TAI for interpolating the sun position
      tdtoff =  (32.184d0 + 19.0d0)/86400.d0

c Define pi and the conversion for degrees to radians
      pi= 4.D0*DATAN(1.D0) 
      dtr = pi/180.d0  
         
c Read the command line
      iarg = iclarg(1,tfname)
      if (iarg.le. 0 )  call report_stat('FATAL',module,'orbits/yawtab'
     .                     ,' ','Missing argument for input t-file',0)
      iarg = iclarg(2,yfile)
      if (iarg.le. 0 )  call report_stat('FATAL',module,'orbits/yawtab'
     .                     ,' ','Missing argument for output y-file',0) 
c**  replace this with read of session.info or x-file
      iarg = iclarg(3,buf4)
      if (iarg.le. 0 )  call report_stat('FATAL',module,'orbits/yawtab'
     .              ,' ','Missing argument for observation interval',0)
      read(buf4,*) xfile_int            
c PT980227: use two intervals - the calculation interval and the output
c           y-file interval (= x-file interval). 30 seconds is probably 
c           adequate for the calculations but if the x-file interval is
c           less then use it for calcs and output     
      if(xfile_int.lt.30)then
        calc_int = xfile_int
        interval_factor = 1.d0
      else
        calc_int = 30  
        interval_factor = xfile_int*1.d0/(calc_int*1.d0)
      endif            
      write(message,'(a,i4,a)')' Yaw Table interval       : ',xfile_int
     .                            ,' seconds'
      call report_stat('STATUS',module,'orbits/yawtab',' ',message,0)
      write(message,'(a,i4,a)')' Yaw calculation interval : ',calc_int
     .                            ,' seconds'  
      call report_stat('STATUS',module,'orbits/yawtab',' ',message,0)


c Define the unit numbers 
      iut=15
      iyawb=17
      iyawp=18 
c     /soltab. (old scheme)
      isun=35                        
c     n-body ephemeris (new scheme)
      ibody=36 
      iprnt=0
      iscrn=6

c Open the t-file
      call lowers(tfname)
      call topens(tfname,'old',iut,ioerr)
      write(message,'(a,a16)' ) ' Ephemeris (T-) File      : ',tfname
      call report_stat('STATUS',module,'orbits/yawtab',' ',message,0)


c Open the output y-file  
      open(iyawb,file=yfile,form='unformatted'
     .     ,access='sequential',status='unknown',iostat=ioerr)
      if (ioerr .ne. 0) then
         call report_stat('FATAL',module,'orbits/yawtab',yfile,
     .     'error opening output yaw table:',ioerr)
      endif 

c Read the t-file header to determine start/stop times of tfile
       frame = 'UNKWN'
      call thdred ( iut,iscrn,iprnt,nsat,gnss,itsat,satnam
     .            , jdb,tb,jdf,tf,sdelt,nepcht,jde,te
     .            , nics,tsatic,nintrs, icsnam
     .            , precmod,nutmod,gravmod,frame,srpmod,eradmod
     .            , antradmod )  

                               
c Get the start, stop, and midpoint epochs for the y-file header and print-file name
* MOD TAH 210202: Add more range (+-sdelt; tablular spacing)
* MOD TAH 210205: Make times consistent with the start of the t-file
*     from arc and consistent with model test.
C     tjdb = jdb*1.d0+(tb+1.5*3600.d0)/86400.d0 
C     tjdf = jdf*1.d0+(tf-1.5*3600.d0)/86400.d0 
*     Code from arc, so make the same. This should generate start
*     stop times that matchs the data start and stop times (in arc.in).
*     Assume dmod(trun0,delt) is zero (Tabular matches data start 
*     and stop (could be an issue with non-standard processing).
*        trun0= trun0 - dmod(trun0,delt) - 3600.d0 - 6*delt
*        trunf= trunf - dmod(trunf,delt) + 3600.d0 + 6*delt  
*     Actual start and stop time are computed below so make this
*     calcualtion consistent.  These values should match the x-file
*     start and stop time exactcly.  Added additional 30 second )calc_int)
*     buffer each end so first and last tabular point should be
*     one 30-sec epoch before and after data.
      tjdb = jdb*1.d0+(tb+3600.d0+6*sdelt-calc_int)/86400.d0 
      tjdf = jdf*1.d0+(tf-3600.d0-5*sdelt+calc_int)/86400.d0  
      tjdmid = (tjdb+tjdf)/2.d0

c Name and open the print file
      call dayjul(int(tjdmid),yr,doy)
      yr2 = mod(yr,100)           
      write(ayr2,'(i2.2)') yr2       
      write(adoy,'(i3.3)') doy
      print_file = 'yawtab.out.'//ayr2//adoy
      open(iyawp,file=print_file,form='formatted',status='unknown'
     .    , iostat=ioerr)

c Open and read the read headers of the ephemeris files
      if( fcheck('nbody') ) then
c       don't need Venus and Jupiter  
        lbody = 0 
c       don't need velocities
        ivel = 0 
c       GPST or UTC ok for determining what records to read into storage
        fjd =  dfloat(jde) + te/86400.d0 
        call ephred(ibody,fjd,lbody,ivel) 
      else                 
        fjd =  dfloat(jde) + te/86400.d0 
        if(debug) print *,'aft ephdrd fjd,fjdbsn,  '
     .                 ,  fjd,fjdbsn  
        call ephdrd(fjd)  
        call evrtcf 
      endif

c Usable part of tfile is 1.5 hours less than this at either end to allow
c for an 11-point interpolator.          
      jdstart = jdb
      tstart = tb
      jdend = jdf
      tend = tf
* MOD TAH 210202: Increased time: Changed 1.5 to 1.0)
* MOD TAH 210205: This seems to duplicate start and stop time calcautions
*     above, so again change to make consistent with the way arc caculates
*     start time from data start and end times.
*     See comments on tjdb tjdf caculation about 1 30-sec (calc_int)  epoch buffer.
C     call timinc(jdstart,tstart,1.0*3600.d0)
      call timinc(jdstart,tstart,3600.d0+6*sdelt-calc_int)
      start_time = jdstart +  tstart/86400.d0 
* MOD TAH 210205: Same change as above t amke consisent with arc (extend
*     one extra tabular point to make sure we could past end of data).
C     call timinc(jdend,tend,-1.0*3600.d0)
      call timinc(jdend,tend,-(3600.d0+5*sdelt)+calc_int)
      end_time = jdend + tend/86400.d0      
cd     print *,'start_time end_time ',start_time,end_time      

c       nepoch =  number of epochs (using calc_int for calculations) which 
c                 cover the span of usable time in t-file
c       nepoch_tab =  number of epochs (using input step size ) which cover 
c                  the span of usable time in t-file
* MOD TAH 210205: Use nint here to make sure we go to nearest inteval.
C     nepoch = (end_time - start_time)*86400.d0/(calc_int*1.d0) 
C     nepoch_tab = (end_time - start_time)*86400.d0/(xfile_int*1.d0) 
      nepoch = nint((end_time - start_time)*86400.d0/calc_int) 
      nepoch_tab = nint((end_time - start_time)*86400.d0/xfile_int) 
 
                           
c  Read svnav.dat to get the the SV antenna/body-type, maximum yaw rates
c  and yaw bias for the satellites at this epoch
c  MOD: Add one day to date for call to svnav_read because the start of the t-
c       file may occur before the start time of a PRN record. MAF (2020-04-14, MIT)
c     call dayjul(jdb,yr,doy)
      call dayjul(jdb+1,yr,doy)
      call monday(doy,month,day,yr)
      hr = dint(tb/3600.d0)
      min = dint(mod(tb/3600.d0,1.d0)*60.d0)
      if(debug) print *,'Fr svnav.dat prn svn antbody yawbias yawrate:'
cd      print *,'nsat itsat ',nsat,itsat
      do isat = 1,nsat                
cd          print *,'isat,itsat(isat) ',isat,itsat(isat)
* MOD TAH 190702: Added antpwrx to snav_read call (different from common version)
          call svnav_read( -1,yr,doy,hr,min,gnss,itsat(isat)
     .                   , isvn(isat),frqchn,svantbody(isat),sbmassx
     .                   , aybias(isat),yrate(isat), antpwrx
     .                   , svnstart,svnstop )
          if(debug) print *,isat,itsat(isat),isvn(isat),svantbody(isat)
     .           ,aybias(isat),yrate(isat)                      
c         Yaw currently coded only for GPS, Glonass, Beidou, and Galileo
          if( svantbody(isat)(1:5).ne.'BLOCK' .and.
     .        svantbody(isat)(1:5).ne.'GLONA' .and.
     .        svantbody(isat)(1:5).ne.'BEIDO' .and.
     .        svantbody(isat)(1:5).ne.'GALIL')  then
            write(message,'(a,a20,a,i2)')  'Yaw not coded for '
     .                     , svantbody(isat),' PRN ',itsat(isat)
            call report_stat('WARNING',module,'orbits/yawtab',' '
     .                       ,message,0)   
          endif
      enddo     
      write(message,'(a,200i3)') 'PRN nos. in channels selected: '
     .                        ,  (itsat(i),i=1,nsat)
      call report_stat('STATUS',module,'orbits/yawtab',' ',message,0)

c Write the binary y-file header
c      set the version number (970 = GAMIT release 9.70, Jan 1998)
       nversn = 970            
c      new version with epoch flags (1051 = GAMIT releast 10.51, Jan 2013)
       nversn = 1051     
c      new version with SV body type rather than block # in header
       nversn = 1061 
       write(iyawb) nversn,tfname,tjdb,tjdf,xfile_int,nepoch_tab,nsat
     .             ,(itsat(i),isvn(i),svantbody(i),i=1,nsat)  
cd       print *,'tb+ tf- ',tb+1.5*3600.d0,tf-1.5*3600.d0

c Write the print-file header 
      write(message,60) IYEAR,IMONTH,IDAY,IHR,IMN,ISEC,UNAME
      write(iyawp,'(a,a16,a,f16.7,a,f16.7,a,i4,a,i6)')    
     .  '  Input T-file: ',tfname,'  Start PJD:'
     .  ,tjdb,'  Stop PJD:', tjdf,'  Interval:',xfile_int
     .  ,'  No. epochs:',nepoch_tab
       write(iyawp,'(a,f16.7)') 
     .   'Max yaw rate, yaw-bias, beta and sun-long at mid-point epoch '
     .     ,tjdmid
      jds = jdstart
      ts = tstart
 
c PT990614: need to set some variables to zero prior to calling get_yaw
c           so the HP doesn't get upset!
      iendf = 0
      jlast = 0
      ji0 = 0
      iy1 = 0
      iy2 = 0    
                  
c Initialize the event counters and bias-warning flag
      do i=1,maxsat
       ievent(i) = 0                     
       bias_warning(i) = .false.
      enddo
                          

c ** Loop over all epochs ***
* MOD TAH 210202: FAKE it an add extra epoch to avoid model error.

      do iepoch = 1,nepoch
         
        if(iepoch.eq.1)then   
          do j = 1,nsat
            first_call(j) = .true.
          enddo
        elseif(iepoch.eq.2)then
          do j = 1,nsat
            first_call(j) = .false.
          enddo
          call timinc(jds,ts,calc_int*1.d0)
        else
          call timinc(jds,ts,calc_int*1.d0)
        endif
c       write status every 500 epochs
        if( mod(iepoch,500).eq.0 ) then
         write(message,'(a,i5)') 'Epoch',iepoch
         call report_stat('STATUS',module,'orbits/yawtab',' ',message,0)
        endif
* MOD TAH 210202: Report last epoch
        if( iepoch.eq.nepoch ) then
          write(message,'(a,i5)') 'Last Epoch',iepoch
          call report_stat('STATUS',module,'orbits/yawtab',' ',
     .           message,0)
        endif

                             
        fjd = jds+ts/86400.d0+tdtoff
        if( fcheck('nbody') ) then                
cc          print *,'YAWTAB calling ephred iepoch,sun ',iepoch,sun  
          call ephtrp(fjd,3,0,sun) 
        else          
          call solred(1,fjd,sun)  
          if(debug) print *,'YAWTAB calling solred iepoch,fjd,sun '
     .         ,iepoch,fjd,sun                    
        endif
c       the ephemeris file has coords of earth rel to sun.; we want sun rel. to earth 
cd        print *,'yawtab aft solred '
cd        if(iepoch.gt.1 ) stop 1
        do kk = 1,6
           sun(kk) =  -sun(kk) 
        enddo                           
        satarg= timdif(jds,ts,jdb,tb)               
        if(debug) then 
          print *,'aft solrad jds ts sun ',jds,ts,sun
          print *,'jds ts jdb tb sararg ',jds,ts,jdb,tb,satarg
        endif

        
c Although the time tag for the calculations in kouba_yaw will be
c seconds from the beginning of the t-file to avoid week-boundary
c discontinuities, compute also GPS week and SOW in order to compare 
c with Kouba's test case
        call jd_to_gps(jds,ts,week,sow)
                        
        if( debug ) then 
          print *,' '
           print *,'YAWTAB iepoch jds ts sow satarg '
     .             ,iepoch,jds,ts,sow,satarg 
        endif

c Loop over satellites and compute the yaw angle  

        do isat = 1,nsat                                                

          if( itsat(isat).eq.iprndb .and. iepoch.lt.2 ) then
            debug =.true.
          else
            debug = .false.
          endif 

          if( debug ) print *,' ' 
          if( debug ) print *,'YAWTAB isat isv svantbody '
     .           ,isat,isvn(isat),svantbody(isat)
cd        if( itsat(isat).eq.6 ) debug = .true.  
cd        if( iepoch.gt.1 ) stop 0   


c         get the satellite coordinates     
cd        print*,' satarg isat sdelt yy11 ',satarg,isat,sdelt,yy(1,1,1)
          call gsatel( 2,satarg,iut,isat,svec,ytrp,yy
     .               , nsat,sdelt,nintrs,ji0,jil,iy1,iy2,jlast
     .               , jnow,iendf,nepcht)
cd        print *,'From gsatel svec: ',svec


c         get the beta angle and sv body unit vectors from the library routine; 
c         we recompute beta later but kouba_yaw needs the body unit vectors
          if(debug) print *,' calling yaw_attit svec sun  ', svec,sun
          call yaw_attit( svec,sun,ideal_yaw,beta,sc_x_t,sc_y_t
     .                  , sc_z,sc_W,svantbody(isat) )   
          if(debug.and.itsat(isat).eq.iprndb.and.iepoch.lt.2 ) 
     .        print *,'From yaw_attit beta SC x y z W :'
     .                    ,beta,sc_x_t,sc_y_t,sc_z,sc_W 

c  Get the orbital angles
          call svsun_angles(svec,sun,betar,sunlon,svbcos)  
c           betar is the angle between the sun vector and the orbit normal (radians)
c           betadg is the conventional angle between the SV plane and the sun direction (deg)
c           For kouba_yaw I've kept the original input variable but for the new routines
c           I've used the conventional definition (betatdg) since that's what's used
c           internally to the routines. 
            betadg = betar/dtr -90.d0  

          sunlon = sunlon/dtr
cd         print *,'TEST SUN_ANGLES betar betadg svbcos sunlon'
cd     .        , betar,betadg,svbcos,sunlon
cd          if (iepoch.gt.1 ) stop 3
cd          print *,'YAWTAB aft get_yaw iepoch jds ts isat  yaw_angle'
cd     .          ,iepoch,jds,ts,isat,yaw_angle(isat)
cd          print *,'jdb tb ',jdb,tb 
          if(debug.and.iepoch.gt.3 ) stop 1

c         Kouba wants the vectors in m, not km
          do i=1,3
            xsv(i) = svec(i)*1.d3
            vsvc(i) = svec(i+3)*1.d3
            xsun(i) = sun(i)*1.d3
          enddo 
          if( debug ) then 
             print *,'YAWTAB iepoch satarg isat svn svantbody yrate'
     .     ,         iepoch,satarg,isat,isvn(isat),svantbody(isat)
     .     ,        yrate(isat)
             print *,'  betar svbcos xsv sc_x_t '
     .               ,betar,svbcos,xsv,sc_x_t
             print *,' svantbody yrate beta ideal_yaw,xhat '
     .           ,isvn(isat),svantbody(isat),yrate,beta,ideal_yaw,sc_x_t
          endif

          if( gnss.eq.'G' ) then 
           if( .not.newcode ) then 
             call kouba_yaw( week,sow,satarg,isat,itsat(isat),isvn(isat)
     .                 , svantbody(isat),yrate(isat),ybias(isat),xsv
     .                 , sc_x_t,vsvc, betar, svbcos
     .                 , yaw_angle(isat),ievent(isat),svbcos_start(isat)
c                      this for debug
     .                 , iepoch )     
            else                         
            call kouba_gps(satarg,svantbody(isat),aybias(isat)
     .                    ,yrate(isat)
     .                    ,xsv,vsvc,sc_x_t,betadg,svbcos
     .                    ,week,sow,isat,itsat(isat),isvn(isat)
     .                    ,yaw_angle(isat),ievent(isat),ybias(isat)
c                         this for debug 
     .                   ,xsun,iepoch ) 
            endif

          elseif( gnss.eq.'R' ) then
            call kouba_glonass( satarg,svantbody(isat),yrate(isat)
     .                        , xsv,vsvc,sc_x_t,betadg,svbcos      
     .                        , week,sow,isat,itsat(isat),isvn(isat)
     .                        , yaw_angle(isat),ievent(isat),ybias(isat)     
c                             this for debug
     .                        , iepoch ) 

          elseif( gnss.eq.'C' ) then 
            if( wang ) then     
              call wang_beidou
            else
              call kouba_beidou( 
     .                   svantbody(isat),vsvc,sc_x_t,betadg,svbcos        
     .                 , week,sow,isat,itsat(isat),isvn(isat)
     .                 , yaw_angle(isat),ievent(isat),ybias(isat)
c                       this for debug only
     .                 , iepoch )
            endif
          elseif( gnss.eq.'E' ) then     
            call kouba_galileo( 
     .            isat,yrate(isat),betadg,svbcos,xsv,vsvc,xsun
     .          , satarg,week,sow,itsat(isat),isvn(isat),svantbody(isat)
     .          , yaw_angle(isat),ievent(isat)
c                 this for debug only 
     .            ,iepoch )         
cd              if( itsat(isat).eq.12) 
cd     .           print *,'YAWTAB iepoch isat ievent yaw_angle '
cd     .           , iepoch,isat,itsat(isat),ievent(isat),yaw_angle(isat)

          elseif( gnss.eq.'I') then                              
            call irnss_yaw( vsvc,sc_x_t,betadg,yaw_angle ) 
            ievent(isat) = 0 
          endif

          if(debug)  print *,'YAWTAB after calls, isat yaw_angle '
     .                      , isat,yaw_angle(isat)
          debug = .false.

c If mid-span, write the orbital information into yawtab.out
          if( iepoch.eq.nepcht/2 ) then
            write(iyawp
     .         ,'(a,i2.2,a,i3,1x,a,f10.4,3f10.2)')
     .         ' PRN ',itsat(isat),'  SVN ',isvn(isat)
     .         ,  svantbody(isat),yrate(isat),ybias(isat)
     .         ,  betadg,sunlon 
          endif

c       end loop over satellites
        enddo

c Write the epoch time, yaw-angles, and event flags to the y-file 
        if(mod(iepoch*1.d0,interval_factor).eq.0)then 
          etime = jds*1.d0 + ts/86400.d0
          hr = int((etime-jds*1.d0)*24.d0)
          min = int(mod(((etime-jds*1.d0)*24.d0),1.d0)*60.d0)
          sec=nint(mod((mod(((etime-jds*1.)*24.d0),1.d0)*60.d0),1.d0)
     .                                                        *60.d0) 
cd        write(*,'(3i5,32f9.2)') hr,min,sec,(yaw_angle(i),i=1,nsat)  
cd        stop 1                                                   
c       version 970
c         write(iyawb) jds*1.d0 + ts/86400.d0,(yaw_angle(i),i=1,nsat)
c        version 1051  
c        for at least the IIA eclipse period, the angle can get more than 360 degrees, so:
         do i=1,nsat
           if(yaw_angle(i).gt.360.d0) yaw_angle(i) = yaw_angle(i)-360.d0
         enddo
         write(iyawb) jds*1.d0 + ts/86400.d0
     .              ,(yaw_angle(i),ievent(i),i=1,nsat)
        endif
                  
c***** end of loop over epochs
      enddo

      write(message,'(a,a)')' Created file: ',yfile
      call report_stat('STATUS',module,'orbits/yawtab',' ',message,0)
      write(message,'(a)')' Normal stop in YAWTAB'
      call report_stat('STATUS',module,'orbits/yawtab',' ',message,0)
      close(iyawb)
      end

c*** temporary conversion from jd,sod to gpsw sow

      subroutine jd_to_gps(jd,sod,week,sow)
                                    
      implicit none

      integer*4 jd,week,yr,doy,hr,min,itflag
      real*8 sod,sow,sec,utcoff
             
      call dayjul(jd,yr,doy)
      call ds2hms(yr,doy,sod,hr,min,sec)            
      itflag = -4
      call timcon(itflag,week,sow,yr,doy,hr,min,sec,utcoff)
      return
      end

