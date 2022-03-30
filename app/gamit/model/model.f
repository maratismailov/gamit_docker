Copyright (c) Massachusetts Institute of Technology and the University of
California, San Diego. 1994; 2005.   All rights reserved.

       Program MODEL
C
C       Written at MIT by R. King in July 1987 using older code from the
C       MIT Planetary Ephemeris Program (PEP) as modified by S. A. Gourevitch,
c       R. Abbot, and Y. Bock (July 1981 - June 1987).  Modifications since
c       January 1987 are recorded in subroutine MVERSN.

C       MODEL computes values for time delay and delay rate, from
C       which it forms the theoretical phase observable.  It also
C       computes partial derivatives of phase with respect to station
C       position, and orbital, clock,and atmospheric parameters.

c       Converted to work in GPST (rather than UTC) July 13, 1994.  X-,
c       and C-file times will be converted to GPST if written originally
c       in UTC.  T-file start and stop times will likewise be converted.
c       The argument for reading the satellite coordinates from the T-file
c       is the difference between the T-file start time and the current epoch,
c       so doesn't depend explicitly on the time used.

C       Input files (unit numbers assigned in subroutine OPEN):
C
C         Tabular ephemeris (T file)
C         Data (X or C) or simulation-control (S) file
C         Station coordinate file (L file)
C         Geodetic datum file
C         Print output (P file)   
c         Antenna and receiver information (station.info)
c         Session span (session.info) - needed only for simulations
c         Receiver antenna phase-center model (file antmod.dat) 
C         Weather table (W file) - optional
C         WVR data table (Z file) - optional 
C         Ocean tide file (file otide.) - optional                  
c         Receiver antenna phase-center model (file antmod.dat) 
c         Satellite antenna offsets (file svant.dat) - optional   
C         UT1 file (file ut1.)         
C         Pole position file  (file pole.)
C         Nutation file (file nutabl.)     
c         (These last three may be omitted if using an Earth-fixed frame, but
c          this option is not yet fully tested)
c
      Implicit none

      include '../../libraries/includes/const_param.h'
      include '../includes/dimpar.h'    
      include '../includes/units.h'  
      include '../includes/global.h'
      include '../includes/errflg.h'
      include '../includes/model.h' 
      include '../includes/arc.h'
                             
      logical  bad, metred, first_iter, fstchan
     .      , fstelev, eop, fcheck,  dcb_found
     .      , dcb_warning(maxsat)
c
      CHARACTER*1 FLET,TRASH
      CHARACTER*4 none
      character*16 uname
      character*22 buf22,date
      character*80 scratch_dir
      character*256 message,tmp_cfile_name,cmd
c 
      integer*4  NEPCHS(MAXSAT)
      integer*4 i,j,k,iarray
     .,  ji0, jil,iy1,iy2, jlast,jnow,iendf,iter

     ., isat,ierr,ncount
     ., iseed,trimlen,idcb
     ., uid,pid,getpid,getuid
     ., iyear,imonth,iday,ihr,imn,isec,ihnsec,kk
c
c     epoch counter
      integer*4 iepochs
c     pass counter
      integer*4 ipass
c     satellite counter
      integer*4 ichan
c     counter for good channel per epoch
      integer*4 ngdchn
c     shadow integer variables
c       note shortened spelling to avoid conflict with differently dimensioned variables in arc.h
      integer*4 neclps(maxsat)
      real*8 eclps_start(maxsat,30),eclps_end(maxsat,30) 

c     function
      integer*4 nblen  
      
c
      REAL*4 noisel1,noisel2,noisep1,noisep2,rand5_unif,rand5_gauss,dum4
      REAL*8 EVEC(6,2),R1(3),R2(3),DR(3),RHAT(3,MAXSAT)
     ., OLDDEL(MAXSAT)
     ., SVEC(MAXYT2),YTRP(5,2,MAXYT2,MAXSAT),YY(10,MAXYTP,MAXSAT)
     ., sitepart(3,3),polepart(3,4)
     ., ut1part(3,2)
      real*8 phctdel(2,maxsat),range1,range2
     ., dot,phase1,phase2,r1leng,r2leng
     ., freq0,atmpart,satarg,fjd
     ., confin,dinter
     ., satcrd(6,2)
     ., pnsmat(3,3),dipoldel(2),iongrpdel(2),ionphsdel(2)
     ., xhat_t(3),yhat_t(3),zhat_t(3)
     ., shadow(maxsat),shadow_old(maxsat)
     ., tshad(maxsat),tshad_old(maxsat)
c variable to estimate max yaw rates
c     ., dyaw_ang_dp                   
     ., svantpart(3,3)
     ., simvec(3,2),simdel(2),simdobs(4) 
     ., rsathat(3,2)  
     ., dcbcorr 
c debug 
c     ., elevd

c
c     All times are kept in Julian day and second of day
c     0    = initial epoch
c     OBS  = observation time tag
c     RCVR = true receive time
c     SEND = true send time
c     LAST = last RCVR time, used for integration
c     C    = Reference time for sat clock polynomial
c     EOP  = Reference time for earth-orientation parameters (center of span)
c
      INTEGER*4 JDRCVR,JDSEND,JDC
      REAL*8 TRCVR,TSEND,TC
c
c     for integral of TRM2
      REAL*8    SUM2(MAXSAT),OLD2(MAXSAT)
      INTEGER*4 ILAST(MAXSAT)
c
c     relative times
      real*8 trmtc
c     function to return the difference between two time tags
      real*8 timdif

c     function for bit-coded flags
      logical kbit
c
c     clock treatment
c     klock   Receiver               Satellite
c     1       0                      poly        [Minimacs]
c     2       poly                   poly        [TI4100s]
c     3       poly + epoch-by-epoch  poly        [Trimbles]
c
c     receiver clock offset in seconds: 1 estimate per sat
      real*8 clksav(maxsat)
c     the estimate we actually use
c
c     0,1,2 order terms in polynomial for satellite (SV) clock offset, and L1 phase correction
c     units are sec,dimless, 1/s, and cycles
      REAL*8 svcepc,svcrat,svcacc,svcl1

c     satellite clock offset and L1 phase correction values written on C-file

c     For debug only
      real*8 sat_l(2), site_l, amag3, r_tmp(3)
c
c     variable to track whether a satellite has been observed before
      logical start(maxsat)
                                        

      data none /'none'/

      real*8 sitrad,satrad,relcon,reldel,gme

* MOD TAH 000312: Original relativity correction (Earth graviational
*     bending term).  The new correction added is the gravitional effect
*     on the satellite oscillator.
      real*8 reldeo

c  yaw/shadow variables
      real*8 yatt(maxsat),sec(2)
      integer*4 yr(2),mo(2),day(2),hr(2),min(2),ievent(maxsat)

c  quantities for Sun and Moon
      real*8 sun(6),moon(6) 
c     for positions only, ispeed=1 for old solred/lunred; ivel=0 for new ephtrp
      integer*4 ivel/0/, ispeed/1/

c  local station unit vectors - the third can be formed by a cross product
      real*8 naxis(3),uaxis(3)
      data naxis/1.d0,0.d0,0.d0/
      data uaxis/0.d0,0.d0,1.d0/

c  station local unit vectors rotated to inertial space
      real*8 nhat_neu_i(3),uhat_neu_i(3)
          
      logical debug/.false./

c  Assign terminal and screen (other units assigned in open.f)
      iterm = 5
      iscrn = 6
                                    
c  Initialize the counter for c-file header text
      ntext = 0 

c  Skip MODEL if a previous step has failed
      if( fcheck('GAMIT.fatal') )
     .  call report_stat('FATAL','MODEL','model',' '
     .                  ,'GAMIT.fatal exists: MODEL not executed',0)

c  Read the batch file
      call read_batch( scratch_dir,trash ) 
                      
c  Open the input and output files and write the run headers to the screen and p-file
 
      call OPEN( trash,scratch_dir )       

c  Set up input controls, read and write file headers

      call SETUP 
      if(debug) print *,'MODEL aft setup jd0 t0 jdend tend  '
     .   ,jd0,t0,jdend,tend

c  Initialize the t-file interpolation pointers
           
      iendf = 0
      ji0 = 0
      jil = 0
      iy1 =0
      iy2 = 0
      jlast = 0
      jnow = 0
                     

c  Initialize the other variables and arrays 
           
      twopi = 2.d0* pi
      rclock0 = 0.0d0
      do k=1,maxsat
        olddel(k)= 0.1d0
        sum2(k) = 0.0d0
        old2(k) = 0.0d0
        nepchs(k) = 0
        ilast(k) = 0
        phctdel(1,k)=0.0d0
        phctdel(2,k)=0.0d0 
        tshad_old(k)=0.0d0
        do j=1,maxytp
          do i=1,10
            yy(i,j,k)=0.d0
          enddo
        enddo
      enddo
      dinter = inter
c     random-number seed for simulation
c     temporary to get stations to repeat: use x coordinate in mm
      iseed = idint(dabs(evec0(1,1)*1.d3)) 
      dum4 = rand5_unif(iseed)
cd      print *,'MODEL: iseed dum4 ',iseed,dum4 
      do i=1,4
        simdobs(i) = 0.d0
      enddo                        
      do i=1,maxsat
        dcb_warning(i) = .false.
      enddo
      first_iter = .true.
      do i = 1,maxsat
        shadow(i) = 1.d0
        shadow_old(i) = 1.d0
        neclps(i) = 0
        start(i) = .true.
        do j = 1,30
          eclps_start(i,j) = 0
          eclps_end(i,j) = 0
        enddo
      enddo
c     Define variable to correct GPST to AT for interpolating the sun position
      tdtoff =  (32.184d0 + 19.0d0)/86400.d0
       

c  Set up planetary ephemeris files for interpolation

      if( fcheck('nbody') ) then
c       no need for Venus and Jupiter, or for Sun and Moon velocities                          
        lbody = 0         
        fjd = jd0+t0/86400.d0+tdtoff
cd        print *,'jd0 t0 fjd ',jd0,t0,fjd 
        call ephred(ibody,fjd,lbody,ivel) 
      else   
cd         print *,'MODEL units ',iterm,iscrn,iprnt,ibody,isun,ilun
cd     .                       ,inut,iut1,ipole                  
cd         print *,'MODEL lunits ',isesfo,istnfo,ihi,ipcv,iuu,iudcb
cd     .        ,iuc,iud,iuf,iui,iuj,iobs,iuw,iuy,iuz,iueqrn,iepcv
cd     .        , iul,iug,iut             
cd             
cd         print *,'MODEL bef EPHDRD frame ',frame
cd        print *,'   jd0 t0  ',jd0,t0   
        fjd = jd0+t0/86400.d0+tdtoff
        call ephdrd(fjd )
c       read the constant Everett coefficients into storage
        call evrtcf 
      endif 

                       
c ******** Main Loops of Program ************************************

      call report_stat('STATUS','MODEL','model',' '
     .                ,'Begin processing',0)
                       
      do 900 iepoch=1,nepoch   
cd      if(iepoch.eq.2 ) stop
             
cd    if(iepoch.eq.6 ) stop
c     Initialize omc (not sure if this is necessary)
      omc(:,:) = 0.d0

c     Read the observation data file or get the epoch for simulation
           
      call obsred  
                     
c     For the Ashtech codeless receiver, add the addhoc clock correction
               
      if( rcvcod.eq.'ASHL12' .and..not.fixash )  
     .    call timinc( jdobs,tobs,.000045d0 )  

c     If any valid observations at this epoch, update the satellite attitudes
c     and station coordinates
                           
      if( msat.ne. 0 ) then   
        if(iuy.gt.0) call satatt( iuy, jdobs*1.d0+tobs/86400.d0 ,ntsat
     .                         ,  yatt, ievent, iepoch )
        call update_coords
      endif      
      if(debug) print *,'MODEL aft satatt ntsat iepoch yatt(1-4) '
     .       ,            ntsat,iepoch,yatt(1),yatt(2),yatt(3),yatt(4)

c*     write(*,'(a,3f15.7)') '   evec0 ',(evec0(1,1),i=1,3)  
cd      if( iepoch.gt.2 ) stop

c     only need to read the met file once per epoch
      metred = .false.

c -------- Do two passes to average receiver clock correction --------

      DO 710 IPASS=1,2

c -------- Begin Loop on Data Channels (Satellites) ------------------

      fstchan = .true.
      fstelev = .true.
      eop = .true.
      ngdchn = 0
c 
      DO 700 ICHAN=1,NCHAN

c     Don't bother if no data, data deleted in cleaning, sv missing from t-file,
c     or the smaller segment in a session with a coordinate shift

      if( sv_missing(ichan)) ier(ichan) = ignosv 
      if( iepoch.lt.nepoch_use(1).or. iepoch.gt.nepoch_use(2) ) 
     .       ier(ichan) = igmodl
      if( ier(ichan).eq.ignone .or. ier(ichan).eq.igchop .or.
     .    ier(ichan).eq.ignosv .or. ier(ichan).eq.igmodl ) goto 700   

c     count epochs per channel over session, and good channels for each epoch
c
      if (ipass .eq. 1) then
         nepchs(ichan) = nepchs(ichan) + 1
         ngdchn = ngdchn + 1
      endif                 

c     get the differential code bias to be applied to  C1 or P2'(see sb obsmod)    
      dcbcorr = 0.d0  
      if( pcncod.eq.'P'.or.pcncod.eq.'C' ) then
        idcb = 0
        if( kbit(data_flag(ichan),28) )  idcb = 1
        if( kbit(data_flag(ichan),29) )  idcb = 2
        if( idcb.gt.0 ) then  
          dcb_found = .false.
          do i=1,numdcb
            if( prndcb(i).eq.ischan(ichan) ) then   
              dcbcorr = dcb(i)
              dcb_found = .true.
            endif
          enddo
          if( .not.dcb_found .and..not.dcb_warning(i)) then
            write(message,'(a,i2)') 'DCB correction not found for PRN '
     .          ,ischan(ichan)
            call report_stat('WARNING','MODEL','model',' ',message,0)
            dcbcorr = 0.d0    
            dcb_warning(i) = .true.
          endif
        endif  
       endif

c     initialize delay, delay rate, and elevation angle
c
      if (ipass .eq. 1) then
         delay(1,ichan)= olddel(ichan)
         delay(2,ichan)= olddel(ichan)
         drate(ichan)= 0.d0
         elev(ichan)=  0.d0
         elvdot(ichan)= 0.d0
         azim(ichan) = 0.d0
         azmdot(ichan)= 0.d0
         atmdel(ichan)= 0.d0
         okmet=  0
C         wetvar= 0.d0
         clksav(ichan)= 0.d0
      endif
      iter= 0

c      Compute the receiver clock corrections from polynomial
                               
      call poly_clk( ipass )                    
c-------------- come here on light iteration --------------------------

200      continue

c        Iteration of send time 

c        reset jdsend to jdobs before adding delays
         jdsend= jdobs
         tsend = tobs    
         if(debug) print *
     .     ,'MODEL 1 iter ipass ichan jdsen  tsent rclock delay '
     .           , iter,ipass,ichan,jdsend,tsend,rclock,delay(1,ichan)
         call timinc( jdsend,tsend,-rclock-delay(1,ichan) )  
         if(debug) print *,'MODEL    jdsend  tsend ', jdsend,tsend

c        Compute the receiver clock correction from satellite clocks,
c        pseudo-range, and the theoretical delay, if epoch clock required.
         call epoch_clk( ipass,iter,jdsend,tsend,ichan,none )

c        Compute the receive time, correcting for receiver clock offset.

         jdrcvr = jdobs
         trcvr= tobs
         call timinc( jdrcvr,trcvr,-rclock )

c        Compute the transmit time, correcting for receiver clock offset

         jdsend= jdobs
         tsend = tobs     
         if(debug) print *
     .      ,'MODEL 2 ipass ichan jdsen  tsent rclock delay '
     .                , ipass,ichan,jdsend,tsend,rclock,delay(1,ichan)
         call timinc( jdsend,tsend,-rclock-delay(1,ichan) )
         if(debug) print *,'MODEL    jdsend  tsend ', jdsend,tsend


c        Compute the precise Sun and Moon positions.....
c        Need only call it for the first iteration. Time passed is AT not GPST

         if (iter.eq.0) then
           if( fcheck('nbody') )  then   
             fjd = jdsend+tsend/86400.d0+tdtoff
             call ephtrp(fjd,3,ivel,sun) 
             call ephtrp(fjd,10,ivel,moon)
           else           
             fjd = jdsend+tsend/86400.d0+tdtoff
cd             print *,'calling SOLRED fjd ',fjd                  
             call solred(ispeed,fjd,sun)
             call lunred(ispeed,fjd,moon) 
c            planetary ephemeres returns coords of earth rel to sun.; we want sun rel. to earth
           endif
          do kk = 1,6
              sun(kk) =  -sun(kk)
           enddo
         endif     
cd         print *,'sun moon ',sun,moon

c        Determine Site Coordinates and Partials

c        evec(6,2) is the site coordinate and velocity vector
c        only need to call sitcor once pre receiver clock iteration pass.
         if (fstchan) then  
           call sitcor ( jdrcvr,trcvr,sun,moon,evec,sitepart,pnsmat
     .                 , eop,polepart,ut1part,simvec,ipass )
           fstchan= .false. 
         endif  

c        Compute the satellite position and velocity coordinates
c        isat is index in itsat (t-file array)

cd         if( ichan.lt.0 .or. ichan.gt.30 )
cd     .       print *,'iepoch ichan ntsat ',iepoch,ichan,ntsat
         isat= iarray( ischan(ichan),itsat,ntsat )
         satarg= timdif(jdsend,tsend,jdbt,tbt)
         call gsatel( 2,satarg,iut,isat,svec,ytrp,yy
     .             , ntsat,sdelt,nintrs,ji0,jil,iy1,iy2,jlast
     .             , jnow,iendf,nepcht)     
 
c        Store the SV position and velocity for yaw and SV antenna PCV calculations
         do j=1,2
           do i = 1,6
             satcrd(i,j) = svec(i)
           enddo 
           sat_l(j) = amag3(satcrd(1,j))       
           do i=1,3
             rsathat(i,j) = satcrd(i,j)/sat_l(j)
           enddo
         enddo
         if(debug) print *,'tsend satcrd ',tsend,satcrd 
          

c        Compute amount of shadowing of the satellite by the earth

         shadow_old(ichan) = shadow(ichan)
         tshad(ichan) = jdsend + tsend/86400.d0
c        Check if last eclipsed epoch was over 2 hour ago. If so then the
c        satellite set while in last eclipse. Therefore set shadow_old(ichan) = 1.d0
         if((tshad(ichan) - tshad_old(ichan)).gt.7200.d0/86400.d0) then
           shadow_old(ichan) = 1.d0 
         endif
         shadow(ichan) = 1.0d0
         call shadow1(satcrd(1,1),sun,shadow(ichan)) 
c        Check whether sat is in eclipse, if so store the time
         if (shadow(ichan).lt.1.d0) then
           tshad_old(ichan) = tshad(ichan)
         endif
c        Determine the first eclipsed epoch of eclipsing satellites
         if ( shadow(ichan) .lt. 1.d0 .and.
     .        shadow_old(ichan) .ge. 1.d0 ) then
           neclps(ichan) = neclps(ichan) + 1 
           eclps_start(ichan,neclps(ichan)) = tshad(ichan)
c        Determine last eclipsed epoch of eclipsing satellites
         elseif ( shadow(ichan) .lt. 1.d0 ) then
            eclps_end(ichan,neclps(ichan)) = tshad(ichan)
         endif            


c        Compute the yaw angle, the satellite unit vectors and correct the unit
c        vectors for any yaw error. Routines from JPL (IGS mail #0591)
                            
         if(debug) print *,'MODEL calling svbody_coords isat yatt '
     .      ,isat,yatt(isat)
         if(debug)  print *,'  satcrd ',satcrd 
         call svbody_coords( iuy,satcrd(1,1),sun,xhat_t,yhat_t,zhat_t
     .          , itsat(isat),svantbody(ichan),yatt(isat),ievent(isat)
c            for debug 
     .          , iepoch)
         if(debug) print *,'aft SVBODY_COORDS satcrd xhat yhat zhat '
     .       ,satcrd,xhat_t,yhat_t,zhat_t 
c         Correct the satellite vector for the antenna phase center offset

         call svant( satcrd,ischan(ichan)
     .             , xhat_t,yhat_t,zhat_t,svantdx(1,1,ichan),svantpart )
         if(debug) then 
           print *,'i iprn satcrd ',ichan,ischan(ichan),satcrd
           print *,'xhat yhat zhat ',xhat_t,yhat_t,zhat_t     
           print *,'svantpart ',svantpart
           print *,'SV PCO (L1/L2) satcrd ',
     .              ((svantdx(i,j,ichan),i=1,3),j=1,2),satcrd
         endif 

c        Compute the delay and observation unit vector

         do i=1,3
           r1(i)=satcrd(i,1)-evec(i,1)
           r2(i)=satcrd(i,2)-evec(i,2)
c*       for debugging:
           r_tmp(i) = evec(i,1)
         enddo                            
         r1leng=dsqrt(dot(r1,r1))
         r2leng=dsqrt(dot(r2,r2))
c        Debug: compute magnitude of sat vectors (sat vector done above)
         site_l = amag3(evec)
         if(debug) print*,' sat and site vector magnitudes '
     .      , sat_l, site_l
         delay(1,ichan)=r1leng/vel_light_km
         delay(2,ichan)=r2leng/vel_light_km
         do  i=1,3
            rhat(i,ichan)= r1(i)/r1leng
         enddo  
         if(debug) print *,'delay ',delay(1,ichan),delay(2,ichan)      


c     Decide whether to iterate the delay calculation

         confin=delay(1,ichan)-olddel(ichan)
cd         print *,'olddel confin ',olddel(ichan),confin
         olddel(ichan)=delay(1,ichan)
         iter= iter + 1
         if( iter.gt.7 ) goto 991
         if(dabs(confin).gt.1.d-9) go to 200
c        save the receiver clock correction to use later in avclck when computing
c        clock offset
         clksav(ichan) = rclock                             
cd         print *,'ichan olddel rclock ',ichan,olddel(ichan),rclock 


c -------------- End of the Light-time Iteration --------------------
                           
                                            
c       Get the satellite velocity and partials
c
         call gsatel( 5,satarg,iut,isat,svec,ytrp,yy
     .              , ntsat,sdelt,nintrs,ji0,jil,iy1,iy2,jlast
     .              , jnow,iendf,nepcht )

c       Compute the delay rate

        do i=1,3
           dr(i)=svec(i+3)-evec(i+3,1)
        enddo
        drate(ichan)=dot(dr,r1)/(vel_light_km*r1leng)
  

        if (ipass .eq. 2) then

c       Compute the delay guess for the next epoch

            olddel(ichan)=delay(1,ichan)+drate(ichan)*dinter

c       Compute the atmospheric corrections to the delay

c           compute the station-satellite elevation angle and azimuth from north
            call az_elev( fstelev,naxis,uaxis,pnsmat,evec0,r1,twopi
     .                  , ichan,azim,elev,nhat_neu_i,uhat_neu_i )
c    rwk 070910: Comment out until logic restored in az_elev
c            fstelev = .false.                   
                                             
c  The minimum elevation angle for processing is now controlled by AUTCLN except
c  with LC_HELP it is possible to set a higher cutoff in SOLVE.  
c  Two cases where an observation is skipped are when the PCV table has a 
c  higher cutoff than the one requested in autcln, and in simulation mode   
c*** debug
c               elevd = elev(ichan)*360.d0/twopi 
c***
            if( elev(ichan)*360.d0/twopi .lt.pcvminelev  ) then 
               ier(ichan) = igmodl     
               goto 700
            endif
            if( simulation ) then       
               if( elev(ichan)*360.d0/twopi .lt.elvcut ) then 
                  ier(ichan) = igloel
                  goto 700
               endif
            endif 
cd         print *,'MODEL: iepoch,ichan elev elevd elvcut '
cd     .      , iepoch,ichan,elev(ichan),elev(ichan)*360.d0/twopi,elvcut

                 
            call get_atmdel( jdobs,tobs,ischan(ichan),azim(ichan)
     .                 ,elev(ichan),atmdel(ichan),atmpart )    

cd       write(*,'(a,3f20.15)')'atdel,part,elev',atmdel(ichan),atmpart
cd     .          ,elev(ichan)*180.d0/3.1415926d0
cd      print*,' drymap and wetmap ',drymap,' ',wetmap
cd      print*,'elevation angle was ',elev(ichan)*180.d0/3.1415926d0
cd      print*,'ichan,ischan(ichan) and tobs =',ichan,ischan(ichan),tobs
cd       stop
               
c      Compute the higher-order ionospheric corrections
cd      Print*, 'Model/model ionsrc: ',ionsrc,'jdobs: ',jdobs,'tobs:'
cd     .             , tobs,'latr',latr,'lonr',lonr,'height',height
cd     .             , 'sitecd ',sitecd,'  prn',ischan(ichan)
cd      Print*,'Model/modelsatcrd:',satcrd,'evec0(1,1)',evec0(1,1),'evec0'
cd     .             , evec0,'rhat(1,ichan)',rhat(1,ichan),'elev'
cd     .             , elev(ichan),'azim',azim(ichan)
cd      Print*,'Model/model pnsmat'
cd      Print*,pnsmat
        call iondel( ichan,ischan(ischan),satcrd(1,1),pnsmat
     .             , rhat(1,ichan),iongrpdel,ionphsdel ) 

cd       print*, 'MODEL\model - returned from iondel'
cd       print *,'MODEL iongrpdel ionphsdel ',iongrpdel,ionphsdel
C EJP debug 
cd     Print*, 'MODEL\model: evec0:', evec0(1,1), 'satcrd',satcrd

c      Compute antenna / transmitter orientation dependent phase
c      corrections for Right Circularly Polarized electro magnetic waves
c         Reference Wu J.T., et al Manuscripta Geogetica (1993) 18:pp91-98...

            call dipole_comp( evec0,pnsmat,r1,r1leng
     .                      , fL1(ichan),fL2(ichan),ichan
     .                      , xhat_t,yhat_t,dipoldel,iepoch )
                                            
c      Compute the General Relativistic time delay due to the Earth

            sitrad = dsqrt(dot(evec(1,1),evec(1,1)))
            satrad = dsqrt(dot(svec,svec))
            gme = 398600.5d0
            relcon = (2.d0*gme/vel_light_km**3)
* MOD TAH 000301: Name of original relativity correction reldeo and
*           added computation of the other component of the model.  Both
*           are summed and used in the phase calculation.
            reldeo = relcon
     .              *log((sitrad+satrad+r1leng)/(sitrad+satrad-r1leng))
*           Gravitional change effects on the satellite oscillators.  Needed
*           for one-clock estimates (add also the bending term)
            reldel = 2.d0*
     .              (svec(1)*svec(4)+svec(2)*svec(5)+svec(3)*svec(6))/
     .              vel_light_km**2 + reldeo
cd        print *,'reldel ',reldel
      
c      Compute direction-dependent phase-center variations for ground and SV antennas
c           get the nadir angle for SV antenna  
            nadang(ichan) = dasin(site_l/sat_l(1))*dsin(twopi/4.d0 
     .                       - elev(ichan)) 
cd            print *,'MODEL svantbody svantmod nadang'
cd     .           ,svantbody(ichan),svantmod(ichan),nadang(ichan)
            call phasecc( ichan,yatt,phctdel )      
            first_iter = .false.  
cd            if( iepoch.eq.894.or.iepoch.eq.1080 )   
cd       print *,'Aft PHASECC ichan iprn iblk '
cd     .        ,'phctdel ',ichan,ischan(ichan),svantbody(ichan)
cd     .        ,           phctdel(1,ichan),phctdel(2,ichan)
 
c      For simulated observations, compute the difference in geometric delays
c      (L1 and L2) due to a station displacement

            if( simulation ) then     
              do j=1,2                                    
                do i=1,3
                  r1(i) = satcrd(i,j) - simvec(i,j)  
                enddo
                simdel(j)=dsqrt(dot(r1,r1))/vel_light_km-delay(j,ichan)
              enddo
            endif                     


c-----Compute the modelled carrier beat phase and pseudorange
    
cd        if( ichan.eq.2.and.(iepoch.eq.894.or.iepoch.eq.1080 ) )  
cd     .    print *,'OBSMOD ',dipoldel,phase1,phase2,range1,range2
             call obsmod( iscrn,iuj,ischan,nchan,ichan,jd0,t0
     .      , jfiln,none,sum2,ilast,iepoch,fL1(ichan),fL2(ichan)
     .      , vel_light,jdrcvr,trcvr,delay,atmdel,phctdel,reldel
     .      , dipoldel,iongrpdel,ionphsdel,idcb
     .      , dcbcorr,svcepc,svcrat,svcacc,svcl1,clkrat,clkacc,clkcub
     .      , inter,jdc,tc,phase1,phase2,range1,range2 
     .      , simulation,simdel,simdobs )         
                                            

cd         print *,'obs ',(obs(i,ichan),i=1,4)
cd         print *,'reldel,dipodel,phase1,phase1,range1,range2 '
cd     .     ,reldel,dipoldel,phase1,phase2,range1,range2

c           Save the SV clock offset and the phase-correction terms for AUTCLN
            svclock(1,ichan)=svcepc
            svclock(2,ichan)=svcl1
            
c           Finally, form observed minus calculated, in cycles.
            if ( simulation ) then  
              obs(1,ichan) = phase1
              obs(2,ichan) = phase2
c             get noise: assume zero-mean with sigma read from input
              noisel1 = rand5_gauss(noise(1),0.)
              noisel2 = rand5_gauss(noise(2),0.)   
cd              print *,'noisel2 noisel2 ',noisel1,noisel2
              omc(1,ichan) = noisel1 + simdobs(1)
              omc(2,ichan) = noisel2 + simdobs(2)
            else         
              omc(1,ichan)= obs(1,ichan) - phase1
              omc(2,ichan)= obs(2,ichan) - phase2
            endif

c           Observed minus calculated pseudoranges in CYCLES!
c           use actual transmitted frequency
            trmtc = timdif( jdrcvr,trcvr,jdc,tc )

c           use a phase measurement for PR which doesn't include phase windup or 
c           receiver antenna phase centre variations ie range1 and range2 (PT 950322)  
            if (simulation ) then
              obs(3,ichan) = range1
              obs(4,ichan) = range2 
              noisep1 = rand5_gauss(noise(3),0.)
              noisep2 = rand5_gauss(noise(4),0.)   
              omc(3,ichan) = noisep1 + simdobs(3)
              omc(4,ichan) = noisep2 + simdobs(4)
            else    
              omc(3,ichan)= fL1(ichan)*(1.0d0 + svcrat 
     .                      + 2.0d0*svcacc*trmtc)
     .                      *obs(3,ichan)/(vel_light) - range1
              omc(4,ichan)= fL2(ichan)*(1.0d0 + svcrat 
     .                      + 2.0d0*svcacc*trmtc)
     .                      *obs(4,ichan)/(vel_light) - range2
            endif                              

c  This also leads to incorrect widelanes as the phase windup and rec.
c  antenna phace centre modelling are included in phase1 and phase2
c            OMC(3,ICHAN)= FREQL1*(1.0d0 + SVCRAT + 2.0d0*SVCACC*TRMTC)
c     .                    *OBS(3,ICHAN)/(vel_light) - PHASE1
c            OMC(4,ICHAN)= FREQL2*(1.0d0 + SVCRAT + 2.0d0*SVCACC*TRMTC)
c     .                    *OBS(4,ICHAN)/(vel_light) - PHASE2

c           this leads to incorrect wide lanes kurt 910304
c            OMC(3,ICHAN)= FREQL1*(1.d0 + SUM2(ICHAN) + PHSCOR)
c     .                    *OBS(3,ICHAN)/(vel_light) - PHASE1
c            OMC(4,ICHAN)= FREQL2*(1.d0 + SUM2(ICHAN) + PHSCOR)
c     .                    *OBS(4,ICHAN)/(vel_light) - PHASE2

C           Compute the Partial Derivatives
            call partl( ichan,rhat(1,ichan),drate(ichan)
     .               , svec,fL1(ichan),norbpart
     .               , npart,sitepart,atmpart,polepart,ut1part,svantpart
     .               , tmpart,iepoch)

         endif
c----    End of IF on second pass (after clock calculation) -----------

C 
 700     continue
C ---    End of Loop on Receiver Channels (Satellites) -----------

c        Average the receiver clocks and check for outliers (skip for simulation)

         if (ipass .eq. 1 .and. ngdchn.gt.0 .and. .not.simulation ) then
            call avclck ( iprnt
     .            ,clksav,ischan,nchan,ier,rclock0,iepoch,bad )
            if (bad) then    
c              --rwk:  no longer compute or use the polynomial; data are unweighted
c               TR = (JDOBS-JD0)*86400.0D0 + TOBS - T0
c               RCLOCK0 = CLKEPC + CLKRAT*TR + CLKACC*TR*TR
               write(iprnt,*)'No valid data for RCLOCK at epoch'
     .         ,' ',iepoch 
               write(message,'(a,a,i4,a)')'No valid data for RCLOCK at '  
     .       , 'epoch: ',iepoch
c     .       , 'epoch: ',iepoch,' using the polynomial rclock estimate.'
               call report_stat('WARNING','MODEL','model',' ',message,0)
            endif
         endif

 710     continue

C ----   End of Loop on Passes: RCLOCK0 is determined -----------

C        Write the Data onto the Output C-File
              
c** ----- store loading N E U (km) values in 'save' until c-file modified  
c**       C-file modified rwk 100827
c           if( ntatml.gt.0) then
c             nsave = 3
c             do i=1,3
c              save(i) = dble(atmlod(i))
c             enddo 
c           else
             nsave = 1
             save(1) = 0.d0
c           endif
c**--------

c----- Write the observation records to the C-file

         call cfout                            

c         save epoch counter for status printout
          iepochs = iepoch
cd          if(iepoch.eq.9 ) stop
900      continue
C ****    End loop on Epochs *******************************************
         
c        Move /tmp/cfile to proper location in local processing directory if necessary.

         if ( scratch_dir .ne. 'NONE' ) then
           close(unit=iuc)
c          Recreate the scratch c-file name (use uid and pid to make it unique)
           pid = getpid()
           uid = getuid()
           write(tmp_cfile_name,'(a,a1,a,a1,i6.6,i5.5)')
     .        scratch_dir(1:trimlen(scratch_dir)),'/',
     .        cfiln(1:trimlen(cfiln)),'.', pid,uid
c           tmp_cfile_name = scratch_dir(1:trimlen(scratch_dir))
c     .                      //'/'//cfiln
c          Move the scratch cfile to its correct location
           cmd = '\mv '//tmp_cfile_name(1:trimlen(tmp_cfile_name))
     .           //' ./'//cfiln
           call system(cmd)
          endif

c        Print out how much data we saw

         write (iprnt,'(/,a)') 'PRN  Ngood'
         ncount = 0
         do i = 1,nchan
            write (iprnt,'(1x,i2.2,1x,i5)') ischan(i),nepchs(i)
            ncount = ncount + nepchs(i)
         enddo
         if( ncount.gt.0 ) then
            write(message,'(i5,a)') ncount,' valid observations'
            call report_stat('STATUS','MODEL','model',' ',message,0)
         else
             call report_stat('WARNING','MODEL','model',' '
     .                       ,'No valid data on input X- or C-file',0)
         endif

c        Print out eclipse summary

         do 930 i = 1,nchan
            if ( neclps(i) .ne. 0 ) then
              do 940 j = 1,neclps(i)
c               check to see that the eclipse end time is greater than eclipse
c               start time. If sat. set in eclipse without rising again then this
c               will not be true so assign the time of the last occurance of the
c               sat. as eclipse end time
                if(eclps_end(i,j).le.eclps_start(i,j))then
                  eclps_end(i,j) = tshad_old(i)
                endif
c               convert time in julian day to yr,mo,day,hr,min,sec
                call jul2hms(eclps_start(i,j),yr(1),mo(1),day(1)
     .                       ,hr(1),min(1),sec(1))
                call jul2hms(eclps_end(i,j),yr(2),mo(2),day(2)
     .                       ,hr(2),min(2),sec(2))
                do k=1,2
                  min(k) = nint(min(k)+sec(k)/60.d0)
                  if(min(k).eq.60) then
                     hr(k) = hr(k) + 1
                     min(k) = 0
                   endif
                enddo  
                write(message,920)ischan(i),yr(1),mo(1),day(1),hr(1)
     .                         ,min(1),yr(2),mo(2),day(2),hr(2),min(2) 
                call report_stat('STATUS','MODEL','model',' ',message,0)
                write(iprnt,920)ischan(i),yr(1),mo(1),day(1),hr(1)
     .                         ,min(1),yr(2),mo(2),day(2),hr(2),min(2)
 920            format(1x,'PRN ',i3,' is seen eclipsing from ',i4
     .          ,2(i3),2x,i2,':',i2.2,'  to  ',i4,2(i3),2x,i2,':',i2.2)
 940          continue
            endif
 930     continue
                                          
c         Close input X or C-file 

         flet=obfiln(1:1)
         call uppers(flet)
         call uppers(trash)
         if(flet.eq.'X') then
           close(unit=12)
         elseif(flet.eq.'C') then
           if(trash.eq.'Y') then
             close(unit=12,status='delete')
           elseif(trash.eq.'N') then
             close(unit=12)
           endif
         endif
                                 
         write(message,'(a,a4,a,i5,a)') 'Site ',sitecd,
     .     ' Normal stop in MODEL after ',iepochs,' epochs'
         call report_stat('STATUS','MODEL','model',' ',message,0)
         STOP

 991     call report_stat('FATAL','MODEL','model',' ',
     .   'Too many iterations. Number of light-time iterations > 7',0)

         END






