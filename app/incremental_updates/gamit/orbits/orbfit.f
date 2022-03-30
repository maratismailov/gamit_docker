      Program ORBFIT
c
c Purpose: To fit (in a least squares sense) tabular points from one or more imported 
c          orbit files (sp3, converted to a T-file), using partials and a priori 
c          tabular points from an ARC integrated tfile.
c
c S. McClusky simon@chandler.mit.edu SEPT 1994 - APRIL 1995
c Modified by R. King (rwk@chandler.mit.edu) to accept multiple imported T-files
c Modified by R. King to allow differencing without estimation (nparam=0) October 2017
c
c Input:  (1) Standard ARC integrated GAMIT T-File with partials.
c         (2) T-Files created from NGS ephemeris using ngstot
c
c Output:  Adjustments to initial state vector, rms of fit of orbit components,
c          and optionally plot files containing postfit residuals in x,y,z or 
c          radial, along and cross track.
c
c Notes on parameterization and current limits: 
c
c    1. The dimensions have been set to accommodate 1 reference T-file and imported 
c       files for 3 days (maxtfil=4 in ../includes/orbits.h), each with 15-minute 
c       tabular points (96/day) (maxepc = 3*96 + 28, where the extra 28 epochs account
c       for additional  tabular points needed for interpolation). 
* MOD TAH 200327: Increased to 9-days; max maxtfil=10.
*       Multi-day fits are strange becaase the first and last 1.5 hrs of each are not
*       used because of the way t-files are read (the tfiles are not linked so we can't
*       cross day boundaries).
c        
c    2. Partials are coded for 19 satellite parameters (6 ICs + 13 radiation-pressure)
c       plus 13 global parameters (maxglb in ../includes/orbits.h) (3+3+1 inertial 
c       translations, rotations, and scale, and 3+3 terrestrial rotations and rates.
c
c    3. Original units (km, time-sec, arc-sec) maintained except in printout,
c       where we use meters (m), milliarcseconds (mas), and parts per billion (ppb).

      implicit none

      include '../includes/dimpar.h'
      include '../includes/orbits.h'
      include 'orbfit.h'
                                 
c     Common variables, for orbfit only, in ./orbfit.h:

c       dimension limits: maxsat and maxorb defined in ../includes/dimpar.h
c                         maxtfil and maxglb defined in ../includes/orbits.h
c     integer*4   maxepc,mxoprm
c     parameter   (maxepc=(maxtfil-1)*96+28)
c     parameter   (mxoprm=maxglb+maxorb*maxsat)
c      
c       Unit numbers
c     integer*4   lucmd,lut(maxtfil),luplt(maxsat),lurms,lufit,lug,lusvs
c    .           ,iscrn,iprnt,inut,iut1,ipole
c     common/units/lucmd,lut,luplt,lurms,lufit,lug,lusvs,iscrn,iprnt
c    .            ,inut,iut1,ipole
c
c       T-file header information 
c     integer*4   nics(maxtfil),nintrs(maxtfil),nepcht(maxtfil)
c    .        ,   nprn(maxtfil),itsat(maxsat,maxtfil)
c    .        ,   jdb(maxtfil),jdstp(maxtfil),jde
c     real*8      tb(maxtfil),tstp(maxtfil),te,delt(maxtfil)
c    .        ,   satics(maxorb,maxsat)   
c     character*1 gnss
c     character*4 icsnam(maxorb)
c     character*5 precmod,frame,srpmod,eradmod,antradmod
c     character*16 satnam(maxsat,maxtfil)
c     common /tfhdrs/tb,tstp,te,delt,satics,nics,nintrs,nepcht,nprn
c    .             ,itsat,jdb,jdstp,jde,satnam,icsnam
c    .             ,precmod,frame,srpmod,nutmod,gravmod
c    .             ,eradmod,antradmod
c
c       SVs and epochs actually used
c     integer*4 nsat,isat(maxsat),nepoch,iepstart,iepstop 
c     common /svepochs/ nsat,isat,nepoch,iepstart,iepstop
c
c       Reference vector, O-Cs and partials for every epoch (save to calculate residuals)
c     real*8      ref_pos_vel(6,maxsat,maxepc),omc(3,maxsat,maxepc)
c    .           ,part(3,mxoprm,maxsat,maxepc),drac(3,maxsat,maxepc)
c    .           ,prefit_sum2
c     integer*4 nobs(maxsat)
c     logical   lobs(maxsat,maxepc)
c     common/omcpart/ ref_pos_vel,omc,part,drac,prefit_sum2,nobs,lobs
c                 
c       Parameter values, adjustments, and normal equations
c     integer*4   nparam,islot(mxoprm) 
c     real*8      apr_prm(mxoprm),adjust(mxoprm)
c     real*8      amat(mxoprm,mxoprm),bvec(mxoprm)
c     character*30 prmnam(mxoprm)
c     common/nrmcom/nparam,islot,amat,bvec,apr_prm,adjust,prmnam

c     Local variables

c       File names
      character*80 cmdnam,tfnam(maxtfil),rmsnam,fitnam,gname,pltnam
 
c       Vectors and interpolation indices for each T-file
c     indices saved from read-to-read
      integer*4 ji0(maxtfil),jil(maxtfil),jnow(maxtfil),jlast(maxtfil)
     .        , iendf(maxtfil),iy1(maxtfil),iy2(maxtfil)
c     values from T-file saved from read-to-read
      real*8   yy(10,maxytp,maxsat,maxtfil)
     .         ,ytrp(5,2,maxyt2,maxsat,maxtfil)   
c     postion/velocity/partials vectors calculated each call of gsatel
      real*8    rvec(maxyt2)
             
c       T-file header information not needed in common
      integer*4 jde2
      real*8 satics2(maxorb,maxsat),te2
      character*4 icsnam2(maxorb)
      character*5 precmod2,frame2,srpmod2,nutmod2,gravmod2
     .          , eradmod2,antradmod2
 

c       Other variables
      integer*4   iepoch,jsat,i,j,k,m,ii,jj
      integer*4   ntfiles,iarg,jd,itf,is,iarray
      integer*4   iop,iter
                  
      real*8      start,end,t,trun,timdif,startsec,stopsec
            
      character*14 ic_time
      character*120 version
      character*80 arg
      character*256 message
 
      logical   cmdargs,iterate,debug/.false./

                              
c  Initialization
               
c     logical unit numbers for I/O
      iscrn = 6
      iprnt = 0 
      lucmd = 8      
c     allow for 20 T-files (ref + 19 days)
      do i=1,maxtfil
         lut(i) = 10+i
      enddo 
c**  this causes a compiler warning on DEC/OSF
c      if( maxtfil.gt.20 ) call report_stat('FATAL','ORBFIT'
c     .   ,'orbits/orbfit',' ','Too many T-files for unit assignments',0)
c     units 31-80 reserved for plots (up to 50 SVs)
      do i=1,maxsat
         luplt(i) = 30 + i  
      enddo
      if( maxsat.gt.50 ) call report_stat('FATAL','ORBFIT'
     .   ,'orbits/orbfit',' ','MAXSAT too large for unit assignments',0) 
      lurms = 81
      lufit = 82
      lug   = 83
      lusvs = 84
      inut  = 90
      iut1  = 91
      ipole = 92

c     filenames
      do i=1,maxtfil
        tfnam(i) = ' '
      enddo
      rmsnam = ' '
      fitnam = ' '
      gname = ' '
      pltnam= '  '   

c     interpolation indices
      do i=1,maxtfil  
        ji0(i) = 0
        jil(i) = 0
        jnow(i)= 0
        jlast(i)=0
        iendf(i)=0
        iy1(i)  =0
        iy2(i)  =0
      enddo  

c     other variables
      iop = 0
      do i=1,maxsat
        ibad_sat(i) = 0
      enddo
      ic_time = ' '  
      frame =  'UNKWN'
      frame2 = 'UNKWN'

c  Get version number and write out status line

      call oversn(version)
      write(iscrn,'(a)')' '
      write(message,'(a,a120)') 'Started ORBFIT ',version
      call report_stat('STATUS','ORBFIT','orbits/orbfit',' ',message,0)

       
c  Read the command line 

      call rcpar(1,cmdnam) 
      if( cmdnam.eq." " ) then 
        write(iscrn,10)
10      format(/,'######################################################
     .######',/,'                    Program ORBFIT         '
     .  ,/,'Purpose: Obtain a T-file that best fits 1 or more imported '
     .            ,'SP3 files'
     .  ,/,'Usage: orbfit cmd-file rms-file plt-option ref-tfile tfile1'
     .  ,   ' (tfile2 ..')
c    .  ,/,'  where cmd-file is the input command-file '
c    .  ,/,'         rms-file is the file to write statistics'
c    .  ,/,'         plt-option indicates what quantities to plot'    
c    .  ,/,'            0 to skip plots'
c    .  ,/,'            1 for  dx,dy,dz'
c    .  ,/,'            2 for  dr,da,dc'
c    .  ,/,'            3 for  vx,vy,vz
c    .  ,/,'            4 for  dN,dE,dU and SV ground track'
c    .  ,/,'         ref-tfile is the name of the reference T-file'
c    .  ,/,'         tfile1,.. are the T-files for SP3 files to fit'
c    .  ,/,'------------------------------------------------------------'
c    .  ,'---',/)  
        write(iscrn,11)
   11   format(/,' orbfit.cmd file format is given below:',/,
     .  '---------------------------------------',/,
     .  ' trans:  0 0 0              ',/,
     .  ' i_rot:  0 0 0              ',/,
     .  ' t_rot:  1 1 1              ',/,
     .  ' t_rat:  0 0 0              ',/,
     .  ' scale:  0                  ',/,
     .  ' pos:    0 0 0              ',/,
     .  ' vel:    0 0 0              ',/,
     .  ' srad:   0 0 0 0 0 0 0 0 0 0 0 0 0 ',/,
     .  ' exclude:                   ',/,
     .  ' max_fit_tol:                   ',/,
     .  '---------------------------------------',/,
     .  'Copy text between dashed lines to orbfit.cmd file',/,
     .  'Estimate parameters by selecting 1',/,
     .  'Fix parameters by selecting 0',/)
        stop
      endif
      call rcpar(2,rmsnam)
      call rcpar(3,arg)
      read(arg,'(i1)') iop 
      call rcpar(4,tfnam(1))
      call rcpar(5,tfnam(2))
      if ( rmsnam.eq." " .or. arg.eq." " .or. tfnam(1).eq." "
     .   .or. tfnam(2).eq." " ) call report_stat('FATAL','ORBFIT'
     .   ,'orbits/orbfit',' ','Too few arguments in command line',0) 
      cmdargs = .true.
      ntfiles = 1
      iarg = 5  
      cmdargs = .true.
      do while (cmdargs)  
        iarg = iarg + 1
        ntfiles = ntfiles + 1  
        arg = " " 
        call rcpar(iarg,arg)
        if( arg.ne." ") then
          tfnam(iarg-3) = arg
        else
          cmdargs = .false.
        endif
      enddo

c  Open the files

      call open_orbfiles( ntfiles,cmdnam,rmsnam,pltnam,tfnam )
 
c  Iterate orbit fitting until misfit tolerance is satisfied

      iterate = .true.
      iter = 1
      do while (iterate)
       write(message,'(a,i2)') 'Running iteration: ',iter
       call report_stat('STATUS','ORBFIT','orbits/orbfit',' ',message,0)
       
c     interpolation indices
       do i=1,maxtfil  
         ji0(i) = 0
         jil(i) = 0
         jnow(i)= 0
         jlast(i)=0
         iendf(i)=0
         iy1(i)  =0
         iy2(i)  =0
       enddo  

c  Read in the header of the reference T-file, and get starting/ending times
c                               
       call thdred( 
     .     lut(1),iscrn,iprnt,nprn(1),gnss,itsat(1,1),satnamt(1,1)
     .   , jdb(1),tb(1),jdstp(1),tstp(1),delt(1),nepcht(1),jde,te
     .   , nics(1),satics,nintrs(1),icsnam
     .   , precmod,nutmod,gravmod,frame,srpmod,eradmod,antradmod )   
       if(debug) then
         write(*,'(a,2(i8,f10.1),i5)')  'Ref start stop epochs '
     .         ,jdb(1),tb(1),jdstp(1),tstp(1)	     	
         write(*,'(a,i3,128(1x,i2.2))') 
     .        'SVs ',nprn(1),(itsat(i,1),i=1,nprn(1)) 
       endif

       if( nepcht(1).gt.maxepc ) then
         write(message,'(a,i3,a,i3,a)') 
     .      'Number of epochs on reference T-file (',nepcht(1)
     .      ,') exceeds maxepc (',maxepc,')'
        call report_stat('FATAL','ORBFIT','orbits/orbfit',' ',message,0)
       endif 
       start = jdb(1) + tb(1)/86400.d0
       end   = jdstp(1) + tstp(1)/86400.d0 
cd       print *,'Ref start stop dt ep '
cd     .    ,jdb(1),tb(1),jdstp(1),tstp(1),delt(1),nepcht(1)


c  Read the header information of the secondary t-files

       do i=2,ntfiles    
         call thdred( 
     .       lut(i),iscrn,iprnt,nprn(i),gnss,itsat(1,i),satnamt(1,i)
     .     , jdb(i),tb(i),jdstp(i),tstp(i),delt(i),nepcht(i),jde2,te2
     .     , nics(i),satics2,nintrs(i),icsnam2
     .     , precmod2,nutmod2,gravmod2,frame2,srpmod2,eradmod2
     .     , antradmod2 )    
         if(debug) then
            write(*,'(a,2i3)') 'Secondary t-file i nsat '
     .          ,i,nprn(i)
            write(*,'(a,2(i8,f10.1),i4)') 'start stop epochs '
     .          ,jdb(i),tb(i),jdstp(i),tstp(i)          
           write(*,'(a,i3,128(1x,i2.2))') 
     .        'SVs ',nprn(i),(itsat(j,i),j=1,nprn(i)) 
         endif  

c       check the frame info for consistency
         if(precmod2.ne.precmod .or. frame2.ne.frame ) then
           write(message,'(a,a16,a,4(1x,a5))') 'T-file ',tfnam(i)(1:16)
     .           ,' does not match frame of reference T-file: '
     .           ,  frame,frame2,precmod,precmod2
           call report_stat('WARNING','ORBFIT','orbits/orbfit',' '
     .                  ,message,0)
         endif                     
c       issue a warning if a satellite is missing
         do k=1,nprn(1)
           if( iarray( itsat(k,1),itsat(1,i),nprn(i) ) .eq.0 ) then
              write(message,'(a,a16,a,i2)') 'T-file ',tfnam(i)(1:16)
     .          ,' missing PRN ',itsat(k,1)
              call report_stat('WARNING','ORBFIT','orbits/orbfit',' '
     .                     ,message,0)
           endif
         enddo
       enddo

c  Read input command file and get the parameters to be estimated

       call read_input  
      
      
c  Initialize the data flags number of observations for each satellite
             
       do i=1,nsat
         nobs(i) = 0
         do j=1,nepcht(1)
            lobs(i,j) = .false.
         enddo
       enddo                            

c  Initialize the o-c's and partial derivatives
                   
       if(debug) print *,'nparam ',nparam 
       do m=1,nepcht(1)  
          do k=1,maxsat
             do j=1,3
               omc(j,k,m) = 0.d0
               if( nparam.gt.0 ) then 
                 do i=1,nparam
                   part(j,i,k,m) = 0.d0
                 enddo
               endif
             enddo
          enddo
       enddo
       
c   Initialize the normal equations
   
       prefit_sum2 = 0.d0
       if( nparam.gt.0 ) then
         do i=1,nparam
           adjust(i) = 0.d0
           bvec(i) = 0.d0
           do ii = 1,nparam
             amat(i,ii) = 0.d0   
           enddo
         enddo   
       endif

c  Read through the epochs of the reference T-file, difference the coordinates
c  with those from the second t-file, and if nparam.gt.0, form the normal 
c  equations. We cannot interpolate within 5 epochs of ends of T-file, so skip 
c  these epochs
               
c     set start,stop epochs of reference T-file for interpolation
       iepstart = 6
       iepstop = nepcht(1) -5     
       nepoch = iepstop - iepstart + 1
cd       print *,'iepstart iepstop nepoch ',iepstart,iepstop,nepoch
       call report_stat('STATUS','ORBFIT','orbits/orbfit',' '
     .         , 'Reading T-files and forming normal equations',0)
       do iepoch = iepstart,iepstop  

         trun = dble(iepoch-1) * delt(1)
         jd = jdb(1)
         t = tb(1)
         call timinc(jd,t,trun)       
        if(debug) print *,'Ref t-file epoch jd t trun ',iepoch,jd,t,trun
c       interpolate to get position, velocity, and partials  
         do  jsat=1,nsat
           jj = iarray(isat(jsat),itsat(1,1),nprn(1))
c         interpolate to get position, velocity, and partials   
           call gsatel(4,trun,lut(1),jj,rvec,ytrp(1,1,1,1,1),yy(1,1,1,1)
     .               ,nprn(1),delt(1),nintrs(1),ji0(1),jil(1)
     .               ,iy1(1),iy2(1),jlast(1),jnow(1),iendf(1),nepcht(1)) 
          if(debug) print *,'jsat rvec1-7 ',jsat,(rvec(i),i=1,7)
c         calculate the partial derivatives
          if( nparam.gt.0 ) call partl( jd,t,iepoch,jsat,rvec )
c         save reference vector for plotting
           do i=1,6
             ref_pos_vel(i,jsat,iepoch) = rvec(i)
           enddo    
         enddo

c       see if data available at this time for each t-file
         do itf = 2, ntfiles
c**   need to change these limits to represent stop.start of original sp3
c**    with no overlap between files              
           startsec = 5.d0*delt(itf)
           stopsec = dble((nepcht(itf)-6))*delt(itf)
          if( timdif(jd,t,jdb(itf),tb(itf)).ge.startsec  .and. 
     .        timdif(jd,t,jdb(itf),tb(itf)).le.stopsec ) then
c          if so, read the vector, form O-C, and increment the nrm.eqs for each SV
            trun = timdif( jd,t,jdb(itf),tb(itf) )
            if(debug) print *,' startsec stopsec ',startsec,stopsec
            do jsat=1,nsat 
              is = iarray( isat(jsat),itsat(1,itf),nprn(itf) ) 
               if(debug) print *,'jd t iepoch isat jsat '
     .                         ,jd,t,iepoch,isat,jsat 
              if( is.gt.0 ) then                                  
                call gsatel(1,trun,lut(itf),is,rvec,ytrp(1,1,1,1,itf)
     .             ,yy(1,1,1,itf)
     .             ,nprn(itf),delt(itf),nintrs(itf),ji0(itf),jil(itf)
     .             ,iy1(itf),iy2(itf),jlast(itf),jnow(itf),iendf(itf)
     .             ,nepcht(itf) )  
                do i=1,3
                 omc(i,jsat,iepoch) = rvec(i)-ref_pos_vel(i,jsat,iepoch)
                enddo                         
                if(debug) write(*,'(a,2i4,2(3f15.6),3f10.4)') 
     .                 'is jsat rvec ref_pos_vel omc(m) '
     .                , is,jsat,(rvec(i),i=1,3)
     .                , (ref_pos_vel(i,jsat,iepoch),i=1,3)
     .                , (omc(i,jsat,iepoch)*1.d3,i=1,3)
cd                if(debug) print *,'part(1,2,jsat,iepoch) '
cd     .             ,part(1,2,jsat,iepoch)
                nobs(jsat) = nobs(jsat) + 3
                if(debug) print *,'iepoch,jsat,nobs '
     .             ,iepoch,jsat,nobs(jsat)
                lobs(jsat,iepoch) =.true.  
                if( nparam.gt.0 ) call norminc( iepoch,jsat )
              endif            
cd              stop 1 
            enddo   
          endif
c        end of loop on T-files
         enddo 
cd         if(iepoch.gt.20) then
cd           print *,'DEBUG stop iepoch 20 '
cd           stop 
cd         endif  
c      end of loop on epochs
       enddo
      
cd       print *,'AGAIN nepoch delt(1) ',nepoch,delt(1)
         
c  Solve the normal equations for the estimated ICs and radiation-pressure parameters

       if( nparam.gt.0 ) call norm_solve 
  
cd       print *,'AGAIN 2 nepoch delt(1) ',nepoch,delt(1)
     
c  Compute the statistics and write out the rms, fit, and g files
           
cd       print *,'Calling write_summary ntfiles tfnam ',ntfiles,tfnam
cd       print *,'lurms lufit ',lurms,lufit
       call write_summary ( ntfiles,tfnam )        
cd       print *,'After write_summary nepoch delt lurms lufit '
cd     .        ,lurms,lufit,delt(1)
      
c  See if we need to iterate the solution
       if ( ibad_sat(iter) .gt. 0 ) then 
c  Rewind T-files.
         do i = 1, ntfiles
           rewind(lut(i))
         enddo  
  	 iter = iter + 1
       else
         iterate = .false.
       endif
             
c End of iteration loop
      enddo

c  Plot the residuals

      if( iop.gt.0 )  call plt_postfit(iop,pltnam )

c  That's all folks

      call report_stat('STATUS','ORBFIT','orbits/orbfit',' '
     .                , 'Normal end of ORBFIT',0)
      stop
      end
      

      Subroutine open_orbfiles( ntfiles,cmdnam,rmsnam,pltnam,tfnam )


      implicit none   

      include '../includes/dimpar.h'
      include '../includes/orbits.h'
               
c       File names
      character*80 cmdnam,tfnam(maxtfil),rmsnam,fitnam,gname,pltnam
     .           , svsnam,ut1fil,polfil,nutfil
        
c       Unit numbers
      integer*4   lucmd,lut(maxtfil),luplt(maxsat),lurms,lufit,lug,lusvs
     .           ,iscrn,iprnt,inut,iut1,ipole
      common/units/lucmd,lut,luplt,lurms,lufit,lug,lusvs,iscrn,iprnt
     .            ,inut,iut1,ipole

      integer*4 ntfiles,ioerr,iwhere,index,nblen,i

      logical fcheck

c  Create file names for the fit-file, g-file and plot-files

        iwhere = index(rmsnam,'.')
c       cover case where user leaves off .
        if (iwhere.eq.0 ) then
          iwhere = nblen(rmsnam) + 1
          rmsnam(iwhere:iwhere) = "."   
          rmsnam(iwhere+1:iwhere+3) = "rms"
        endif
        fitnam = rmsnam(1:iwhere)//"fit"
        pltnam = rmsnam(1:iwhere)
        iwhere = index(tfnam(1),'.')
        gname = tfnam(1)(1:iwhere+3)//".fit"
        gname(1:1)='g'
        svsnam = "svs_"//tfnam(1)(2:iwhere+3)//"_fit.apr"

c open the command file
      open(unit=lucmd,file=cmdnam,status='old',iostat=ioerr)
      if( ioerr.ne.0 ) call report_stat('FATAL','ORBFIT','orbits/orbfit'
     .    ,cmdnam,'Error opening command file: ',ioerr)

c open the rms  file
      open(unit=lurms,file=rmsnam,status='unknown',iostat=ioerr)
      if (ioerr.ne. 0) call report_stat('FATAL','ORBFIT','orbits/orbfit'
     .    ,rmsnam,'Error opening the rms file: ',ioerr)
  
c open the SVS apr file
      open(unit=lusvs,file=svsnam,status='unknown',iostat=ioerr)
      if (ioerr.ne. 0) call report_stat('FATAL','ORBFIT','orbits/orbfit'
     .    ,svsnam,'Error opening output svs.apr file: ',ioerr)
                          
c open the fit-file
      open(unit=lufit,file=fitnam,status='unknown',iostat=ioerr)
      if (ioerr.ne. 0) call report_stat('FATAL','ORBFIT','orbits/orbfit'
     .    ,fitnam,'Error opening output fit file: ',ioerr)

c open the g-file
      open(unit=lug,file=gname,status='unknown',iostat=ioerr)
      if (ioerr.ne. 0) call report_stat('FATAL','ORBFIT','orbits/orbfit'
     .    ,gname,'Error opening output fit file: ',ioerr)

c open the t-files
            
      do i=1,ntfiles
         call topens ( tfnam(i),'old',lut(i),ioerr) 
         call report_stat('STATUS','ORBFIT','orbits/orbfit',tfnam(i)
     .              ,'Opened T-file: ',ioerr)
      enddo
        
c open the ut1, pole, and nutation files (WARNING, not FATAL since may 
c     not be needed if no terrestrial rotation parameters estimated
                          
      ut1fil = 'ut1.'
      open (unit=iut1,file=ut1fil,status='old',iostat=ioerr)
      if (ioerr .ne. 0) then
       call report_stat('WARNING','ORBFIT','orbits/orbfit',ut1fil,
     .    'Error opening ut1 table: ',ioerr)
      endif   
      polfil = 'pole.'
      open (unit=ipole,file=polfil,status='old',iostat=ioerr)
      if (ioerr .ne. 0) then
       call report_stat('WARNING','ORBFIT','orbits/orbfit',polfil,
     .    'Error opening pole table: ',ioerr)
      endif     
c     open the nutation table only if needed (no 'nbody' ephemeris)
      if( .not.fcheck('nbody') ) then
        nutfil = 'nutabl.'
        open (unit=inut,file=nutfil,status='old',iostat=ioerr)
        if (ioerr .ne. 0) then
         call report_stat('WARNING','ORBFIT','orbits/orbfit',nutfil,
     .      'Error opening nutation table: ',ioerr)
        endif 
      endif

      return
      end
