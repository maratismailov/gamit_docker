      Subroutine clkera
     .  ( kfile, sname, nojumpfit, jumpfit, clkfit, it0, t00, lprint
     .  , istats, epoch, rate, accel, cubic )
c
c     read clock log tables to determine epoch, rate, and acceleration
c
c     Input:
c       kfile              clock log file from MAKEX or MAKEK
c       sname              site name - for a check 
c       nojumpfit          order of polynomial for clock fit with no jumps allowed
c                             L=linear, Q=quadratic, C=cubic
c       jumpfit            order of polynomial for clock fit with jumps allowed (L, Q, C)
c       clkfit             choice of polynomial for output clock fit (for I-file in FIXDRV)
c                             N=no-jump, J=jump, blank=pick from residuals, 
c                             I=determine interactively (PLOTK call)
c       it0(1),(2),(3)     month, day and year of start time (GPST)
c       t00(1),(2),(3)     hour, minute and second of start time (GPST) 
c       lprint             unit for print file (=9 for FIXDRV, =0 for PLOTK)

c     Output   
c       istats             = 0, then clock epoch, rate, acc in table
c                            1, not used (but used in clkprm)
c                            2, then no clock info in table
c       epoch              epoch offset in seconds
c       rate               rate offset in seconds per second
c       accel              acceleration offset in 1/sec
c       cubic              cubic term in 1/sec/sec   
c

c     The approach here is to find jumps by doing a linear, cubic, or quadratic fit 
c     to each segment.  A segment is defined as the points between jumps.  Once a 
c     jump is identified, it is modeled as a step function,  whose magnitude constitutes 
c     a new parameter.  

c     The current version of this routine, allowing jumps in the clock values,
c     was written by K. Feigl (circa 1990).  Modified by R. King 970522
c     Modified for version 2 k-files R. King 160901

      implicit none
        
c     array dimensions  

c     maximum number of models to fit (must change printout if > 2)
      integer maxfit
      parameter (maxfit=2)
c     maximum number of poly coefficients
      integer maxp
      parameter (maxp=4)
c     maximum number of jumps  (240 allows every 6 minutes for 24 hrs)
      integer maxj
      parameter (maxj=240)
c     maximum number of parameters
      integer maxm
      parameter (maxm=maxj+maxp)
     
c     number of terms in polynomials actually fit (order+1)

c     single polynomial (no jumps)
      integer mall
c     each segment 
      integer mseg

c     tolerance for jump detection
      real*8 jumptol
      parameter (jumptol=500.d-6)

c     character variables
      character*1 lowerc,nojumpfit,jumpfit,clkfit,pikfit
      character*4 sname,site
      character*9 poly1,poly2
      character*16 kfile
      character*32 buff32
      character*51 kfmt 
      character*105 kheader 
      character*128 aline(13),prog_name
      character*256 message  

c     k-file version
      real*4 kversion

c     unused variables from the k-file
      character*1 gnss 
      integer*4 iwkn
      real*8 sow,pr,svclk

c     order of polynomial for no-jumps and jump, respectively    
c     should correspond to values of 'mall' and 'mseg' above
c      parameter (poly1='Linear   ',poly2='Cubic    ') 
                     
c     sampling interval
      real*8 dt
      parameter (dt=120.d0)
c     time tags
      real*8 tj(maxj),tmtj
c     time of initial epoch of X- or C-file
      integer it0(3)
c     clock offset
      real*8 clkoff
c     time read from k-file
      integer*4 month,ihr,min,iyr,idyoyr
      real*8 sec
c     values to return for S-file
      real*8 epoch,rate,accel,cubic
c     local values
      real*8 epoch1(maxfit),rate1(maxfit),accel1(maxfit),cubic1(maxfit)
     .     , rcond(maxfit)
      real*8 t00(3)
      real*8 resid0,resid1,resid2,rmax(maxfit),rms(maxfit),fit1,fit2
      real*8 resid0_last
      real*8 epoch0,rate0,acc0,cub0
      integer ifit,njump,mm,m1,len,rcpar

c     times: data, ref. epoch for poly., first data point
      real*8 tsec,tsec0
      integer*4 jd,jd0
c     time when clock is first synchronized
      real*8 tsynch
c     seconds between current entry and initial epoch of X- or C-file
      real*8 tmt0
c     time of first k-file entries-save to skip on second pass if not synched
      real*8 tmt0_first

c     for Allen variance
      real*8 alsum,allen,allentol,oldoff,oldtmt,offset
      integer nallen

c     needed for LINPACK
c     for complete fit with jumps
      real*8 a1(maxm,maxm),a2(maxm,maxm),b1(maxm),b2(maxm)
      real*8 x1(maxm),scale(maxm),determ(maxm)
      integer info
c     for segment fit
      real*8 aa1(maxp,maxp),bb1(maxp),aa2(maxp,maxp)
      real*8 bb2(maxp),x2(maxp)
      real*8 rcond0,dterm2(2)

c     partial derivatives
      real*8 pp(maxm),pp1(maxp)

      integer ioerr,i,j,n,istats,nprn,iday,lprint
      integer nseg,ijump,ncall
      integer*4 julday

c**   debug
      integer nn
      logical first,synchd,rlook,lask,site_check,first_epoch,skip_first
      character*256 line 


      save rlook,ncall

      data ncall /0/

c     get calling program name and m-file name for report_stat
      len = rcpar(0,prog_name) 

c     set the order of the polynomial and the summary headers for no-jump 
c     and jump-allowed fits (the order for jump-allowed fits is the same 
c     as for each segment)   
      if( lowerc(nojumpfit).eq.'l' ) then
         mall = 2
         poly1 = 'Linear---'
      elseif (lowerc(nojumpfit).eq.'q' ) then
         mall = 3
         poly1 = 'Quadratic'
      elseif (lowerc(nojumpfit).eq.'c' ) then
         mall = 4
         poly1 = 'Cubic----'
      else     
        call report_stat('FATAL',prog_name,'lib/clkera',' '
     .    , 'Invalid no-jumps polynomial ',0)
      endif
      if( lowerc(jumpfit).eq.'l' ) then
          mseg = 2
          poly2 = 'Linear---'
      elseif( lowerc(jumpfit).eq.'q' ) then  
          mseg = 3
          poly2 = 'Quadratic'
      elseif( lowerc(jumpfit).eq.'c' ) then
          mseg = 4
          poly2 = 'Cubic----'
      else   
        call report_stat('FATAL',prog_name,'lib/clkera',' '
     .    , 'Invalid with-jumps polynomial ',0)
      endif


c     question for interactive use (PLOTK)
      if (ncall .eq. 0 .and. lowerc(clkfit) .eq. 'i') then
         print '(2a,$)','Do you wish to see the residuals for the '
     .,  'station clock polynomial fit?'
         rlook = lask()
      else
         rlook = .false.
      endif
      ncall = ncall + 1

c     default - no clock information
      istats=2

c     check site code only once
      site_check = .true.

c     logical to allow a second pass
      skip_first = .false.

c     Begin initializing variables for clock fit
c     If the fit has failed the first time through, come here a second time
c     and this time skip the first point in the k-file, which may be bad
   1  continue
c
c     initial values
      epoch0 = 0.d0
      rate0 = 0.d0
      acc0 = 0.d0
      cub0 = 0.d0
      tsynch = 0.d0
      tmt0_first = 0.d0
      if ( mall .eq. 2 ) poly1 = 'Linear---'
      if ( mall .eq. 3 ) poly1 = 'Quadratic'
      if ( mall .eq. 4 ) poly1 = 'Cubic----'
      if ( mseg .eq. 2 ) poly2 = 'Linear---'
      if ( mseg .eq. 3 ) poly2 = 'Quadratic'
      if ( mseg .eq. 4 ) poly2 = 'Cubic----'

c     number of jumps
      njump = 0

c     clock is not synched at start.
      synchd = .false.

c     set logical for first_epoch in k-file - to be used for skipping if problems
      first_epoch = .true.

c     initialize sum for allen variance
      oldoff = 0.0d0
      alsum  = 0.0d0
      allen  = 1.0d0
      offset = 1.0d0
      nallen = 0

c     time of jumps
      do i=1,maxj
         tj(i) = 99999.0d0
      enddo
c
      do i=1,maxfit
         epoch1(i) = 0.0d0
         rate1(i)  = 0.0d0
         accel1(i) = 0.0d0
         cubic1(i) = 0.0d0
      enddo
c
      epoch = 0.d0
      rate  = 0.d0
      accel = 0.d0

c     initialize normal equations
c       full data set
      n = 0
      do j=1,maxm
         x1(j) = 0.0d0
         b1(j) = 0.0d0
         pp(j) = 0.0d0
         do i=1,maxm
            a1(i,j) = 0.0d0
            a2(i,j) = 0.0d0 
          enddo
      enddo 
c        between-jumps segment
      nseg = 0
      do j=1,maxp
         bb1(j) = 0.0d0
         do i=1,maxp
            aa1(i,j) = 0.0d0
            aa2(i,j) = 0.0d0
         enddo
      enddo

c     first point in segment gets special treatment
      first = .true.
      resid0 = 0.0d0
      resid0_last = 0.0d0
                                                   
c     open clock log (k-file) and determine the version
      call lowers(kfile)    
      open (unit=20,file=kfile,status='old',iostat=ioerr)
      if  (ioerr .eq. 0) then 
        if( lprint.gt. 0 ) then 
          write(lprint,'(a,a16)') 
     .      'Estimating receiver clock poly from ', kfile
        else    
          write(*,'(a,a16)') 
     .      'Estimating receiver clock poly from ', kfile
         endif
      else
        call report_stat('WARNING',prog_name,'lib/clkera',kfile
     .          , 'Could not open K-file: ',0)
        return
      endif
      read(20,'(a)',iostat=ioerr) kheader                                       
      if( ioerr.ne.0) call report_stat('FATAL',prog_name
     .     ,'lib/clkera',kfile,'Error reading k-file header',ioerr)
      if( kheader(1:3).eq.'Ver' ) then
        read(kheader(5:8),'(f4.0)',iostat=ioerr) kversion
        if(ioerr.ne.0) call report_stat('FATAL',prog_name
     .     ,'lib/clkera',kfile,'Error reading k-file version',ioerr)
c       skip the column-header line, then read the format
        read(20,'(/,a)',iostat=ioerr) kfmt 
      else
        kversion = 1.0 
      endif                           
cd      print *,'kversion ',kversion
cd      print *,'kfmt ',kfmt 

c     compute start time in Julian day and sec of day
      jd0 = julday(it0(1),it0(2),it0(3))
      tsec0 = t00(1)*3600.0d0 + t00(2)*60.d0 + t00(3)

c     initialize time difference from initial epoch
      tmt0=0.d0                                

c------------------------------------------------------------------------------------

c     begin loop on records of k-file

c     come here to read one record
  200 continue
      oldoff = clkoff
      oldtmt = tmt0     
      if( kversion.lt.2.0 ) then       
        read(20,100,end=997,iostat=ioerr)
     .    site,nprn,iyr,idyoyr,ihr,min,sec,clkoff
  100   format(1x,a4,2x,i2,2x,i4,1x,i3,1x,i2,1x,i2,f9.4,30x,f12.8) 
        call fix_y2k(iyr)
cd      print *, site,nprn,iyr,idyoyr,ihr,min,sec,clkoff
c       a blank line indicates end of file, too
        if (site .eq. '   ') goto 997 
        if (ioerr .ne. 0) 
     .    call report_stat('WARNING',prog_name,'lib/clkera',kfile
     .          , 'Error reading Version 1 K-file for values: ',0)
cd      print *, site,nprn,iyr,idyoyr,ihr,min,sec,clkoff
c       check if site codes match at first point (now warning only)
        if( site_check ) then
          call uppers(site)
          call uppers(sname)
          if(site(1:4).ne.sname(1:4))  then
             write(message,'(2a,a4)')
     .         '  Site code does not match K-file site, ',sname
             call report_stat('WARNING',prog_name,'lib/clkera',' '
     .                       , message,0)
          endif
          site_check = .false.
        endif
      else       
        read(20,'(a)',iostat=ioerr) line
        if( ioerr.eq.0.and.line(1:1).ne.' ' ) then
cd          print *,'line ',line 
          read(line,kfmt,iostat=ioerr) iyr,idyoyr,ihr,min,sec
     .                        , iwkn,sow,gnss,nprn,pr,svclk,clkoff
cd           print *,iyr,idyoyr,ihr,min,sec,iwkn,sow,gnss,nprn,pr,svclk
cd     .             ,clkoff
        elseif( ioerr.eq.-1 ) then
          goto 997
        else
          call report_stat('FATAL',prog_name,'lib/clkera',kfile
     .          ,'Error reading Version 2 K-file for values: ',ioerr)
        endif
      endif   

c     convert K-file UTC day of year to GPST month day year
      call utc2gps( idyoyr,iyr,ihr,min,sec )
      call monday(idyoyr,month,iday,iyr)

c     compute data  time in Julian day and sec of day, store if first
      jd = julday(month,iday,iyr)
      tsec =  ihr*3600.0d0 + min*60.0d0 + sec

c     compute time difference in seconds
      tmt0 = (jd-jd0)*86400.d0 + tsec - tsec0
      if( first_epoch) then
         tmt0_first = tmt0
         first_epoch = .false.
      endif

c     Fix for problem of bad first epoch.   If we failed to synch on the
c     first pass through, we're now trying a second time, in which we
c     skip the first epoch.
      if( skip_first .and. (tmt0 .le. (tmt0_first+10.) ) ) goto 200

c     compute 2-pt Allen variance using pairs of points separated by at least 10 s
      if (dabs(tmt0-oldtmt) .gt. 10.d0) then
         nallen = nallen + 1
         if (nallen .gt. 1) then
c            alsum = alsum + (clkoff - oldoff)*(clkoff - oldoff)
            alsum = alsum + ((clkoff - oldoff)/(tmt0-oldtmt))**2
            allen = dsqrt(alsum/(2.0d0 * (nallen-1))) 
            offset = dsqrt((allen**2)*2)*(tmt0-oldtmt)
         endif
      endif
                      
c     If the Allen variance is greater than a value corresponding to the jump 
c     tolerance, assume that we don't yet have a valid rate estimation
c        
c         2-pt Allen variance = jumptol/dt/sqrt(2)
c
c         so., eg., if jumptol = 500 microseconds and dt = 120 s
c                      allentol = 2.95d-6
      allentol = jumptol/dt/sqrt(2.)
cd      print *,'allentol allen ',allentol,allen
      if (.not. synchd) then
         if (dabs(allen) .le. allentol ) then 
cd           print*,'synchd by allen',dabs(allen),dabs(offset),tmt0
            synchd = .true.
            tsynch = tmt0  
         elseif (0.0005d0 .ge. dabs(offset) ) then
cd           print*,'synchd by offset',dabs(offset),dabs(allen),tmt0
            synchd = .true.
            tsynch = tmt0  
         endif
      endif

      if (synchd) then
         n = n + 1  
c        nseg counts only non-synchronous points for performing cubic fit (to avoid bad fits)
         if( dabs(tmt0-oldtmt) .gt. 10.d0) nseg = nseg + 1
cd        print *,'synched tmt0 oldtmt nseg ',tmt0,oldtmt,nseg
      else
         goto 200
      endif  
c     --if not synched to read another record


c     time since last jump
      if (njump .eq. 0) then
         tmtj = tmt0
      else
         tmtj = tmt0 - tj(njump)
      endif       
cd      print *,'njump tmtj first ',njump,tmtj,first

      if (.not. first) then
c        form current residual from segment fit
        resid0 = clkoff
     .           -epoch0
     .           -rate0*tmtj
     .           -acc0*tmtj*tmtj
     .           -cub0*tmtj*tmtj*tmtj
cd         print *, 'Poly',epoch0, rate0, acc0, cub0, tmtj        
      else
         resid0 = 0.0d0
      endif
cd        print *,'tmt0 tmtj resid0 resid0_last '
cd     .  ,tmt0,tmtj,resid0,resid0_last
                
c**RWK 010725: Check for gross outliers (jumps > 5 ms) to catch TI 4100 problems
c**RWK 040115: Later receivers (e.g. TRMSSE circa 1997) can have clock jumps of over 100 ms,
c**             followed by continued steady drift; the 5 ms test prevents inserting a jump.
c**             Change this to 300 ms and see if most data are handled correctly.
c**rwk 040115      if( dabs(resid0-resid0_last) .gt. 5.d-3 ) then
       if( dabs(resid0-resid0_last) .gt. 300.d-3 ) then
cd         print *,'TEST Outlier tmtj resid0 resid0_last '
cd     .           ,tmtj,resid0,resid0_last
         write(message,'(a,1pd12.4)')
     .       'Outlier in k-file, residual =',resid0
         call report_stat('WARNING',prog_name,'lib/clkera',kfile
     .         ,message,0 )
         n = n -1
         if( dabs(tmt0-oldtmt) .gt. 10.d0) nseg = nseg - 1 
         goto 200
      endif
c
c     If the change in the residual is large, assume current point to 
c     be first point of next jump.  Start a new segment.
      if (dabs(resid0-resid0_last) .gt. jumptol) then  
         if (njump .lt. maxj) then 
cd           print*,'resid0, jumptol ',resid0,jumptol
            njump = njump + 1
         else              
           write(message,'(a,i4)')
     .         'Number of clcok jumps exceeds maxj = ',maxj
           call report_stat('FATAL',prog_name,'lib/clkera'
     .                     ,kfile,message,0) 
         endif
c        record the time of the jump
         tj(njump) = tmt0              
cd         print *,'JUMP njump tmt0 resid0 resid0_last '
cd     .         , njump,tmt0,resid0,resid0_last
         tmtj = 0.d0
c        zero for next segment estimate
         nseg = 0
         resid0 = 0.d0
         first = .true.
c        zero segment normal equation
cd         print *,'zeroing segment neqs'
         do j=1,maxp
            bb1(j) = 0.0d0
            do i=1,maxp
               aa1(i,j) = 0.0d0
            enddo
         enddo
      endif
      resid0_last = resid0

c      if (dabs(resid0) .gt. jumptol) then
c         buff16 = 'Jump'
c      else
c         buff16 = '    '
c      endif
c      write (6,'(5(1pe20.10),1x,a)') tmt0,clkoff,rate0,
c     .         tmtj,resid0,buff16
 

c     Increment the normal equations 

c     For bookkeeping convenience we include in the full-span
c     normal equations polynomial terms up the maximum allowed 
c     (maxp=4 to accommodate a cubic), but then replace the diagonals
c     by 1.0 and the off-diagonals by 0. for terms we want to fix.
c     For each segment, inverted at every point to allow for predictions,
c     we fill out only the terms needed (mseg=1,2,3, or 4).  

c     partials: epoch, rate, acceleration, step function
c               for the full span and segment. 
      pp(1)=1.d0  
      pp(2) = tmt0
      pp(3) = tmt0*tmt0
      pp(4) = tmt0*tmt0*tmt0
      pp1(1) = 1.d0
      if( mseg.ge.2 ) pp1(2) = tmtj
      if( mseg.ge.3 ) pp1(3) = tmtj*tmtj
      if( mseg.ge.4 ) pp1(4) = tmtj*tmtj*tmtj
      do i = 1,njump
         if (tmt0 .lt. tj(i)) then
            pp(maxp+i)=0.0d0
         else
            pp(maxp+i)=1.0d0
         endif
      enddo 

      mm = maxp + njump
   
cd      print *,'mall mm pp ',mall,(pp(i),i=1,10)
cd      print *,'mseg mm pp1 ',mseg,(pp1(i),i=1,10)

c     fill normal equations for all points read thus far
c     upper triangle of left hand side
      do j = 1,mm
         do i = 1,j
            a1(i,j) = a1(i,j) + pp(i)*pp(j)
         enddo
      enddo
c     right hand side
      do i = 1,mm
         b1(i) = b1(i) + pp(i)*clkoff
      enddo

c     fill normal equations for a polynomial fit to this segment
c     left hand side
      do j = 1,mseg
         do i= 1,mseg
            aa1(i,j) = aa1(i,j) + pp1(i)*pp1(j)
         enddo
      enddo
c     right hand side
      do i = 1,mseg
         bb1(i) = bb1(i) + pp1(i)*clkoff
      enddo

c     solve for segment rate only if we have enough data
cd      print *,'nseg mseg tmtj ',nseg,mseg,tmtj
c        -this test formerly nseg.gt.mseg+2 : tradeoff between robustness and missing jumps
      if( nseg.gt.mseg .and.  tmtj .gt. 300.d0) then
c        save to avoid destruction
         do j=1,mseg
            bb2(j) = bb1(j)
            do i=1,mseg
               aa2(i,j) = aa1(i,j)
            enddo
         enddo
c        scale matrix
         do i = 1,mseg
            scale(i) = 1.0d0/dsqrt(aa1(i,i))
         enddo
         do j = 1,mseg
            do i = 1,mseg
               aa1(i,j) = scale(i)*scale(j)*aa1(i,j)
            enddo
         enddo 
cd        print *,'solving segment neqs mseg ',mseg

c        invert segment normal matrix AA1 using linpack (blas 1) 
c        routines assuming that the matrix is positive definite 
c        (and thus symmetric). DPODI returns inverse in upper triangle, 
c        so the full matrix must be filled in

         call dpoco(aa1,maxp,mseg,rcond0,x2,info)
         if (info.eq.0 .and. rcond0 .ge. 1.0d-16) then
            call dpodi(aa1,maxp,mseg,dterm2,11)
c           fill in lower half and rescale
            do i = 1,mseg
               do j = 1, i
                  aa1(j,i) = scale(i)*scale(j)*aa1(j,i)
                  aa1(i,j) = aa1(j,i)
               enddo
            enddo
c           multiply inv(AA1) * BB1 = X2
            call dgemv ('N', mseg, mseg
     .                 , 1.d0, aa1, maxp
     .                 , bb1, 1
     .                 , 0.d0, x2, 1)
            epoch0 = x2(1)
            rate0  = x2(2)
	    if( mseg.gt.2 ) then
               acc0   = x2(3)
	    endif
	    if( mseg.gt.3 ) then
               cub0   = x2(4)
	    endif
c           raise the flag that the solution is OK to use
c           for predictions
cd           print*,'setting first = false '
            first = .false.
         else
c**            write (6,*) 'Matrix for segment fit is ill conditioned.',
c**     .      ' RCOND0 = ',rcond0
         endif

c        restore segment normal equations
         do j=1,mseg
            bb1(j) = bb2(j)
            do i=1,mseg
               aa1(i,j) = aa2(i,j)
            enddo
         enddo
      else
cd         print *,'no solution tmtj nseg ',tmtj,nseg
      endif

c     read next data point
      go to 200     

c---------------------------------------------------------------------

c     come here on end of K-file
  997 continue    


c     now solve for two models (or add code for more) 
       
c        ifit = 1   :  linear fit (mall = 2), no jumps
c        ifit = 2   :  cubic  fit (mseg = 4), with jumps

c     solve the normal equations only if there are 3 or more points

      if (n .gt. 2) then  
c     --ELSE for this IF after 'goto 40', ENDIF before '998'

c        copy A1 into A2 to avoid destroying it
c        copy B1 into B2 to avoid destroying it
         do j=1,maxm
            b2(j) = b1(j)
            do i=1,maxm
               a2(i,j) = a1(i,j)  
            enddo
         enddo

         do ifit = 1,2
            if (ifit .eq. 1) then
               mm = mall             
            else
               mm = maxp + njump 
            endif
c           scale matrix
            do i = 1,mm
               scale(i) = 1.0d0/dsqrt(a1(i,i))
            enddo
            do j = 1,mm
               do i = 1,mm
                  a1(i,j) = scale(i)*scale(j)*a1(i,j)
               enddo
            enddo     
c           For second fit, fix polynomial parameters not requested  
            if( ifit.eq.2 ) then
               if (mseg.lt.maxp ) then
                 do i=1,maxp 
                   m1 = maxp
                   if( i.gt.mseg ) m1 = mm
                   do j=mseg+1,m1
                     if( i.eq.j ) a1(i,j) = 1.d0
                     if( i.ne.j)  a1(i,j) = 0.d0
                   enddo
                 enddo   
                 m1 = maxp
                 do i=mseg+1,m1
                   b1(i) = 0.d0
                 enddo
               endif
            endif
cd           print *,'test ifit mm solving neqs ',ifit,mm  
c            do i=1,8
c              write(*,'(8d14.5)') (a1(i,j),j=1,8)
c            enddo
c            write(*,'(8d14.5)') (b1(i),i=1,8)
c           invert normal matrix A1 using linpack (blas 1) routines
c           assuming that the matrix is positive definite (and thus
c           symmetric). DPODI returns inverse in upper triangle, so
c           the full matrix must be filled in

            call dpoco(a1,maxm,mm,rcond(ifit),x1,info)

c            if (info.eq.0 .and. rcond(ifit)+1.0d0 .ne. 1.0d0) then
            if (info.eq.0 .and. rcond(ifit) .ge. 1.0d-16) then

               call dpodi(a1,maxm,mm,determ,11)

c              fill in lower half and rescale
               do i = 1,mm
                  do j = 1, i
                     a1(j,i) = scale(i)*scale(j)*a1(j,i)
                     a1(i,j) = a1(j,i)
                  enddo
               enddo

c              multiply inv(A1) * B1 = X1
               call dgemv ('N', mm, mm
     .                    , 1.d0, a1, maxm
     .                    , b1, 1
     .                    , 0.d0, x1, 1)

c              save the estimates
               epoch1(ifit) = x1(1)
               rate1(ifit)  = x1(2)
               if (ifit .eq. 2) then
                  accel1(ifit) = x1(3)
                  cubic1(ifit) = x1(4)
               else
                  accel1(ifit) = 0.0d0
                  cubic1(ifit) = 0.0d0
               endif
            else
             if (info .ne. 0) then
c             a1 is not positive definite.   This case
c             only seems to arise due to degeneracies.
c             There shouldn't be any such cases for this problem
              write (message,'(a,i1,a,i1,a)') 'Matrix for fit ',ifit
     .             ,' of order ',mm,' is not positive definite'      
              call report_stat('WARNING',prog_name,'lib/clkera',' ',
     .             message,0)
             endif
c            if (rcond(ifit)+1.0d0 .eq. 1.0d0) then
               if (rcond(ifit) .lt. 1.0d-16) then  
                 write (message,'(a,i1,a,i1,a)') 'Matrix for fit ',ifit
     .                 ,' of order ',mm,' is ill conditioned' 
                  call report_stat('WARNING',prog_name,'lib/clkera',
     .                 ' ',message,0)
                endif
                epoch1(ifit) = 0.0d0
                rate1(ifit)  = 0.0d0
                accel1(ifit) = 0.0d0
                cubic1(ifit) = 0.0d0
                do i =1,njump
                  x1(maxp+i) = 0.0d0
                enddo
             endif

c           reload matrices
            do j=1,maxm
               b1(j) = b2(j)
               do i=1,maxm
                  a1(i,j) = a2(i,j) 
               enddo
            enddo

         enddo
c        end do for loop on ifit


c        Compute the residuals by re-reading the data
         rewind 20             
c        skip the headers         
         if( kversion.ge.2.0 ) read(20,'(///)') 
          

         if (rlook) then
            write(*,'(a)')
            write(*,'(a,1x,a9,a,4x,a9,a)')
     .         '    Time    Sat  Observed offset  '
     .          ,poly1,' no jumps  ',poly2,' with jumps'
            write(*,'(2a)')
     .   '     sec    (prn)    (micros)    Fit         Residual'
     .        ,'    Fit        Residual'
            write(*,'(1x)')
         endif

         ijump = 1 
         do i =1,2
           rmax(i) = 0.d0
           rms(i) = 0.0d0
         enddo
          
         nn = 0  
           
   40    continue   
c        --come here to read each record
         if( kversion.eq.1.0 ) then 
           read(20,100,end=998,iostat=ioerr) 
     .        site,nprn,iyr,idyoyr,ihr,min,sec,clkoff 
           if( iyr.lt.1900 ) iyr = iyr + 1900
           if( ioerr.ne.0 ) 
     .       call report_stat('WARNING',prog_name,'lib/clkera',kfile
     .            , 'Error reading K-file for residuals: ',0)
           nn = nn + 1
c          blank line indicates end of file, too
           if (site .eq. '   ') goto 998
         elseif( kversion.ge.2.0 ) then  
           read(20,'(a)',iostat=ioerr) line
           if( ioerr.eq.0.and.line(1:1).ne.' ' ) then
cd             print *,'line ',line 
             read(line,kfmt,iostat=ioerr) iyr,idyoyr,ihr,min,sec
     .                        , iwkn,sow,gnss,nprn,pr,svclk,clkoff
cd             print *,iyr,idyoyr,ihr,min,sec,iwkn,sow,gnss,nprn,pr,svclk
cd     .              ,clkoff
             nn = nn + 1 
             elseif( ioerr.eq.-1 ) then
               goto 998
             else
               call report_stat('FATAL',prog_name,'lib/clkera',kfile
     .          ,'Error reading Version 2 K-file for values: ',ioerr)
           endif
         endif   

c        convert K-file UTC day of year to GPST month day year
         call utc2gps( idyoyr,iyr,ihr,min,sec )
         call monday(idyoyr,month,iday,iyr)

c        compute data  time in julian day and sec of day
         jd   =  julday(month,iday,iyr)
         tsec =  ihr*3600.0d0 + min*60.0d0 + sec

c        compute time difference in seconds
         tmt0 = (jd-jd0)*86400.d0 + tsec - tsec0

c        compute residuals seconds  
         fit1   = epoch1(1)
     .          + rate1(1)*tmt0
     .          + accel1(1)*tmt0*tmt0
     .          + cubic1(1)*tmt0*tmt0*tmt0
         resid1 = clkoff-fit1
         fit2   = epoch1(2)
     .          + rate1(2)*tmt0
     .          + accel1(2)*tmt0*tmt0
     .          + cubic1(2)*tmt0*tmt0*tmt0
         do i = 1,njump
cd            print *,'fit2 tmt0 i tj ',fit2,tmt0,i,tj(i)
            if (.not. (tmt0 .lt. tj(i))) then
               fit2 = fit2 + x1(maxp+i)
            endif
cd            print *,'fit2 maxp i x1 ',fit2,maxp,i,x1(maxp+i)
         enddo
         resid2 = clkoff - fit2

c        display in microseconds with notation if jump or not used
         if (tmt0 .lt. tsynch) then
            write (buff32,'(1x,a)') 'Not synched.  (value ignored)'
            resid1 = 0.d0
            resid2 = 0.d0
         else
            buff32 = ' '
         endif
cd         print *,'At display tmt0 tj(ijump) ijump ',tmt0,tj(ijump),ijump
         if (tmt0-tj(ijump) .ge. 0.d0 .and.
     .       tmt0-tj(ijump) .lt. 10.d0) then
            write (buff32,'(1x,a,f10.2)') 'Jump ',x1(maxp+ijump)*1.d6
            ijump = ijump + 1
         endif          
         if (dabs(resid2).gt.5.d-3 ) then  
            write (buff32,'(1x,31a)') 'Outlier gt 5 ms (value ignored)'  
            resid1 = 0.d0
            resid2 = 0.d0
         endif
         if (rlook) then
          write(6,'(f11.2,2x,i2,5(2x,f10.2),a)') tmt0,nprn,clkoff*1.0d6
     .         ,fit1*1.d6,resid1*1.d6,fit2*1.d6,resid2*1.0d6,buff32
         endif
         rmax(1) = max(dabs(resid1),rmax(1))
         rmax(2) = max(dabs(resid2),rmax(2))
         rms(1) = rms(1) + resid1*resid1
         rms(2) = rms(2) + resid2*resid2
         go to 40  
c        --got read another record
                                  

      else

c        there are not enough data  

c        first try skipping to first point in the k-file
         if( .not.skip_first ) then   
            write(message,'(2a,d10.2)') 
     .            'Fewer than 3 synched data points in K-file: '
     .           ,'try again, skipping first epoch. Allen variance = '
     .           ,dabs(allen)
            call report_stat('WARNING',prog_name,'lib/clkera'
     .                   ,kfile,message,0) 
            skip_first = .true.
            rewind 20
            goto 1
         else   
            write(message,'(2a,d10.2)') 
     .            'Fewer than 3 synched data points in K-file: clock '
     .           ,'coefficients set = 0. in I-file. Allen variance = '
     .           ,dabs(allen)
            call report_stat('WARNING',prog_name,'lib/clkera'
     .                   ,kfile,message,0) 
            do i =1,2
              rmax(i) = 0.0d0
              epoch1(i) = 0.d0
              rate1(i) = 0.d0
              accel1(i) = 0.d0
              cubic1(i) = 0.0d0
              rcond(i) = 0.d0
            enddo
            write(message,'(2a)') 
     .            'Try remaking k-file using a k-file sampling interval' 
     .           ,' equal to your x-file sampling interval: ' 
            call report_stat('WARNING',prog_name,'lib/clkera',kfile
     .                   ,message,0) 
         endif
      endif

c     come here on end of file or blank line
 998  continue
            
      rms(1)=dsqrt(rms(1)/max(n-1,1))
      rms(2)=dsqrt(rms(2)/max(n-1,1))

c     come here to display statistics again
 999  continue                                                          
      write(aline(1),'(a,f5.2,a10)') 
     .       'Receiver clock 2-pt allen variance: '
     .      ,dabs(allen)*1000000.d0,' micro secs'       
      write(aline(2),'(1x)') 
      write(aline(3),'(19x,a,6x,a9,3x,a9)') 'Clock parameters:'
     .        ,poly1,poly2
      write(aline(4),'(42x,a)') 'no jumps    with jumps'   
      write(aline(5),130) kfile,' epoch (sec):     ',epoch1(1),epoch1(2)
      write(aline(6),130) kfile,' rate  (s/s):     ',rate1(1),rate1(2)     
      write(aline(7),130) kfile,' accel (1/s):     ',accel1(1),accel1(2)   
      write(aline(8),130) kfile,' cubic (1/s/s):   ',cubic1(1),cubic1(2)  
      write(aline(9),130) kfile,' 1/condition:     ',rcond(1),rcond(2)      
      write(aline(10),132) kfile,' jumps      :     ',njump        
      write(aline(11),132) kfile,' synched N  :     ',n       
      write(aline(12),133)    kfile,' max res (micros):',
     .                rmax(1)*1.d6,rmax(2)*1.d6  
      write(aline(13),133)    kfile,' rms res (micros):',
     .                rms(1)*1.d6,rms(2)*1.d6  
  130 format(a16,1x,a19,3x,2(1pe12.3))
  132 format(a16,1x,a19,3x,12x,i12)
  133 format(a16,1x,a19,3x,2(f12.3))
      do i=1,13
        if( lprint.gt.0 ) then 
           write(lprint,'(a)') aline(i)
        else   
           write(*,'(a)') aline(i)
        endif
      enddo      

c     if estimation failed, use other fit.

      if (dabs(epoch1(2)) .gt. 0.0d0) then
         if (rms(1) .gt. 1.d-6) then
            pikfit = 'j'
            if (rms(2) .gt. 100.d-6 ) then
              if( rms(2).lt.1.d-3) then
                call report_stat('WARNING',prog_name,'lib/clkera',' '
     .      , 'Clock residuals > 100 microseconds - possible problem',0)
              else
               call report_stat('WARNING',prog_name,'lib/clkera',kfile
     .      , 'Clock residuals >  millisecond - edit the K-file:',0)
              endif
               call report_stat('WARNING',prog_name,'lib/clkera',kfile
     .      , '--See Section 5.2 of GAMIT manual ',0)
            endif
         else
            pikfit = 'n'
         endif
      else
         pikfit = 'n'
      endif
            
      if( lprint.gt.0) write(lprint,'(1x)')
      if( lprint.le.0) write(*,'(1x)')
      if (lowerc(clkfit) .eq. 'i') then
         write (6,1001) pikfit
 1001    format  ('Choose fit without (N) or with (J) jumps.',/
     .         'I recommend ',a1,'. Do you agree?')
         if (.not.lask()) then
            if (lowerc(pikfit) .eq. 'j') then
                pikfit = 'n'
            else
                pikfit = 'j'
            endif
         endif
      else if (lowerc(clkfit).eq. 'n') then
         pikfit = 'n'
      else if (lowerc(clkfit).eq. 'j') then
         pikfit = 'j'
      endif

      if(lowerc(pikfit) .eq. lowerc('n')) then
         epoch = epoch1(1)
         rate  = rate1(1)
         accel = accel1(1)
         cubic = cubic1(1)  
         write (aline(1),'(a,a9,a)') 'CLKERA: using ',poly1
     .         ,' fit with no jumps removed' 
         if( lprint.gt.0 ) then
           write(lprint,'(a)') aline(1)
         else
           write(*,'(a)') aline(1)
         endif
      else if(lowerc(pikfit) .eq. lowerc('j')) then
         epoch = epoch1(2)
         rate  = rate1(2)
         accel = accel1(2)
         cubic = cubic1(2)
         write (aline(1),'(a,a9,a)') 'CLKERA: using ',poly2
     .         ,' fit with jumps removed '
         if( lprint.gt.0 ) then
           write (lprint,'(a)') aline(1)
         else
           write(*,'(a)') aline(1)
         endif
      else  
         call report_stat('FATAL',prog_name,'lib/clkera',' '
     .    , 'Invalid clock polynomial selection ',0)
         goto 999
      endif

c     got clock information from table
      istats=0
      close(unit=20)
      return
      end




