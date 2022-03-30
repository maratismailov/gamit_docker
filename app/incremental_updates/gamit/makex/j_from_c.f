      subroutine j_from_c (afmt,nprns,iprns,debug,iprndb)

c     Write a J-file based on phase residuals from C-files
c     at sites with atomic oscillators
c     The model fits three terms from 3 data points as described by
c     Feigl, 1991 (thesis) and Feigl et al. [GRL 18, 1289, 1991]

c     Feigl 1990; Dong/King June 1992

      implicit none

      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/errflg.h'

      character*60     afmt
      character*80 afile(maxsit),wildcard,pickfn

c     Order of polynomial and number of points --must change also
c     in AVGCLK, RATE_EST, and ALLANV
      integer mm,nn
      parameter ( mm=3, nn=3 )

c     Quantities from the C-files
c       the L1, O minus C phase data in cycles
      real*8 pl1(maxepc,maxsat,maxsit)
c       the theoretical delay in seconds (not used)
c**      real*8 tau(maxepc,maxsat,maxsit)
c       the accompaning error flag
      integer*4 ierrf (maxepc,maxsat,maxsit)
c       and the accompanying time tag (tsend-t0 in seconds)
      real*8 ttag (maxepc)
c       receiver clock values from C-file
      real*8 clkepc,clkrat,clkacc

c     Reference (Broadcast) satellite clock quantities
c     units should be s,s/s,1/s respectively
      real*8 svcepc,svcrat,svcacc,valid

c     Start epoch, current epochs
      integer*4 jd,jd0,jdc0,jdj0
      real*8    sod,sod0,sodc0,sodj0

c     Times from initial epoch (C-file epochs, center epoch, J-file epoch)
      real*8 t(maxepc,maxsit),ti,tc

c     Quantities for estimation
c       partial derivatives
      real*8 pp(mm,nn)
c       estimated phase polynomial coefficients
      real*8 c(mm,maxsit),c1,c2
c       estimated satellite oscillator terms
      real*8 a0,a1,a2

c     number and logical flag for good files at each epoch
      integer ngoodf
      logical lgoodf(maxsit)

c     count of good estimates and outliers in the polynomial estimates
      integer*4 goodest(maxsat),outliers(maxsit,maxsat)

c     satellite to debug
      integer iprndb

      logical          debug,first

      integer*4        itflag,iwkn,iyr,idoy,ihr,min,npnts
     .               , i,j,k
     .               , iyrstart,idoystart
     .               , ihrstart,minstart,iprns(maxsat),iprn,nprns
     .               , ioerr,nepoch,idumb,icall,inter
     .               , iepoch,isite,nsites,len

      real*8           sec,sow,utcoff,sodgps,sodm10
     .               , freql1,svdt,dumb

c*      integer ngoods(maxsat)
c*     real*8 allan(maxsat),deriv(maxepc),pi

c     message variable for routine report_stat
      character*256 message

c     L1 frequency in Hz
      DATA FREQL1/1.57542D9/

c       Initialization

c     the format for user questions:
    1 format (/,1x,a,1x)
c     first record
      first = .true.
c*      pi = 4.0d0*datan(1.0d0)
      do 10 i=1,maxsat
        iprns(i)=0
c*       ngoods(i) = 0
c*       allan(i) = 0.d0
   10 continue

c       Select the C-files
                 
      wildcard = 'c*.???'
      write (6,1)
     .'Choose one or more C-files from stations with atomic standards'
      call pickfns(wildcard,afile,nsites)

c       Open the input J-file

      print *,' '
      print *,' Choose as input reference the J-file used by MODEL'
     .       ,' to produce the C-files'
      wildcard =  'j*.???' 
      len = 6
      fclock = pickfn ( wildcard,len )
* MOD TAH 980904: set fclock to 1:len to get blanks at end
      fclock = fclock(1:len)
      if( fclock.eq.fsvclk ) then
         call report_stat('FATAL','MAKEJ','j_from_c',fclock,
     .   'Reference J-file name same as output J-file: ',0)
      endif
      open (unit=uclock,file=fclock,status='old',iostat=ioerr)
      if (ioerr .eq. 0) then
         call report_stat('STATUS','MAKEJ','j_from_c',fclock,
     .   'Opened Reference J-file: ',ioerr)
      else
         call report_stat('FATAL','MAKEJ','j_from_c',fclock,
     .   'Error opening Reference J-file: ',ioerr)
      endif

c       Loop over the C-files, reading in the observed-calculated L1 values

      do 20 isite = 1, nsites

        fxfile = afile(isite)
c       comment this since the C-file is opened in JREADC.  Perhaps reverse this
c       later and keep the open here, removing it from JREADC.  rwk 980316
c       call copens (fxfile,'old',uxfile,ioerr)
c       if (ioerr .eq. 0) then
c          call report_stat('STATUS','MAKEJ','j_from_c',fxfile,
c    .     'Opened C-file: ',ioerr)
c       else
c          call report_stat('FATAL','MAKEJ','j_from_c',fxfile,
c    .     'Error opening C-file: ',ioerr)
c       endif
        call jreadc ( uxfile, fxfile(1:16), nepoch, inter, nprns, iprns
     .              , jdc0  , sodc0       , pl1,       ttag
     .              , ierrf , clkepc, clkrat,       clkacc, isite )
c       jreadc changes the sign of the phase O-C so that it is now in
c       Doppler convention, consistent with the equations in Feigl [1991]
        if (dabs(clkepc) .lt. 1.d-9) then
           write (message,*) 'Warning: the C-file has receiver clock ',
     .    'corrections < 1.d-9 s. This is OK if its a MiniMac. '
           call report_stat('WARNING','MAKEJ','j_from_c',fxfile,
     .     message,0)
        endif
        call report_stat('STATUS','MAKEJ','j_from_c',fxfile,
     .  'Finished reading C-file: ',ioerr)
        close(uxfile)

c      set the nominal starting epoch to the closest even GPS second
c       --commented code is for C-files in UTC, no longer allowed --rwk 960517
       if( isite.eq.1 ) then
c*           utcoff = taiutc(jdc0) - 19.d0
c*           sodgps = sodc0 + utcoff
             sodgps = sodc0
c          add 5 seconds and then mod(10) to account for offsets w/in
c          5 seconds of even 10-s marks
           sodm10 = dmod(sodgps + 5.d0,10.d0)
           sodgps = sodgps - sodm10 + 5.d0
c*           sod0   = sodgps - utcoff
           sod0 = sodgps
           jd0 = jdc0
           if( debug ) print *,'Site 1  jdc0, sodc0, jd0, sod0 :'
     .                        , jdc0,sodc0,jd0,sod0
       endif

c      compute the time tags from the nominal  starting epoch
       if( jdc0.ne.jd0 ) then 
         write(message,'(a,i5)') 
     .      'Starting epoch does not match for isite: ',isite
         call report_stat('FATAL','MAKEJ','j_from_c',' ',message,0)
       else
         do i= 1,nepoch
           t(i,isite) = sodc0 - sod0 + ttag(i)
         enddo
       endif

 20   continue

c       Read the reference values for all satellites into storage

        rewind(uclock)
        icall = 0
        call readj ( uclock,iprns,nprns,idumb,idumb,dumb,dumb
     .             , icall,idumb,dumb,svcepc,svcrat,svcacc,valid)
        call report_stat('STATUS','MAKEJ','j_from_c',fclock,
     .  'Finished reading J-file: ',ioerr)

c       Setup to estimate a third-order polynpolynomial from the phase residuals

c     0th order we do not care about
c     1st order term is slope and thus frequency
c     2nd order term has units of 1/s, but we don't trust it.
c     number of data points to use
c      nn = 3  - set in data statement
c     number of parameters to estimate
c      mm = 3 - set in data statement

c       Loop over all the satellites

        do 200 iprn=1,nprns
         
          write(message,'(a,i3)') 'Estimating freqency for PRN '
     .                            ,iprns(iprn)
          call report_stat('STATUS','MAKEJ','j_from_c',' ',message,0)
          goodest(iprn) = 0

c        Loop over all the epochs

          do 100 iepoch = 1+nn/2, nepoch-nn/2

c           time of validity is middle epoch, but at nominal (even) value
            ti = dble(iepoch-1)*inter

c           Get the broadcast offset term from the stored input J-file values

               jd = jd0
               sod = sod0
               call timinc (jd,sod,ti)
               icall = 1
               call readj ( uclock, iprns, nprns, iprns(iprn)
     .                    , jd,     sod
     .                    , svdt,   icall,  jdj0,  sodj0
     .                    , svcepc, svcrat, svcacc,valid )
c              the J-file epoch may be different from that of the polynomial fit
c              tc = (jdj0,sodj0) - (jd0,sod0)
               tc = 86400.0d0*dble(jdj0-jd0) + sodj0-sod0

c           Loop over all C-files, estimating the coefficients from each

            ngoodf = 0
            do 50 isite = 1, nsites
              lgoodf(isite) = .false.

c           Compute the partial derivatives for this epoch

            j = 0
            do k=iepoch-nn/2,iepoch+nn/2
               j = j + 1
               pp(1,j) = 1.d0
               pp(2,j) = t(k,isite) - ti
               pp(3,j) = 2.0d0*(t(k,isite) - ti)**2
            enddo

c             estimate the polynomial coefficients from the C-file residuals
              call rate_est( iepoch, iprn, isite, pl1, ierrf, pp
     .                     , c, npnts, debug )
c             use only values with three good points
              if( npnts.ge.3 )  then
                  ngoodf = ngoodf + 1
                  lgoodf(isite) = .true.
              endif
 50         continue

c           Average the estimates from each of the C-files, discarding (but counting) outliers

            if ( ngoodf.gt.0 )
     .         call avgclk ( nsites, iprn, ngoodf, lgoodf, outliers, c
     .                     , c1 ,c2, debug )
            if ( ngoodf.gt.0 ) then
               goodest(iprn) = goodest(iprn) + 1
            else
               c1 = 0.d0
               c2 = 0.d0
               write(message,'(a,i3,a,i5,a)')'Estimate failed for PRN ',
     .         iprns(iprn), ' at epoch ',iepoch,
     .         '; using reference values'
               call report_stat('WARNING','MAKEJ','j_from_c',' ',
     .         message,0)
            endif

c           Add back the reference values

c           we use the clock correction from the broadcast message
c           as the zeroth order poly coefficient valid at epoch iepoch
            a0 = svcepc

c           Now we should have a good estimate for the first order term
c           The extra term is because the two polynomials are not referenced
c           to the same time.  Eqn (9):
            a1 = c1/freql1 + svcrat + 2.d0*svcacc*(ti-tc)
            a2 = c2/freql1 + svcacc

c           print some stuff
            if (debug .and. iprns(iprn) .eq. iprndb) then
              print *,' J_FROM_C : iepoch, jd, sod, iprn, ti, tc: '
     .               ,  iepoch, jd, sod, iprn, ti, tc
              print *,'  svcepc, svcrat, svcacc : ',svcepc,svcrat,svcacc
              print *,'  a0    , a1    , a2     : ',a0,a1,a2
              print *,' '
            endif

c           Write a record of the J-file

c            get the year and day of year
             call dayjul (jd,iyr,idoy)
c            convert UTC sod into hh:mm:ss
             call ds2hms (iyr,idoy,sod,ihr,min,sec)
c            the J-file format also requires GPS time (week number + s.o.w.)
             itflag= -2
             call timcon(itflag,iwkn,sow,iyr,idoy,ihr,min,sec,utcoff)
             write( usvclk,afmt ) iyr, idoy, ihr, min, sec, iwkn, sow
     .                          , iprns(iprn), a0, a1 ,a2

c              Save the first epoch
               if (first) then
                  iyrstart=iyr
                  idoystart=idoy
                  ihrstart=ihr
                  minstart=min
                  first=.false.
               endif

c--------end of loop on epochs
 100     continue

c        Compute the Allan variance of the estimates for information only

c*        if( lallanv )
      call allanv

c-----end loop on satellites
 200  continue

c     Print a summary

      write(message,410) nprns,iyrstart,idoystart,ihrstart,minstart
     .                , iyr,idoy,ihr,min
  410 format('J-File written for ',i2,' satellites',
     .' Start: ',i2,2x,i4,1x,2i3,' Stop : ',i2,2x,i4,1x,2i3)
      call report_stat('STATUS','MAKEJ','j_from_c',' ',message,0)

c     write(uscren,420)
c 420 format(//,1x,'PRN  N    sqrt(allan variance)')
c     do 430 i = 1,nprns
c        if (ngoods(i) .gt. 1) then
c           allan(i) = dsqrt(allan(i)/(2.d0*ngoods(i)))
c           write (uscren,'(1x,i2.2,1x,i4,1x,1pe16.4)')
c    .      iprns(i),ngoods(i),allan(i)
c        endif
c 430  continue

      write(uscren,440) (iprns(i),i=1,nprns)
* MOD TAH 200618: Updated 32I to 50I to allow for 35 Beidou satellites
  440 format(//,1x,'Valid estimates  PRN:',50i8)
      write(uscren,450) (goodest(i),i=1,nprns)
  450 format(/,22x,50i8)

      write(uscren,460) (iprns(i),i=1,nprns)
  460 format(//,1x,'Outliers   PRN:',50i8)
      write(uscren,'(/)')
      do i=1,nsites
        write(uscren,470) afile(i)(1:12),(outliers(i,j),j=1,nprns)
  470   format(a12,3x,50i8)
      enddo

      return
      end
