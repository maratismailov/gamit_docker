      subroutine editor (asites,nobs,nsat,ncfls,imode,llist)
c
c     Interactive editor for GAMIT C-files.
c
c     Kurt Feigl MIT June 1989
c
c     Modified by Shimon Wdowinski, April 1991
c
c     based on WAVED by Greg Beroza
c
c     Although written for the Apollo GPR routines, this routine is designed
c     to be machine independent.  As such, there are "wrappers" around
c     the machine-specific graphics calls.  These all have names starting
c     with the letter G and may be found in the file GSUBS.FTN.
c     They include:
c        GMOUSE     wait for input from mouse
c        GMVCUR     move the cursor
c        GSETUP     initialize the graphics package
c        GEND       termininate the graphics package
c        GMSG       write a message in a given part of the screen
c        GDOT       make a dot at a given point with the given radius
c        GLINE      draw a line between 2 points
c        GTSIZE     how big is the text?
c        GCLIP      set clipping area
c
c     It is menu driven with several levels of menus.
c     The variable icomm tells which menu box was picked.
c     The variable imenu tells you which level of menu is active.
c
c     The L1, L2, P1 and P2 time series (one-way, singly-differenced, or
c     doubly-differenced) are pulled from the common block arrays by
c     the time series.  These series are then combined into "combined
c     observables" (LG, LC, WL, etc.) by the subroutine combo.  This
c     has the advantage of only requiring a trip into the big 3-D data
c     arrays at the beginning and end of the work on a given series.
c     This way, the user can see the effect of his changes on all observables
c     before he decides to save them.
c
c
c     standards for portability:
c     ANSI 77 standard FORTRAN, except where noted by %%
c     no enddo

      implicit none
c
c     standard GAMIT array dimensions
      include '../includes/dimpar.h'

c     some non GAMIT declarations
      include '../includes/makex.h'

c     common block declarations
      include '../includes/cview.h'

c     standard GAMIT error flags
      include '../includes/errflg.h'

c     include machine dependent things:
      include '../includes/macdep.h'

c     FUNCTIONS:
c     functions to sort out error codes
c        returns .true. for observations good enough to be used in a solution
         logical lgood
c        returns .true. for marginal observations which might potentially useful
         logical lmarg
c        returns .true. for observations with a bias flag
c        logical lbias
c        returns .true. if a file exists
         logical fcheck
c        returns the nearest half-integer
c        real*8 rnhalf
c        returns the nearest quarter-integer
c        real*8 rnqrtr
c        returns .true. if a data series is available.
         logical getser
c        evaluates coefficients from maximum entropy estimate
         real*8 evlmem

c     pre-fit = 0 post-fit = 1
      integer imode 

c     new function to compute the abs value of an I*2
      integer*2 iabs2

c     list parameters
      integer    llist,ilist

c     4 character site codes
      character*4 asites(maxsit)

c     month codes
      character*3 amonth(12)

c     file name to dump hardcopy
      character*32 pfile

c     Code for observatable quantity (L1, etc.)
      character*2 obcode(-3:ntypes)
c     Code for units of the observable quantity (cycles, etc.)
      character*6 obunit(-3:ntypes)
c     list of commands (menu items)
      character*16 amenus(nmenus)
c     labels sat1,sat2,site1,site2,poly,obs
      character*12 alabel(6)
c     some file names
      character*96 amsg,title,shfile

c     these flags allow the loop to spin all the time TRUE TO:
c     redraw menus, etc
      logical        redraw
C     re-plot time series
      logical        replot
C     do it once
      logical        doonce
C     change the time series
      logical        newser
C     quit this routine
      logical        quitit
C     find next set of observables
      logical        seekob
C     blank out menu boxes
      logical        newbox
C     color screen
      logical        lcolor
C     time series to file
      logical        penplt
C     plot marginal points
      logical        domarg
C     plot hidden biases
      logical        hide_b
C     show the whole series
      logical        allpts
c     stack up plots
      logical        lstack
c     show movie
      logical        lmovie,stopmv
c     span the available X axis
      logical        lspan
c     have switched the chan/site index with find
      logical        lfound
c     flags occurance of a cycle-slip
      logical        slip,nextslip,lslip,simple

c     types of graphs
c     igraph = 1     time series
c     igraph = 2     sky map
c     igraph = 3     spectral estimate
c     igraph = 4     Allan standard deviation vs. period
      integer        igraph


c     can an observable be formed?
      logical        lcando(-3:ntypes)
c     Has clipping occurred?
      logical        lclipt
c     Useful
      logical        ltemp


c%%   values which must be integer*2 for the GSUBS routines
      integer*2 ix1,iy1,ix2,iy2,pos(2),ipos(2)
     .,         nx,ny,iradi,ibmx0,ibmy0
     .,         lbx,rbx,nxcurs

c     PLOTTING WINDOWS
c     In the convention used here, a window is defined by the elements of array w
c     w(1) x coordinate of u.l. corner
c     w(2) y coordinate of u.l. corner (down from u.l. corner of screen)
c     w(3) size of window in x dimension
c     w(4) size of window in y dimension
c     thus the corners are
c     w(1),w(2)          w(1)+w(3),w(2)
c     w(1),w(2)+w(4)     w(1)+w(3),w(2)+w(4)
C     clipping window
      integer*2         mclip(4)
C     time series only windows
      integer*2         mplot(maxwin,4)
C     station record window
      integer*2         mtext(4)
C     time scale window
      integer*2         mtaxi(4)
C     y-scale window
      integer*2         myaxi(maxwin,4)
C     mwarn window
      integer*2         mwarn(4)
C     mhelp window
      integer*2         mhelp(4)
C     working rectangle
      integer*2         mrect(4)
C     upper menu panel nmenus/2 + 1 - nmenus
      integer*2         menuu(4)
C     lower menu panel 1-nmenus/2
      integer*2         menul(4)
C     choices
      integer*2         menus(nmenus,4)

C     origin
      integer*2         ix0,iy0
C     origin for SkyMaps
      integer*2         ix0sm,iy0sm
C     x,y sizes of fonts
      integer*2         csize(2)
C     constant 2

c%%   this system status variable may need to be integer*2 on some systems
      integer*4 igerr

c     INDICES INTO THE TIME SERIES ARRAYS
C     plot from here to there
      integer         ii0,ii1,imarg
c     temporary limits
      integer         jj0,jj1
C     right and left bracket indices
      integer         iibl,iibr,iibl0
C     plot from here to there (buffered)
      integer         ii0old,ii1old
C     for extrapolation
      integer         iibeg,iiend
C     for slip search
      integer         istart,ilast
C     good obs
      integer         iig01,iig11,iig02,iig12
c     marginal observations
      integer         iim01,iim11,iim02,iim12

C     number of brackets on screen
      integer         nbrack
C     number of good points
      integer         ngood
C     number of points fit in PATCH
      integer         nmlfit
C     number of good or marginal points
      integer         nmarg
C     number of good points L1
      integer         ngood1
C     number of good points L2
      integer         ngood2
C     number of save requests (need 2)
      integer         nsave
C     number of bias-eliminationsave requests (need 2)
      integer         nelim
C     number of patch sequances in the multi patch process
      integer         npatch
C     number of slips in one sequance of multi patch process
      integer         nslip
C     number of time series to plot
      integer         nwins
C     pointer to window
      integer         iwin
C     number of the command
      integer         icomm
C     menu level
      integer         imenu
C     menu level (old buffer)
      integer         imenub
c     pointer to choice of SEEK
      integer         iseek
c     number of stacks so far
      integer         nstack
c     number of times a window has clipped
      integer         nclip(maxwin)

C     pointers to satellite channels (i is current, j is buffered, k is more buffered)
      integer         ichan1,ichan2,jchan1,jchan2,kchan1,kchan2
C     pointers to stations (i is current, j is buffered, k is more buffered)
      integer         isite1,isite2,jsite1,jsite2,ksite1,ksite2
C     number of observations in C-files
      integer         nobs
C     number of satellites
      integer         nsat
C     number of C-files (= # stations)
      integer         ncfls

c     weights for polynomial fit
      real*8  wgt(maxepc)

c     residuals from polynomial fit
      real*8  res(maxepc)

c     working vectors for polynomial fit
      real*8  xft(maxepc),yft(maxepc),fit(maxepc)

c     what we are plotting
      real*8  plt(maxepc)
      integer*4 kpl(maxepc)

c     problematic bias array and its counter
      integer       iprob,nprob,pbias(maxepc)

c     file unit number for plot file
      integer luout

c     DATA TYPES
c          -2 AZ azimuth from station to sat (CW from North)
c          -1 EL elevation angle in degrees
c           1 L1 phase
c           2 LC phase
c           3 LG phase
c           4 L2 phase
c           5 P1 pseudorange (cycles)
c           6 P2 pseudorange (cycles)
c           7 WL Widelane bias number
c           8 W* ModifiedWidelane bias number
c           9 PG P2 - gear*P1
c          10 BG LG - PG
c          11 I1 observable I1 = -g/(1-g*g) PG
c          12 I2 observable I2 = -1/(1-g*g) PG
c          13 M1 multipath delay in L1
c          14 M2 multipath delay in L2
c          15 MC multipath delay in LC
c          16 MG multipath delay in LG
c          17 CL Clock offset at station in microseconds
c          18 DT data - 4 1-WAY set of data
c     The following variables can take on any of the above values
C     what we are editing
      integer         l12e
C     what we are currently plotting
      integer         l12p
C     what is on the screen
      integer         l12v(maxwin)

C     plotting polynomial
      integer ipoly,jpoly


c     silly plotting variables
      real*8 dt
c            scales: X, Y, Radial and Character
     .,      xs,ys,rs,cs,ttt,resmin
     .,      ybeg,yend,tlen,yslp,yadd,ymult,yok0


c     statistics and extrema
      real*8 ytemp,ymin,ymax,ymid,yavg,yssd,yold,yjump,yjump0
      real*8 ytot,yflag0,yflag1,ystep
c     buffered values (for stacks)
      real*8 ymin2(maxwin),ymax2(maxwin),epsilon(4)
c     maximum entropy cofficients
      real*8 pm
      integer mcof
c     number of periods in Allan calculation
      integer nper

c     0 if screen dump is OK
      integer idump

      integer i,j,k,i0,it,nl,itemp,iok0,kfit,jchan,jsite,ntemp
      integer istep,ijump

      character*1 key,key0
      character*100 star100

      write (star100,2)
 2    format (100('*'))

      amonth(1) =  'JAN'
      amonth(2) =  'FEB'
      amonth(3) =  'MAR'
      amonth(4) =  'APR'
      amonth(5) =  'MAY'
      amonth(6) =  'JUN'
      amonth(7) =  'JUL'
      amonth(8) =  'AUG'
      amonth(9) =  'SEP'
      amonth(10)=  'OCT'
      amonth(11)=  'NOV'
      amonth(12)=  'DEC'

      obcode(-2)= 'AZ'
      obcode(-1)= 'EL'
      obcode(0) = '..'
      obcode(1) = 'L1'
      obcode(2) = 'LG'
      obcode(3) = 'LC'
      obcode(4) = 'L2'
      obcode(5) = 'P1'
      obcode(6) = 'P2'
      obcode(7) = 'WL'
      obcode(8) = 'W*'
      obcode(9) = 'PG'
      obcode(10)= 'BG'
      obcode(11)= 'I1'
      obcode(12)= 'I2'
      obcode(13)= 'M1'
      obcode(14)= 'M2'
      obcode(15)= 'MC'
      obcode(16)= 'MG'
      obcode(17)= 'CL'
      obcode(18)= 'DT'

      obunit(-2)= 'degree'
      obunit(-1)= 'degree'
      obunit(0) = '..'
      obunit(1) = 'cycles'
      obunit(2) = 'cycles'
      obunit(3) = 'cycles'
      obunit(4) = 'cycles'
      obunit(5) = 'cycles'
      obunit(6) = 'cycles'
      obunit(7) = 'cycles'
      obunit(8) = 'cycles'
      obunit(9) = 'cycles'
      obunit(10)= 'cycles'
      obunit(11)= 'cycles'
      obunit(12)= 'cycles'
      obunit(13)= 'cycles'
      obunit(14)= 'cycles'
      obunit(15)= 'cycles'
      obunit(16)= 'cycles'
      obunit(17)= 'micros'
      obunit(18)= 'points'

c     number of brackets
      nbrack = 0

c     width of cursor
      nxcurs = 16

c     a dumb constant

c     unit number for plot file
      luout = 20

c     initialize the graphics package
      call gsetup (ibmx0,ibmy0,nx,ny,lcolor,csize,ipos)

c     radius of dots for data points
      iradi = 2
      if (nx .gt. 1000) iradi = 2 * iradi

c     initialize plotting values
      yadd = 0.0d0
      iok0 = 0
      yok0 = 0.d0

c     set initial position for the list location
      ilist = 0

c%%?  delete shell script file for plotting commands
      shfile = 'cview.plots.cmd'
      if (fcheck(shfile)) then
         open (unit=1,file=shfile,status='unknown')
         close (unit=1,status='delete')
      endif

c     initial setting
      redraw = .true.
      replot = .false.
      quitit = .false.
      newbox = .true.
      newser = .false.
      seekob = .false.
      penplt = .false.
      domarg = .false.
      hide_b = .true.
      allpts = .false.
      lstack = .false.
      lmovie = .false.
      stopmv = .false.
      lspan  = .true.
      lfound = .false.
      slip   = .false.
      nextslip = .false.
      lslip = .false.
      nelim  = 0
      npatch = 0
      igraph = 1
      iseek  = 1
      imenu  = -1
      imenub = 0

      if (nsat .gt. 1) then
         ichan1 = 2
         ichan2 = 1
      else
         ichan1 = 1
         ichan2 = 0
      endif
      if (ncfls .gt. 1) then
         isite1 = 2
         isite2 = 1
      else
         isite1 = 1
         isite2 = 0
      endif

c     start with no polynomial fit
      ipoly   = 0

c     start with zero stack
      nstack = 0

c     initialize extrema
      do 35 i =1,maxwin
         nclip(i) = 0
         ymin2(i) =-1.0d8
         ymax2(i) = 1.0d8
 35   continue
      yflag0 = 0.d0
      yflag1 = 0.d0
      yjump0 = 0.d0
      yjump  = 0.d0
      yold   = 0.d0

c     start with L1,L2,LC,LG if this is a dual channel receiver
      if (iabs2(lambds(isite1,ichan1,2)) .ne. 0) then
         l12v(1) = 1
         l12v(2) = 4
         l12v(3) = 3
         l12v(4) = 2
         l12v(5) = 7
         nwins = 5
c        if we are only looking at one-ways, show us CL, too.
c         if (ncfls .eq. 1) then
c            l12v(5) = ntypes
c            nwins = 5
c         else
c            l12v(5) = 0
c            nwins = 4
c         endif
      else
         l12v(1) = 1
         l12v(2) = ntypes
         l12v(3) = 0
         l12v(4) = 0
         l12v(5) = 0
         nwins = 2
      endif

c     plot the whole series
      ii0 = 1
      ii1 = nobs

c     time interval
      dt = dble(inter)

c     **** START LOOP TO PLOT *****
      do 9000 while (.not. quitit)
      doonce = .true.

C     *** START LOOP TO FIND DATA ***
      do 1000 while (doonce .or. seekob)

c        *** look for next series with data
         if (seekob) then
c           buffer channel and site indexes
            ksite1 = isite1
            ksite2 = isite2
            kchan1 = ichan1
            kchan2 = ichan2
            call nextobserv(
     .           ichan1,ichan2,isite1,isite2,nsat,ncfls,
     .           iseek,llist,ilist,mwarn,seekob)
            iibl = 1
            nbrack = 0
            if (.not. seekob) then
               newser = .false.
               redraw = .true.
               newbox = .false.
               lmovie = .false.
               stopmv = .false.
               imenu = -1
               if (nextslip) then
                  nextslip = .false.
                  imenu = 0
               endif
               call gmsg (mwarn,'Last series!')
            endif
         endif


c        *** Get the observations from the big arrays
         if (newser) then
c           buffer channel and site indexes
            jsite1 = isite1
            jsite2 = isite2
            jchan1 = ichan1
            jchan2 = ichan2

c           assume all observables can be formed
            lcando(-2)= .true.
            lcando(-1)= .true.
            lcando(0) = .true.
            lcando(1) = .true.
            lcando(2) = .true.
            lcando(3) = .true.
            lcando(4) = .true.
            lcando(5) = .true.
            lcando(6) = .true.
            lcando(7) = .true.
            lcando(8) = .true.
            lcando(9) = .true.
            lcando(10)= .true.
            lcando(11)= .true.
            lcando(12)= .true.
            lcando(13)= .false.
            lcando(14)= .false.
            lcando(15)= .false.
            lcando(16)= .false.
            lcando(17)= .true.
            lcando(18)= .true.

c           get pseudoranges
c           note that kww will get overwritten
            it = 5
            lcando(it) = getser (ichan1,ichan2,isite1,isite2,it,imode,
     .              domarg,pc1,kww,nobs,ngood1,iig01,iig11,iim01,iim11)

            it = 6
            lcando(it) = getser (ichan1,ichan2,isite1,isite2,it,imode,
     .              domarg,pc2,kww,nobs,ngood1,iig02,iig12,iim02,iim12)

c           get azimuth
            it = -2
            lcando(it) = getser (ichan1,ichan2,isite1,isite2,it,imode,
     .              domarg,azw,kww,nobs,ngood1,iig01,iig11,iim01,iim11)

c           get elevation angle
            it = -1
            lcando(it) = getser (ichan1,ichan2,isite1,isite2,it,imode,
     .              domarg,elw,kww,nobs,ngood1,iig01,iig11,iim01,iim11)

c           get station clocks
            it = 16
            ltemp = getser (ichan1,ichan2,isite1,isite2,it,imode,
     .              domarg,clw,kww,nobs,ngood1,iig01,iig11,iim01,iim11)

c           get multipath delay for L1
c           it = 12
c           lcando(it) = getser (ichan1,ichan2,isite1,isite2,it,imode,
c    .              domarg,rw1,kcc,nobs,ngood1,iig01,iig11,iim01,iim11)

c           get multipath delay for L2
c           it = 13
c           lcando(it) = getser (ichan1,ichan2,isite1,isite2,it,imode,
c    .              domarg,rw2,kcc,nobs,ngood1,iig01,iig11,iim01,iim11)

c           get L1, which is always needed
            it = 1
            lcando(it) = getser (ichan1,ichan2,isite1,isite2,it,imode,
     .              domarg,cl1,kcc,nobs,ngood1,iig01,iig11,iim01,iim11)
c           load working vectors
            call vcopy (kcc,cl1,kww,wl1,nobs)

c           get L2 if it is there
            it = 4
            lcando(it) = getser (ichan1,ichan2,isite1,isite2,it,imode,
     .              domarg,cl2,kcc,nobs,ngood2,iig02,iig12,iim02,iim12)
            call vcopy (kcc,cl2,kww,wl2,nobs)
                
c RWK 151230: Save in common the index for the first satellite to use for
c           frequency quantities, now in satellite-indexed arrays.  The
c           frequency quantities are the same for all satellites for most 
c           GNSS, but for Glonass FDMA, they will be different.  This kluge 
c           of using the frequency for the first satellite will allow plotting
c           of one-way and between-site differences correctly.  Between-satellite
c           and double differences cannot be correctly formed.
            jsat1 = ichan1     

            ngood = ngood1 + ngood2

c           continue only if we have data
            if (ngood .gt. 0) then
               seekob = .false.
               if (.not. lmovie) then
                  newser = .false.
               endif
               replot = .true.
               redraw = .true.
               newbox = .true.
               if (.not.penplt .and.
     .             .not.lstack .and.
     .             .not.lmovie) then
                  imenu = 0
               endif
            else
               if (seekob) then
                  call gmsg (mwarn,'Looking for data.')
               else
                  call gmsg (mwarn,'No data here.')
                  newser = .false.
                  imenu = -1
               endif
               redraw = .true.
               replot = .false.
               newbox = .false.
            endif
         endif

c        *** REDRAW MENU ***
         if (redraw) then
c           number of save requests
            nsave  = 0

c           number of brackets
            nbrack = 0

c           count the number of plots to make
            nwins = 0
            do 1005 i = 1, maxwin
               if (l12v(i) .ne. 0) nwins = nwins + 1
 1005       continue

c           set up window sizes
            call editor_windos
     .      (mtext,mwarn,mhelp,menuu,menul,menus,mclip,
     .      myaxi,mplot,mtaxi,
     .      ibmx0,ibmy0,nx,ny,nxcurs,csize,nwins)

c           set up labels for buttons
            call editor_labels
     .      (ichan1,ichan2,isite1,isite2,ipoly,
     .      nsat,ncfls,asites,alabel)

c           set up menu choices
            call editor_menus
     .      (amenus,menus,menul,menuu,mwarn,mhelp,
     .      alabel,imenu,imenub,iseek,obcode,
     .      l12v,ipoly,yadd,amonth,pfile,title,nextslip,igraph,
     .      newbox,domarg,hide_b,lstack,lmovie,lspan,allpts,csize)

            redraw = .false.
         endif
c        *** END LOOP TO FIND DATA ***
         doonce = .false.
 1000 continue

c     *** PLOT DATA ***
      if (replot) then

c        choose data span
         ii0old = ii0
         ii1old = ii1
         if (allpts) then
            ii0 = 1
            ii1 = nobs
         else if (lspan) then
            if (domarg) then
               ii0 = min(iim01,iim02)
               ii1 = max(iim11,iim12)
            else
               ii0 = min(iig01,iig02)
               ii1 = max(iig11,iig12)
            endif
            imarg = (ii1-ii0)/20
            ii0 = max(ii0-imarg,1)
            ii1 = min(ii1+imarg,nobs)
         endif

c         nl = ii1 - ii0
c         if (nl .lt. 5) then
c            call gmsg(mwarn,'Too few points.')
c            ii1 = ii1 + 2
c            ii0 = ii0 - 2
c         endif
c        don't go out of bounds
         if (ii1 .gt. nobs) then
            call gmsg(mwarn,'Reached end of series.')
            ii1 = nobs
         endif
         if (ii1 .lt. 1) then
            call gmsg(mwarn,'At start of series.')
            ii1 = mplot(1,4)/(2*iradi+1)
         endif
         if (ii0 .lt. 1) then
            call gmsg(mwarn,'At start of series.')
            ii0 = 1
         endif
         if (ii0 .gt. nobs) then
            call gmsg(mwarn,'Reached end of series.')
            ii0 = nobs
         endif

c        clean out the mclip window
C        call gmsg (mclip,' ')

c        set the clipping window
c        call gclip (mclip,.true.)

c        title line
         call gmsg (mtext,title)

c        *** LOOP OVER PLOTTING WINDOWS ***
         do 3000 j = 1,nwins
c           what type of combination observable to form?
            l12p = l12v(j)

c           load mrect with this window
            do 3003 i = 1,4
               mrect(i) = mplot(j,i)
 3003       continue

c           set clipping window
            call gclip(mrect,.true.)

c           clean out the plotting window if not stacking things
            if ((.not. lstack) .or. (nstack .eq. 0)) then
               call gmsg(mrect,'   ')
            endif

c           load plt vector with combined observable or elevation
c           form the observable, get extrema, and count good points
            if(l12p .ne. 18) call combo (nobs,ii0,ii1,l12p,domarg,
     .          kpl,plt,ymin,ymax,yavg,yssd,yjump,ystep,
     .          ngood,nmarg,ijump,istep,lcando)

c           **** set a flag to stop movie ***
            if (l12p .eq. 3) then
               yflag0 = yflag1
               yflag1 = yssd
               yjump0 = yjump
            endif

c           DERIVATIVE
c           Do not perform on elevation or azimuth
            if ((ipoly.eq. -1 .or.  ipoly .ge. 10) .and.
     .          (l12p .ne. -1 .and. l12p  .ne. -2)) then
c              find first good point
               do 3012 i=ii0,ii1
                  if (lgood(kpl(i)).or.(domarg.and.lmarg(kpl(i)))) then
                     iok0 = i
                     yok0 = plt(i)
                     kpl(i) = ignone
                     goto 3013
                  endif
 3012          continue
 3013          continue

               ymax = -1.0d20
               ymin =  1.0d20
               do 3015 i=iok0+1,ii1
                  if (lgood(kpl(i)).or.(domarg.and.lmarg(kpl(i)))) then
                     ttt = tx(i) - tx(iok0)
                     if (dabs(ttt) .gt. 1e-10) then
                        ytemp = (plt(i)-yok0)/ttt
                     else
                        ytemp = 0.0d0
                     endif
c                    buffer the most recent value
                     yok0 = plt(i)
                     plt(i) = ytemp
                     ymax = dmax1(ymax,plt(i))
                     ymin = dmin1(ymin,plt(i))
                     iok0 = i
                  endif
 3015          continue
            endif


c           POLYNOMIAL FIT
c           do this before smoothing
c           Do not perform polynomial fit on elevation or azimuth.
            if ((ipoly .eq. -2 .or. ipoly .gt. 0)
     .           .and. l12p .ne. -1 .and. l12p .ne. -2) then
c              load working arrays
               k = 0
               do 3005 i=ii0,ii1
                  if (lgood(kpl(i)).or.(domarg.and.lmarg(kpl(i)))) then
                     k = k + 1
                     xft(k) = tx(i) - tx(ii0)
                     yft(k) = plt(i)
c                    uniform weighting
                     wgt(k) = 1.0d0
                  endif
 3005          continue
               ngood = k

c              perform the polynomial fit
               kfit = 2
               if (ipoly .eq. -2) then
c                 do second order fit prior to SMOOTH
                  jpoly = 2
               else
                  jpoly = mod(ipoly,10)
               endif
               call polyft
     .            (ngood,yft,xft,0.0d0,wgt,jpoly,fit,res,yssd,kfit)
               if (kfit .ge. 0) then
c                 unload the arrays, refigure extrema, too
                  k = 0
                  ymax = -1.0d20
                  ymin =  1.0d20
                  do 3007 i=ii0,ii1
                     if (lgood(kpl(i)).or.
     .                  (domarg.and.lmarg(kpl(i)))) then
                        k = k+1
                        plt(i) = res(k)
                        ymax = dmax1(ymax,plt(i))
                        ymin = dmin1(ymin,plt(i))
                     else
                        plt(i) = 0.d0
                     endif
 3007             continue
               else
                  call gmsg (mwarn,'Poly fit failed.')
                  call gmsg (mrect,
     .            'Polynomial fit failed. Click ALL or SPAN')
                  replot = .false.
               endif
            endif


c           SMOOTH
c           Do not perform on elevation or azimuth.
            if (ipoly.eq.-2 .and. l12p.ne.-1 .and. l12p.ne.-2) then
c              load working arrays
               k = 0
               ngood = 0
               do 3105 i=ii0,ii1
                  k = k + 1
                  if (lgood(kpl(i)).or.(domarg.and.lmarg(kpl(i)))) then
                     yft(k) = plt(i)
                     ngood = ngood+1
                  else
                     yft(k) = 0.d0
                  endif
 3105          continue

               if (ngood .gt. 50) then
c                 perform the smoothing over a window of about 1/6 number of points
                  ttt = ngood/6.d0
                  ttt = dmax1(ttt,10.d0)
                  ttt = dmin1(ttt,50.d0)
                  call smooft(yft,k,ttt)
c                 unload the arrays, refigure extrema, too
                  ymax = -1.0d20
                  ymin =  1.0d20
                  k = 0
                  do 3107 i=ii0,ii1
                     k = k+1
                     if (lgood(kpl(i)).or.
     .                  (domarg.and.lmarg(kpl(i)))) then
                        plt(i) = plt(i) - yft(k)
                        ymax = dmax1(ymax,plt(i))
                        ymin = dmin1(ymin,plt(i))
                     endif
 3107             continue
               else
                  call gmsg (mwarn,'SMOOTH failed.')
                  call gmsg (mrect,
     .            'Not enough points to smooth. Click SPAN or ALL.')
                  replot = .false.
               endif
            endif


c           SPECTRUM
c           estimate the power spectrum using maximum entropy
            if (igraph .eq. 3) then
c              if data are not good, then set to zero
               ngood = 0
               k = 0
               do 3205 i=ii0,ii1
                  k = k + 1
                  if (lgood(kpl(i)).or.(domarg.and.lmarg(kpl(i)))) then
                     ngood = ngood + 1
                     plt(k) = plt(i)
                  else
                     plt(k) = 0.d0
                  endif
 3205          continue
               mcof  =  20
               if (ngood .gt. 4 * mcof) then
c                 use the polynomial arrays for work space
c                 ltemp returned .true. if OK
                  call memcof (plt,k,mcof,pm,xft,yft,wgt,res,ltemp)
                  if (ltemp) then
c                    use tx for frequency, plt for power
                     ymin = +1.d20
                     ymax = -1.d20
                     do 3207 i = 1,maxepc
c                       go from 0 to the Nyquist frequency.
                        tx(i) =  0.5d0 * dble(i)/dble(maxepc)
                        plt(i) = evlmem(tx(i),xft,mcof,pm)
c                       take log
                        plt(i) = dlog10(plt(i))
                        ymin = dmin1(plt(i),ymin)
                        ymax = dmax1(plt(i),ymax)
c                       convert from cycles/epoch to cycles/sec
                        tx(i) = tx(i)/dt
 3207                continue
                  else
                     call gmsg (mrect, 'SPECTRUM failed')
                     replot = .false.
                  endif
               else
                  call gmsg (mrect,
     .            'Not enough points for SPECTRUM. Click ALL')
                  replot = .false.
               endif
            endif

c           ALLAN
c           estimate the Allan standard deviation at several periods
            if (igraph .eq. 4) then
               if (l12p.eq.1.or.l12p.eq.2.or.l12p.eq.3.or.l12p.eq.4)then
c                 if data are not good, then set to obscene value (-9.d99)
                  k = 0
                  do 3305 i=ii0,ii1
                     k = k + 1
                     xft(k) = tx(i)
                     if(lgood(kpl(i)).or.(domarg.and.lmarg(kpl(i))))then
                        ngood = ngood + 1
                        plt(k) = plt(i)
                     else
                        plt(k) = -9.d99
                     endif
 3305             continue

                  if (ngood .gt. 0) then
c                    use the polynomial arrays for work space
c                    nper returned 0 on error
c                    Select nominal frequency
crwk 151229: It appears that the frequency use here is crude, so for Glonass, don't
c            worry about the diffences among SVs
                     if (l12p .eq. 4) then
                        ytemp = fL2(1)
                     else
                        ytemp = fL1(jsat1)
                     endif
                     call allan (ytemp,plt,xft,yft,k,nper)

                     if (nper .gt. 0) then
c                       use tx for frequency, plt for power
                        ymin = +1.d20
                        ymax = -1.d20
                        do 3307 i = 1,nper
c                          X-axis is sample period
                           tx(i) =  i * dt
                           plt(i) = yft(i)
c                          take log
                           plt(i) = dlog10(plt(i))
                           ymin = dmin1(plt(i),ymin)
                           ymax = dmax1(plt(i),ymax)
 3307                   continue
                     else
                        call gmsg (mrect, 'ALLAN failed')
                        replot = .false.
                     endif
                  else
                     call gmsg (mrect,
     .               'Not enough points for ALLAN. Click ALL')
                     replot = .false.
                  endif
               else
                  call gmsg (mrect,
     .            'ALLAN not defined for this frequency')
                   replot = .false.
               endif
            endif

c           Announce it if the observable cannot be formed.
            if (.not. lcando(l12p)) then
               amsg = 'Cannot form '//obcode(l12p)//'.'
               if (lstack) then
                  call gmsg (mwarn,amsg)
               else
                  do 3019 i = 1,4
                     mrect(i) = mplot(j,i)
 3019             continue
                  call gmsg (mrect,amsg)
               endif
            endif

c           set a minimum resolution, which depends on the observable.
            if (igraph .eq. 3 .or. igraph .eq. 4) then
c              no minimum for spectrum
               resmin = 0.0d0
            else if (l12p .eq. 16) then
c              .01 microseconds for CL
               resmin = 0.01d0
            else if (lstack) then
               resmin = 1.0d0
            else
c              .25 cycles for everybody else
ckf910630      resmin = 0.25d0
c              .50 cycles is the popular favorite.
C MOD TAH 931011: Change resmin to smaller value if yscale is small
               if( ymax-ymin.lt.0.2d0 ) then
                   resmin = 0.01d0
               else
                   resmin = 0.50d0
               end if
            endif

c           resolution should be tighter for derivatives
            if (ipoly .eq. -1 .or. ipoly .ge. 10) resmin = resmin/dt

            if (dabs(ymax - ymin) .lt. resmin) then
               ymax = ymax + resmin/2.0d0
               ymin = ymin - resmin/2.0d0
            endif

c           set X scale
            if (igraph .eq. 3) then
c              all the way to Nyquist frequency for spectrum
               tlen = tx(maxepc) - tx(1)
            else if (igraph .eq. 4) then
c              only to 1/8 of total possible for Allan
               tlen = tx(maxepc/8) - tx(1)
            else
               tlen = tx(ii1) - tx(ii0)
            endif

            if (dabs(tlen) .gt. 1.0d-10 .and. nmarg .gt. 0) then
               xs = (mplot(j,3)-3*iradi)/tlen
c              buffer the limiting indexes
               ii0old = ii0
               ii1old = ii1
            else
               call gmsg (mwarn,'No data in range.')
               lcando(l12p) = .false.
            endif


c           BEGIN PLOTTING
            if (replot .and. lcando(l12p)) then
c              x origin
               ix0 = mplot(j,1) + iradi + 1

c              initial y scale and origin
c              select the midpoint
c              remove midpoint for everything
c              but elevation angle and azimuth and power spectrurm
c              plotting elevation or azimuth on a skymap should
c              just give a trace, without changing symbol size
               if (l12p .eq. 18) then
c                 data plot
                  ymax =  1.0d0
                  ymin = -1.0d0
                  ymid = 0.0d0
                  ys = mplot(j,4)/(1.1d0*(ymax-ymin))
                  cs = 40.0d0/nwins/(ymax-ymin)
                  iy0 = mplot(j,2) + mplot(j,4)/2
               else if (l12p .eq. -2) then
c                 azimuth plot
                  ys   = mplot(j,4)/(740.0d0)
                  ymax =  370.0d0
                  ymin = -370.0d0
                  ymid = 0.0d0
                  iy0 = mplot(j,2) + mplot(j,4)
                  if (igraph .eq. 2) cs = 0.0d0
               else if (l12p .eq. -1) then
c                 elevation angle plot
c                 can only absolutly plot oneway elevation angles.
c                 will plot elev angle differences for single and double diffs!!!!!
                  ymid =  0.0d0
                    if( ymax .lt. 0 .and. ymin .lt. 0 ) then
                      iy0 = mplot(j,2)
                      ys   = mplot(j,4)/(1.10d0*dabs(ymin))
                    elseif( ymax .gt. 0 .and. ymin .lt. 0 ) then
                      iy0 = mplot(j,2) + mplot(j,4)/((dabs(ymin)+ymax)
     .                      /ymax)
                      ys   = mplot(j,4)/(1.10d0*(dabs(ymin)+ymax))
                    elseif( ymax .gt. 0 .and. ymin .gt. 0 ) then
                      iy0 = mplot(j,2) + mplot(j,4)
                      ys   = mplot(j,4)/(1.10d0*dabs(ymax))
                    endif
cold              ys   = mplot(j,4)/(200.0d0)
cold              ymax =  100.0d0
cold              ymin = -100.0d0
cold              ymid = 0.0d0
cold              iy0 = mplot(j,2) + mplot(j,4)
                  if (igraph .eq. 2) cs = 0.0d0
               else
c                 any other kind of plot
                  ys = mplot(j,4)/(1.1d0*(ymax-ymin))
                  cs = 40.0d0/nwins/(ymax-ymin)
                  ymid = (ymin+ymax)/2.0d0
                  iy0 = mplot(j,2) + mplot(j,4)/2
               endif

               if (igraph .eq. 2) then
c                 set origin and scales for skymap
                  rs = mplot(j,4)/2.1d0
                  ix0sm = mplot(j,1) + mplot(j,3)/2
                  iy0sm = mplot(j,2) + mplot(j,4)/2
               endif

c              Keep the scale the same for stacking plots.
               if (nstack .eq. 0) then
                  nclip(j) = 0
                  ymin2(j) = ymin
                  ymax2(j) = ymax
               endif

c              Tell us if clipping occurs.
               if (lstack) then
                  iradi = 1
                  if (l12p .gt. 0) then
                     ytemp = dabs(ymax2(j)-ymin2(j))
                     if (dabs(ymax-ymin) .gt. ytemp) then
                        nclip(j) = nclip(j) + 1
                        mrect(3) = (nclip(j) + 1) * csize(1)
                        mrect(4) = 1.5 * csize(2)
                        mrect(1) = mplot(j,1)+mplot(j,3)-mrect(3)
                        mrect(2) = mplot(j,2)
                        if (nclip(j) .gt. 99) then
                           call gmsg(mwarn,'Over 100 clips!')
                           nclip(j) = 101
                        else
                           call gmsg (mrect,star100(1:nclip(j)))
                        endif
                        lclipt = .true.
                     else
                        lclipt = .false.
                     endif
                     ymin =  ymin2(j)
                     ymax =  ymax2(j)
                     ys = mplot(j,4)/(1.1d0*(ymax-ymin))
                  endif
               else
                  lclipt = .false.
                  iradi = 2
               endif

c              Set clipping window
               do 1040 i = 1,4
                  mrect(i) = mplot(j,i)
 1040          continue
               call gclip(mrect,.true.)

c              can save a little bit of time here by skipping
c              over the empty epochs
               if (allpts) then
                  call bound (ichan1,ichan2,isite1,isite2,jj0,jj1)
               else
                  jj0 = ii0
                  jj1 = ii1
               endif

               if (igraph .eq. 1) then
c                 plot time series
                  if (l12p .ne. 18) then
                    call editor_plotts (ii0,ii1,jj0,jj1,ix0,iy0,xs,ys,
     .              ymid,ymin,ymax,iradi,plt,kpl,domarg,hide_b,lclipt,
     .              isite1,isite2)
                  else
                    call data_yaxis
     .                (j,myaxi,ymin,ymax,yssd,ymid,ix0,iy0,xs,ys,
     .                csize,obcode(l12v(j)),obunit(l12v(j)),ipoly,
     .                lstack,igraph,alabel)
                    call editor_plotdt (ii0,ii1,jj0,jj1,ix0,iy0,xs,ys,
     .                ymid,ymin,ymax,iradi,plt,kpl,domarg,hide_b,lclipt,
     .                isite1,isite2,ichan1,ichan2)
                  endif
               else if (igraph .eq. 2) then
c                 plot sky map
                  call editor_plotsm
     .            (jj0,jj1,ix0sm,iy0sm,rs,cs,ymid,dt,plt,kpl,domarg)
               else if (igraph .eq. 3) then
c                 plot power spectrum
                  call editor_plotps (1,maxepc,ix0,iy0,xs,ys,ymid,
     .            ymin,ymax,iradi,plt,kpl,domarg,lclipt)
               else if (igraph .eq. 4) then
c                 plot Allan variance
                  call editor_plotal
     .            (1,nper,ix0,iy0,xs,ys,ymin,ymid,ymax,plt)
               else
                  call report_stat('STATUS','CVIEW','editor',' ',
     .            'Error, bogus igraph',0)
               endif

c              dump time series (with midpoint removed) to pfile
               if (penplt) then
c                 observable type to file name
                  pfile(27:29) = '.'//obcode(l12p)

                  call editor_penplt
     .            (luout,pfile,shfile,title,obcode(l12p),ii0,ii1,
     .             kpl,plt,ymin,ymax,ymid,yssd,
     .             j,nwins,mwarn,penplt,igraph)
               endif

c              draw a Y-axis and label it with observable type and
c              the standard deviation
               if ((nstack .eq. 0) .and. (l12p .ne. 18)) then
* MOD TAH 031223: Added passing of yavg into editor_yaxis so that it
*                 can be output on the window (add before yssd argument).
                  call editor_yaxis
     .            (j,myaxi,ymin,ymax,yavg,yssd,ymid,ix0,iy0,xs,ys,
     .            csize,obcode(l12v(j)),obunit(l12v(j)),ipoly,
     .            lstack,igraph)
               endif

c              line to separate the plots
               call gclip(mrect,.false.)
               ix1 = myaxi(j,1)
               ix2 = mplot(j,1) + mplot(j,3)
               iy1 = mplot(j,2)
               iy2 = mplot(j,2)
               call gline (ix1,iy1,ix2,iy2,igerr)
               if (j .eq. nwins) then
                  iy1 = mplot(j,2) + mplot(j,4)
                  iy2 = mplot(j,2) + mplot(j,4)
               endif

c              seem to need this
               lbx = mplot(j,1)
               rbx = mplot(j,1)+mplot(j,3)

            endif
c           end of loop over multiple time series plots
 3000    continue

c        draw a horizontal X-axis (time)
         call editor_xaxis(mtaxi,ii0,ii1,ix0,xs,csize,igraph)

c        increment stack counter
         if (lstack) nstack = nstack + 1

c        turn off clipping
         call gclip (mclip,.false.)

c        *** END OF PLOTTING ***
         replot = .false.
      endif

c**********************************************************************
c**********************************************************************
c
c     Now we are done plotting, so wait for input from mouse,
c     by using the menus at file ED_MNU0.FTN.
c
c**********************************************************************
c**********************************************************************

      call ed_menus (
     .         lmovie,stopmv,lstack,lspan,lfound,lcolor,quitit,allpts,
     .         seekob,newser,replot,redraw,penplt,newbox,domarg,hide_b,
     .         nbrack,nobs,nsave,nstack,ncfls,nsat,nmlfit,ngood,ntemp,
     .         npatch,nwins,iiend,iibeg,nelim,
     .         isite1,isite2,ichan1,ichan2,
     .         jsite1,jsite2,jchan1,jchan2,
     .         ksite1,ksite2,kchan1,kchan2,
     .         jsite,jchan,iwin,nl,luout,
     .         key,key0,pos,icomm,imenu,imenub,idump,ipoly,igraph,itemp,
     .         ymin2,ymax2,yadd,yend,ybeg,yslp,ymult,yold,yjump0,slip,
     .         nextslip,llist,ilist,iseek,
     .         ibmx0,ibmy0,nx,ny,l12v,l12e,l12p,ii0,ii1,iibl,iibr,iibl0,
     .         mplot,mrect,mtaxi,menus,mclip,csize,mwarn,mhelp,amsg,
     .         igerr,amenus,lbx,rbx,ix1,ix2,iy1,iy2,ix0,i0,dt,xs,ttt,
     .         plt,lslip,nslip,istart,ilast,epsilon,ytot,yflag0,yflag1,
     .         iprob,nprob,pbias,simple,imode,lcando)


c     **************************************************************
c     **** bottom of while loop
c     **************************************************************
 9000 continue

      call gmsg (mwarn,'Bye-bye!')

c     put the cursor back to avoid confusing user
      call gmvcur (ipos, igerr)

c     terminate the gracphics package
      call gend (igerr)

      return
      end
c


c**********************************************************************

