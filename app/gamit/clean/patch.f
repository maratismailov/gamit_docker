
CTITLE PATCH

      subroutine patch(iibl,iibr,ii0,ii1,mwarn,site1,chan1,
     .                idd,nobs,replot,nmfitl,yadd,confid,doit)

*     This routine will try to patch cycle slips at both L1 and L2 carrier
*     frequencies.  The point which is fixed is determined by the left most
*     brace on the current plot (point iibl), and the cycle slip will be
*     fixed to point iibr (if there is no second brace then to the end of
*     data is fixed.  The algorithm used projects linear trends from the
*     TnumfitU points before the slip and the numfit points after the slip,
*     at both L1 and L2.  The estimate is the number of cylce slips is made
*     by minimizing the discontinity in Lc and Lg at the midpoint of the
*     last and fisrt points in the two projected lines.  The relative
*     weighting of the Lc ang Lg discontinuity is based on the projection of
*     the rms scatters in L1 and L2 projected into the Lc and Lg combinations.

* PARAMETERS
*
*       mxfitl  - Maximum number of data to use in fiting the linear
*               - trends.  Currently this is also the number of values
*               - to use. MUST BE GREATER THAN OR EQUAL TO 3.

      integer mxfitl

*                                 ! This is a guess as to RbestS value
      parameter ( mxfitl = 15 )
*
* INCLUDE FILES
*                                         ! Not sure if I need this
*                                         ! GAMIT array dimensions
      include '../includes/dimpar.h'
      include '../includes/makex.h'
*                                         ! CVIEW common
      include '../includes/cview.h'
*                                         ! SOLVE error flags (may not be needed)
      include '../includes/errflg.h'

* VARIABLES FROM COMMON.
*     WL1(maxepc)     ! LI carrier phase in cylces (modified in this
*                     ! routine.
*     WL2(maxepc)     ! L2 carrier phase in cylces (modified here)
*     kww(maxepc)     ! Lookup entries for data quality (copied)
*     lambds(site,chan,1/2)   ! Wavelength factors for site, satellite and
*                         ! L1 or L2.
*     BL1(maxepc)     ! copy of WL1 so we can undo
*     BL2(maxepc)     ! Copy of WL2 so we can undo.
*     kbb(maxepc)     ! copy of the point to quality flag look up
*                     ! array indices. (used for undo)


* PASSED Variables.
*
*       iibl, iibr  - Point numbers of left and right braces on
*                   - plot (if no right brace then iibr=number of
*                   - observations.
*       ii0, ii1        - Start and stop points on current plot.  These
*                   - are not used in current code but may be used
*                   - in the future to define the linear fit range.
*       site1       - site number (used to get wavelength factor)
*       chan1       - Channel number (for wavelength factor)
*       nobs        - Total number of observations (used for coping
*                   - arrays before RfixingS cycle slips.
*       nmfitl      - set the number of points to use in the linear fits
*       doit        - true to perform the correction,
*                   - false to just return the amount in yadd
*                   - added to call list by kurt 910619
*       idd         - double difference indicator

      integer iibl, iibr, ii0, ii1, site1, chan1, idd, nobs

*         mwarn(4)  - Window for warning messages (used by graphics
*                   - routine)

      integer*2 mwarn(4)

*       Replot      - Indicates we should replot.  Only done
*                   - if patch was successful in finding data.

      logical Replot,doit

* LOCAL VARIABLES

*       nmfitl      - Number of values to use in linear fits.
*                   - The full number must be availabel or else
*                   - this routine will exit.
*   indfl(mxfitl)   - Indexes to obervations to be used in the
*                   - data left of the brace.  Assumed to the same
*                   - for L1 and L2. LGOOD is used to decide on
*                   - data to be used.
*   indfr(mxfitl)   - Indices to observations right of the left
*                   - brace.
*   errl,errr       - Error returns from getting the data to the
*                   - left and right of left brace.  Currently has
*                   - the following meanings
*                   -    0 - all OK, enough data found
*                   - -  1 - Not enough data found.  Do not proceed.
*                   - -999 - Invalid direction passed to scan routine.

      integer nmfitl, indfl(mxfitl), indfr(mxfitl), errl,errr

*       ambi(2)     - Ambiguities in cycles for L1 and L2.
*       estl(2)     - Estimate of the L1 and L2 carry phase at the
*                   - midpoint, obtained from data left of the brace
*       estr(2)     - Estimates of the L1 and L3 carry phase at the
*                   - midpoint, obtained from data right of the brace
*       ratel(2)        - Rate of change of carrier phase (by index number)
*                   - for L1 and L2 in the left segment
*       rater(2)        - Rate of change of carrier phase (by index number)
*                   - for L1 and L2 in the right segment
*       dest(2)     - Difference in estimates of L1 and L2 at the
*                   - midpoint.

*       midpt       - Midpoint fracttional epoch number.  Currently
*                   - just the mean of the first points from the left
*                   - and right series.  Could be updated for position
*                   - of minimum variance in the estimate in the
*                   - left and right difference.
*       dcg(2)      - Difference in the estimates of LC and LG from
*                   - from the left and right series.
*       varcg(2)        - variances of the dcg differences.
*       Ltocg(2,2)  - Transformation matrix from L1 and L2 to
*                   - Lc and Lg.  (Here L1 = Lc + Lg, and
*                   - L2 = gLc + Lg/g (g is ratio of L1 freq/L2 freq)
*       confid      - Confidence of slip repair in %.  Based on
*                   - (1-best chi**2/next best chi**2)*100.

*       yadd(2)     - Number of cycles to be added to L! and L2.

      real*8 ambi(2), estl(2), estr(2), ratel(2), rater(2), dest(2),
     .    midpt, dcg(2), varcg(2), Ltocg(2,2), confid, yadd(2)

*            mess       - Message string used with GMSG to send message
*                   - to user about patch.

      character*40 mess

      integer*4 len,rcpar

      character*80 prog_name
 
c     get calling program name to know which version of gmsg to call
      len = rcpar(0,prog_name)

      if (nmfitl .le. 0 .or. nmfitl .gt. mxfitl) nmfitl = mxfitl

*     Scan the available data and get the indices of the data to be used
*     left and right of the left brace point.

      call scanlgood(iibl,kww, nobs, nmfitl, -1, indfl, errl )
      call scanlgood(iibl,kww, nobs, nmfitl, +1, indfr, errr )

*     If no error in getting data then continue
      if( errl.eq.0 .and. errr.eq.0 ) then

*         Get the fractional epoch number of the midpoint
          midpt = (float(indfl(1)) + float(indfr(1)) )/2.d0

*         Now do the linear fits to the segments, returning the estimates
*         at the midpt and the RMS scatters of the data.

          call fitlin(WL1, indfl, nmfitl, estl(1), ratel(1), midpt,
     .                errl)
*                                     ! Fit L2 if L1 OK
          if( errl.eq.0 ) then
              call fitlin(WL2, indfl, nmfitl, estl(2), ratel(2),
     .                     midpt, errl)
          end if
          call fitlin(WL1, indfr, nmfitl, estr(1), rater(1), midpt,
     .                errr)
*                                     ! Fit L2 if L1 OK
          if( errr.eq.0 ) then
              call fitlin(WL2, indfr, nmfitl, estr(2), rater(2),
     .                     midpt, errr)
          end if
      end if

***** Check error again to see if we should continue
      if( errl.eq.0 .and. errr.eq.0 ) then

*         Get the misfit and their variances of L1 and L2 at the midpoint
          dest(1) = estl(1) - estr(1)
          dest(2) = estl(2) - estr(2)

CD         write(*,*) 'Dest ',dest(1), dest(2)

*         Now get the discontinuity in Lc and Lg at the midpoint
          call mislcg( dest, dcg, Ltocg, gear(jsat1))

          call gvarcg ( WL1, WL2, indfl, indfr, estl, estr, ratel,
     .                rater, midpt, ltocg, nmfitl, varcg )


*         Get the ambiquiuties spacing at L1 and L2.
          ambi(1) = abs(lambds(site1, chan1, 1))
          if( ambi(1).gt.0 ) then
              ambi(1) = 1/ambi(1)
          else
              errl = -997
          end if

          ambi(2) = abs(lambds(site1, chan1, 2))
          if( ambi(2).gt.0 ) then
              ambi(2) = 1/ambi(2)
          else
              errr = -997
          end if
      end if

***** continue if no error yet
      if ( errl.eq.0 .and. errr.eq.0 ) then

*         Now get the estimates of the ambiquities at L1 and L2 which
*         minimize the discontinuity in Lc and Lg.
CD         write(*,*) 'VARCG, AMBI ',varcg, ambi
          call estn( dcg, varcg, dest, Ltocg, yadd, confid, ambi, idd)

*         Save the current WL1 and WL2 values for undo (If we wanted to
*         here we could make the following code conditional in high
*         confidence in the results obtained)

*         Only perform the correction if so instructed
          if (doit) then
             call vcopy(kww, WL1, kbb, BL1, nobs)
             call vcopy(kww, WL2, kbb, BL2, nobs)
          endif

          call adyadd(yadd, iibl, iibr, WL1, WL2)
CD         write(*,*) 'YADD', yadd
*         Put out message giving confidence.
          write(mess,100) yadd, confid
 100      format(2(f5.1,1x), ' with ',f4.0,'% confid')
          if( prog_name(1:5).eq.'cview' ) then
              call gmsg(mwarn, mess)
          else
              call gmsgd(mwarn, mess) 
          endif
          replot = .true.

*                                 ! We hit error in finding data
      ELSE

*         write message and get out with replot .false.
c          write(*,*) 'Errscan errors ',errl,errr
          call erscan( mwarn, mess, errl, errr)
          replot = .false.
      END IF

***** Thats all
      return
      end

CTITLE SCANLGOOD

      subroutine scanlgood(iibl, kww, nobs, nmfitl, direct,
     .             indf, err)

*     This routine scans the LGOOD flags forward or backward from iibl
*     and builds up the indices of the NMFITL good points available.

* PASSED VARIABLES

*   iibl        - Start point of scan if direct is 1 or one plus
*               - point is direction is negative.
*   nobs        - Number of observations avaiable.  Used for checking
*               - bounds and dimensioning.
*   kww(nobs)   - Pointers to Lgood flags.
*   nmfitl          - Number of points to use in fit
*   direct      - Direction of scan (+1 for forward, -1 for backward)
*   indf(nmfitl)    - Indices to points to be used.
*   err         - error flag. 0 if OK, -1 if not enough data

      integer iibl, nobs, kww(nobs), nmfitl, direct, indf(nmfitl), err

*       lgood    - function to indicate good data

      logical lgood

* LOCAL VARIABLES

*   Start       - Start index (iibl is +direction, ibbl-1 if
*               - negative direction).
*   i           - Loop counter for scanning
*   j           - Current position in index array

      integer Start, i, j

****  Start, since we step based on direct, make sure we have valid value
      IF( direct.eq.1 .or. direct.eq.-1 ) THEN

*         get the first index to use
          if( direct.eq. -1 ) then
              start = iibl - 1
          else
              start = iibl
          end if

*         Initialize and start to scan.
          i = start
          j = 0
          do 100 while ( j.lt.nmfitl .and. i.ge.1 .and. i.le. nobs)

*             Check data
              if( lgood(kww(i)) ) then
                  j = j + 1
                  indf(j) = i
              end if
              i = i + direct
 100      continue

*         See if we found enough values
          if( j.lt. nmfitl ) then
              err = -1
          else
              err = 0
          end if

*                         ! We have a coding error.
      ELSE

*                         ! This will stop process
          err = -999

      END IF

****  Thats all
      return
      end

CTITLE FITLIN

      subroutine fitlin(WL, indf, nmfitl, est, rate, midpt, err)

*     Routine to take the data in WL pointed to be the INDF array
*     fit a line to it.  The estimate EST if referred to epoch MIDPT.
* MOD TAH 930519: Changed the slope estimation so that over large gaps
*     the slope is effectively set to zero.  This should avoid the problem
*     100's of cycles difference when large gaps are patched thus resulting
*     in much smaller residuals (although the bias flag can not be removed).
*     The implementation is to check if the number of epochs between the
*     first point and the last point is less than the number bewteen the
*     first and the midpoint.  If the gap is larger than the number of
*     data point then an apriori weight is put on the slope by introducing
*     a "fake" observation at the midpoint epoch whose sigma is inversely
*     proportional to the gap/data ratio.

* PASSED VARIABLES

*   nmfitl      - Number of data used in fit
*   indf(nmfitl)    - pointers to the data to be used.
*   err         - Error return: 0 if OK, -998 if det is zero or
*               - less than 0

      integer nmfitl, indf(nmfitl), err

*   WL(*)       - the carrier phase data (either L1 or L2)
*   est         - the estimate of the carrier phase at the
*               - epoch MIDPT
*   rate        - Rate of change of phase (units of time in index
*               - units.)
*   midpt       - epoch of for referring est to.  It is midway
*               - between the left and right data.

      real*8 WL(*), est, rate, midpt

* LOCAL VARIABLE

*   i,j         - Loop counters for forming estimates.
*   iel         - Current element in data arrays

      integer i,j, iel

*   det         - Determinate of the normal equations for the
*               - estimation.
*   norm(2,2)   - Normal equations of estimation.
*   inv(2,2)        - Inverse of normal equations
*   bvec(2)     - B vector for estimation.
*   sol(2)      - offset and slope of the phase wrt epoch number
*   ref         - Reference value for WL (removed from all data
*               - before estimation to avoid rounding errors in
*               - computing RMS)
*   dt          - difference in epoch number between current obs
*               - and midpt.

      real*8 det, norm(2,2), inv(2,2), bvec(2), sol(2), ref, dt

* MOD TAH 930519: New variables for weighting the slope in large
*     gaps
*   gaprat     - Ratio of gap beteen (first and midpoint) and (first
*                and last)
*   wgh        - Weight to be added to slope.

      real*8 gaprat, wgh

****  START, Clear the estimation arrays
      do 110 i = 1,2
          bvec(i) = 0.d0
          do 100 j = 1,2
              norm(i,j) = 0.d0
 100      continue
 110  continue
*
* MOD TAH 930519: Implement the slope weighting if the gap is large.`
*     NOTE: for right hand side data, the data runs backwards in time
*     and therefore indf(1) is the point closest to midpt.
*     See if we should weight the slope.
      if( abs(indf(1)-indf(nmfitl)).lt.abs(indf(1)-midpt) ) then
*         Yes, the gap is "large" so we will weight.
          gaprat = abs( (indf(1)-midpt)/(indf(1)-indf(nmfitl)) )
          wgh = 1.d0/(exp(-(gaprat-1.d0)) + 0.001)**2
*         Add the weight to the normal equations
          norm(2,2) = norm(2,2) + wgh
      end if
* END MOD TAH 930519:

      ref = WL(indf(1))

*     Loop over the data forming the estimation equations
      do 200 i = 1, nmfitl
          dt = indf(i) - midpt
          iel = indf(i)

*               T       T
*         form A  A and A  (WL-ref)
          bvec(1) = bvec(1) + (WL(iel)-ref)
          bvec(2) = bvec(2) + (WL(iel)-ref)*dt

          norm(1,1) = norm(1,1) + 1
          norm(1,2) = norm(1,2) + dt
          norm(2,2) = norm(2,2) + dt*dt
 200  continue

*     Invert normal equations and solve estimate
      norm(2,1) = norm(1,2)

      det = norm(1,1)*norm(2,2) - norm(1,2)*norm(2,1)

*     Check det before continuing
      if( det.gt.0.d0 ) then
          inv(1,1) = norm(2,2)/det
          inv(1,2) = -norm(1,2)/det
          inv(2,1) = inv(1,2)
          inv(2,2) = norm(1,1)/det

*         Now get the estimates
          sol(1) = inv(1,1)*bvec(1) + inv(1,2)*bvec(2)
          sol(2) = inv(2,1)*bvec(1) + inv(2,2)*bvec(2)

****      Finish computation, return full values
          est = sol(1) + ref
          rate = sol(2)
          err = 0
      ELSE

*****     Determinate is zero.... Something is vry wrong
          err = -998

      END IF

***** Thats all
      return
      end

CTITLE MISLCG

      subroutine mislcg( dest,  dcg, Ltocg, gear)

*     Routine to compute the differences in LC and LG accross the boundary
*     of the left and right segments.  We define here that
*     L1 = LC + LG
*     L2 = gear*LC + LG/gear
*     where is the ratio of L2 frequency to L1 frequency


* PASSED VARIABLES

*   dest(2)     - Differences in L1 and L2 computed in sense of
*               - right segment minus left segment
*   dcg(2)      - Differences of LC and LG from the left and right
*               - segmenst. (cycles of L1 phase)
*   ltocg(2,2)  - Transformation from L1 and L2 phase to LC and LG.
*   gear        - Ratio of L2 and L1 frequencies (f2/f1)


      real*8 dest(2), dcg(2), ltocg(2,2), gear

* LOCAL VARIABLES


*   beta        - beta = 1- g**2 (used enough to be defined)


      real*8 beta

***** Start, define ratio of L1 to L2
      beta = (1.d0 - gear*gear)

*     Get the transformation from L1 and L2 to Lc and Lg.  Derived from
*     the equation above.
      Ltocg(1,1) = 1/beta
      Ltocg(1,2) = -gear/beta

      Ltocg(2,1) = -gear*gear/beta
      Ltocg(2,2) = gear/beta

*     Compute the values of Lc and Lg differences
      call getlcg( dest, ltocg, dcg )

***** Thats all
      return
      end

CTITLE ESTN

      subroutine estn( dcg, varcg, dest, Ltocg, yadd, confid, ambi, idd)

*     Routine to estimate the ambiguity in L1 and L2 based on trying to
*     keep continuity in the Lc and Lg observables.  We weight the Lc and
*     Lg discontinuity by


* PASSED VARIABLES

*   dcg(2)      - Differences in Lc and Lg between the left and
*               - right segments
*   varcg(2)        - Variances of Lc and Lg computed from rms scatter
*   dest(2)     - Difference in L1 and L2 at midpoint of left and
*               - right boundaries. (cycles)
*   Ltocg(2,2)  - Transformation matrix which gets us from L1 and L2
*               - to Lc and Lg
*   yadd(2)     - Number of cycles to be added to L1 and L2 to get
*               - as continuous LC and LG as possible.
*   confid      - Estimated confidence in patch.
*   ambi(2)     - Ambiguities at L1 and L2 in cycles.
*   wambi(2)    - Working ambiguities at L1 and L2 in cycles.


      real*8 dcg(2), varcg(2), dest(2), Ltocg(2,2), yadd(2), confid,
     .    ambi(2),wambi(2),tmp(2)

* LOCAL VARIABLES

*       maxchi  - Maximum number of values of chi^2 to save

      integer maxchi

      parameter ( maxchi = 3 )

*   Nguess(2)   - First guess of cycles skips at L1 and L2 in
*               - units of the ambiquity spacing at the two
*               - frequencies
*   i,j         - Loop counters for testing trial numbers of cycle
*               - slips.
*   bestij(2, maxchi)   - Values of i and j ambiquiuities which generate
*               - the smallest value of chi squared.

      integer Nguess(2), i,j, idd, bestij(2, maxchi)

*   chi(maxchi) - Best values of chi**2 found, ranked in ascending
*               - order.  Values here go with bestij values.
*   Loff(2)     - trial values of the offsets between L1 and L2.
*   dcgtry(2)   - Values of the Lc and Lg discontinuity with the
*               - current values of the ambiguity
*   chitry      - Value of chi^2 for current trial values of
*               - ambiguities.


      real*8 chi(maxchi), Loff(2), dcgtry(2), chitry

****  Load the working ambiguities - regular or 0.5 for double-difference
*                                             which allows 0.5 patching

C***rwk : comment out this code (pending removal) since we don't want to
c         PATCH half-cycles (too dangerous), only MOVE them occasionally.
c     if ((idd .ne. 0) .and.
c    .   (dabs(dest(1)) .le. 0.8) .and.
c    .   (dabs(dest(2)) .le. 0.8)) then
c           wambi(1) = 0.5d0
c           wambi(2) = 0.5d0
c     else
            wambi(1) = ambi(1)
            wambi(2) = ambi(2)
c     endif

****  Get the first guess at the cycle slips from the L1 and L2 differences

      tmp(1) = dest(1)/wambi(1)
      tmp(2) = dest(2)/wambi(2)
      if(dabs(tmp(1)) .le. 200000000.) then
           Nguess(1) = idnint(tmp(1))
      else
           Nguess(1) = 200000000
           if(tmp(1) .lt. 0.0) Nguess(1) = -200000000
      endif

      if(dabs(tmp(2)) .le. 200000000.) then
           Nguess(2) = idnint(tmp(2))
      else
           Nguess(2) = 200000000
           if(tmp(2) .lt. 0.0) Nguess(2) = -200000000
      endif

*     Initialize the chi^2 values to be saved
      do 100 i = 1, maxchi
          chi(i) = 1.d20
 100  continue

***** Put some limits on the variance estimates of LC and LG so that one
*     doesnot totally dominate the solution
      varcg(1) = max(varcg(1), 0.005d0)
      varcg(2) = max(varcg(2), 0.005d0)
CD     write(*,*) ' dcg, varcg ', dcg, varcg

****  Now loop over the search space of +-2 about first guess

      do 210 i = Nguess(1)-2, Nguess(1)+2
          do 200 j = Nguess(2)-2, Nguess(2)+2

*             Load up the ambiguity array
              Loff(1) = i*wambi(1)
              Loff(2) = j*wambi(2)

*             Get the errors in discontinuity in Lc and Lg
              dcgtry(1) = dcg(1) -
     .                    (Ltocg(1,1)*Loff(1) + Ltocg(1,2)*Loff(2))
              dcgtry(2) = dcg(2) -
     .                    (Ltocg(2,1)*Loff(1) + Ltocg(2,2)*Loff(2))

*             Compute Chi^2 of for this trial
              chitry = dcgtry(1)**2/varcg(1) + dcgtry(2)**2/varcg(2)

*             Now rank the trial chi and save indices
              call rankc( i, j, chitry, bestij, chi, maxchi)

 200      continue
 210  continue

****  Now see how we did, return the best values and give confidence
*     based on chi ratios.

      yadd(1) = bestij(1,1)*wambi(1)
      yadd(2) = bestij(2,1)*wambi(2)
CD     write(*,*) 'Bestij, chi',bestij, chi

      confid = (1 - chi(1)/chi(2))*100.d0

****  Thats all
      return
      end

CTITLE RANKC

      subroutine rankc( offL1, offL2, chitry, bestij, chi, maxchi)

*     Routine to rank chi^2 in increasing order.  Current value is
*     tested against values in table and inserted at the appropriate
*     place.  The values of the ambiguity at L1 and L2 are saved in
*     bestij array.,

* PASSED VARIABLES

*   maxchi          - Maximum number of chi^2 values to save
*   offL1, offL2        - Current values of the ambiquities at L1
*                   - and L2
*   bestij(2,maxchi)    - Saved values of offL1, and offL2.

      integer maxchi, offL1, offL2, bestij(2,maxchi)

*   chi(maxchi)     - Saved values of chi^2 (in ascending order)
*   chitry          - current value of chi^2 to be tested.


      real*8 chi(maxchi), chitry

* LOCAL VARIABLES

*   i,j             - Loop counters for scanning table

      integer i,j

*   bigger          - Inticates chitry is bigger than table
*                   - entry.

      logical bigger

****  Find value that we are less than

      i = 0
      bigger = .true.
      do 150 while ( bigger .and. i.lt. maxchi)
          i = i + 1
*                                     ! found place to insert
          if( chitry.lt.chi(i) ) then

*             Move up current values
              do j = maxchi-1, i, -1
                  chi(j+1) = chi(j)
                  bestij(1,j+1) = bestij(1,j)
                  bestij(2,j+1) = bestij(2,j)
              end do

*             Now add new values
              chi(i) = chitry
              bestij(1,i) = offL1
              bestij(2,i) = offL2

              bigger = .false.
          end if
 150  continue

****  Thats all, values inserted
      return
      end

CTITLE ADYADD

      subroutine adyadd(yadd, iibl, iibr, WL1, WL2)

* PASSED VARIABLES

*   iibl        - Left point to start at in adding ambiquity
*   iibr        - Right point to stop at.

      integer iibl, iibr

*   yadd(2)     - The number cycles to add at L1 and L2 (units of
*               - cycles at each frequency, but may be half cycle
*               - for codeless recievers.
*   WL1(*)      - L1 carrier phase values
*   WL2(*)      - L2 carrier phase values

      real*8 yadd(2), WL1(*), WL2(*)

* LOCAL VARIABLES

*   i           - Loop counter for doing the addition

      integer i

****  Simply add the values
      do 100 i = iibl, iibr
          WL1(i) = WL1(i) + yadd(1)
          WL2(i) = WL2(i) + yadd(2)
 100  continue

***** Thats all
      return
      end

CTITLE GVARCG

      subroutine gvarcg ( WL1,WL2, indfl, indfr, estl, estr, ratel,
     .                rater, midpt, ltocg, nmfitl, varcg )

*     Routine to compute the rms scatters of Lc and Lg on in the left
*     and right segments of data.  The routine estimates of the variances
*     of the Lc and Lg differences accross the boundary.


* PASSED VARIABLES

*   nmfitl          - NUmber of data fit on each side.
*   indfl(nmfitl)   - Indices pointing to the good data in the left
*                   - segment.
*   indfr(nmfitl)   - Indices pointing to the good data in the right
*                   - segment

      integer nmfitl, indfl(nmfitl), indfr(nmfitl)

*   WL1(*)          - L1 carrier phase data (cycles)
*   WL2(*)          - L2 carrier phase data (cycles)
*   estl(2)         - Estimates of offset in left segment for L1 and
*                   - L2.
*   estr(2)         - Estimates of offset in right segment for L1 and
*                   - L2.
*   ratel(2)            - Carrier rate at L1 and L2 in left segment
*                   - in units of index.
*   rater(2)            - Carrier rate at L1 and L2 for right segment.

*   midpt           - epoch of the midpoint between the two segment.
*   ltocg(2,2)      - Trandformation from L1 and L2 to Lc and Lg.

*   varcg(2)            - Estimated variances of the data in the two
*                   - segments for Lc anf Lg.  Used for weighting the
*                   - discontinuity.

      real*8 WL1(*), WL2(*), estl(2), estr(2), ratel(2), rater(2),
     .    midpt, ltocg(2,2), varcg(2)

* LOCAL VARIABLES

*   i           - Loop counter for samplinf data in line fit
*   iel         - index in Wl1 or WL2 of the residual being computed.

      integer i, iel

*   dt          - Time difference between current data and midpt (units
*               - of index)
*   resL12(2)   - residuals of L1 and L2 to fit to line
*   resLcg(2)   - computed Lc and Lg residuals
*   chi(2)      - summation variable for lc and Lg residuals squared.


      real*8 dt, resL12(2), resLcg(2), chi(2)

****  Start, Loop over the data in each segment, compute Lc and Lg residuals
*     to line and compute the rms^2 of these residuals

      chi(1) = 0.d0
      chi(2) = 0.d0

      do 100 i = 1, nmfitl
          dt = indfr(i) - midpt
          iel = indfr(i)

          resL12(1) = WL1(iel) - (estr(1) + dt*rater(1))
          resL12(2) = WL2(iel) - (estr(2) + dt*rater(2))

          call getlcg( resl12, ltocg, reslcg)

          chi(1) = chi(1) + reslcg(1)**2
          chi(2) = chi(2) + reslcg(2)**2

 100  continue

      do 150 i = 1, nmfitl
          dt = indfl(i) - midpt
          iel = indfl(i)

          resL12(1) = WL1(iel) - (estl(1) + dt*ratel(1))
          resL12(2) = WL2(iel) - (estl(2) + dt*ratel(2))

          call getlcg( resl12, ltocg, reslcg)

          chi(1) = chi(1) + reslcg(1)**2
          chi(2) = chi(2) + reslcg(2)**2

 150  continue

****  Now finishup computation of RMS^2 values
      varcg(1) = chi(1)/(2*(nmfitl-2))
      varcg(2) = chi(2)/(2*(nmfitl-2))

****  Thats all
      return
      end

CTITLE GETLCG

      subroutine getlcg( resL12, ltocg, resLcg )

*     Mutliplies ResL12 by Ltocg to compute the residuals at Lc ang Lg.

* PASSED VARIABLES

*   resL12(2)       - residuals at L1 and L2
*   resLcg(2)       - Computed residuals at Lc and Lg
*   Ltocg(2,2)      - Transformation between the two types of
*                   - observables

      real*8 resL12(2), resLcg(2), Ltocg(2,2)

* LOCAL VARIABLES

*       i,j         - Loop counters for doing mulitplication

      integer i,j

****  This is just a matrix mulitply

      do 110 i = 1, 2
          resLcg(i) = 0.d0
          do 100 j = 1,2
              resLcg(i) = resLcg(i) + ltocg(i,j)*resL12(j)
 100      continue
 110  continue

****  Thats all
      return
      end

CTITLE ERSCAN

      subroutine erscan( mwarn, mess, errl, errr)

*     This routine will output any error messages from the patch
*     routine in the the warning window of the graphics package.

* PASSED VARIABLES

*   mwarn(4)        - The warning window

      integer*2 mwarn(4)

*   errl        - Error number from left segment of data
*   errr        - error number from righ segment of data.
*               - Current values are:
*               -    0 - All OK
*               -   -1 - Not enough data
*               - -997 - frequency not correct at L1 or L2
*               - -998 - Det zero in fiting line to data.
*               - -999 - Invalid direction passed to scan

      integer errl, errr

*   mess        - Message to be output to screen


      character*(*) mess

* LOCAL VARIABLES

*   pos         - Start position in message string for right
*               - segment message.

      integer pos  

      integer*4 len,rcpar

      character*80 prog_name
 
c     get calling program name to know which version of gmsg to call
      len = rcpar(0,prog_name)



****  See what message we have.
      mess = ' '
      if( errl.eq.-1 ) then
          mess = 'Not enough data in left segment'
      else if (errl.eq.-997 ) then
          mess = 'No L1 frequnecy'
      else if (errl.eq.-998 ) then
          mess = 'Line fit singular in left segmt'
      else if (errl.eq.-999 ) then
          mess = 'Invalid direction in left segmt'
      end if

*     Now concatinate the right segment error (if I knew correct call I
*     would just get length of string.
      if( errl.eq.0 ) then
          pos = 1
      else
          pos = 31
      end if

      if( errr.eq.-1 ) then
          mess(pos:) = 'Not enough data in right segment'
      else if( errr.eq.-997 ) then
          mess(pos:) = 'No L2 frequency'
      else if (errr.eq.-998 ) then
          mess(pos:) = 'Line fit singular in right segmt'
      else if (errr.eq.-999 ) then
          mess(pos:) = 'Invalid direction in right segmt'
      end if

****  write message 
      if( prog_name(1:5).eq.'cview' ) then
          call gmsg(mwarn, mess)
      else
          call gmsgd(mwarn, mess) 
      endif

****  Thats all
      return
      end


C******************************************************************************
C******************************************************************************
C******************************************************************************
C******************************************************************************
CTITLE PATCH_LGWL

      subroutine patch_lgwl(iibl,iibr, ii0, ii1, mwarn, site1, chan1,
     .                nobs,replot, nmfitl, yadd, confid, doit)


      integer mxfitl

      parameter ( mxfitl = 15 )

      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
      include '../includes/errflg.h'


      integer iibl, iibr, ii0, ii1, site1, chan1, nobs
      integer*2 mwarn(4)

      logical Replot,doit


      integer nmfitl, indfl(mxfitl), indfr(mxfitl), errl,errr, i


      real*8  estl(2), estr(2), ratel(2), rater(2), dest(2)
     .     ,   midpt, dcg(2), Ltocg(2,2), confid, yadd(2), rnhalf

      logical         slipwl
      real*8          dwl,dy1,dy2
      character*96    amsg


      character*40 mess

       
      integer*4 len,rcpar

      character*80 prog_name
 
c     get calling program name to know which version of gmsg to call
      len = rcpar(0,prog_name)


      if (nmfitl .le. 0 .or. nmfitl .gt. mxfitl) nmfitl = mxfitl

*     Scan the available data and get the indices of the data to be used
*     left and right of the left brace point.

      call scanlgood(iibl,kww, nobs, nmfitl, -1, indfl, errl )
      call scanlgood(iibl,kww, nobs, nmfitl, +1, indfr, errr )

*     If no error in getting data then continue
      if( errl.eq.0 .and. errr.eq.0 ) then

*         Get the fractional epoch number of the midpoint
          midpt = (float(indfl(1)) + float(indfr(1)) )/2.d0

*         Now do the linear fits to the segments, returning the estimates
*         at the midpt and the RMS scatters of the data.

          call fitlin(WL1, indfl, nmfitl, estl(1), ratel(1), midpt,
     .                errl)
*                                     ! Fit L2 if L1 OK
          if( errl.eq.0 ) then
              call fitlin(WL2, indfl, nmfitl, estl(2), ratel(2),
     .                     midpt, errl)
          end if
          call fitlin(WL1, indfr, nmfitl, estr(1), rater(1), midpt,
     .                errr)
*                                     ! Fit L2 if L1 OK
          if( errr.eq.0 ) then
              call fitlin(WL2, indfr, nmfitl, estr(2), rater(2),
     .                     midpt, errr)
          end if
      end if

***** Check error again to see if we should continue
      if( errl.eq.0 .and. errr.eq.0 ) then

*         Get the misfit and their variances of L1 and L2 at the midpoint
          dest(1) = estl(1) - estr(1)
          dest(2) = estl(2) - estr(2)


*         Now get the discontinuity in Lc and Lg at the midpoint
          call mislcg( dest, dcg, Ltocg, gear(jsat1))

*         Now get the discontinuity in WL
          call check_wl (ii0,ii1,iibl,slipwl,dwl)

c        if user clicks in window for LG, PATCH according
c            dWL = L2 - L1      = dwl
c            dLG = L2 - gL1     = dcg(1)
CD        write(*,*) 'Dest ',dest(1), dest(2), dcg, dwl

* MOD TAH 930128: Haved the number of cycles in patch.  Should
*         looked at whey this is needed.
C         dy1 = (dcg(2)-dwl)/(1.0d0-gear(jsat1))
          dy1 = (dcg(2)-dwl)/(1.0d0-gear(jsat1))/2
          dy2 = dwl+dy1

c         For lamba = +/- 2 cycle slip may be be a half-integer.
c         Otherwise, it should be a whole integer.

* MOD TAH 930519: Changed the rnhalf call to only L2 so that half cycles
*         will not be added to L1 (To be complete should check wavelength
*         factor at L1 but will assume here that it is always full
*         wavelength for wide-lane patches.

          dy1 = dnint(dy1)

          if (abs(lambds(site1,chan1,2)) .eq. 2) then
c           half-cycle -slip
            dy2 = rnhalf(dy2)
          else
c            cycle slip must be an integer
             dy2 = dnint(dy2)
          endif

          yadd(1) = dy1
          yadd(2) = dy2

c         add yadd to all points
          do i=iibl,iibr
             wl1(i) = wl1(i) + dy1
             wl2(i) = wl2(i) + dy2
          enddo

          write (amsg,100) site1,chan1,dy1,dy2
 100      format('LG-WL: ',2(i2,1x),2(f8.1,1x),f4.0)  
          if( prog_name(1:5).eq.'cview' ) then
              call gmsg(mwarn, mess)
          else
              call gmsgd(mwarn, mess) 
          endif

      ELSE

*         write message and get out with replot .false.
c          write(*,*) 'Errscan errors ',errl,errr
          call erscan( mwarn, mess, errl, errr)
          replot = .false.
      END IF

***** Thats all
      return
      end

