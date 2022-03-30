      subroutine check_wl (
     .                 ii0,ii1,islip,slipwl,dwl)

      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
      include '../includes/errflg.h'

      integer         ii0,ii1,islip
      integer         i,i1,n1,n2,istart,iend
      real*8          ywl0,ywl1,ywli,dy
      real*8          sum1,sum2,dsum,dwl
      logical         lgood,slipwl,good

* MOD TAH 930128: Added code so that widelane patching will
*     work for sparse sampled data.
* VARIABLES
*  inext - The largest sampling interval for the possible
*          two sites in the data. (Value obtained from jnext array
*     include '../includes/cview.h'

      i1 = islip - 1
      good = .false.
      slipwl = .false.
      sum1 = 0.d0
      sum2 = 0.d0
      dwl  = 0.d0
      n1 = 0
      n2 = 0

* MOD TAH 930128: Added comment
*     Start at the one epoch before the point to patch
*     (defined by islip) and search bak until we find
*     a good data point.  (Don't go back past the
*     start of the data).
      do 100 while ((.not. good) .and. (i1 .gt. ii0))
           if(lgood(kww(i1))) then
                good = .true.
           else
                i1 = i1 - 1
           endif
  100 continue

*     Get the sampling interval to make sure that we
*     go enough epochs to get some data.  We don't
*     the site information so we will choose the
*     maximuim number to make sure we are covered. We
*     don't number of cfiles either at this point
*     so we will loop over maximum number.
      do i = 1, ncvsit
         inext = max(inext,jnext(i))
      end do

*     This should not be needed but just in case
      if( inext.eq.0 ) inext = 8

*
* MOD TAH 930128: Added comment
*     If we found data before the patch point, then
*     contiue on.
      if (good) then

*          Get the wide lane values at the slip and
*          at the first good point before this.
           ywl0 = wl2(islip)-wl1(islip) +
     .           facwl(jsat1)*(pc1(islip)+pc2(islip))
           ywl1 = wl2(i1)-wl1(i1) +
     .           facwl(jsat1)*(pc1(i1)+pc2(i1))

           dy = dabs(ywl0 - ywl1)

           if(dy .gt. 0.85) then

*              Now scan out far enough so that we get
*              some data.  (Use inext to cover 4 four
*              data points at least). (Multiplied the
*              original 10 points by 10*inext
               istart = max(ii0,i1-10*inext)
               iend   = min(ii1,islip+10*inext)

               do 200 i=istart,i1
                   if(lgood(kww(i))) then
                       ywli = wl2(i)-wl1(i) +
     .                               facwl(jsat1)*(pc1(i)+pc2(i))

* MOD TAH 930128: Changed the limit acceptable to 2.0 cycles
*                      from 5 since we can now be scanning
*                      a longer interval of data.
                       if (dabs(ywl1-ywli) .lt. 2.0) then
                            sum1 = sum1 + ywli
                            n1 = n1 + 1
                       endif
                   endif
  200          continue

               do 300 i=islip,iend
                   if(lgood(kww(i))) then
                       ywli = wl2(i)-wl1(i) +
     .                               facwl(jsat1)*(pc1(i)+pc2(i))
* MOD TAH 930128: Changed the limit acceptable to 2.0 cycles
*                      from 5 since we can now be scanning
*                      a longer interval of data.
                       if (dabs(ywl0-ywli) .lt. 2.0) then
                            sum2 = sum2 + ywli
                            n2 = n2 + 1
                       endif
                   endif
  300          continue

               if ((n1 .ge. 4) .and. (n2 .ge. 4)) then

                   sum1 = sum1/dble(n1)
                   sum2 = sum2/dble(n2)

                   dsum = dabs(sum1-sum2)
                   if(dsum .gt. 0.8) then
                         slipwl = .true.
                         dwl = sum1-sum2
                   endif
               else
* MOD TAH 930128: Place holder for warning.
*                  Warn user that there is not enough data.
*                  Can't call this since mwarn is not avaiable.
C                  call gmsg (mwarn,'**Not Enough WL data**')
               endif
           endif

      endif

      return
      end
c************************************************************************


