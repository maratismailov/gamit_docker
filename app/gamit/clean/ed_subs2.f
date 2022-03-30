c**********************************************************************
      subroutine editor_plotts (ii0,ii1,jj0,jj1,ix0,iy0,xs,ys,
     .ymid,ymin,ymax,iradi,plt,kpl,domarg,hide_b,lclipt,
     .isite1,isite2 )

c     modified to plot vertical bar on last point kurt June 91
c
c     actually plot the time series
c
c     standard GAMIT array dimensions
      include '../includes/dimpar.h'

c     some non GAMIT declarations
      include '../includes/makex.h'

c     common block declarations
      include '../includes/cview.h'

c     standard GAMIT error flags
      include '../includes/errflg.h'
c     FUNCTIONS:
c     functions to sort out error codes
c        returns .true. for observations good enough to be used in a solution
         logical lgood
c        returns .true. for marginal observations which might potentially useful
         logical lmarg
c        returns .true. for observations with a bias flag
         logical lbias,llbias
c        returns .true. for low elevation observations which might potentially useful
         logical lloel
c

      logical  domarg,hide_b

c     Has clipping occurred?
      logical lclipt

c%%   values which must be integer*2 for the GSUBS package
      integer*2 iradi,ix1,iy1,ix2,iy2,ix0,iy0
c     values for diagonal bias flags
      integer*2 ix11,iy11,ix12,iy12, ix21,iy21,ix22,iy22

c%%   this system status variable may need to be integer*2 on some systems
      integer*4 igerr

c     indexes into the time series
      integer  jj0,jj1,iii

c     actual start and end
      integer ii0,ii1

c     scales and such
      real*8 xs,ys,ymid,ymin,ymax,ytt,ttt
      integer i,iinext,ilast

c     what we are plotting
      real*8  plt(maxepc)
      integer*4 kpl(maxepc)

c     For diagonal bias flags
      character*4 bias_count,char_bias,ctmp
      common /biascount/bias_count(maxepc)

      integer*4 isite1, isite2
      integer*4 ioffst 

      ytt = 0.d0 
      ilast = 0
      ioffst = 0

* MOD TAH 930128: Fixed data plotting problems when
*     data not recorded at the sampling interval and does
*     not start on epoch number 1.
* New variables
*  ioffst -- Gives the number of epochs between the
*            current start epoch (jj0) and when the first
*            data point should appear minus -1.  (The -1
*            is needed becuase of the incrementation in the
*            loop below where inext is incremented before it
*            is checked.
* WARNING: Since the channel is not passed to this routine
*          the first channel is used to get the ioffst value
*          If there is no data in the first channel then
*          ioffst will not be correct.  Solution is to pass
*          ichan1 and ichan2 into routine.

c     plot de-meaned data
c     filled circles at good data points
c     open circles at marginal data points
c     line segments connecting marginal points
c     vertical line for points with new bias parameter

c     OK, do a dumb thing.
c     When plotting with STACK (lstack = .true.) some HUGE
c     values are possible.  So, have 2 version of the routine
c     depending.

c     finds number of epoches between two data points
* MOD TAH 930128: Changed index calculation so we will
*     start on a point with data when data not recorded at
*     sampling interval in the Xfile.

      inext = jnext(isite1)
*     Scan up the channels until we find a non-3000
*     value. (Just in case chan 1 has no data).  We
*     don't know the number of satellites so use
*     maximum number (maxsat)
      do i = 1, maxsat
***          if( kk0(i,isite1).ne.3000 ) then
***Modified  M.B. Oral Wed May 12 20:27:55 EDT 1993
          if( kk0(i,isite1).ne.maxepc ) then
              ioffst = mod(int(kk0(i,isite1)),inext)
          end if
      end do

      if(isite2 .gt. 0) then
          if(jnext(isite2) .gt. inext) then
             inext = jnext(isite2)
*            Again scan to make sure we get a value
*            value for ioffst. Solved if ichan passed.
             do i = 1,maxsat
***                 if( kk0(i,isite2).ne.3000 ) then
***Modified  M.B. Oral Wed May 12 20:27:55 EDT 1993
                 if( kk0(i,isite2).ne.maxepc ) then
                     ioffst = mod(int(kk0(i,isite2)),inext)
                 end if
             end do
          end if
      endif

c     Fast, dangerous way.
      if (.not. lclipt) then
          iinext = mod(jj0-ioffst,inext)-1

* MOD TAH 930128: Change the calculation of the last
*         point to make sure we get all points.
C         ilast  = jj1 - mod(jj1-1,inext)
* MOD TAH 930129: set ilast to be last point, the
*         loop below only goes to the next to last
*         point becuase the last point is done
*         by itself (after the end of the 1070 loop).

          ilast  = jj1 - mod(jj1-ioffst,inext)

* MOD TAH 930128: Comment added but no changes made
*         to ilast calculation.
*         Strictly the computation of ilast should be
*         changed but is seems at the moment that we
*         may just go a little further than we should
*         so this is probably not a problem

         do 1070 i = jj0,ilast-inext
            iinext = iinext + 1
            if (iinext .eq. inext) iinext = 0
            if ((iinext .eq. 0) .and.
     .      (lgood(kpl(i)).or.(domarg.and.lmarg(kpl(i))))) then
               ytt = ys * (plt(i) - ymid)
               ttt = xs * (tx(i) - tx(ii0))
               ix1 = int(ix0 + ttt)
               ix11 = int(ix0 + ttt)
               iy1 = int(iy0 - ytt)
               if (lgood(kpl(i))) then
                  call gdot(ix1,iy1,iradi,igerr)
                  if (lgood(kpl(i+inext))) then
                     ytt = ys * (plt(i+inext) - ymid)
                     ttt = xs * (tx(i+inext) - tx(ii0))
                     ix2 = int(ix0 + ttt)
                     iy2 = int(iy0 - ytt)
                     call gline (ix11,iy1,ix2,iy2,igerr)
                  endif

                  llbias = lbias(kpl(i))
                  iii = max(i-1,1)
                  if(llbias .or. lgood(kpl(iii))) then
                       char_bias = bias_count(i)
                  else
CCC                    PUSH BIAS
                       char_bias = "0000"
                       do while (.not.lgood(kpl(iii)) .and.
     .                                          (iii.gt.ii0))
                           ctmp = bias_count(iii)
                           if( ctmp .ne. "0000") then
                              if( ctmp(1:1) .eq. "1")
     .                                char_bias(1:1) = "1"
                              if( ctmp(2:2) .eq. "1")
     .                                char_bias(2:2) = "1"
                              if( ctmp(3:3) .eq. "1")
     .                                char_bias(3:3) = "1"
                              if( ctmp(4:4) .eq. "1")
     .                                char_bias(4:4) = "1"
                           endif
                           iii = iii - 1
                       enddo
                       if( char_bias .ne. "0000") then
                           ytt = ys * (plt(i) - ymid)
                           ttt = xs * (tx(i) - tx(ii0))
                           ix1 = int(ix0 + ttt)
                           iy1 = int(iy0 - ytt)
                           ix11 = int(ix1 -  4*iradi)
                           iy11 = int(iy1 +  4*iradi)
                           call  gline(ix1,iy1,ix11,iy11,igerr)
                           iy11 = int(iy1 -  4*iradi)
                           call  gline(ix1,iy1,ix11,iy11,igerr)
                           iy1  = int(iy1 +  4*iradi)
                           call  gline(ix11,iy1,ix11,iy11,igerr)
C                           char_bias = "0000"

                           llbias = .true.
                       endif
                  endif

                  if (llbias) then
                     ytt = ys * (plt(i) - ymid)
                     ttt = xs * (tx(i) - tx(ii0))
                     ix1 = int(ix0 + ttt)
                     ix2 = int(ix0 + ttt)
                     iy1 = int(iy0 - ytt - 8*iradi)
                     iy2 = int(iy0 - ytt + 8*iradi)
                     call gline (ix1,iy1,ix2,iy2,igerr)


c* indicating the one-way
c*
c*        21       22                     Bottom: Sat 1     Top   : Sat 2
c*            \ /
c*             |                          Left  : Site 1    Right : Site 2
c*             |
c*             |
c*            / \
c*
c*       11        12

c                channel:1     site:1
                   if( char_bias(1:1) .eq. "1") then
                      ix11 = int(ix2 -  4*iradi)
                      iy11 = int(iy2 +  4*iradi)
                      call  gline(ix2,iy2,ix11,iy11,igerr)
                   endif
c                channel:1     site:2
                   if( char_bias(2:2) .eq. "1")  then
                      ix12 = int(ix2 +  4*iradi)
                      iy12 = int(iy2 +  4*iradi)
                      call gline(ix2,iy2,ix12,iy12,igerr)
                   endif
c                channel:2     site:1
                   if( char_bias(3:3) .eq. "1") then
                      ix21 = int(ix1 -  4*iradi)
                      iy21 = int(iy1 -  4*iradi)
                      call gline(ix1,iy1,ix21,iy21,igerr)
                   endif
c                channel:2     site:2
                    if( char_bias(4:4) .eq. "1") then
                       ix22 = int(ix1 +  4*iradi)
                       iy22 = int(iy1 -  4*iradi)
                       call gline(ix1,iy1,ix22,iy22,igerr)
                    endif

                  endif
               else if (lloel(kpl(i))) then
                  call grect(ix1,iy1,iradi,igerr)
               else if (lmarg(kpl(i))) then
                  call gcirc(ix1,iy1,iradi,igerr)
               endif
            endif
            if (.not. hide_b) then
               char_bias = bias_count(i)
               if(char_bias.ne."0000")   then
                   ytt = ys * (ymin - ymid)
                   ttt = xs * (tx(i) - tx(ii0))
                   ix1 = int(ix0 + ttt)
                   iy1 = int(iy0 - ytt)
                   iy2 = int(iy0 - ytt - 8*iradi)
                   call gline (ix1,iy1,ix1,iy2,igerr)
                   call gcirc(ix1,iy1,iradi,igerr)
               endif
            endif
 1070    continue

c        symbol for for final point
         if (lgood(kpl(ilast)).or.
     .              (domarg.and.lmarg(kpl(ilast)))) then
            ytt = plt(ilast) - ymid
            ttt = tx(ilast)-tx(ii0)
            ix1 = ix0 + int(xs*ttt)
            iy1 = iy0 - int(ys*ytt)
            if (lgood(kpl(ilast))) then
               call gdot(ix1,iy1,iradi,igerr)
            else if (lloel(kpl(ilast))) then
               call grect(ix1,iy1,iradi,igerr)
            else if (lmarg(kpl(ilast))) then
               call gcirc(ix1,iy1,iradi,igerr)
            endif
c           vertical line for bias flag
            if (lbias(kpl(ilast))) then
               ytt = ys * (plt(ilast) - ymid)
               ttt = xs * (tx(ilast) - tx(ii0))
               ix1 = int(ix0 + ttt)
               ix2 = int(ix0 + ttt)
               iy1 = int(iy0 - ytt - 8*iradi)
               iy2 = int(iy0 - ytt + 8*iradi)
               call gline (ix1,iy1,ix2,iy2,igerr)

c**       modifications for diagonal bias flags
               char_bias = bias_count(i)
      if(char_bias.ne."0000")   then
channel:1     site:1
         if( char_bias(1:1) .eq. "1") then
            ix11 = int(ix2 -  4*iradi)
            iy11 = int(iy2 +  4*iradi)
            call  gline(ix2,iy2,ix11,iy11,igerr)
         endif
channel:1     site:2
         if( char_bias(2:2) .eq. "1")  then
            ix12 = int(ix2 +  4*iradi)
            iy12 = int(iy2 +  4*iradi)
            call gline(ix2,iy2,ix12,iy12,igerr)
         endif
channel:2     site:1
         if( char_bias(3:3) .eq. "1") then
            ix21 = int(ix1 -  4*iradi)
            iy21 = int(iy1 -  4*iradi)
            call gline(ix1,iy1,ix21,iy21,igerr)
         endif
channel:2     site:2
         if( char_bias(4:4) .eq. "1") then
            ix22 = int(ix1 +  4*iradi)
            iy22 = int(iy1 -  4*iradi)
            call gline(ix1,iy1,ix22,iy22,igerr)
         endif
      endif
c**   end modifications for diagonal bias flags
            endif
         endif
c     slow, conservative way
      else
         do 2070 i = jj0,ilast-inext
            if (lgood(kpl(i)).or.(domarg.and.lmarg(kpl(i)))) then
               ytt = ys * (plt(i) - ymid)
               if (dabs(ytt) .lt. dabs(ymax-ymin)) then
                  ttt = xs * (tx(i) - tx(ii0))
                  ix1 = int(ix0 + ttt)
                  iy1 = int(iy0 - ytt)
                  if (lgood(kpl(i))) then
                     call gdot(ix1,iy1,iradi,igerr)
                     if (lgood(kpl(i+1))) then
                        ytt = ys * (plt(i+1) - ymid)
                        if (dabs(ytt) .lt. dabs(ymax-ymin)) then
                           ttt = xs * (tx(i+1) - tx(ii0))
                           ix2 = int(ix0 + ttt)
                           iy2 = int(iy0 - ytt)
                           call gline (ix1,iy1,ix2,iy2,igerr)
                        endif
                     endif
                     if (lbias(kpl(i))) then
                        ytt = ys * (plt(i) - ymid)
                        if (dabs(ytt) .lt. dabs(ymax-ymin)) then
                           ttt = xs * (tx(i) - tx(ii0))
                           ix1 = int(ix0 + ttt)
                           ix2 = int(ix0 + ttt)
                           iy1 = int(iy0 - ytt - 8*iradi)
                           iy2 = int(iy0 - ytt + 8*iradi)
                           call gline (ix1,iy1,ix2,iy2,igerr)
                        endif
                     endif
                  else if (lloel(kpl(i))) then
                     call grect(ix1,iy1,iradi,igerr)
                  else if (lmarg(kpl(i))) then
                     call gcirc(ix1,iy1,iradi,igerr)
                  endif
               endif
            endif
 2070    continue

c        symbol for for final point
         if (lgood(kpl(ilast)).or.
     .           (domarg.and.lmarg(kpl(ilast)))) then
            if (dabs(ytt) .lt. dabs(ymax-ymin)) then
               ytt = plt(ilast) - ymid
               ttt = tx(ilast)-tx(ii0)
               ix1 = ix0 + int(xs*ttt)
               iy1 = iy0 - int(ys*ytt)
               if (lgood(kpl(ilast))) then
                  call gdot(ix1,iy1,iradi,igerr)
               else if (lloel(kpl(ilast))) then
                  call grect(ix1,iy1,iradi,igerr)
               else if (lmarg(kpl(ilast))) then
                  call gcirc(ix1,iy1,iradi,igerr)
               endif

c              vertical line for bias flag
               if (lbias(kpl(ilast))) then
                  ytt = ys * (plt(ilast) - ymid)
                  ttt = xs * (tx(ilast) - tx(ii0))
                  ix1 = int(ix0 + ttt)
                  ix2 = int(ix0 + ttt)
                  iy1 = int(iy0 - ytt - 8*iradi)
                  iy2 = int(iy0 - ytt + 8*iradi)
                  call gline (ix1,iy1,ix2,iy2,igerr)
               endif
            endif
         endif
      endif

      return
      end
c**********************************************************************
      subroutine data_yaxis (iwin,myaxi,ymin,ymax,yssd,ymid,
     .ix0,iy0,xs,ys,csize,yaxlab,yaxuni,ipoly,lstack,igraph,
     .alabel)

c     draw the Y-axis and label it

      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'

      integer        igraph
      character*2 yaxlab
      character*6 yaxuni
      real*8 ymin,ymax,ymid,xs,ys,yssd
      character*80 amsg
      logical lstack
      character*12 alabel(6)

      integer iwin,i
      integer*4 igerr,ipoly
      integer*2 csize(2),myaxi(maxwin,4),mrect(4)
      integer*2 ix1,iy1,ix2,iy2,ix0,iy0

c     set clipping in the window and clean it
      do 1082 i= 1,4
         mrect(i) = myaxi(iwin,i)
 1082 continue
      call gclip (mrect, .true.)
      call gmsg  (mrect,' ')

      mrect(3) = int(32* csize(1))
      mrect(4) = int(1.1 * csize(2))

      mrect(1) = myaxi(iwin,1)
      mrect(2) = myaxi(iwin,2) + 1*myaxi(iwin,4)/16
      write (amsg,2045) alabel(1),alabel(3)
      call gmsg (mrect,amsg)

      mrect(1) = myaxi(iwin,1)
      mrect(2) = myaxi(iwin,2) + 5*myaxi(iwin,4)/16
      write (amsg,2045) alabel(1),alabel(4)
      call gmsg (mrect,amsg)

      mrect(1) = myaxi(iwin,1)
      mrect(2) = myaxi(iwin,2) + 9*myaxi(iwin,4)/16
      write (amsg,2045) alabel(2),alabel(3)
      call gmsg (mrect,amsg)

      mrect(1) = myaxi(iwin,1)
      mrect(2) = myaxi(iwin,2) + 13*myaxi(iwin,4)/16
      write (amsg,2045) alabel(2),alabel(4)
      call gmsg (mrect,amsg)

 2045 format (a8,1x,a8)

c     draw the y-axis
      ix1 = myaxi(iwin,1) + myaxi(iwin,3) - 1
      ix2 = myaxi(iwin,1) + myaxi(iwin,3) - 1
      iy1 = myaxi(iwin,2)
      iy2 = myaxi(iwin,2) + myaxi(iwin,4)
      call gline (ix1,iy1,ix2,iy2,igerr)

      return
      end

c**********************************************************************
      subroutine editor_plotdt (ii0,ii1,jj0,jj1,ix0,iy0,xs,ys,
     .ymid,ymin,ymax,iradi,plt,kpl,domarg,hide_b,lclipt,
     .isite1,isite2,ichan1,ichan2)

c     modified to plot vertical bar on last point kurt June 91
c
c     actually plot the time series
c
c     standard GAMIT array dimensions
      include '../includes/dimpar.h'

c     some non GAMIT declarations
      include '../includes/makex.h'

c     common block declarations
      include '../includes/cview.h'

c     standard GAMIT error flags
      include '../includes/errflg.h'
c     FUNCTIONS:
c     functions to sort out error codes
c        returns .true. for observations good enough to be used in a solution
         logical lgood
c        returns .true. for marginal observations which might potentially useful
         logical lmarg
c        returns .true. for observations with a bias flag
         logical lbias
c
c     standard error flag decoders
      character*4 errcod

      logical  domarg,hide_b

c     Has clipping occurred?
      logical lclipt

c%%   values which must be integer*2 for the GSUBS package
      integer*2 iradi,ix1,iy1,iy11,iy12,ix0,iy0

c%%   this system status variable may need to be integer*2 on some systems
      integer*4 igerr

c     indexes into the time series
      integer  jj0,jj1

c     actual start and end
      integer ii0,ii1

c     scales and such
      real*8 xs,ys,ymid,ymin,ymax,ytt,ttt,yyy

      integer*4 kpl(maxepc)
      real*8 plt(maxepc)

      integer*4 isite1,isite2,ichan1,ichan2,ksite,kchan,icom,i

      yyy = 1.25*ymax

      do 800 icom = 1,4

          ksite = isite1
          kchan = ichan1
          if((icom .eq. 3) .or. (icom .eq. 4)) kchan = ichan2
          if((icom .eq. 2) .or. (icom .eq. 4)) ksite = isite2

          yyy = yyy - 0.5
          ytt = ys * yyy * 1.1
          iy1 = int(iy0 - ytt)
          iy11 = int(iy0 - ytt - 4*iradi)
          iy12 = int(iy0 - ytt + 4*iradi)

          if((ksite .ne. 0) .and. (kchan .ne. 0)) then

               do 400 i = ii0,ii1


                  if (lgood(ierr(i,kchan,ksite))) then
                      ttt = xs * (tx(i) - tx(ii0))
                      ix1 = int(ix0 + ttt)
                      call gdot(ix1,iy1,iradi,igerr)
                      if (lbias(ierr(i,kchan,ksite)))
     .                      call gline (ix1,iy11,ix1,iy12,igerr)
                  else if (domarg.and.lmarg(ierr(i,kchan,ksite))
     .                 .and.errcod(ierr(i,kchan,ksite)).ne.'loel') then
                      ttt = xs * (tx(i) - tx(ii0))
                      ix1 = int(ix0 + ttt)
                      call gcirc(ix1,iy1,iradi,igerr)
                  else if (domarg.and.lmarg(ierr(i,kchan,ksite))
     .                 .and.errcod(ierr(i,kchan,ksite)).eq.'loel') then
                      ttt = xs * (tx(i) - tx(ii0))
                      ix1 = int(ix0 + ttt)
                      call grect(ix1,iy1,iradi,igerr)
                  endif

  400          continue

          endif

  800 continue

      return
      end
c**********************************************************************
