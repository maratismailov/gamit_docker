c**********************************************************************
      subroutine editor_windos
     .    (mtext,mwarn,mhelp,menuu,menul,menus,mclip,myaxi,mplot,mtaxi,
     .     ixbm0,iybm0,nx,ny,nxcurs,csize,nwins)

c     set up windows for editor

c     standard GAMIT array dimensions
      include '../includes/dimpar.h'

c     some non GAMIT declarations
      include '../includes/makex.h'

      include '../includes/cview.h'

c     plotting windows
c     in the apollo convention a window is defined by the elements of array w
c     w(1) x coordinate of u.l. corner
c     w(2) y coordinate of u.l. corner (down from u.l. corner of screen)
c     w(3) size of window in x dimension
c     w(4) size of window in y dimension
c     thus the corners are
c     w(1),w(2)          w(1)+w(3),w(2)
c     w(1),w(2)+w(4)     w(1)+w(3),w(2)+w(4)
      integer    nwins,i

      integer*2
C     clipping window
     .          mclip(4)
C     time series only windows
     .,         mplot(maxwin,4)
C     station record window
     .,         mtext(4)
C     time scale window
     .,         mtaxi(4)
C     y-scale window
     .,         myaxi(maxwin,4)
C     warning window
     .,         mwarn(4)
C     help window
     .,         mhelp(4)
C     upper menu panel nmenus/2 + 1 - nmenus
     .,         menuu(4)
C     lower menu panel 1-nmenus/2
     .,         menul(4)
C     choices
     .,         menus(nmenus,4)

c%%   values which must be integer*2 for the Apollo GPR package
      integer*2 nx,ny,nxcurs,ixbm0,iybm0,csize(2)

c     widths
      mtext(3) = nx
      mwarn(3) = 3*nx/10
      mhelp(3) = mwarn(3)
      menuu(3) = nx - mwarn(3)
      menul(3) = nx - mhelp(3)
      mclip(3) = nx - nxcurs
      do 310 i=1,nwins
         myaxi(i,3) = 26 * csize(1)
         mplot(i,3) = mclip(3) - myaxi(1,3)
 310  continue
      mtaxi(3) = mplot(1,3)

c     heights
      mtext(4) = 2 * csize(2)
      mtaxi(4) = 3 * csize(2)
      menul(4) = csize(2) + int(12)
      menuu(4) = csize(2) + int(12)
      mclip(4)  = ny - menul(4) - menuu(4) - 1
      do 320 i=1,nwins
         mplot(i,4) = (mclip(4) - mtext(4) - mtaxi(4))/nwins
         myaxi(i,4) = (mclip(4) - mtext(4) - mtaxi(4))/nwins
 320  continue
      mwarn(4) = menuu(4)
      mhelp(4) = menul(4)

c     x positions
      menuu(1) = ixbm0
      menul(1) = ixbm0
      mtext(1) = ixbm0
      mclip(1) = ixbm0
      do 330 i=1,nwins
         myaxi(i,1) = ixbm0
         mplot(i,1) = myaxi(i,1) + myaxi(i,3)
 330  continue
      mtaxi(1) = myaxi(1,1) + myaxi(1,3)
      mwarn(1) = nx - mwarn(3)
      mhelp(1) = nx - mhelp(3)

c     y positions
      menuu(2) = iybm0
      mtext(2) = menuu(2) + menuu(4)
      mclip(2) = mtext(2)
      mtaxi(2) = mclip(2) + mclip(4) - mtaxi(4)
      myaxi(1,2) = mtext(2) + mtext(4)
      mplot(1,2) = mtext(2) + mtext(4)
      do 340 i=2,nwins
         myaxi(i,2) = myaxi(i-1,2) + myaxi(i-1,4)
         mplot(i,2) = mplot(i-1,2) + mplot(i-1,4)
 340  continue
      menul(2) = mtaxi(2) + mtaxi(4)
      mwarn(2) = menuu(2)
      mhelp(2) = menul(2)

CD     print 360,'mtext',0,(mtext(j),j=1,4)
CD     print 360,'menuu',0,(menuu(j),j=1,4)
CD     print 360,'menul',0,(menul(j),j=1,4)
CD     print 360,'Mclip',0,(mclip(i),i=1,4)
CD     print 360,'mtaxi',0,(mtaxi(j),j=1,4)
CD     print 360,'myaxi',0,(myaxi(1,j),j=1,4)
CD     print 360,'mwarn',0,(mwarn(i),i=1,4)
CD     do 350 i = 1,nwins
CD        print 360,'mplot',i,(mplot(i,j),j=1,4)
CD350  continue
 360  format (a,1x,5(i5,1x))

      return
      end

c**********************************************************************
      subroutine editor_menus
     .(amenus,menus,menul,menuu,mwarn,mhelp,
     .alabel,imenu,imenub,iseek,obcode,
     .l12v,ipoly,yadd,amonth,pfile,title,nextslip,igraph,
     .newbox,domarg,hide_b,lstack,lmovie,lspan,allpts,csize)

c     set up the menus for the editor


c     standard GAMIT array dimensions
      include '../includes/dimpar.h'

c     some non GAMIT declarations
      include '../includes/makex.h'

      include '../includes/cview.h'

      integer imenu,imenub,iseek

c     file name to dump hardcopy
      character*32 pfile

c     month codes
      character*3 amonth(12)

      character*2 obcode(-3:ntypes),char2
      character*96 amsg,title
      character*16 amenus(nmenus)
c     labels sat1,sat2,site1,site2,poly,obs
      character*12 alabel(6)

      integer*2 size(2),menus(nmenus,4),mrect(4),csize(2),mwarn(4)
     .,mhelp(4),menul(4),menuu(4)
      integer*4 igerr

      logical domarg,newbox,lstack,lmovie,lspan,allpts,nextslip,hide_b
      real*8 yadd

C     plotting polynomial
      integer   ipoly

c     types of graphs
c     igraph = 1     time series
c     igraph = 2     sky map
c     igraph = 3     spectral estimate
      integer        igraph


C     data type on the screen
      integer    l12v(maxwin)
      integer   i,j

c     title for plot
      write (title,155)
     .      iit0(2),amonth(iit0(1)),mod(iit0(3),100),
     .      alabel(1),alabel(2),alabel(3),alabel(4),
     .      alabel(5)
 155  format (i2,1x,a3,1x,i2,1x,a8,1x,a8,1x,a7,1x,a7,1x,a7)

c     file name of file in which to dump time series for later plotting
c     name of form 26MAY89.SIT1.SIT2.01.02.00.LC
c                                          ^^
c                                          P1 for Poly 1
c                                          T1 for Deriv 1 + Poly
c                                          SM for Smooth
      if (ipoly .eq. -1) then
         char2 = 'D1'
      else if (ipoly .eq. -2) then
         char2 = 'SM'
      else if (ipoly .gt. 9) then
         write(char2,'(''T'',i1.1)') mod(ipoly,10)
      else
         write(char2,'(''P'',i1.1)') ipoly
      endif
      write (pfile,157)
     .      iit0(2),amonth(iit0(1)),mod(iit0(3),100),
     .      alabel(1)(7:8),alabel(2)(7:8),
     .      alabel(3)(4:7),alabel(4)(4:7),char2
 157  format(i2.2,a3,i2.2,'.',a2,'.',a2,'.',a4,'.',a4,'.',a2)

c     the individual menu choices
c     menu level -1 is the initial, top level menu
      if (imenu .eq. -1) then
         amenus(1) = 'SEEK'
         amenus(2) = 'PLOT'
         amenus(3) = 'LIST'
         amenus(4) = '1-WAY'
         amenus(5) = alabel(1)
         amenus(6) = alabel(2)
         amenus(7) = alabel(3)
         amenus(8) = alabel(4)
         amenus(9) = 'FILE'

         if (iseek .eq. 1) then
            amenus(10) = 'SEEK-FORW'
         else if (iseek .eq. 2) then
            amenus(10) = 'SEEK-BACK'
         else if (iseek .eq. 3) then
            amenus(10) = 'SEEK-LIST'
         else
            amenus(10) = '  ERROR!'
         endif
        
c        Comment this out to avoid generating a compiler warning.
c        Statement added to cview.h in includes to reinstate it
c        if the number of menus exceeds 10.
c         if( nmenus/2 .gt.11 ) then
c           do i = 11,nmenus/2
c              amenus(i) = ' '
c           enddo
c         endif

c        observable code: L1, LC etc.
         do 133 i = 1, maxwin
            if (l12v(i) .lt. -2)     l12v(i) = -2
            if (l12v(i) .gt. ntypes) l12v(i) = ntypes
            write (amenus(nmenus/2+i),'(a2)') obcode(l12v(i))
 133     continue

c        polynomial
         amenus(nmenus/2+maxwin+1) = alabel(5)

c        type of graph
         if (igraph .eq. 1) then
            amenus(nmenus/2+maxwin+2) = 'TIME SER'
         else if (igraph .eq. 2) then
            amenus(nmenus/2+maxwin+2) = ' SKY MAP'
         else if (igraph .eq. 3) then
            amenus(nmenus/2+maxwin+2) = 'SPECTRUM'
         else if (igraph .eq. 4) then
            amenus(nmenus/2+maxwin+2) = ' ALLAN  '
         else
            amenus(nmenus/2+maxwin+2) = '  ERROR!'
         endif

c        stack
         if (lstack) then
            amenus(nmenus/2+maxwin+3) = 'STACK [*]'
         else
            amenus(nmenus/2+maxwin+3) = 'STACK [ ]'
         endif

c        movie
         if (lmovie) then
            amenus(nmenus/2+maxwin+4) = 'MOVIE [*]'
         else
            amenus(nmenus/2+maxwin+4) = 'MOVIE [ ]'
         endif

         amenus(nmenus) = 'STOP'

c     menu level to manipulate the plot
      else if (imenu .eq. 0) then
         if (igraph .eq. 1) then
           amenus(1) = 'SAVE'
         else
           amenus(1) = ' '
         endif

         amenus(2) = '---ABORT---'

         if (lspan) then
            amenus(3) = 'SPAN[*]'
         else
            amenus(3) = 'SPAN[ ]'
         endif
         if (allpts) then
            amenus(4) = 'ALL[*]'
         else
            amenus(4) = 'ALL[ ]'
         endif
         if (igraph .eq. 1 .or. igraph .eq. 2) then
            amenus(5) = '<< T >>'
            amenus(6) = '>> T <<'
         else
            amenus(5) = ' '
            amenus(6) = ' '
         endif
         if (domarg) then
            amenus(7) = 'MARG[*]'
         else
            amenus(7) = 'MARG[ ]'
         endif
         if (hide_b) then
            amenus(8) = 'HIDE[ ]'
         else
            amenus(8) = 'HIDE[*]'
         endif
         if (igraph .eq. 1 .or. igraph .eq. 2) then
            amenus(9) = '<<'
            amenus(10) = '>>'
         else
            amenus(9) = ' '
            amenus(10) = ' '
         endif
         if (igraph .ne. 1) then
            amenus(11)= ' '
            amenus(12)= ' '
            amenus(13)= ' '
            amenus(14)= ' '
            amenus(15)= ' '
            amenus(16)= ' '
            amenus(17)= ' '
            amenus(18)= ' '
            amenus(19)= ' '
            amenus(20)= ' '
         else
            if (nextslip) then
                amenus(11) = 'SLIP [*]'
            else
                amenus(11) = 'SLIP [ ]'
            endif
            amenus(12)= 'MOVE'
            amenus(13)= 'PATCH'
            amenus(14)= 'FIND'
            amenus(15)= 'SAVE'
            amenus(16)= 'UNWT'
            amenus(17)= 'REWT'
            amenus(18)= 'BIAS -/?/+'
            amenus(19)= 'UNDO'
            amenus(20)= 'ELIM'
         endif
c     menu level to do a manual move
      else if (imenu .eq. 12) then
         amenus(11) = '0.25'
         amenus(12) = '0.5'
         amenus(13) = '1'
         amenus(14) = '2'
         amenus(15) = '5'
         amenus(16) = '10'
         amenus(17) = '50'
         amenus(18) = '100'
         amenus(19) = '500'
         amenus(20) = '1000'
         amenus(1)= 'CANCEL (M)'
         amenus(2)= 'PERFORM (R)'
         do 167 i = 3,10
            amenus(i) = ' '
 167     continue
      endif

CD     do 141 i=1,nmenus
CD        print *,amenus(i)
CD141  continue

      if (newbox) then
c        make both menu panels white
         call gbox (menul,igerr)
         call gbox (menuu,igerr)
         if (imenu .ne. imenub) then
ckf         if (imenu .eq. 0) then
ckf            call gmsg (mwarn,'Click mouse on menu or data.')
ckf            else if (imenu .eq. -1) then
            if (imenu .eq. -1) then
               call gmsg (mwarn,
     .           'Mouse: 1) down 2) default 3) up.')
ckf         else if (imenu .eq. 13) then
ckf            call gmsg (mwarn,'Click on cycle slip in LG.')
            else if (imenu .eq. 12) then
c              show the user the move amount
               write (amsg,4265) yadd
 4265          format (' Add ',f20.2)
               call gmsg (mwarn, amsg)
            endif
         endif
         newbox = .false.
      endif


c     make boxes for each command
      do 410 i=1,nmenus
c        how big is the text?
         call gtsize (amenus(i),size,igerr)

c        width
         if (size(1) .gt. 0) then
            menus(i,3) = int(size(1)+6)
         else
            menus(i,3) = int(size(1))
         endif

c        height
         if (size(2) .gt. 0) then
            menus(i,4) = int(size(2)+6)
         else
            menus(i,4) = int(size(2))
         endif

c        x position
         if (i.eq.1) then
            menus(i,1) = menul(1) + csize(1)
         else if (i.eq.nmenus/2 + 1) then
            menus(i,1) = menuu(1) + csize(1)
         else if (i.le.nmenus/2) then
            menus(i,1) = menul(1) + csize(1)
            do 420 j=1,i-1
               menus(i,1) = menus(i,1)+menus(j,3)+2*csize(1)
 420        continue
         else if (i.le.nmenus) then
            menus(i,1) = menul(1) + csize(1)
            do 430 j=nmenus/2+1,i-1
               menus(i,1) = menus(i,1)+menus(j,3)+2*csize(1)
 430        continue
         endif

c        y position
         if (i.le.nmenus/2) then
            menus(i,2) = menul(2)+int((menul(4)-menus(i,4))/2)
         else
            menus(i,2) = menuu(2)+int((menuu(4)-menus(i,4))/2)
         endif

c        write the command in the box
         do 405 j = 1,4
            mrect(j) = menus(i,j)
 405     continue
         call gmsg (mrect,amenus(i))
CD        print 12,'menus',i,(menus(i,j),j=1,4)
 12      format (1x,a5,5(1x,i5))
 410  continue


      return
      end

c**********************************************************************
      subroutine editor_labels
     .   (ichan1,ichan2,isite1,isite2,ipoly,
     .    nsat,ncfls,asites,alabel)

c     set up labels for buttons and plot

c     standard GAMIT array dimensions
      include '../includes/dimpar.h'

c     some non GAMIT declarations
      include '../includes/makex.h'

      include '../includes/cview.h'

      integer ichan1,ichan2,isite1,isite2,ipoly,ncfls,nsat
      character*4 asites(maxsit)
      character*12 alabel(6)

c     sat 1
 113  format (i2,1x,a3,i2.2)
      if (ichan1 .gt. 0 .and. ichan1 .le. nsat) then
         write (alabel(1),113) ichan1,'PRN',isprn(ichan1)
      else
         write (alabel(1),113) ichan1,'...',0
      endif

c     sat 2
      if (ichan2 .gt. 0 .and. ichan2 .le. nsat) then
         write (alabel(2),113) ichan2,'PRN',isprn(ichan2)
      else
         write (alabel(2),113) ichan2,'...',0
      endif

c     site 1
 120  format (i2,1x,a4)
      if (isite1 .gt. 0 .and. isite1 .le. ncfls) then
         write (alabel(3),120) isite1,asites(isite1)
      else
         write (alabel(3),120) isite1,'....'
      endif
c     upper case
      call uppers(alabel(3))

c     site 2
      if (isite2 .gt. 0 .and. isite2 .le. ncfls) then
         write (alabel(4),120) isite2,asites(isite2)
      else
         write (alabel(4),120) isite2,'....'
      endif
c     upper case
      call uppers(alabel(4))

c     polynomial
      if (ipoly .eq. -1) then
         write (alabel(5),'(a7)') 'DERIV 1'
      else if (ipoly .eq. -2) then
         write (alabel(5),'(a7)') ' SMOOTH'
      else if (ipoly .eq. 0) then
         write (alabel(5),'(a7)') 'NO POLY'
      else if (ipoly .ge. 1 .and. ipoly .le.9) then
         write (alabel(5),'(a4,2x,i1.1)') 'POLY',ipoly
      else if (ipoly .ge. 10 .and. ipoly .le. maxply) then
         write (alabel(5),'(a4,2x,i1.1)') 'D1+P',mod(ipoly,10)
      else
         write (alabel(5),'(a7)') 'ERROR ?'
      endif

      return
      end


c**********************************************************************
      subroutine editor_plotsm
     .(ii0,ii1,ix0sm,iy0sm,rs,cs,ymid,dt,plt,kpl,domarg)
c
c     plot residuals on skymap
c
c     plus signs for positive residuals
c     open circles for negative residuals
c
c     if ys is zero, then we just get a sky track
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

      logical  domarg

c%%   values which must be integer*2 for the GSUBS package.
      integer*2 ix,iy,ix0sm,iy0sm,iradib

c%%   this system status variable may need to be integer*2 on some systems
      integer*4 igerr

c     indexes into the time series
      integer  ii0,ii1

c     scales and such
      real*8 cs,ymid,rs,rho,pihalf,dt
      integer i

c     what we are plotting
      real*8  plt(maxepc)
      integer*4 kpl(maxepc)

      pihalf = 2.0d0 * datan(1.0d0)

c     draw a circle around it
      iradib = int(rs)
      call gcirc(ix0sm,iy0sm,iradib,igerr)

c     this is an equal-area projection to space the symbols better
      do 1070 i = ii0,ii1
         if (lgood(kpl(i)).or.(domarg.and.lmarg(kpl(i)))) then
            rho = rs * dsin((pihalf - elw(i))/2.0d0)
            ix = ix0sm + int(rho*dsin(azw(i)))
            iy = iy0sm - int(rho*dcos(azw(i)))
            iradib = int(dabs(cs*(plt(i)-ymid)))
            if (plt(i) .gt. ymid) then
               call gcirc(ix,iy, iradib, igerr)
            else
               call gplus(ix,iy, iradib, igerr)
            endif
         endif
 1070 continue


      return
      end
c**********************************************************************
      subroutine editor_plotal
     .(ii0,ii1,ix0,iy0,xs,ys,ymin,ymid,ymax,plt)
c
c     plot the Allan standard deviation as a function of sample period
c
c     standard GAMIT array dimensions
      include '../includes/dimpar.h'

c     some non GAMIT declarations
      include '../includes/makex.h'

c     common block declarations
      include '../includes/cview.h'

c     standard GAMIT error flags
      include '../includes/errflg.h'

c     what we are plotting
      real*8 plt(maxepc)


c%%   values which must be integer*2 for the GSUBS package
      integer*2 ix1,iy1,ix2,iy2,ix0,iy0

c%%   this system status variable may need to be integer*2 on some systems
      integer*4 igerr

c     indexes into the time series
      integer  ii0,ii1

c     scales and such
      real*8 xs,ys,ymid,ymin,ymax,ytt,ttt
      integer i

      do 1070 i = ii0,ii1-1
         ytt = ys * (plt(i) - ymid)
         ttt = xs * (tx(i) - tx(ii0))
         ix1 = int(ix0 + ttt)
         iy1 = int(iy0 - ytt)
         ytt = ys * (plt(i+1) - ymid)
         ttt = xs * (tx(i+1) - tx(ii0))
         ix2 = int(ix0 + ttt)
         iy2 = int(iy0 - ytt)
         call gline (ix1,iy1,ix2,iy2,igerr)
 1070 continue

      return
      end

c**********************************************************************
      subroutine editor_plotps
     .(ii0,ii1,ix0,iy0,xs,ys,ymid,ymin,ymax,iradi,plt,kpl,domarg,lclipt)
c
c     plot the power spectrum
c
c     standard GAMIT array dimensions
      include '../includes/dimpar.h'

c     some non GAMIT declarations
      include '../includes/makex.h'

c     common block declarations
      include '../includes/cview.h'

c     standard GAMIT error flags
      include '../includes/errflg.h'

c     Has clipping occurred?
      logical lclipt

c     plot marginal points
      logical domarg

c%%   values which must be integer*2 for the GSUBS package
      integer*2 iradi,ix1,iy1,ix2,iy2,ix0,iy0

c%%   this system status variable may need to be integer*2 on some systems
      integer*4 igerr

c     indexes into the time series
      integer  ii0,ii1

c     scales and such
      real*8 xs,ys,ymid,ymin,ymax,ytt,ttt
      integer i

c     what we are plotting
      real*8  plt(maxepc)
      integer*4 kpl(maxepc)


      do 1070 i = ii0,ii1-1
         ytt = ys * (plt(i) - ymid)
         ttt = xs * (tx(i) - tx(ii0))
         ix1 = int(ix0 + ttt)
         iy1 = int(iy0 - ytt)
         ytt = ys * (plt(i+1) - ymid)
         ttt = xs * (tx(i+1) - tx(ii0))
         ix2 = int(ix0 + ttt)
         iy2 = int(iy0 - ytt)
         call gline (ix1,iy1,ix2,iy2,igerr)
 1070 continue

      return
      end
c**********************************************************************
* MOD TAH 031223: Added yavg being passed into routine so that it
*     can be output on the window above the s = line
      subroutine editor_yaxis (iwin,myaxi,ymin,ymax,yavg,yssd,ymid,
     .ix0,iy0,xs,ys,csize,yaxlab,yaxuni,ipoly,lstack,igraph)

c     draw the Y-axis and label it

      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'

c     types of graphs
c     igraph = 1     time series
c     igraph = 2     sky map
c     igraph = 3     spectral estimate
      integer        igraph

      character*2 yaxlab
      character*6 yaxuni
      real*8 ymin,ymax,ymid,xs,ys,ynice1(2),ynice2(2),ytick,ytemp
     .,      yavg, yssd,xtemp,y0
      character*80 amsg
      logical lstack

      integer iwin,i,icount
      integer*4 igerr,ipoly
      integer*2 csize(2),myaxi(maxwin,4),mrect(4)
      integer*2 ix1,iy1,ix2,iy2,ix0,iy0

c     set clipping in the window and clean it
      do 1082 i= 1,4
         mrect(i) = myaxi(iwin,i)
 1082 continue
      call gclip (mrect, .true.)
      call gmsg  (mrect,' ')

c     find a nice division
c     Keep in mind: that for stacked plots, ymid floats with
c     the data, but ymin and ymax are fixed.

      if (igraph .eq. 3 .or. igraph .eq. 4) then
         ynice1(1) = ymin
         ynice1(2) = ymax
         y0 = ymid
      else
         ynice1(1) = ymin - ymid
         ynice1(2) = ymax - ymid
         y0 = 0.0d0
      endif

      call nicer (ynice1,ynice2,ytick)
 1083 continue
      if (ytick*ys .lt. 1.3*csize(2)) then
         ytick = 2.0d0 * ytick
         goto 1083
      endif

c     draw the y-axis
      ix1 = myaxi(iwin,1) + myaxi(iwin,3) - 1
      ix2 = myaxi(iwin,1) + myaxi(iwin,3) - 1
      iy1 = myaxi(iwin,2)
      iy2 = myaxi(iwin,2) + myaxi(iwin,4)
      call gline (ix1,iy1,ix2,iy2,igerr)

c     do tick marks on Y-axis, labelled with time in left half of myaxi
c     Modifications by pjm, nov 2008 as non-integer do loops no longer 
c     supportrd.

c      do 1150 ytemp = ynice2(1),ynice2(2),ytick
      icount = 1
      do
         ytemp = ynice2(1) + (icount-1)*ytick
         icount = icount + 1
         if(ytemp .gt.ynice2(2)) exit
         xtemp  = myaxi(iwin,1)+myaxi(iwin,3)-1
         ix1 = xtemp - csize(1)
         ix2 = xtemp
         iy1 = iy0 - ys * (ytemp - y0)
         iy2 = iy0 - ys * (ytemp - y0)
         call gline (ix1,iy1,ix2,iy2,igerr)
         mrect(3) = 9 * csize(1)
         mrect(4) = csize(2)
         mrect(1) = myaxi(iwin,1) + myaxi(iwin,3)-mrect(3)-2*csize(1)-1
         mrect(2) = int(iy0 - ys * (ytemp-y0))
         if (dabs(ytemp) .le. 1.0d-5) then
            write (amsg,'(i8)') 0
         elseif (dabs(ytemp) .lt. 0.1d0) then
            write (amsg,'(1pe8.1)') ytemp
         elseif (dabs(ytemp) .lt. 1.0d3) then
            write (amsg,'(f8.1)') ytemp
         else
            write (amsg,'(1pe8.1)') ytemp
         endif
         call gmsg (mrect,amsg)
c 1150 continue
      end do
   

c     label at left edge of Y scale telling us observable
c     units are per second if this is a derivative
c     units are log power if this is power spectrum
c     units are dimensionless if this is allan
      if (igraph .eq. 3) then
         write (amsg,202) yaxlab
  202    format (a2,' ','log10 power')
      else if (igraph .eq. 4) then
         write (amsg,203) yaxlab
  203    format (a2,' ','log10')
      else if (ipoly .eq. -1 .or. ipoly .ge. 10) then
         write (amsg,201) yaxlab,yaxuni
  201    format (a2,' ',a6,'/s')
      else
         write (amsg,200) yaxlab,yaxuni
  200    format (a2,' (',a6,')')
      endif
      mrect(3) = int(14* csize(1))
      mrect(4) = int(1.1 * csize(2))
      mrect(1) = myaxi(iwin,1)
      mrect(2) = myaxi(iwin,2) + myaxi(iwin,4)/4
      call gmsg (mrect,amsg)

c     label at left edge of Y scale telling us midpoint
      if (igraph .eq. 1 .or. igraph .eq. 2) then
         if (dabs(ymid) .lt. 1.0d3 .and. dabs(ymid) .gt. 0.1) then
            write (amsg,2043) ymid
 2043       format ('0@ ',f9.2)
         else
            write (amsg,2045) ymid
 2045       format ('0@ ',1pe9.2)
         endif
         mrect(3) = int(16 * csize(1))
         mrect(4) = int(1.1 * csize(2))
         mrect(1) = myaxi(iwin,1)
         mrect(2) = myaxi(iwin,2) + myaxi(iwin,4)/2
         call gmsg (mrect,amsg)
      endif

c     label at Y scale with sample standard deviation
      if (igraph .eq. 1 .or. igraph .eq. 2) then
* MOD TAH 031223: Added output of the average value as well
         if (abs(yavg) .lt. 1.0d3) then
            write (amsg,1040) yavg
 1040       format ('Av= ',f8.2)
         else
            write (amsg,1042) yavg
 1042       format ('Av= ',1pe8.2)
         endif
         mrect(3) = int(14 * csize(1))
         mrect(4) = int(1.1 * csize(2))
         mrect(1) = myaxi(iwin,1)
         mrect(2) = myaxi(iwin,2) + 5*myaxi(iwin,4)/8
         call gmsg (mrect,amsg)

         if (yssd .lt. 1.0d3) then
            write (amsg,1043) yssd
 1043       format ('s = ',f8.2)
         else
            write (amsg,1045) yssd
 1045       format ('s = ',1pe8.2)
         endif
         mrect(3) = int(14 * csize(1))
         mrect(4) = int(1.1 * csize(2))
         mrect(1) = myaxi(iwin,1)
         mrect(2) = myaxi(iwin,2) + 3*myaxi(iwin,4)/4
         call gmsg (mrect,amsg)
      endif

      return
      end

c**********************************************************************
      subroutine editor_xaxis(mtaxi,ii0,ii1,ix0,xs,csize,igraph)

c     Draw a horizontal X axis
c     for time series and skymaps, scale is in HH:MM
c     for power spectra, it is in Hz
c     for allan variance, it is in seconds

c     Plot from goes from index ii0 to ii1

      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'

c     types of graphs
c     igraph = 1     time series
c     igraph = 2     sky map
c     igraph = 3     spectral estimate
      integer        igraph

      real*8 xs,xnice1(2),xnice2(2),xtick,ttt,tt
      integer*2 mtaxi(4),ix1,iy1,ix2,iy2,csize(2),mrect(4),ix0
      integer ii0,ii1
      integer ihr,min
      integer*4 igerr, icount
      character*80 amsg

c     set the clipping window
      call gclip (mtaxi, .true.)

c     clean out the window
      call gmsg (mtaxi, ' ')

c     draw the t-axis
      ix1 = mtaxi(1)
      ix2 = mtaxi(1) + mtaxi(3)
      iy1 = mtaxi(2) + mtaxi(4)/2
      iy2 = mtaxi(2) + mtaxi(4)/2
      call gline (ix1,iy1,ix2,iy2,igerr)

      if (igraph .eq. 3) then
c        full spectrum
         xnice1(1) =  tx(1)
         xnice1(2) =  tx(maxepc)
      else if (igraph .eq. 4) then
c        last good period
         xnice1(1) =  tx(1)
         xnice1(2) =  tx(maxepc/8)
      else
c        start and end time in minutes
         xnice1(1)=  60.0d0*tt00(1) + tt00(2) + (tt00(3)+tx(ii0))/60.0d0
         xnice1(2)=  60.0d0*tt00(1) + tt00(2) + (tt00(3)+tx(ii1))/60.0d0
      endif

c     find a nice division
      call nicer (xnice1,xnice2,xtick)
c     xtick is returned with zero on disaster
      if (xtick .le. 0.0d0) then
         call gmsg (mtaxi,'NICER failed.')
         return
      endif

      if (xnice2(1) .lt. xnice1(1)) xnice2(1) = xnice2(1) + xtick
      if (xnice2(2) .gt. xnice1(2)) xnice2(2) = xnice2(2) - xtick

c     do tick marks on t-axis
c     labels in lower half, ticks go up
c      do 3010 ttt = xnice2(1),xnice2(2),xtick
      icount = 1
      do
         ttt = xnice2(1) + (icount-1)*xtick
         icount = icount + 1
         if(ttt .gt. xnice2(2)) exit
         tt = ttt - xnice1(1)
         if (igraph .eq. 3 .or. igraph .eq. 4) then
c           frequency
            ix1 = ix0 + int(xs*tt)
            ix2 = ix0 + int(xs*tt)
         else
c           time in minutes since first epoch we are plotting
            ix1 = ix0 + int(60*xs*tt)
            ix2 = ix0 + int(60*xs*tt)
         endif

         iy1 = mtaxi(2) + mtaxi(4)/3
         iy2 = mtaxi(2)

         if (ix1 .lt. mtaxi(1)+mtaxi(3) .and.
     .       ix1 .gt. mtaxi(1)) then
             call gline (ix1,iy1,ix2,iy2,igerr)

            if (igraph .eq. 3 .or. igraph .eq. 4) then
c              label scale with frequency or
               mrect(3) = int(11 * csize(1))
               mrect(4) = int(1.1*csize(2))
               mrect(1) = ix1 - mrect(3)/2
               mrect(2) = mtaxi(2) + mtaxi(4) - int(1.1*csize(2))
               if (mrect(1)+mrect(3) .lt. mtaxi(1)+mtaxi(3) .and.
     .             mrect(1)          .gt. mtaxi(1)         ) then
                  if (igraph .eq. 3) then
                     write (amsg,3003) ttt
 3003                format (1pe8.1,'Hz')
                  else
                     write (amsg,3005) ttt
 3005                format (1f8.0,'s')
                  endif
                  call gmsg (mrect,amsg)
               endif
            else
c              label scale with HH:MM
               ihr = int(ttt)/60
c              minutes
               min = nint(ttt - 60.0d0 * ihr)
               if (min .eq. 60) then
                  min = 0
                  ihr = ihr + 1
               endif
               if (ihr .gt. 24) ihr = ihr - 24

               mrect(3) = int(6 * csize(1))
               mrect(4) = int(1.1*csize(2))
               mrect(1) = ix1 - mrect(3)/2
               mrect(2) = mtaxi(2) + mtaxi(4) - int(1.1*csize(2))
               if (mrect(1)+mrect(3) .lt. mtaxi(1)+mtaxi(3) .and.
     .             mrect(1)          .gt. mtaxi(1)         ) then
                  write (amsg,3009) ihr,min
 3009             format (i2.2,':',i2.2)
                  call gmsg (mrect,amsg)
               endif
            endif
         endif
c3010 continue
      end do
c
      return
      end
c**********************************************************************
      subroutine editor_penplt
     .(luout,pfile,shfile,title,obkode,ii0,ii1,kpl,plt,
c              or period and Allan SD
     .ymin,ymax,ymid,yssd,
     .iwin,nwins,mwarn,penplt,igraph,nper)


c     write time series to file and prepare plotting commands
c
c%%   This whole routine is bound to be system dependent.
c%%   It will also depend on the plotting package you choose.

      character*2 obkode
      character*32 pfile
      character*96 shfile,title
      real*8 raddeg,ymin,ymax,ymid,yssd
      integer luout,i
      integer iwin,nwins
      logical penplt
      integer nblen
      integer nper

c     types of graphs
c     igraph = 1     time series
c     igraph = 2     sky map
c     igraph = 3     spectral estimate
      integer        igraph

      integer*2 mwarn(4)
      integer*4 igerr

c     error function
      logical lgood

c     standard GAMIT array dimensions
      include '../includes/dimpar.h'

c     some non GAMIT declarations
      include '../includes/makex.h'

c     common block declarations
      include '../includes/cview.h'

c     standard GAMIT error flags
      include '../includes/errflg.h'

c     indexes into the time series
      integer ii0,ii1
c     what we are plotting
      real*8  plt(maxepc)
      integer*4 kpl(maxepc)

      integer iplot,ipage
      save iplot,ipage
      data iplot /0/
      data ipage /1/


c     radians to degrees
      raddeg = 45.0d0/datan(1.0d0)

      igerr = -2
      open (file    = pfile,
     .      unit    = luout,
     .      status  = 'unknown',
     .      err     = 75,
     .      iostat  = igerr,
     .      form    = 'formatted')
   75 continue

      call ljust (96,title)

      iplot = iplot + 1

      if (igerr .eq. 0) then
C        write (*,*) 'Opened plot file ',pfile
         write (luout,'(a)') obkode//' '//title(1:nblen(title))
         if (igraph .eq. 3) then
            write (luout,100)
  100       format (1x,'FREQ (Hz)    Log10 (Power)')
         else if (igraph .eq. 4) then
            write (luout,110)
  110       format (1x,'period (s)    Log10(period)  Log10 (Allan SD)')
         else
            write (luout,120)   (ymax-ymin)/2.0d0,ymid,yssd
  120       format (1x,'SPAN= ',1PE20.12,
     .              1x,'YMID= ',1PE20.12,
     .              1x,'SDEV= ',1PE20.12)
         endif

         if (igraph .eq. 3) then
c           power spectrum
            do 1079 i = 1,maxepc
c              frequency in Hz and log10 power
               write (luout,1077) tx(i),plt(i)
 1077          format (2(1pe12.4,1x))
 1079       continue
         else if (igraph .eq. 4) then
c           allan variance
            do 1179 i = 1,nper
c              period (s),log10(per),log10(Allan SD)
               write (luout,1177) tx(i),log10(tx(i)),plt(i)
 1177          format (3(1pe12.4,1x))
 1179       continue
         else
            do 1279 i = ii0,ii1
c              elevation is like a latitude
c              azimuth is like a longitude
c              time in minutes
c              residuals with midpoint removed
               if (lgood(kpl(i))) then
                  write (luout,1277) i,plt(i)-ymid ,
     .            elw(i)*raddeg,
     .            azw(i)*raddeg,
     .            (tx(i)-tx(ii0))/60.0d0,
     .            plt(i)-ymid
 1277             format (i5,1pe20.12,2(f11.5,1x),2(1pe20.12,1x))
               endif
 1279       continue
         endif
         close (luout)
      else
         call gmsg (mwarn,'Could not open plot file.')
         print *,'Could not open plot file ',pfile
         call ferror (6,igerr)
      endif

c     write the list of plotting commands to a shell script file
      igerr = -2
      open (file    = shfile,
     .      unit    = luout,
c**  .      status  = 'append',
c**  status = 'append' is not ANSI standard, unfortunately
     .      status  = 'unknown',
     .      access  = 'append',
     .      err     = 1175,
     .      iostat  = igerr,
     .      form    = 'formatted')
 1175 continue

      if (igerr .eq. 0) then
CD        write (*,*) 'Opened shell file ',shfile

         if (iplot .eq. 1) then
            write (luout,1210) ipage,ipage
 1210       format ('if existf myplot',i2.2,'.gmr then dlf myplot'
     .,     i2.2,'.gmr endif')
         endif

         if (igraph .eq. 2) then
            write (luout,1220) pfile(1:nblen(pfile)),ipage
 1220       format ('skymap ',a,' 3 myplot',i2.2)
         else
            write (luout,1250) pfile(1:nblen(pfile)),ipage
 1250       format ('xyplotc ',a,' 3 myplot',i2.2)
         endif

         if (iwin.ge.nwins) then
            write (luout,1260) ipage,iplot
 1260       format ('stackplot myplot',i2.2,1x,i1)
            write (luout,1270) ipage
 1270       format ('penvec myplot',i2.2)
            ipage = ipage+1
            iplot = 0
            penplt = .false.
         endif
         close (luout)
      endif

      return
      end
