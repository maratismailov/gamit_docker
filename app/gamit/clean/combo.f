      subroutine combo
     .   (nobs,ii0,ii1,l12p,domarg,kpl,plt,
     .    ymin,ymax,yavg,yssd,yjump,ystep,
     .    ngood,nmarg,ijump,istep,lcando)
c
c     Form a combined observable from data in common block
c
c     Note: there is no test to see if the observable can be formed.
c     e.g. Will form WL without P2.
c
c     input:
c        nobs    number of points in the vector
c        ii0,ii1 start and end points
c        l12p
c          -2 AZ azimuth from station to sat (CW from North)
c          -1 EL elevation angle in degrees
c           1 L1 phase
c           2 LC phase
c           3 LG phase
c           4 L2 phase
c           5 P1 pseudorange (cycles)
c           6 P2 pseudorange (cycles)
c           7 WL Widelane bias number
c           8 W* Modified Widelane bias number
c           9 PG P2 - gear*P1
c          10 BG LG - PG
c          11 I1 observable I1 = -g/(1-g*g) PG
c          12 I2 observable I2 = -1/(1-g*g) PG
c          13 M1 multipath delay in L1
c          14 M2 multipath delay in L2
c          15 MC multipath delay in LC
c          16 MG multipath delay in LG
c          17 CL Clock offset at station in microseconds
c
c        wl1,wl2 L1,L2 phases
c        pc1,pc2 P1,P2 pseudoranges in cycles
c        kww,kww error flags corresponding to 2 channels
c        lcando        .true. if the observable can be formed
c     output
c        plt         combined observable
c        kpl         error flags for it
c        ijump       epoch at which largest jump occurs
c        istep       epoch at which largest step occurs
c        ymin, ymax  its extrema
c        yavg        mean value of good points
c        yssd        sample variance of points
c        yjump       maximum jump in series
c                    (a jump is an UNflagged offset from one
c                     epoch to the next good epoch. If the jump
c                     occurs between epoch 1 and 2, ijump = 2)
c                     By definition, there is NO bias flag at 2.
c        ystep       maximum step in series
c                    (a step is a FLAGGED offset from one
c                     epoch to the next good epoch. If the step
c                     occurs between epoch 1 and 2, istep = 2).
c                     By definition, there is a bias flag at 2.
c        ngood       the number of good ones
c        nmarg       number of marginal or good points

      implicit none

      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
      include '../includes/errflg.h'

c     functions to sort out error codes
      logical lgood,lmarg,lbias,debug

      integer nobs,ngood,l12p,ii0,ii1,i,nmarg,istep,ijump,iprev
      real*8 ymin,ymax,yavg,ylg,ypg,yssd,raddeg,pi,ylast,yflat
     .      ,yjump,ystep
      logical domarg

c     what we are plotting
      real*8  plt(maxepc)
      integer*4 kpl(maxepc)

c     can an observable be formed?
      logical        lcando(-3:ntypes)

c     variables for prog_name
      character*80 prog_name
      integer*4 len,rcpar
  
      data debug/.false./


c     get the main program
      len = rcpar(0,prog_name)

c     dummy statement for unused calling argument
      if( debug ) print *,'nobs ',nobs

c     to keep things fast, try to minimize the number of if statements
c     within the do loops

      ymin =  1.0d25
      ymax = -1.0d25
      yavg =  0.0d0
      ngood = 0
      nmarg = 0

c     radians to degrees conversion factor
      raddeg = 180.0d0 / (4.0d0 * datan(1.0d0))

c     decide if the observable can be formed
c     L1 and L2 needed for LC,LG,WL,PG,BG,I1,I2
      if (.not. (lcando(1).and.lcando(4))) then
         lcando(2) = .false.
         lcando(3) = .false.
         lcando(7) = .false.
         lcando(8) = .false.
         lcando(9) = .false.
         lcando(10) = .false.
         lcando(11) = .false.
      endif

c     P1 and P2 needed for WL
      if (.not. (lcando(5).and.lcando(6))) then
         lcando(7) = .false.
         lcando(8) = .false.
      endif

      if (l12p .eq. -2) then
c        azimuth (degrees clockwise from North)
         do 3 i = ii0, ii1
            kpl(i) = kww(i)
            if (lgood(kpl(i)) .or. (domarg.and.lmarg(kpl(i)))) then
               plt(i) = azw(i) * raddeg
               ymax = dmax1(ymax,plt(i))
               ymin = dmin1(ymin,plt(i))
               nmarg = nmarg + 1
               if (lgood(kpl(i))) then
                  ngood = ngood+1
                  yavg = yavg + plt(i)
               endif
            endif
 3       continue
      else if (l12p .eq. -1) then
c        elevation angle in degrees
         do 5 i = ii0, ii1
            kpl(i) = kww(i)
            if (lgood(kpl(i)) .or. (domarg.and.lmarg(kpl(i)))) then
               plt(i) = elw(i) * raddeg
               ymax = dmax1(ymax,plt(i))
               ymin = dmin1(ymin,plt(i))
               nmarg = nmarg + 1
               if (lgood(kpl(i))) then
                  ngood = ngood+1
                  yavg = yavg + plt(i)
               endif
            endif
 5       continue
      else if (l12p .eq. 1) then
c        L1 phase
         do 10 i = ii0, ii1
            kpl(i) = kww(i)
            if (lgood(kpl(i)) .or. (domarg.and.lmarg(kpl(i)))) then
               plt(i) = wl1(i)
               ymax = dmax1(ymax,plt(i))
               ymin = dmin1(ymin,plt(i))
               nmarg = nmarg + 1
               if (lgood(kpl(i))) then
                  ngood = ngood+1
                  yavg = yavg + plt(i)
               endif
            endif
 10      continue
      else if (l12p .eq. 2) then
c        form LG = L2 - gear(jsat1)*L1
         do 20 i = ii0, ii1
            kpl(i) = kww(i)
            if (lgood(kpl(i)) .or. (domarg.and.lmarg(kpl(i)))) then
               plt(i) = wl2(i) - gear(jsat1)*wl1(i)
               ymax = dmax1(ymax,plt(i))
               ymin = dmin1(ymin,plt(i))
               nmarg = nmarg + 1
               if (lgood(kpl(i))) then
                  ngood = ngood+1
                  yavg = yavg + plt(i)
               endif
            endif
 20      continue
      else if (l12p .eq. 3) then
         do 30 i = ii0, ii1
c        form LC
            kpl(i) = kww(i)
            if (lgood(kpl(i)) .or. (domarg.and.lmarg(kpl(i)))) then
               plt(i) = wl1(i)-faclr(jsat1)*(wl2(i)-gear(jsat1)*wl1(i))
               ymax = dmax1(ymax,plt(i))
               ymin = dmin1(ymin,plt(i))
               nmarg = nmarg + 1
               if (lgood(kpl(i))) then
                  ngood = ngood+1
                  yavg = yavg + plt(i)
               endif
            endif
 30      continue
      else if (l12p .eq. 4) then
c        L2 phase
         do 40 i = ii0, ii1
            kpl(i) = kww(i)
            if (lgood(kpl(i)) .or. (domarg.and.lmarg(kpl(i)))) then
               plt(i) = wl2(i)
               ymax = dmax1(ymax,plt(i))
               ymin = dmin1(ymin,plt(i))
               nmarg = nmarg + 1
               if (lgood(kpl(i))) then
                  ngood = ngood+1
                  yavg = yavg + plt(i)
               endif
            endif
 40      continue
      else if (l12p .eq. 5) then
c        P1 pseudorange (cycles)
         do 50 i = ii0, ii1
            kpl(i) = kww(i)
            if (lgood(kpl(i)) .or. (domarg.and.lmarg(kpl(i)))) then
               plt(i) = pc1(i)
               ymax = dmax1(ymax,plt(i))
               ymin = dmin1(ymin,plt(i))
               nmarg = nmarg + 1
               if (lgood(kpl(i))) then
                  ngood = ngood+1
                  yavg = yavg + plt(i)
               endif
            endif
 50      continue
      else if (l12p .eq. 6) then
c        P2 pseudorange (cycles)
         do 60 i = ii0, ii1
            kpl(i) = kww(i)
            if (lgood(kpl(i)) .or. (domarg.and.lmarg(kpl(i)))) then
               plt(i) = pc2(i)
               ymax = dmax1(ymax,plt(i))
               ymin = dmin1(ymin,plt(i))
               nmarg = nmarg + 1
               if (lgood(kpl(i))) then
                  ngood = ngood+1
                  yavg = yavg + plt(i)
               endif
            endif
 60      continue
      else if (l12p .eq. 7) then
c        wide lane observable WL
         iprev = ii0
         yflat = 0.d0
         do 70 i = ii0, ii1
            kpl(i) = kww(i)
            if (lgood(kpl(i)) .or. (domarg.and.lmarg(kpl(i)))) then
               plt(i)=wl2(i)-wl1(i)+facwl(jsat1)*(pc1(i)+pc2(i))
               ymax = dmax1(ymax,plt(i))
               ymin = dmin1(ymin,plt(i))
               nmarg = nmarg + 1
               if (lgood(kpl(i))) then
                  ngood = ngood+1
                  yavg = yavg + plt(i)
               endif
            endif
 70      continue
      else if (l12p .eq. 8) then
c        wide lane observable WL
         iprev = ii0
         yflat = 0.d0
         do 80 i = ii0, ii1
            kpl(i) = kww(i)
            if (lgood(kpl(i)) .or. (domarg.and.lmarg(kpl(i)))) then
               plt(i)=wl2(i)-wl1(i)+facwl(jsat1)*(pc1(i)+pc2(i))

ccc            flatening the data
               if((i-iprev) .gt. 50)
     .                  yflat = anint(plt(iprev)-plt(i))
c     .                  yflat = yflat + anint(plt(iprev)-plt(i))
               plt(i) = plt(i) + yflat
               iprev = i

               ymax = dmax1(ymax,plt(i))
               ymin = dmin1(ymin,plt(i))
               nmarg = nmarg + 1
               if (lgood(kpl(i))) then
                  ngood = ngood+1
                  yavg = yavg + plt(i)
               endif
            endif
 80      continue
      else if (l12p .eq. 9) then
c        PG observable PG = P2 - gear(jsat1)*P1
         do 90 i = ii0, ii1
            kpl(i) = kww(i)
            if (lgood(kpl(i)) .or. (domarg.and.lmarg(kpl(i)))) then
               plt(i)=pc2(i)-gear(jsat1)*pc1(i)
               ymax = dmax1(ymax,plt(i))
               ymin = dmin1(ymin,plt(i))
               nmarg = nmarg + 1
               if (lgood(kpl(i))) then
                  ngood = ngood+1
                  yavg = yavg + plt(i)
               endif
            endif
 90      continue
      else if (l12p .eq. 10) then
c        BG observable BG = B2 - gear(jsat1)*B1
c                         = LG + PG
         do 100 i = ii0, ii1
            kpl(i) = kww(i)
            if (lgood(kpl(i)) .or. (domarg.and.lmarg(kpl(i)))) then
               ylg =  wl2(i) - gear(jsat1)*wl1(i)
               ypg =  pc2(i)- gear(jsat1)*pc1(i)
               plt(i) = ylg + ypg
               ymax = dmax1(ymax,plt(i))
               ymin = dmin1(ymin,plt(i))
               nmarg = nmarg + 1
               if (lgood(kpl(i))) then
                  ngood = ngood+1
                  yavg = yavg + plt(i)
               endif
            endif
 100      continue
      else if (l12p .eq. 11) then
c        I1 observable I1 = -g/(1-g*g) PG
c
         do 110 i = ii0, ii1
            kpl(i) = kww(i)
            if (lgood(kpl(i)) .or. (domarg.and.lmarg(kpl(i)))) then
               plt(i)= -gear(jsat1)*(pc2(i)-gear(jsat1)*pc1(i))/
     .                   (1.d0-gear(jsat1)*gear(jsat1))
               ymax = dmax1(ymax,plt(i))
               ymin = dmin1(ymin,plt(i))
               nmarg = nmarg + 1
               if (lgood(kpl(i))) then
                  ngood = ngood+1
                  yavg = yavg + plt(i)
               endif
            endif
 110     continue
      else if (l12p .eq. 12) then
c        I2 observable I2 = -1/(1-g*g) PG
c
         pi = 4.0d0 * datan(1.0d0)
         do 120 i = ii0, ii1
            kpl(i) = kww(i)
            if (lgood(kpl(i)) .or. (domarg.and.lmarg(kpl(i)))) then
               plt(i)= -(pc2(i)-gear(jsat1)*pc1(i))/
     .                         (1.d0-gear(jsat1)*gear(jsat1))
               ymax = dmax1(ymax,plt(i))
               ymin = dmin1(ymin,plt(i))
               nmarg = nmarg + 1
               if (lgood(kpl(i))) then
                  ngood = ngood+1
                  yavg = yavg + plt(i)
               endif
            endif
 120     continue
      else if (l12p .eq. 13) then
c        L1 multipath delay
c        disabled to same time and space for now kurt 900307
         do 130 i = ii0, ii1
            kpl(i) = kww(i)
            if (lgood(kpl(i)) .or. (domarg.and.lmarg(kpl(i)))) then
ckf            plt(i) = rw1(i)
               plt(i) = 0.0d0
               ymax = dmax1(ymax,plt(i))
               ymin = dmin1(ymin,plt(i))
               nmarg = nmarg + 1
               if (lgood(kpl(i))) then
                  ngood = ngood+1
                  yavg = yavg + plt(i)
               endif
            endif
 130     continue
      else if (l12p .eq. 14) then
c        L2 multipath delay
c        disabled to same time and space for now kurt 900307
         do 140 i = ii0, ii1
            kpl(i) = kww(i)
            if (lgood(kpl(i)) .or. (domarg.and.lmarg(kpl(i)))) then
ckf            plt(i) = rw2(i)
               plt(i) = 0.0d0
               ymax = dmax1(ymax,plt(i))
               ymin = dmin1(ymin,plt(i))
               nmarg = nmarg + 1
               if (lgood(kpl(i))) then
                  ngood = ngood+1
                  yavg = yavg + plt(i)
               endif
            endif
 140     continue
      else if (l12p .eq. 15) then
c        LC multipath delay
c        disabled to same time and space for now kurt 900307
         do 150 i = ii0, ii1
            kpl(i) = kww(i)
            if (lgood(kpl(i)) .or. (domarg.and.lmarg(kpl(i)))) then
ckf             plt(i) = rw1(i)-faclr(jsat1)*(rw2(i)-gear(jsat1)*rw1(i))
               plt(i) = 0.0d0
               ymax = dmax1(ymax,plt(i))
               ymin = dmin1(ymin,plt(i))
               nmarg = nmarg + 1
               if (lgood(kpl(i))) then
                  ngood = ngood+1
                  yavg = yavg + plt(i)
               endif
            endif
 150     continue
      else if (l12p .eq. 16) then
c        form MG = M2 - gear(jsat1)*M1
c        disabled to same time and space for now kurt 900307
         do 160 i = ii0, ii1
            kpl(i) = kww(i)
            kpl(i) = 0.0d0
            if (lgood(kpl(i)) .or. (domarg.and.lmarg(kpl(i)))) then
ckf            plt(i) = rw2(i) - gear(jsat1)*rw1(i)
               plt(i) = 0.0d0
               ymax = dmax1(ymax,plt(i))
               ymin = dmin1(ymin,plt(i))
               nmarg = nmarg + 1
               if (lgood(kpl(i))) then
                  ngood = ngood+1
                  yavg = yavg + plt(i)
               endif
            endif
 160     continue
      else if (l12p .eq. 17) then
c        form clock offset in microseconds
         do 170 i = ii0, ii1
            kpl(i) = kww(i)
            if (lgood(kpl(i)) .or. (domarg.and.lmarg(kpl(i)))) then
               plt(i) = dble(clw(i)) * 1.0d6
               ymax = dmax1(ymax,plt(i))
               ymin = dmin1(ymin,plt(i))
               nmarg = nmarg + 1
               if (lgood(kpl(i))) then
                  ngood = ngood+1
                  yavg = yavg + plt(i)
               endif
            endif
 170     continue
      else
         call report_stat('WARNING',prog_name,'combo',l12p,
     .   'Unknown data type: ',0)
      endif

      if (ngood .gt. 0) then
         yavg = yavg/ngood
      else
         yavg = 0.0d0
      endif

c     loop through again to get sample standard deviation
c     Also find largest jump and largest step
      ijump = 0
      istep = 0
      yssd =  0.0d0
      yjump = 0.0d0
      ystep = 0.0d0
      ylast = plt(ii0)
      do 200 i = ii0,ii1
         if (lgood(kpl(i)) .or. (domarg.and.lmarg(kpl(i)))) then
            yssd = yssd + (plt(i)-yavg)**2
            if (.not. lbias(kpl(i))) then
               if (dabs(plt(i)-ylast) .gt. dabs(yjump)) then
                  yjump = plt(i) - ylast
                  ijump = i
               endif
            else if (lbias(kpl(i))) then
               if (dabs(plt(i)-ylast) .gt. dabs(ystep)) then
                  ystep = plt(i) - ylast
                  istep = i
               endif
            endif
            ylast = plt(i)
         endif
 200  continue

      if (ngood. gt. 1) then
          yssd = dsqrt(yssd/dble(ngood-1))
      else
          yssd = 0.0d0
      endif

CD     if (ngood .gt. 0) then
CD        do 900 i=ii0,ii1
CD           print *,'COMBO: ',i,kpl(i),plt(i)
CD900     continue
CD     endif
CD     print *,'COMBO: ngood,nmarg ',ngood,nmarg
CD     print *,'COMBO: ymin,ymax,yavg,yssd ',ymin,ymax,yavg,yssd
CD     print *,'COMBO: yjump ',yjump


      return
      end

