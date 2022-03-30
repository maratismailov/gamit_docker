      Subroutine SHADOW1 ( satcrd,sun,lambda )
c
c Purpose: Computation of fraction (lambda) of solar disk seen by spacecraft
c          Beebe, King, Reasonberg, Preston: June 1971
c          McClusky 1995
c
c PARAMETERS:
c         IN: sun     : coordinates of the sun wrt earth CM    R*8(6)
c             satcrd  : coordinates of SV wrt earth CM         R*8(6)
c
c        OUT: lambda  : shadow factor                          R*8
c
c       USED: sbcor   : coordinates SV wrt earth
c             bcor    : coordinates SV wrt sun
c             ccor    : coordinates earth wrt sun
c             rb      : dist SV to sun
c             rc      : dist earth to sun
c
c SUBROUTINES CALLED:  timinc, timdif
c
c CREATED:  MAR 1971               LAST MODIFIED: 16th MAR 1995
c
c AUTHOR: rwk, modified by S McClusky.
c
c COPYRIGHT: DEPARTMENT OF EARTH AND PLANETRY SCIENCES
c            M.I.T. 1995
c
      implicit none
c
      include '../includes/dimpar.h'
c
      integer i
      integer j

      real*8 lambda,satcrd(6),sun(6)
      real*8 sunrad,pcrad,pi,rb,rc
      real*8 sbcor(3),bcor(3),ccor(6)
      real*8 twopi,vlight
      real*8 dot,amag3,tlag,pcvel(3)
      real*8 rsbx,sbcorx(3),area3,area1,area2,ari,rs,rp,sep,phi
      real*8 r1,r2,thet,hgt,sepp(3),ubcor(3)
c             
      data twopi/6.283185307179586d0/ 
      data vlight/299792.458d0/
 
      sunrad=6.96D+05
      pcrad=6378.145D0
      lambda=1.0d0
      pi=twopi/2.d0
         
c
c shadow computation - geometric model
c lambda = 1 - no shadow
c lambda = 0 - no sunlight
c 0 .lt. lambda .lt. 1 - partial shadow
c
      do j = 1,3
        ccor(j) = -sun(j)
        ccor(j+3) = -sun(j+3)
        sbcor(j) = satcrd(j)
        bcor(j) = sbcor(j) - sun(j)
      enddo
c      print *, 'CCOR ',(ccor(i),i=1,3)
c      print *, 'SBCOR ',(sbcor(i),i=1,3)
c      print *, 'BCOR ',(bcor(i),i=1,3)
c
      rb = amag3(bcor)
      rc = amag3(ccor)

c      print*, 'SV-SUN ',rb,' EARTH-SUN ',rc
c      print *, 'UBCOR ',(ubcor(i),i=1,3)
c
c no consideration is given to the change of direction associated with partial shadow.
c
      if(rb.gt.rc) then
c
c        adjust sbcor for movement of earth since light passed it
        do 363 i=1,3
          ubcor(i)=bcor(i)/rb
          pcvel(i)=ccor(i+3)/86400.d0
  363   continue
c        tlag is the light propagation time (seconds) from the earth
c        to the satellite along the sun-satellite line of sight.
c        sbcorx is the position of the s/c wrt the earth -
c        satellite at current time / earth earlier by tlag

c       print*,' ubcor = ',ubcor

        tlag = dot(sbcor,ubcor)/vlight
c       print*,' tlag =',tlag
c       print*,' pvel = ',pcvel
        do 364 i=1,3
          sbcorx(i)=sbcor(i)+tlag*pcvel(i)
  364   continue
c
        call cross(sbcorx,ubcor,sepp)
c        rsbx is the projection of sbcorx along bcor
        rsbx=dot(sbcorx,ubcor)
c
c        rs, rp are apparent (from satellite) radii of sun and earth
c        sep is apparent separation of their centers
        rs=sunrad/rb
        rp=pcrad/rsbx
        sep=dsqrt(sepp(1)**2+sepp(2)**2+sepp(3)**2)/rsbx
cd      write (6,55555) rs,rp,sep
cd 55555   format (1x,'shadow...',3(d12.5,1x))
        if(rs+rp.le.sep)then
         go to 373
        elseif(rp-rs.ge.sep)then
          go to 372
        elseif(sep.le.rs-rp)then
          go to 361
        else
        endif
c
c        set r1 = smaller disc, r2 = larger
        if(rs.gt.rp)go to 371
        r1=rs
        r2=rp
        go to 370
  371   continue
        r1=rp
        r2=rs
  370   continue
c
c        phi = 1/2 angle subtended in disc 1 by arc of intersection
        phi = dacos((r1*r1+sep*sep-r2*r2)/(2.0d0*r1*sep))
        if (phi.lt.0.d0) phi = pi + phi
        if(r2/r1.gt.5.0d0)go to 365
c        thet = 1/2 angle subtended in disc 2 by arc of intersection
c        hgt  = 1/2 linear distance between ends of arc of intersection
        hgt=r1*dsin(phi)
        thet=dasin(hgt/r2)
        area2=sep*hgt
        area3=thet*r2**2
        go to 366
  365   continue
c        one disc much bigger - treat boundary as a straight line
        hgt=dsqrt(r1**2-(sep-r2)**2)
        area2=hgt*(sep-r2)
        area3=0.0d0
  366   continue
        area1=(pi-phi)*r1**2
c       ari = area of non-overlapped portion of small disc
        ari=area1+area2-area3
        area1=pi*rs**2
        if(rs.gt.rp)go to 362
c
c        sun is small disc
        lambda=ari/area1
        go to 373
  362   continue
c       earth is small disc
        area2=pi*rp**2
        lambda=(area1+ari-area2)/area1
        go to 373
  361   continue
c       earth lies within sun's disc - what fraction of sun's disk is blocked
        lambda=(rs**2-rp**2)/rs**2
        go to 373
  372   continue
c        sun is completely eclipsed by planet
        lambda=0.0d0
      endif

  373 continue
c
cd      write(6,101)lambda
cd 101   format(1x,'shadow... lambda= ',f7.3)
c
      return
      end
