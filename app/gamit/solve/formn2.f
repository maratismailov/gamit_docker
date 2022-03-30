c
      subroutine formn2(gearf,iphilv,iones,isame,key,bl1,bl2)
c
c     form second part of normal matrix (stackable)
c
      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'

      real*8 gearf,gsq,g1,g2
      integer*4 iones,isame,key

      integer iphilv(maxobs)
      real*8 bl1(maxobs),bl2(maxobs)
c*     this now 'gear' in solve.h' gearf may be 'gear' or 0 
c*      data gears/0.779220779220779d+00/

      gsq = gearf*gearf
c
c form n22 portion of normal matrix
      if(l2flag.ge.-1) then
      g1 = 1.d0-gearf
      g1 = g1*g1
      g2 = 1.d0+gearf
      g2 = g2*g2
      elseif(l2flag.eq.-2) then
      g1=1.d0
      g2=0.d0
      endif    
      call nm2233(iphilv,ipntb1,ibcnt1,g1,g2,isame,key)
c
c form n21 portion of normal matrix
      if(l2flag.ge.-1) then
      g1 = 1.d0-gsq
      g1 = (g1*g1)/(1.d0+gearf)
      g2 = 2.d0*gear*(1.d0+gearf)
      elseif(l2flag.eq.-2) then
      g1=gear
      g2=0.d0
      endif
      call nm2131(lpart,iones,iphilv,ipntb1,ibcnt1,g1,g2,key)
c
c form u2 portion of normal equations
      if(l2flag.eq.-2) then
       g1=1.d0
       g2=0.d0
      endif
      call nu23(iones,iphilv,ipntb1,ibcnt1,bl1,bl2,g1,g2,key)   
c
c skip if l1 or l2
      if(l2flag.le.1) go to 100
c
c form n33 portion of normal matrix
      g1 = gsq
      g2 = 1.d0
      call nm2233(iphilv,ipntb2,ibcnt2,g1,g2,isame,0)   
c
c form n31 portion of normal matrix
      g1 = -gearf*(1.d0-gsq)
      g2 = 2.d0*gearf
      call nm2131(lpart,iones,iphilv,ipntb2,ibcnt2,g1,g2,0)
c
c form u3 portion of normal equations
      call nu23(iones,iphilv,ipntb2,ibcnt2,bl1,bl2,g1,g2,0)   
c
c form n32 portion of normal matrix
      g1 = -gearf*(1.d0-gearf)
      g2 = (1.d0+gearf)
      call nm32(iphilv,ibcnt1,g1,g2,isame) 
c
 100  continue
c
      return
      end


