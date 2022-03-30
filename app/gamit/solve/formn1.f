c
      Subroutine FORMN1(gearf,iphi,iphilv,iones,iflg,key,al1,al2)
c
c     form first part of normal matrix (not stackable)
c
      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'

      real*8 gearf,gsq,g1,g2
      integer*4 iones,iflg,key,ipos,i,j
      integer iphilv(maxobs),iphi(maxsit,maxsat)
      real*8 al1(maxobs),al2(maxobs)

      logical debug /.false./


c     this now 'gear' in solve.h; gearf may be 'gear' or 0 
c      data gears/0.779220779220779d+00/
c                      
      gsq = gearf*gearf
c
c     form n11 portion of normal equations epoch by epoch
      if(l2flag.ge.-1) then
         g1 = 1.d0-gsq
         g1 = g1*g1
         g2 = 4.d0*gsq
      elseif(l2flag.eq.-2) then
         g1=gear
         g2=0.d0
      endif             
      call nmu11(lpart,iones,al1,al2,iflg,g1,g2,key)          
c
c indicate live and dead phases (0-dead, 1-live)
      ipos = 0
      do 50 i = 1,nsite
         do 60 j = 1,nsat
         ipos = ipos+1
         iphilv(ipos) = 0
c        iphilv is ordered by station, like biases
         if(iphi(i,j).ne.0) iphilv(ipos) = 1
         if(debug) print *,'FORMN1 kepoch isite isat ipos iphilv(ipos) '
     .      , kepoch,i,j,ipos,iphilv(ipos) 
 60      continue
 50   continue    

c
      return
      end


