      subroutine bexpec(ii,adn,sclerr,sigma,smax,bepc,key)
c
C     Calculate expectation bias value according to JPL.
C            -DND-  880309
c
c     key = 1 :  just calculate probability of nearest integer.

      integer*4 key,ii,j
      real*8 adn,sclerr,sigma,smax,bepc,factor,c,d1,d2,a1,a2,bint
      real*8 trun,s1
c

      factor = 1.0d0
      c = sclerr*sigma*factor
      bint = dint(adn+.5d0*dsign(1.d0,adn))
      s1 = 0.5d0/(c**2)
      smax = dexp(-(bint-adn)**2*s1)
      if (key.eq.1) goto 100
      trun = 1.d-5*smax
      c = smax
      bepc = bint*smax
      do 430 j = 1,50
         a1 = bint+dble(j)
         a2 = bint-dble(j)
         d1 = dexp(-(a1-adn)**2*s1)
         d2 = dexp(-(a2-adn)**2*s1)
         if (d1.lt.trun.and.d2.lt.trun) goto 440
         c = c+d1+d2
         bepc = bepc+a1*d1+a2*d2
 430  continue
 440  bepc = bepc/c   
c** rwk 070328:  'ii' is passed only for this debug, but to avoid a compiler warning,
c                 use it in a dummy way here
       j= ii 
c      write (6,340) ii,j,adn,bepc
c      write (10,340) ii,j,adn,bepc
c 340  format(1x,'index,truncate,initial,expectation =',2i4,2x,2f8.2)
 100  continue
      return
      end

