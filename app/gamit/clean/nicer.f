      subroutine nicer(xin, xout, xtick)
c     calls only standard library routines
c     routine for scaling intervals and providing tick marks for plotter axes.
c     between 7 and 15 ticks are made, which is suitable for 10in axes.
c     input parameters
c     xin(1),xin(2)  extremes of variable  x  in its own units
c     output parameters
c     xout(1),xout(2)  adjusted extremes, made to be round numbers.  note new
c     interval always covers old one.
c     xtick  distance between tick marks in  x  units (not inches).  this
c     number is always a round number.
      real*8 divs(4),xin(2),xout(2),e,plus,bias,xtick,units,diff
      integer index
      data e/1.0e-7/,divs/.1,.2,.5,1.0/
c
      xout(1)=xin(1)
      xout(2)=xin(2)
c     Handle the awful case Kurt 910507
c     if (xout(2).eq.xout(1)) xout(2)=1.0 + 1.1*xout(2)
      diff = xout(2)-xout(1)
      if (diff .le. 1.0d-99) then
         xout(2)=1.0 + 1.1*xout(2)
         diff = xout(2)-xout(1)
      endif
      plus=1000.0+dlog10(diff)
      index=1.4969 + 2.857*dmod(plus,1.0D0)
      if (index .lt. 1 .or. index .gt. 4) then
c        error case
         xtick = 0.
      else
         units=divs(index)*10.0**(int(plus)-1000)
         bias=(1.+e)*units*dint(1.+dmax1(dabs(xout(1)),
     .                          dabs(xout(2)))/units)
         xout(1)=xout(1) - dmod(xout(1)+bias,units)
         xout(2)=xout(2) - dmod(xout(2)-bias,units)
         if (dabs(xout(1)/units) .le. .01) xout(1)=0.0
         if (dabs(xout(2)/units) .le. .01) xout(2)=0.0
         xtick=units
      endif
      return
      end

