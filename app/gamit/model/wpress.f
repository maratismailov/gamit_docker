      subroutine wpress(iflag,rh,e,t)
c
c     If iflag =  1  Gets partial pressure of water vapor from relative humidity
c     if iflag = -1  Get relative humidity from partial pressure of water vapor
c
c     INPUT:
c       RH     Relative humidity (0-1)
c       E      Partial pressure of water vapor (mb)
c       T      Temperature, deg C
c                    
      integer*4 iflag
      real*8 rh, t, e
c                    
      if( iflag.eq.1 ) then
        e = rh * 6.11D+00
     .         * 10.0D+00 ** (7.5D+00 * t / (t + 2.373D+02))
      elseif( iflag.eq.-1 ) then
        rh = e / ( 6.11D+00
     .            * 10.0D+00 ** (7.5D+00 * t / (t + 2.373D+02)))
      endif
      end
