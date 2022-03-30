      Logical function wl_fixed(half,adjwl )
                 
c       wavelength for bias (=1 or 2 for full or half) 
      integer*4 half  

c       adjustment for associated WL 
      real*8 adjwl

c       local
      real*8 bl,bl1,bl2,bld 
    
      bl1=adjwl
      bl=bl1*dble(half)
      bl2=dint(bl+.5d0*dsign(1.d0,bl))/dble(half)
      bld=dabs(bl1-bl2)    
c       arbitrary tolerance to determine if bias has been fixed to an
c       integer or half-integer
      if (bld.le.0.0001) then
         wl_fixed = .true.
      else
         wl_fixed = .false.
      endif                                
c      print *,'FUN half adjwl bl1 bl bl2 bld  wl_fixed '
c     .       ,      half, adjwl,bl1,bl,bl2,bld,wl_fixed
      return
      end


