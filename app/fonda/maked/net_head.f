      subroutine net_head(inetf)
c
      implicit real*8(a-h,o-z)
      include 'maked.fti'
      integer i2,inetf
      integer nblen
      character*20 txt

c
      i2 = nblen(netfil)
c
      call pline(inetf,125,'=',1)
      write (inetf,'(a1,30x,a44,49x,a1)') '!',
     .  'Network apriori coordinate and velocity file','!'
      call pline(inetf,125,'=',1)
      write (inetf,'(a1,1x,a8)') '*','history:'
      if (i2.gt.0) write (inetf,'(a1,1x,a17,2x,a)') '*',
     .   'transfered from: ',netfil(1:i2)
c     write (inetf,'(a1,1x,a62)') '*',
c    .'question:  which frame is prefered by most users, xyz or neu ?'
      if (iomode(8).eq.1) then
         txt='m/yr'
      else
         txt='mm/yr'
      endif
      write (inetf,'(a1,1x,a62)') '*',
     .'Velocity (u,v,w) are '//txt//' in the NEU system             '
      write (inetf,'(a1)') '*'
      write (inetf,'(a1,2x,a4,4x,a9,6x,a8,8x,a9,
     .   9x,a6,7x,a1,2(7x,a1),4x,a5,3(6x,a2))')
     .   '*','site','full-name','latitude','longitude','height',
     .   'u','v','w','epoch','Sn','Se','Su'
      
c
      return
      end

