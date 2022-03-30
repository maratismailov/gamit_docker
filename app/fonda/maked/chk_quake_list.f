      subroutine chk_quake_list(isit,dt,iek,fac)
c
c     check if the site is in the earthquake site list
c
c     unit:
c        dt  : year
c
      implicit real*8(a-h,o-z)
      include 'maked.fti'

      integer isit,is,i0,i1,i,j,iek
c
c     check earthquake list 
      iek = 0
      fac = 0.0d0
      do 20 i = 1,iq_sit
         is = quake_sit(i)
         if (is.eq.isit) then
            if (i.ge.iq_ind(nquake).or.nquake.eq.1) then
               tq = quake_time(nquake)-rtime
            else
               i0 = 1
               do 10 j = 2,nquake
                  i1 = iq_ind(j)
                  if (i.ge.i0.and.i.lt.i1) then
                     tq = quake_time(j-1)-rtime
                     goto 30
                  endif
                  i0 = i1
 10            continue
            endif
         else 
            goto 20
         endif
 30      if (tq.lt.0.0d0.and.dt.gt.tq) goto 20
         if (tq.gt.0.0d0.and.dt.lt.tq) goto 20
         iek = i
         if (tq.lt.0.0d0.and.dt.lt.tq) fac = -1.0d0
         if (tq.gt.0.0d0.and.dt.gt.tq) fac = 1.0d0
         goto 100
 20   continue

 100  continue
      return
      end
