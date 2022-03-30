      subroutine chk_quake_list(isit,dt,iek,fac,imatch)
c
c     check if the site is in the earthquake site list
c
c     unit:
c        dt  : year
c
      implicit real*8(a-h,o-z)
      include 'solvem.fti'

      integer isit,is,i0,i1,i,j,iek,imatch
      dimension fac(10),iek(10)
c
c     check earthquake list 
c     note: same site could have multiple earthquakes. Dong 11/19/92
      imatch = 0
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
         imatch = imatch+1
         iek(imatch) = i
         fac(imatch) = 0.0d0
         if (tq.lt.0.0d0.and.dt.lt.tq) fac(imatch) = -1.0d0
         if (tq.gt.0.0d0.and.dt.gt.tq) fac(imatch) = 1.0d0
 20   continue

 100  continue
      return
      end
