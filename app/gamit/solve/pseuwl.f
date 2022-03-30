c
       Subroutine PSEUWL( iphi,kkey )

c     calculate pseudo-range wide-lane biases

      implicit none

      include '../includes/dimpar.h' 
      include 'solve.h'

      integer iphi(maxsit,maxsat),kkey,isite
     .     ,  isit,jsit,jsat,ndata

      real*8 dwl,owl,var,var8,dev,vv,small

      parameter (small = 2.5d-9)

      go to (100,200,300,400),kkey
c
c---- kkey = 1: calculate wl and sdv sequentially
c

 100  do 20 isite = 1,nsite
         jsit =  isite
         do 10 jsat = 1,nsat
            if (iphi(isite,jsat).eq.0) goto 10
            dwl = wl1(isite,jsat)
            owl = wl0(jsit,jsat)
            var = vwl(jsit,jsat)
            var8 = 8.0d0 * var
            dev = dabs(dwl - owl)
            if (var.lt.small .or. dev.lt.var8) then
               ndata = idwl(jsit,jsat) + 1
               call avgvar(ndata,dwl,owl,var)
               wl0(jsit,jsat) = owl
               vwl(jsit,jsat) = var
               idwl(jsit,jsat) = ndata
            endif
cd         print *,'PSEUWL kkey jsit jsat wl1 wl0 vwl idwl  '
cd     .                 , kkey,jsit,jsat,wl1(jsit,jsat),wl0(jsit,jsat)
cd     .                 ,vwl(jsit,jsat),idwl(jsit,jsat)
 10      continue
 20   continue
      goto 400

c
c---- kkey = 2: remove first effective epoch value from wl
c               calculate sdv for the mean value(wl)
c
 200  continue
            do 120 isit = 1,nsite
               jsit =  isit
               do 110 jsat = 1,nsat
               ndata = idwl(jsit,jsat)
               if (ndata.le.0) goto 110
               dwl = wl1(isit,jsat)
               wl0(jsit,jsat) = wl0(jsit,jsat) - dwl
               vv = vwl(jsit,jsat)/dsqrt(dble(ndata))
               vwl(jsit,jsat) = vv
cd         print *,'PSEUWL kkey jsit jsat wl1 wl0 vwl idwl  '
cd     .                 , kkey,jsit,jsat,wl1(jsit,jsat),wl0(jsit,jsat)
cd     .                 ,vwl(jsit,jsat),idwl(jsit,jsat)

 110           continue
 120        continue
      goto 400

c---- calculate standard deviation
 300  continue
 400  continue
      return
      end

