      subroutine chk_epi_list(isit,et,iek,imatch)
c
c     check if the site is in the episodic site list
c
c     unit:
c        et  : year
c
      implicit real*8(a-h,o-z)
      include 'solvem.fti'

      integer isit,is,i0,i1,i,j,iek,imatch,ic,i2
c
c     check episodic site list 
      imatch = 0
      i0 = nlive-iaux-jeaux
      do 10 j = 1,nquake
         tq = quake_time(j)
         if (dabs(tq-et).gt.1.0d-2) goto 10
         i1 = iq_ind(j)
         i2 = iq_sit
         if (j.lt.nquake) i2 = iq_ind(j+1)-1
         is = 0
         do 50 i = i1,i2
            if (quake_use(i).le.0) goto 50
            if (quake_sit(i).eq.isit) is = i
 50      continue
         if (is.le.0) goto 10
         ic = 0
         do 40 i1 = 1,is
            if (quake_use(i1).le.0) goto 40
            ic = ic+3
 40      continue
         if (ic.le.0) goto 100
         if (ic.gt.0) goto 30
 10   continue
      goto 100
 30   imatch = imatch+1
      iek = i0+ic-2

 100  continue
      return
      end
