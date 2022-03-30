      subroutine lswght(job)
c
c     put loose constraints on coordinates and velocities
c     job = 1: weight coordinate and velocity
c     job = 2: weight auxiliary parameters
c
      implicit real*8(a-h,o-z)
      include 'solvem.fti'
c
      integer job,is,indx,indv,ia,ic,id
      dimension indx(12),coef(12),indv(12)
c
      if (job.eq.2) goto 100
      f = 1.0d0/finv
      e2 = 2.0d0*f-f*f
      wcoe = 1.0d0/wcoe**2
      wcon = 1.0d0/wcon**2
      wcou = 1.0d0/wcou**2
      wvee = 1.0d0/wvee**2
      wven = 1.0d0/wven**2
      wveu = 1.0d0/wveu**2
      do 50 is = 1,nsit
         sl = dsin(slon(is))
         cl = dcos(slon(is))
         sf = dsin(slat(is))
         cf = dcos(slat(is))
         fac = 1.0d0/(1.0d0-e2*sf**2)
c        insert position weighting submatrix
         indx(1) = (is-1)*6+1
         indx(2) = indx(1)+1
         indx(3) = indx(1)+2
c        insert velocity weighting submatrix
         indv(1) = is*6-2
         indv(2) = indv(1)+1
         indv(3) = indv(1)+2
c        east component
         coef(1) = -sl
         coef(2) = cl
         coef(3) = 0.0d0
         call norms(3,indx,coef,wcoe,omc,2)
         call norms(3,indv,coef,wvee,omc,2)
c        north component
         coef(1) = -sf*cl
         coef(2) = -sf*sl
         coef(3) = cf
         call norms(3,indx,coef,wcon,omc,2)
         call norms(3,indv,coef,wven,omc,2)
c        up component
         coef(1) = cf*cl*fac
         coef(2) = cf*sl*fac
         coef(3) = sf*(1.0d0-e2)*fac
         call norms(3,indx,coef,wcou,omc,2)
         call norms(3,indv,coef,wveu,omc,2)
 50   continue
      goto 200
c
c     loose constrain auxiliary parameters
 100  if (iaux.gt.0) then
         ia = iaux
         wght = (1.0d0/3.6d4)**2
         do 60 is = 1,ia
            ic = nsit*6+iq_sit*3+jaux+is
            id = ic*(ic+1)/2
            anorm(id) = anorm(id)+wght
 60      continue
      endif
      if (jaux.gt.0) then
         wght = 1.0d-4
         do 80 is = 1,jaux
            ic = nsit*6+iq_sit*3+is
            id = ic*(ic+1)/2
            ia = is-(is-1)/6*6
            if (ia.le.3) anorm(id) = anorm(id)+wght
            if (ia.gt.3) anorm(id) = anorm(id)+1.0d0
 80      continue
      endif
      if (iq_sit.gt.0) then
         ia = nsit*6
         wght = 1.0d-4
         do 70 is = 1,iq_sit*3
            ic = nsit*6+is
            id = ic*(ic+1)/2
            anorm(id) = anorm(id)+wght
 70      continue
      endif
       
 200  continue
      return
      end

