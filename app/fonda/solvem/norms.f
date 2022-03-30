      subroutine norms(iterm,indx,coef,wght,omc,mode)
c
c     form submatrix of normal matrix and right term related
c     to the observation with multi sites;  insert the submatrix
c     into global normal matrix
c
c     iterm : size of the submatrix
c             6 for one site coordinates and velocities
c     indx  : index indicator to the global normal matrix
c     omc   : observed - calculated 
c     wght  : weight (currently neglect correlation)
c     coeff : calculated coefficients
c     anorm : global normal matrix
c     bnorm : global right hand term
c     mode  : 1 -- calculate bnorm
c             2 -- not calculate bnorm
c
      implicit real*8(a-h,o-z)
      include 'solvem.fti'

      integer iterm,indx,mode,i1,i2,is1,is2,id
      dimension indx(iterm),coef(iterm)
c
      do 30 i1 = 1,iterm
         is1 = indx(i1)
         do 20 i2 = i1,iterm
            is2 = indx(i2)
            if (is1.ge.is2) id = is1*(is1-1)/2+is2
            if (is2.gt.is1) id = is2*(is2-1)/2+is1
            anorm(id) = anorm(id)+coef(i1)*coef(i2)*wght
 20      continue
         if (mode.ne.1) goto 30
c        right hand term
         bnorm(is1) = bnorm(is1)+coef(i1)*omc*wght
 30   continue
c
 100  continue
      return
      end
c
