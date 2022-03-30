      subroutine add_quake_c(isit1,isit2,iterm,dt,coef,indx)
c
c     add earthquake related coefficients
c
c     unit:
c        dt  : year
c
      implicit real*8(a-h,o-z)
      include 'solvem.fti'

      integer imatch,i,j
      integer isit1,isit2,iterm,iek(10),nbase,indx(30)
      dimension coef(30),fac(10)
      nbase = nsit*6
c
c     check earthquake list 
      call chk_quake_list(isit1,dt,iek,fac,imatch)
      if (imatch.gt.0) then
         do i = 1,imatch
            j = iek(i)
            indx(iterm+1) = nbase+j*3-2
            indx(iterm+2) = nbase+j*3-1
            indx(iterm+3) = nbase+j*3
            quake_use(j) = quake_use(j)+1
            coef(iterm+1) = fac(i)*coef(1)
            coef(iterm+2) = fac(i)*coef(2)
            coef(iterm+3) = fac(i)*coef(3)
            iterm = iterm+3
         enddo
      endif
      call chk_quake_list(isit2,dt,iek,fac,imatch)
      if (imatch.gt.0) then
         do i = 1,imatch
            j = iek(i)
            indx(iterm+1) = nbase+j*3-2
            indx(iterm+2) = nbase+j*3-1
            indx(iterm+3) = nbase+j*3
            quake_use(j) = quake_use(j)+1
            coef(iterm+1) = fac(i)*coef(7)
            coef(iterm+2) = fac(i)*coef(8)
            coef(iterm+3) = fac(i)*coef(9)
            iterm = iterm+3
         enddo
      endif
         
      return
      end
