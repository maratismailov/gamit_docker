      subroutine calc_res(ns,isit1,isit2,isit3,iterm,dt,coef,adj)
c
c     add earthquake related coefficients
c
c     unit:
c        dt  : year
c
      implicit real*8(a-h,o-z)
      include 'solvem.fti'

      integer ns,isit1,isit2,isit3,iterm,iek(10),imatch
      integer nbase,i,im,j,j1
      integer istsav,iauxsav,indsav,iazsav
      dimension coef(30),tempc(30),fac(10)
      dimension indsav(12),coesav(12)
      common/saved/azisav,coesav,obssav,esave,tsave,iazsav,indsav
      common/save2/istsav,iauxsav
c
c     get correction to parameters (tempc)
      call getadj(ns,isit1,isit2,isit3,tempc)

      if (iterm.eq.13)
     .   tempc(13) = bnorm(nlive-iaux+iauxsav)
c
c     check earthquake list 
      if (iomode(3).le.0) goto 30
      nbase = nsit*6
      call chk_quake_list(isit1,dt,iek,fac,imatch)
      if (imatch.gt.0) then
         do j = 1,imatch
            j1 = quake_use(iek(j))
            if (j1.gt.0) then
               do i = 1,3
                  coef(iterm+i) = fac(j)*coef(i)
                  im = map(nbase+(iek(j)-1)*3+i)
                  tempc(iterm+i) = bnorm(im)
               enddo
               iterm = iterm+3
            endif
         enddo
      endif
      if (ns.le.1) goto 30
      call chk_quake_list(isit2,dt,iek,fac,imatch)
      if (imatch.gt.0) then
         do j = 1,imatch
            j1 = quake_use(iek(j))
            if (j1.gt.0) then
               do i = 1,3
                  coef(iterm+i) = fac(j)*coef(6+i)
                  im = map(nbase+(iek(j)-1)*3+i)
                  tempc(iterm+i) = bnorm(im)
               enddo
               iterm = iterm+3
            endif
         enddo
      endif
      if (ns.le.2) goto 30
      call chk_quake_list(isit3,dt,iek,fac,imatch)
      if (imatch.gt.0) then
         do j = 1,imatch
            j1 = quake_use(iek(j))
            if (j1.gt.0) then
               do i = 1,3
                  coef(iterm+i) = fac(j)*coef(12+i)
                  im = map(nbase+(iek(j)-1)*3+i)
                  tempc(iterm+i) = bnorm(im)
               enddo
               iterm = iterm+3
            endif
         enddo
      endif
c     get residual

cmk   Need to call tempc2 for observation residuals
 30   call axb(1,iterm,1,coef,tempc,adj,1,0)
c30   call axb(1,iterm,1,coef,tempc2,adj,1,0)
         
      return
      end
