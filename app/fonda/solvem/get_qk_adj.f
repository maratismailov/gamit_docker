      subroutine get_qk_adj(isit1,iterm,dt,adj)
c
c     get earthquake related adjusts
c
c     unit:
c        dt  : year
c
      implicit real*8(a-h,o-z)
      include 'solvem.fti'

      integer isit1,iterm,iek(10),imatch
      integer nbase,i,im,j,j1
      dimension adj(6),fac(10)
c
      iterm = 0
c
c     check earthquake list 
      if (iomode(3).le.0.or.iq_sit.le.0) goto 30
      nbase = nsit*6
      call chk_quake_list(isit1,dt,iek,fac,imatch)
      if (imatch.gt.0) then
         do i = 1,3
            adj(i) = 0.0d0
         enddo
         do j = 1,imatch
            j1 = quake_use(iek(j))
            if (j1.gt.0) then
               do i = 1,3
                  temp = fac(j)
                  im = map(nbase+(iek(j)-1)*3+i)
                  adj(i) = adj(i)+bnorm(im)*temp
               enddo
            endif
         enddo
      endif
      iterm = imatch
c     
 30   continue
         
      return
      end
