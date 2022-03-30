      subroutine exert_link1(npar,m,indx,coef,obsc,chi,cl)
c
c     exert single constraint with most elements = 0
c     theory:
c       constraint eqn:                      | xo |
c                        L = A X = ( 0  A2 ) | xc |
c       new solution:
c                        X(new) = X(ori) + C(ori) (0 A2)~ W ax
c                        C(new) = C(ori) - C(ori) (0 A2)~ W C(ori) (0 A2) 
c          where      
c                     ax = L - A2 xc
c
      implicit real*8(a-h,o-z)
      include 'solvem.fti'
c
      integer npar,i,i1,i2,j,indx(m)
      integer j1,id1,ij,m
      dimension coef(m)
c      
c     calculate ax
      work = 0.0d0
      do i = 1,m
         i1 = indx(i)
         work = work+coef(i)*bnorm(i1)
         work1 = 0.0d0
         id1 = i1*(i1-1)/2
         do j = 1,m
            i2 = indx(j)
            if (i1.ge.i2) ij = id1+i2
            if (i1.lt.i2) ij = i2*(i2-1)/2+i1 
            work1 = work1+coef(j)*anorm(ij)
         enddo
         scale(i) = work1
      enddo
      ax = obsc - work
c
c     calculate W
      work = 0.0d0
      do i = 1,m
         work = work+coef(i)*scale(i)
      enddo
c     print*, 'work= ',work
      aqa = cl+work
      if (aqa.le.1.0d-8) then
c     if (aqa.le.1.0d-12) then
         print*,' exert_link1: zero diag',indx(1),indx(2),aqa
         aqa = 0.0d0
         goto 100
      else
         aqa = 1.0d0/aqa
      endif
c
c     update chi2
      aq1 = aqa*ax
      chi = chi+aq1*ax
c
c     update x
      do 20 i = 1,npar
         work = 0.0d0
         i1 = i*(i-1)/2
         do j = 1,m
            j1 = indx(j)
            if (i.ge.j1) then
               id1 = i1+j1
            else 
               id1 = j1*(j1-1)/2+i
            endif
            work = work+coef(j)*anorm(id1)
         enddo
         scale(i) = work
         bnorm(i) = bnorm(i)+scale(i)*aq1
 20   continue
c
c     update C
      do i = 1,npar
         i1 = i*(i-1)/2
         do j = 1,i
            j1 = i1+j
            anorm(j1) = anorm(j1)-scale(i)*scale(j)*aqa
         enddo
      enddo
c
 100  continue
      return
      end
c-------------------------------------------------------------------------
      subroutine exert_link(npar,l,m,indl,indm,indx,coef,obsc,chi)
c
c     exert several absolute constraints with most elements = 0
c     theory:
c       constraint eqn:                      | xo |
c                        L = A X = ( 0  A2 ) | xc |
c       new solution:
c                        X(new) = X(ori) + C(ori) (0 A2)~ W ax
c                        C(new) = C(ori) - C(ori) (0 A2)~ W C(ori) (0 A2) 
c          where     
c                     ax = L - A2 xc
c
      implicit real*8(a-h,o-z)
      include 'solvem.fti'
c
      integer npar,i,i1,i2,j,indx(m*l),indm(m),ierr,j2
      integer j1,id1,ij,m,l,il,im,indl(l),m1,icount,jm
      dimension coef(m*l),obsc(l),aqa(36),ax(6),aq1(6)
c      
c     calculate ax
      icount = 0
      do il = 1,l
         work = 0.0d0
         m1 = indl(il)
         do im = 1,m1
            i = icount+im
            i1 = indx(i)
            i2 = indm(i1)
            work = work+coef(i)*bnorm(i2)
         enddo
         ax(il) = obsc(il)-work
         do jm = 1,m
            i2 = indm(jm)
            work1 = 0.0d0
            id1 = i2*(i2-1)/2
            do im = 1,m1
               i = icount+im
               j1 = indx(i)
               i1 = indm(j1)
               if (i2.ge.i1) ij = id1+i1
               if (i2.lt.i1) ij = i1*(i1-1)/2+i2 
               work1 = work1+coef(i)*anorm(ij)
            enddo
            gmvm((il-1)*m+jm) = work1
         enddo
         icount = icount+m1
      enddo
c
c     calculate W
      do il = 1,l
         icount = 0
         do jm = 1,il
            work = 0.0d0
            m1 = indl(jm)
            do im = 1,m1
               i = icount+im
               i1 = indx(i)
               work = work+coef(i)*gmvm((il-1)*m+i1)
            enddo
            i2 = il*(il-1)/2+jm
            aqa(i2) = work
            icount = icount+m1
         enddo
      enddo
      call cholsk(aqa,aq1,1,l,ierr) 
      print*, 'ierr=',ierr,'EXERT_LINK'
      if (ierr.ge.120) goto 100
c
c     update chi2
      call axb(l,l,1,aqa,ax,aq1,2,0)
      call axb(1,l,1,aq1,ax,dchi,1,1)
      chi = chi+dchi
c
c     update x 
      do 20 i = 1,npar
         icount = 0
         work1 = 0.0d0
         i1 = i*(i-1)/2
         do 30 il = 1,l
            work = 0.0d0
            m1 = indl(il)
            do jm = 1,m1
               j = icount+jm
               j2 = indx(j)
               j1 = indm(j2)
               if (i.ge.j1) then
                  id1 = i1+j1
               else 
                  id1 = j1*(j1-1)/2+i
               endif
               work = work+coef(j)*anorm(id1)
            enddo
            gvm((i-1)*l+il) = work
            work1 = work1+work*aq1(il)
            icount = icount+m1
 30      continue
         bnorm(i) = bnorm(i)+work1
 20   continue
c
c     update C
      call latwa(npar,l,gvm,aqa,anorm,gmvm,1,1)
c
 100  continue
      return
      end

