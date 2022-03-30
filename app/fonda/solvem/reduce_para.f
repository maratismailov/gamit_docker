      subroutine reduce_para(isit,idx,chi)
c
c     fix parameters, supress covariance matrix and update solution
c
      implicit real*8(a-h,o-z)
      include 'solvem.fti'
      integer id,icon,npar,isit,j,idx,i,ifirst,i1,ij
      dimension idx(6)
      logical ill
c
      npar = nsit*6+iaux
      if (iomode(3).gt.0) npar = npar+iq_sit*3
      i1 = (isit-1)*6
      do 20 j = 1,6
         if (idx(j).le.0) goto 20
         fix(i1+j) = idx(j)
         id = map(i1+j)
         if (id.le.0) goto 20
         ill = .false.
c        copy Q12 
         ifirst = id*(id-1)/2
         do i = 1,id-1
            icon = ifirst+i
            scale(i) = anorm(icon)
         enddo
         xfix = bnorm(id)
         do i = id+1,nlive
            icon = i*(i-1)/2+id
            scale(i-1) = anorm(icon)
            bnorm(i-1) = bnorm(i)
         enddo
c        copy Q22 and inverse it 
         q22 = anorm(id*(id+1)/2)
         if (q22.gt.1.0d-8) then
            q22 = 1.0d0/q22
         else
            print*,' REDUCE_PARA: ill condition:',
     .         sname(isit),isit,id,q22
            q22 = 0
            ill = .true.
         endif
c        compress normal matrix to get Q11
         call trim_matrix(nlive,anorm,id,1)
c        update pointer
         call update_pointer(id,1,map,npar)
         nlive = nlive-1
         if (ill) goto 20
c        update solution
         do 10 i = 1,nlive
            xd = bnorm(i)
            bnorm(i) = xd-scale(i)*q22*xfix
            xd = bnorm(i)-xd
            ifirst = i*(i-1)/2
            do 30 ij = 1,i
               icon = ifirst+ij
               anorm(icon) = anorm(icon)-scale(i)*scale(ij)*q22
 30         continue
            chi = chi+q22*xd**2
 10      continue
 20   continue
c
 100  continue
      return
      end

