      subroutine model_sln(l,m,indx,xs,ys,zs,slxy)
c
c     get model coordinate solution  --- Dong
c
c     We only use velocity model of 2-D translation + 2-D 
c     rotation.  In this case, (Gm~Gm)=I is still valid.
c     legends: ~ transpose;  ^ inverse
c     dimension: I(m,m); G(m,4); Gm(m,4); H(m,m)
c
      implicit real*8(a-h,o-z)
      include 'solvem.fti'
      integer ip,i,j,ix,l,m,lsit,jsit,j1,i1,i2,i3,ii,mn
      integer irow,icol,in,in1,in2,in3,indx(m),j2
      integer indx1(maxprm*2)
      real*8 fj(3),temp(49)
      real*8 bigg(maxnet*12)
      real*8 tfm(4)
      common /tform/tfm
      equivalence (bigg,cova)
      equivalence (indx1,indx_ele)
 
      print*,' Calculate model coordinate solution ...'
c
c     get G and Gm~
      lsit = m/3
      fac = 1.0d0/dble(lsit)
      fac1 = 1.0d0/dsqrt(slxy)
      ip = l*m
      call zero1d(1,ip,gvm)
      facg= 1.0d0/dsqrt(dble(lsit))
c
c     get Gm
      do 30 j = 1,lsit
         j1 = (j-1)*3+1
         j2 = (j-1)*l*3
         if (l.le.3) then
            do i = 1,3
               i1 = (i-1)*l
               gvm(j2+i1+i) = facg
            enddo
            goto 30
         endif
         jsit = (indx(j1)+5)/6
         jsit = itoj(jsit)
         dxj = (x(jsit)-xs)*fac1
         dyj = (y(jsit)-ys)*fac1
         dzj = (z(jsit)-zs)*fac1
         s1 = dsin(slat(jsit))
         c1 = dcos(slat(jsit))
         s2 = dsin(slon(jsit))
         c2 = dcos(slon(jsit))
         fj(1) = -dyj*s1+dzj*s2*c1
         fj(2) = dxj*s1-dzj*c1*c2
         fj(3) = -dxj*s2*c1+dyj*c2*c1
         do i = 1,3
            i1 = (i-1)*l
            gvm(j2+i1+i) = facg
            gvm(j2+i1+4) = fj(i)
         enddo
 30   continue
c
c     get G
      mn = 0
      do 130 j = 1,nsit
         i1 = map((j-1)*6+4)
         i2 = map((j-1)*6+5)
         i3 = map((j-1)*6+6)
         if (i1.le.0.or.i1.gt.nlive) goto 130
         if (i2.le.0.or.i2.gt.nlive) goto 130
         if (i3.le.0.or.i3.gt.nlive) goto 130
         mn = mn+1
         indx1(mn) = i1
         mn = mn+1
         indx1(mn) = i2
         mn = mn+1
         indx1(mn) = i3
         j2 = (mn-3)*l
         if (l.le.3) then
            do i = 1,3
               ix = (i-1)*l
               bigg(j2+ix+i) = facg
            enddo
            goto 130
         endif
         dxj = (x(j)-xs)*fac1
         dyj = (y(j)-ys)*fac1
         dzj = (z(j)-zs)*fac1
         s1 = dsin(slat(j))
         c1 = dcos(slat(j))
         s2 = dsin(slon(j))
         c2 = dcos(slon(j))
         fj(1) = -dyj*s1+dzj*s2*c1
         fj(2) = dxj*s1-dzj*c1*c2
         fj(3) = -dxj*s2*c1+dyj*c2*c1
         do i = 1,3
            ix = (i-1)*l
            bigg(j2+ix+i) = facg
            bigg(j2+ix+4) = fj(i)
         enddo
 130  continue
c
c     V(model) 
      do i = 1,mn
         i1 = (i-1)*l
         work = 0.0d0
         do j = 1,l
            work = work+bigg(i1+j)*tfm(j)
         enddo
         bnorm(indx1(i)) = bnorm(indx1(i))-work
      enddo
c
c     C(model)
      do 80 irow = 1,l
         in1 = (irow-1)*mn
         do 70 i1 = 1,m
            icol = indx(i1)
            ii = icol*(icol-1)/2
            temp1 = 0.0d0
            do 60 in = 1,m
               in2 = (in-1)*l+irow
               i2 = indx(in)
               if (i2.le.icol) then
                  ix = ii+i2
               else
                  in3 = i2*(i2-1)/2
                  ix = in3+icol
               endif
               temp1 = temp1+anorm(ix)*gvm(in2)
 60         continue
            j1 = 0
            do j2 = 1,mn
               if (icol.eq.indx1(j2)) then
                  j1 = j2
                  indx_row(i1) = j2
                  goto 65
               endif
            enddo
            if (j1.eq.0) print*,' index mismatch for icol ',icol
 65         gmvm(in1+j1) = temp1
 70      continue
 80   continue
c
      do 110 irow = 1,l
         in1 = (irow-1)*mn
         ii = irow*(irow-1)/2
         do 100 i1 = 1,irow
            temp1 = 0.0d0
            do 90 in = 1,m
               in2 = (in-1)*l+i1
               i2 = in1+indx_row(in)
               temp1 = temp1+gmvm(i2)*gvm(in2)
 90         continue
            temp(ii+i1) = temp1
 100     continue
 110  continue

      call get_hcht(nlive,mn,l,bigg,gmvm,temp,scale,anorm,indx1)
c
 200  continue
      return
      end
c-------------------------------------------------------------
