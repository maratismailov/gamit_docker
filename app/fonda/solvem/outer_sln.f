      subroutine outer_sln(l,m,indx,xs,ys,zs,slxy)
c
c     get outer coordinate solution  -- Dong
c
c     1. The transform matrix in this case is no longer symmetric.
c        We have to use full starage for the transform matrix H
c     2. I don't like the rotation approach.  Therefore I use the deducted
c        formula for minimum velosity with azimuth azio.
c     3. In geocentric coordinate and 3-D case, I don't see the 
c        rigorous definition of 'outer' coordinate.  Current definition
c        just fits a regional network and use plate coordinate app-
c        roximation.  Therefore I use approximate formula:
c        constrainted direction by the azimuth angle defined in the 
c        local topocentric coordinate from the network center.
c     4. If we constrain our attention to zero translation and minimizing
c        velocity in one horizontal direction, (Gm~ Gm) = I still
c        valid.
c
      implicit real*8(a-h,o-z)
      include 'solvem.fti'
c
      integer ip,i,j,ix,l,m,lsit,jsit,j1,i1,i2,i3,ii,mn
      integer irow,icol,in,in1,in2,in3,indx(m),j2
      integer indx1(maxprm*2)
      real*8 cjaco(9),fi(6),fj(3),coef(3),temp(49)
      real*8 bigg(maxnet*12)
      equivalence (bigg,cova)
      equivalence (indx1,indx_ele)
c
      print*,' Calculate outer coordinate solution ...'
c
      lsit = m/3
      ip = m*l
      call zero1d(1,ip,gvm)
      call zero1d(1,nsit*l,gmvm)
      call zero1d(1,nsit*l,bigg)

      fac = 1.0d0/dble(lsit)
      fac1 = 1.0d0/dsqrt(slxy)
      facg= 1.0d0/dsqrt(dble(lsit))
      aco = dcos(azio)
      asn = dsin(azio)
c     geodetic coordinate of the network center
      call geoxyz(radius,finv,tx,ty,tz,a1,a2,a3,xs,ys,zs,2,hght)
      call getjac(a2,a1,hght,cjaco,5)
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
         dxj = (x(jsit)-xs)
         dyj = (y(jsit)-ys)
         dzj = (z(jsit)-zs)
         call sph_ca(a1,a2,dxj,dyj,dzj,de,dn,du,2)
         fj(1) = (dn*asn-de*aco)*asn*fac1
         fj(2) = (dn*asn-de*aco)*aco*fac1
         fj(3) = 0.0d0
         call axb(1,3,3,fj,cjaco,coef,1,0)
         do i = 1,3
            i1 = (i-1)*l
            gvm(j2+i1+i) = facg
            gvm(j2+i1+4) = coef(i)
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
         c1 = anorm(i1*(i1+1)/2)
         c2 = anorm(i2*(i2+1)/2)
         c3 = anorm(i3*(i3+1)/2)
         if (c1.le.0.0d0) goto 130
         if (c2.le.0.0d0) goto 130
         if (c3.le.0.0d0) goto 130
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
         dxj = (x(j)-xs)
         dyj = (y(j)-ys)
         dzj = (z(j)-zs)
         call sph_ca(a1,a2,dxj,dyj,dzj,de,dn,du,2)
         fi(1) = (dn*cjaco(1)-de*cjaco(4))*fac1
         fi(2) = (dn*cjaco(2)-de*cjaco(5))*fac1
         fi(3) = (dn*cjaco(3)-de*cjaco(6))*fac1
         do i = 1,3
            ix = (i-1)*l
            bigg(j2+ix+i) = facg
            bigg(j2+ix+4) = fi(i)
         enddo
 130  continue
c
c     V(outer) 
      do 50 i = 1,l
         work = 0.0d0
         do 40 j = 1,m
            jsit = (indx(j)+5)/6
            jsit = itoj(jsit)
            j1 = j-j/3*3
            if (j1.eq.1) vap = vx(jsit)
            if (j1.eq.2) vap = vy(jsit)
            if (j1.eq.0) vap = vz(jsit)
            work = work+gvm((j-1)*l+i)*(bnorm(indx(j))+vap)
 40      continue
         fi(i) = work
 50   continue
      do i = 1,mn
         i1 = (i-1)*l
         work = 0.0d0
         do j = 1,l
            work = work+bigg(i1+j)*fi(j)
         enddo
         bnorm(indx1(i)) = bnorm(indx1(i))-work
      enddo
c
c     C(outer)
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
c
      return
      end

