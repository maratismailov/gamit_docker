      subroutine inner_sln(l,m,indx,xs,ys,zs,slxy)
c
c     get inner coordinate solution  --- Dong
c
c     Unfourtunately, in 3-D case (Gm~Gm)=I is no longer valid
c     We have to calcuculate (I-G(Gm~Gm)^Gm~) in stead of (I-GGm~)
c     But currently we focus our attention to horizontal velocity
c     only, therefore, we only constrain 3-D translation + 2-D 
c     rotation.  In this case, (Gm~Gm)=I is still valid.
c     legends: ~ transpose;  ^ inverse
c     dimension: I(m,m); G(m,4); Gm(m,4); H(m,m)
c
      implicit real*8(a-h,o-z)
      include 'solvem.fti'
      integer ip,i,j,ix,l,m,lsit,jsit,j1,i1,i2,i3,ii,mn
      integer irow,icol,in,in1,in2,in3,indx(m),j2
      integer indx1(maxprm*2)
      real*8 fi(6),fj(3),temp(49)
      real*8 bigg(maxnet*12)
      equivalence (bigg,cova)
      equivalence (indx1,indx_ele)
 
      print*,' Calculate inner coordinate solution ...'
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
c     V(inner) 
      do 50 i = 1,l
         work = 0.0d0
         do 40 j = 1,m
            work = work+gvm((j-1)*l+i)*bnorm(indx(j))
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
c     C(inner)
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
      subroutine axb_indx(l,n,m,a,b,c,indx,modea,modeb,type)
c
c     modea = 0: full matrix for array a
c     modea = 1: symmetric matrix for array a
c     modeb = 0: full matrix for array b
c     modeb = 1: symmetric matrix for array b
c     type  = 1: array a uses index.
c     type  = 2: array b uses index.
c
      integer l,n,m,irow,icol,in1,in2,ix,in3,i2,in
      integer indx(n),modea,modeb,ib,type,ia
      real*8 a(l*(n+modea)/(modea+1)),c(l*m),temp
      real*8 b(m*(n+modeb)/(modeb+1))

      do 80 ia = 1,l
         if (type.eq.1) then
            irow = indx(ia)
         else
            irow = ia
         endif
         if (modea.eq.1) then
            in1 = irow*(irow-1)/2
         else
            in1 = (irow-1)*n
         endif
         do 70 ib = 1,m
            if(type.eq.2) then
               icol = indx(ib)
            else
               icol = ib
            endif
            temp = 0.0d0
            do 60 in = 1,n
               i2 = indx(in)
               if (type.eq.1) then
                  in2 = in1+i2
                  if (modea.eq.1.and.i2.gt.irow)
     .               in2 = i2*(i2-1)/2+irow
               else
                  in2 = in1+in
               endif
               if (type.eq.2) then
                  if (i2.le.icol) then
                     ix = icol*(icol-1)/2+i2
                  else
                     in3 = i2*(i2-1)/2
                     ix = in3+icol
                  endif
               else
                  ix = (in-1)*m+icol
               endif
               temp = temp+a(in2)*b(ix)
 60         continue
            c((ia-1)*m+ib) = temp
 70      continue
 80   continue

      return
      end

c-------------------------------------------------------------
      subroutine get_hcht(m,l,n,avec,bvec,tvec,tmp,cvec,
     .   indx)
c	
c     matrix multiplication C(m,m) = (I - H) T (I - H)~
c
      implicit real*8(a-h,o-z)
      integer ic,l,n,m,in,ii,i1,i2,icol
      integer irow,ir,in1,in2,in3,kn,ic1,ir1
      integer indx(l)
      dimension avec(l*n),bvec(n*l),cvec(m*(m+1)/2)
      dimension tvec(n*(n+1)/2),tmp(n)
  
      do 30 irow = 1,l
         in1 = (irow-1)*n
         do 20 icol = 1,n
            ii = icol*(icol-1)/2
            temp1 = 0.0d0
            do 10 in = 1,n
               in2 = in1+in
               if (in.le.icol) then
                  kn = ii+in
               else
                  in3 = in*(in-1)/2
                  kn = in3+icol
               endif
               temp1 = temp1+tvec(kn)*avec(in2)
 10         continue
            tmp(icol) = temp1
 20   continue

         do 40 ir = irow,l
            ir1 = (ir-1)*n
            temp1 = 0.0d0
            temp2 = 0.0d0
            do 50 ic = 1,n
               ic1 = ir1+ic
               i1 = (ic-1)*l+irow
               temp1 = temp1+tmp(ic)*avec(ic1)
               temp2 = temp2+bvec(i1)*avec(ic1)
               i1 = (ic-1)*l+ir
               temp2 = temp2+bvec(i1)*avec(in1+ic)
 50         continue
            i1 = indx(ir)
            i2 = indx(irow)
            if (i1.ge.i2) ii = i1*(i1-1)/2+i2
            if (i1.lt.i2) ii = i2*(i2-1)/2+i1
            cvec(ii) = cvec(ii)+temp1-temp2
 40      continue

 30   continue

      return
      end
