      subroutine seqsln(n,m,xs,c1,sl2,a2,cl2,chi,tmp1,tmp2,
     .                  tmp3,tmp4,mode)
c
c     general form for sequential least square solution
c     suppose original solution is x(1) with covariance C1,
c     the condition equations are l2 = A2 x  with covariance Cl2,
c     then the updated solutions are:
c        x(1+2) = x(1) + C1 A2~(Cl2 + A2 C1 A2~)^ (l2 - A2 x(1))
c        C(1+2) = C1 - C1 A2~(Cl2 + A2 C1 A2~)^ A2 C1
c        chi(1+2) = chi(1)+dx~C1^dx+(l2-A2x(1+2))~Cl2^(l2-A2x(1+2))
c                 = chi(1) + (l2 - A2 x(1))~ Cl2^ (l2 - A2 x(1+2))
c     where
c        ~ represents transpose
c        ^ represents inverse
c
c     dimension:
c        x:      (n x 1) vector
c        l2:     (m x 1) vector
c        A2:     (m x n) matrix
c        C1:     (n x n) symmetric matrix
c        Cl2:    (m x m) symmetric matrix
c        chi:    scalar
c        x and C1 will be overwritted by x(1+2) and C(1+2)
c        Cl2 will be overwritted by (Cl2 + A2 C1 A2~)^
c        chi will be updated
c        sl2 will become residuals (not finished)
c
      implicit real*8(a-h,o-z)
      dimension xs(n),sl2(m),c1(n*(n+1)/2),a2(m*n),cl2(m*(m+1)/2)
      dimension tmp1(n*m),tmp2(m),tmp3(m),tmp4(n)
      integer n,m,mode,in,im,ik,im0,id,i0,i1,i2,ier
c
c     calculate C1 A2~ (n x m)
      do 30 in = 1,n
         i0 = in*(in-1)/2
         do 20 im = 1,m
            id = (in-1)*m+im
            work = 0.0d0
            im0 = (im-1)*n
            do 10 ik = 1,n
               i1 = i0+ik
               if (ik.gt.in) i1 = ik*(ik-1)/2+in
               i2 = im0+ik
               work = work+c1(i1)*a2(i2)
 10         continue
            tmp1(id) = work
 20      continue
 30   continue
c
c     calculate l2 - A2 x(1)  (m x 1)
      do 80 im = 1,m
         work = 0.0d0
         i0 = (im-1)*n
         do 70 in = 1,n
            work = work+a2(i0+in)*xs(in)
 70      continue
         sl2(im) = sl2(im)-work
 80   continue
c
c     calculate (l2 - A2 x(1))~ Cl2  (1 x m)
      do 50 im = 1,m
         work = 0.0d0
         i0 = im*(im-1)/2
         do 40 ik = 1,m
            id = i0+ik
            if (ik.gt.im) id = ik*(ik-1)/2+im
            work = work+sl2(ik)*cl2(id)
 40      continue
         tmp3(im) = work
 50   continue
         do im = 1,m
           print*,' before: ',(cl2(ik),ik=1,im)
         enddo
c
c     calculate (Cl2 + A2 C1 A2~) (m x m)
      call latwa(m,n,a2,c1,cl2,tmp2,1,2)
         do im = 1,m
           print*,' after: ',(cl2(ik),ik=1,im)
         enddo
c      do 60 im = 1,m
c         i0 = im*(im-1)/2
c         i10 = (im-1)*n
c         do 50 in = 1,n
c            id = i0+in
c            if (in.gt.im) id = in*(in-1)/2+im
c            i20 = (in-1)*m
c            work = 0.0d0
c            do 40 ik = 1,m
c               i1 = i10+ik
c               i2 = i20+ik
c               work = work+a2(i1)*tmp1(i2)
c 40         continue
c            cl2(id) = cl2(id)+work
c 50      continue
c 60   continue
c     
      print*,'before cl2: ',cl2(1)
c     calculate (Cl2 + A2 C1 A2~)^ (m x m)
      call cholsk(cl2,tmp2,1,m,ier)
      print*,'after cl2: ',cl2(1)
c     
c     calculate (Cl2 + A2 C1 A2~)^ (l2 - A2 x(1))  (m x 1)
      do 90 im = 1,m
         work = 0.0d0
         i0 = (im-1)*im/2
         do 100 in = 1,m
            id = i0+in
            if (in.gt.im) id = in*(in-1)/2+im
            work = work+cl2(id)*sl2(in)
 100     continue
         tmp2(im) = work
 90   continue
c     
c     get the updated solution x(1+2) and chi
      do 120 in = 1,n
         work = 0.0d0
         i0 = (in-1)*m
         do 110 im = 1,m
            id = i0+im
            work = work+tmp1(id)*tmp2(im)
 110     continue
         xs(in) = xs(in)+work
         tmp4(in) = work
 120  continue
c
c     update chi
      work1 = 0.0d0
      do 55 im = 1,m
         work = 0.0d0
         i0 = (im-1)*n
         do 60 in = 1,n
            work = work+a2(i0+in)*tmp4(in)
 60      continue
         work1 = work1+tmp3(im)*(sl2(im)-work)
 55   continue
      chi = chi+work1
c     
c     get the updated covariance C(1+2)
      call latwa(n,m,tmp1,cl2,c1,tmp4,1,1)
c
      return
      end

