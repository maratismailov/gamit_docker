      subroutine sln_update(n,c1,sl2,cl2,chi,tmp1,tmp2,mode)
c
c     least square solution with apriori covariance matrix
c     It is equivalent to sequential least square or Kalman filter
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
c
c     In this subroutine, we only use simplified form.  The general form
c     is performed by subroutine seqsln.f.
c     Simplification:
c        A2 = I
c        x(1+2) - x(1) = C1 (Cl2 + C1)^ (l2 - A2 x(1))
c        C(1+2) = C1 - C1 (Cl2 + C1)^ C1
c        chi(1+2) = chi(1) + (l2 - A2 x(1))~ C1^ (x(1+2) - x(1))
c                 = chi(1) + (l2 - A2 x(1))~ (Cl2 + C1)^ (l2 - A2 x(1))
c        Cl2 will be overwritted by C(1+2)
c        chi will be updated
c        sl2 will be updated by x(1+2)
c
      real*8 c1,sl2,cl2,chi,tmp1,work,tmp2
      dimension sl2(n),c1(n*(n+1)/2),cl2(n*(n+1)/2)
      dimension tmp1(n),tmp2(n*(n+1)/2)
      integer n,mode,in,im,id,i0,ier,ip
c
c     calculate (C1 + Cl2) (n x n)
      do 30 in = 1,n
         i0 = in*(in-1)/2
         do 20 im = 1,in
            id = i0+im
            cl2(id) = cl2(id)+c1(id)
 20      continue
 30   continue
      print*,' inverse (C(apr)+C(obs)) ...'
c
c     calculate (C1 + Cl2)^ (n x n)
      call nrmscl(cl2,sl2,tmp1,n,1,0)
      call cholsk(cl2,sl2,1,n,ier)
      call nrmscl(cl2,sl2,tmp1,n,2,0)
c     
c     calculate (Cl2 + C1)^ (l2 - x(1))  (n x 1)
      do 90 im = 1,n
         work = 0.0d0
         i0 = (im-1)*im/2
         do 100 in = 1,n
            id = i0+in
            if (in.gt.im) id = in*(in-1)/2+im
            work = work+cl2(id)*sl2(in)
 100     continue
         tmp1(im) = work
 90   continue
      print*,' update solution and chi2 ...'
c     
c     get the updated solution x(1+2) and chi
      do 120 in = 1,n
         work = 0.0d0
         i0 = (in-1)*in/2
         do 110 im = 1,n
            id = i0+im
            if (im.gt.in) id = im*(im-1)/2+in
            work = work+c1(id)*tmp1(im)
 110     continue
         chi = chi+sl2(in)*tmp1(in)
         sl2(in) = work
 120  continue
      print*,' update covariance matrix ...'
c     
c     get the updated covariance C(1+2)
c     C(1+2) = C1 - C1 (Cl2 + C1)^ C1
      call satwa(n,c1,cl2,tmp2,tmp1,1)
      id = n*(n+1)/2
      do ip = 1,id
         cl2(ip) = c1(ip)-tmp2(ip)
      enddo
c
      return
      end

