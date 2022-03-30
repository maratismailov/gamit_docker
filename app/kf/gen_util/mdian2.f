CTITLE MDIAN2
 
      subroutine mdian2(x,n,xmed)

      implicit none 

*     Routine to return median of vector x with out modifying 
*     x.

      integer*4 n   ! Length of x array to use
      real*8 x(n)   ! Array to be checked
      real*8 xmed   ! Median of values in x(i)

*     Parameters used to guide the search for the median
      real*8 big, afac, amp
      parameter (big=1.d30,afac=1.2d0,amp=1.2d0)

*     Declarations
      real*8 eps, a, ap, am, sum, sumx, xp, xm, xx, dum, aa
      real*8 rng
      integer*4 np, nm, j, iter 

*
* MOD TAH:
      xp=big
      xm=-big
      do j = 1, n
         if( x(j).lt.xp ) xp = x(j)
         if( x(j).gt.xm ) xm = x(j)
      end do
      a = 0.5*(xp+xm)
C     a=0.5*(x(1)+x(n))

      eps=abs(x(n)-x(1))
      ap=big
      am=-big
      iter = 0
1     sum=0.
      sumx=0.
      np=0
      nm=0
      xp=big
      xm=-big
      do 11 j=1,n
        xx=x(j)
        if(xx.ne.a)then
          if(xx.gt.a)then
            np=np+1
            if(xx.lt.xp)xp=xx
          else if(xx.lt.a)then
            nm=nm+1
            if(xx.gt.xm)xm=xx
          endif
          dum=1./(eps+abs(xx-a))
          sum=sum+dum
          sumx=sumx+xx*dum
        endif
11    continue
      rng = xp - xm
      if( rng.eq.0 ) then 
          xmed = xp
          RETURN
      end if
      if(np-nm.ge.2 .and. abs(eps/rng).gt. 1.d-10 )then
        am=a
        aa=xp+max(0.d0,sumx/sum-a)*amp
        if(aa.gt.ap)aa=0.5*(a+ap)
        eps=afac*abs(aa-a)
        a=aa
        if( a.gt.xp ) a = xp
        iter = iter + 1
        if( iter.gt.n*n ) eps = 0
        go to 1
      else if(nm-np.ge.2 .and. abs(eps/rng).gt. 1.d-10 )then
        ap=a
        aa=xm+min(0.d0,sumx/sum-a)*amp
        if(aa.lt.am)aa=0.5*(a+am)
        eps=afac*abs(aa-a)
        a=aa
        if( a.lt.xm ) a = xm
        iter = iter + 1
        if( iter.gt.n*10) eps = 0
        go to 1
      else
        if( abs(eps/rng).gt. 1.d-10 ) then
           if(mod(n,2).eq.0)then
             if(np.eq.nm)then
               xmed=0.5*(xp+xm)
             else if(np.gt.nm)then
               xmed=0.5*(a+xp)
             else
               xmed=0.5*(xm+a)
             endif
           else
             if(np.eq.nm)then
               xmed=a
             else if(np.gt.nm)then
               xmed=xp
             else
               xmed=xm
             endif
           endif
        else
           xmed = a
        end if
      endif
      return
      end
