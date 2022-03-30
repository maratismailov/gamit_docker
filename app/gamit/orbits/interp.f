	subroutine interp(iref,i2nd,nrec,nprn,iprn,tag,orb,delt)
c
c purpose:	interpolate x/y/z values of the second orbit in
c		order to match up the reference orbit time tag
c		the interpolated values will be stored in orb(3,*,*,*)
c		also interpolate x/y/z values of the reference orbit
c		at each time tag + 1sec in order to transform the
c		dx/dy/dz to along track, cross track, and radial
c		the interpolated values will be stored in orb(4,*,*,*)
c
c note:		maxepc is set to the max number of epoches of 7 day at
c		15 min interval
c
      integer*4 maxsat,maxepc
        parameter (maxsat=24,maxepc=673)
	integer nrec(2),nprn(2),iprn(2,85)
      integer*4 iref,i2nd,ii,j,i
	real*8 t(maxepc),x(maxepc),y(maxepc),z(maxepc),tmp(maxepc),as
	real*8 tag(2,maxepc),orb(4,maxsat,maxepc,3),delt(2)
c
c do the interpolation on the second orbit only if time intervals are
c	different, or starting times are different

	if (delt(1).ne.delt(2).or.tag(1,1).ne.tag(2,1)) then

c loop over all satellites in ref orbit
	do i=1,nprn(iref)
	ii=1
c advance to next satellite until both ref. and second satellites match
	do while (iprn(i2nd,ii).ne.iprn(iref,i).and.ii.le.nrec(i2nd))
		ii=ii+1
	enddo
c do natural spline interpolation only when ref. and second sat's match
	if (iprn(iref,i).eq.iprn(i2nd,ii)) then
c prapare the time, x, y, z arrays
	do j=1,nrec(i2nd)
	t(j)=tag(i2nd,j)
	x(j)=orb(i2nd,ii,j,1)
	y(j)=orb(i2nd,ii,j,2)
	z(j)=orb(i2nd,ii,j,3)
	enddo

c do interpolation for x
	call spline(t,x,nrec(i2nd),1.e30,1.e30,tmp)
	do j=1,nrec(iref)
	call splint(t,x,tmp,nrec(i2nd),tag(iref,j),orb(3,i,j,1))
	enddo
c do interpolation for y
	call spline(t,y,nrec(i2nd),1.e30,1.e30,tmp)
	do j=1,nrec(iref)
	call splint(t,y,tmp,nrec(i2nd),tag(iref,j),orb(3,i,j,2))
	enddo
c do interpolation for z
	call spline(t,z,nrec(i2nd),1.e30,1.e30,tmp)
	do j=1,nrec(iref)
	call splint(t,z,tmp,nrec(i2nd),tag(iref,j),orb(3,i,j,3))
	enddo
	endif

c advance to next satellite of reference orbit
	enddo

	endif

c convert 1 milsec of time in terms of day
	as=1./68400/1000
c loop over each satellite
	do i=1,nprn(iref)

c do natural spline interpolation for each epoch+1sec

c prapare the time, x, y, z arrays
	do j=1,nrec(iref)
	t(j)=tag(iref,j)
	x(j)=orb(iref,i,j,1)
	y(j)=orb(iref,i,j,2)
	z(j)=orb(iref,i,j,3)
	enddo
c do interpolation for x
	call spline(t,x,nrec(iref),1.e30,1.e30,tmp)
	do j=1,nrec(iref)
	call splint(t,x,tmp,nrec(iref),tag(iref,j)+as,orb(4,i,j,1))
	enddo
c do interpolation for y
	call spline(t,y,nrec(iref),1.e30,1.e30,tmp)
	do j=1,nrec(iref)
	call splint(t,y,tmp,nrec(iref),tag(iref,j)+as,orb(4,i,j,2))
	enddo
c do interpolation for z
	call spline(t,z,nrec(iref),1.e30,1.e30,tmp)
	do j=1,nrec(iref)
	call splint(t,z,tmp,nrec(iref),tag(iref,j)+as,orb(4,i,j,3))
	enddo
c advance to next satellite of reference orbit
	enddo

	return
	end

      subroutine spline(x,y,n,yp1,ypn,y2)
      implicit real*8 (a-h,o-z)
      integer*4 nmax,i,k,n
      parameter (nmax=336)
      dimension x(n),y(n),y2(n),u(nmax)
      if (yp1.gt..99e30) then
        y2(1)=0.
        u(1)=0.
      else
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     *      /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
11    continue
      if (ypn.gt..99e30) then
        qn=0.
        un=0.
      else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      return
      end

      subroutine splint(xa,ya,y2a,n,x,y)
      implicit real*8 (a-h,o-z)
      integer*4 klo,khi,k,n
      dimension xa(n),ya(n),y2a(n)
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.) pause 'bad xa input.'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+
     *      ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
      return
      end
