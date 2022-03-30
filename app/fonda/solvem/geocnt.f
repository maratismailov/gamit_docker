      subroutine geocnt(xs,ys,zs,slxy,slxz,slyz,so2y,mode)
c
c     get geometric center coordinate and sum of baseline length.
c     All calculations are in the geocentric coordinate.
c     Usually the outer coordinates are related to local geodetic
c     system, the normalized factor(so2y) is an invariant, therefore
c     we can calculate so2y in local geodetic sustem.
c
c     mode = 1: without excluded sites
c
      implicit real*8(a-h,o-z)
      include 'solvem.fti'
c
      integer mode,ntemp,igd,i
      integer itemp(maxnet)
c     
c     sort out the sites which are included in the network
      if (mode.eq.1) then
         ntemp = 0
         do 30 igd = 1,gdsit
            i = itoj(igd)
            itemp(igd) = iexc(i)
            if (itemp(igd).eq.1) ntemp = ntemp+1
 30      continue
      endif
      print*,' get geocenter .. effective sites =',gdsit-ntemp,gdsit
c
c     calculate geometric center's geocentric coordinate
      sumx = 0.0d0
      sumy = 0.0d0
      sumz = 0.0d0
      do 20 igd = 1,gdsit
         i = itoj(igd)
         if (itemp(igd).eq.1) goto 20
         sumx = sumx+x(i)
         sumy = sumy+y(i)
         sumz = sumz+z(i)
 20   continue
      xs = sumx/dble(gdsit-ntemp)
      ys = sumy/dble(gdsit-ntemp)
      zs = sumz/dble(gdsit-ntemp)
c
c     get horizontal geometric center's geodetic coordinate
c     (for outer coordinate constraint)
      call geoxyz(radius,finv,tx,ty,tz,a1,a2,a3,xs,ys,zs,2,hght)
c
c     get normalization factor for inner coordinate
c     get normalization factor for outer coordinate
      so2x = 0.0d0
      so2y = 0.0d0
      sxy = 0.0d0
      slxy = 0.0d0
      slxz = 0.0d0
      slyz = 0.0d0
      do 60 igd = 1,gdsit
         i = itoj(igd)
         if (itemp(igd).eq.1) goto 60
         dx = x(i)-xs
         dy = y(i)-ys
         dz = z(i)-zs
c         slxy = slxy+dx**2+dy**2
         slxz = slxz+dx**2+dz**2
         slyz = slyz+dz**2+dy**2
         call sph_ca(a1,a2,dx,dy,dz,de,dn,du,2)
         so2x = so2x+de**2
         so2y = so2y+dn**2
         sxy = sxy+de*dn
 60   continue
      slxy = so2x+so2y
c
      as = dsin(azio)
      ac = dcos(azio)
      stemp = so2x*ac*ac+so2y*as*as-2.d0*sxy*as*ac
      so2y = stemp
c
 100  continue
      return
      end

