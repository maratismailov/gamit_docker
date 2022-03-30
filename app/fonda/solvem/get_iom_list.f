      subroutine get_iom_list(ifile,isnum,indx,xs,ys,zs,soxy,job)
c
c     get site list for inner/outer/model coordinate solution
c
c     job = 1:   inner coordinate site list
c     job = 2:   outer coordinate site list
c     job = 3:   model coordinate site list
c
      implicit real*8(a-h,o-z)
      include 'solvem.fti'
c
      integer job,indx,is1,ifile,j,i1
      integer ib,isnum,is2,is3,isit,ic,i0
      integer match_name,lift_arg,count_arg
      dimension indx(maxsit*18)
      character*128 line
      real*8 tfm(4)
      common /tform/tfm

c     discard junk job
      if (job.lt.1.or.job.ge.4) goto 70
c
      cut_lvl = 0.3d-2
      goto (100,150,110,200),job
c
c     inner coordinate
 100  continue
      call get_iom_site(ifile,isnum,indx,xs,ys,zs,job)
c
      if (isnum.le.0) goto 200
      call geoxyz(radius,finv,tx,ty,tz,a1,a2,a3,xs,ys,zs,2,hght)
      soxy = 0.0d0
      do 20 j = 1,isnum
         ib = iesit(j)
         is1 = map((ib-1)*6+4)
         is2 = map((ib-1)*6+5)
         is3 = map((ib-1)*6+6)
         if (is1.lt.0.or.is1.gt.nlive) goto 20
         if (is2.lt.0.or.is2.gt.nlive) goto 20
         if (is3.lt.0.or.is3.gt.nlive) goto 20
         c1 = anorm(is1*(is1+1)/2)
         c2 = anorm(is2*(is2+1)/2)
         c3 = anorm(is3*(is3+1)/2)
         if (c1+c2+c3.gt.cut_lvl) goto 20
         dx = x(ib)-xs
         dy = y(ib)-ys
         dz = z(ib)-zs
         call sph_ca(a1,a2,dx,dy,dz,de,dn,du,2)
         soxy = soxy+dn**2+de**2
 20   continue
      goto 200
c  
c     outer coordinate
 150  continue
      call get_iom_site(ifile,isnum,indx,xs,ys,zs,job)
c
      if (isnum.le.0) goto 200
      call geoxyz(radius,finv,tx,ty,tz,a1,a2,a3,xs,ys,zs,2,hght)
      aco = dcos(azio)
      asn = dsin(azio)
      soxy = 0.0d0
      do 220 j = 1,isnum
         ib = iesit(j)
         is1 = map((ib-1)*6+4)
         is2 = map((ib-1)*6+5)
         is3 = map((ib-1)*6+6)
         if (is1.lt.0.or.is1.gt.nlive) goto 220
         if (is2.lt.0.or.is2.gt.nlive) goto 220
         if (is3.lt.0.or.is3.gt.nlive) goto 220
         c1 = anorm(is1*(is1+1)/2)
         c2 = anorm(is2*(is2+1)/2)
         c3 = anorm(is3*(is3+1)/2)
         if (c1+c2+c3.gt.cut_lvl) goto 220
         dx = x(ib)-xs
         dy = y(ib)-ys
         dz = z(ib)-zs
         call sph_ca(a1,a2,dx,dy,dz,de,dn,du,2)
         soxy = soxy+(dn*asn)**2+(de*aco)**2-2.0d0*de*dn*asn*aco
 220  continue
      goto 200
c
 110  rewind (ifile)
      call get_iom_site(ifile,isnum,indx,xs,ys,zs,job)
c
      if (isnum.le.0) goto 200
      call geoxyz(radius,finv,tx,ty,tz,a1,a2,a3,xs,ys,zs,2,hght)
      soxy = 0.0d0
      i1 = 0
      isit = 0
      call zero1d(1,4,tfm) 
      rewind (ifile)
      do 130 j = 1,maxsit
         read (ifile,'(a)',end=50) line
         ic = count_arg(line)
         if (ic.le.0) goto 50
         i0 = lift_arg(line,sitnam,6)
         if (i0.le.0) goto 130
         ib = match_name(nsit,i0,sname,sitnam)
         if (ib.le.0) then
            goto 130
         endif
         is1 = map((ib-1)*6+4)
         is2 = map((ib-1)*6+5)
         is3 = map((ib-1)*6+6)
         if (is1.le.0.or.is1.gt.nlive) goto 130
         if (is2.le.0.or.is2.gt.nlive) goto 130
         if (is3.le.0.or.is3.gt.nlive) goto 130
         c1 = anorm(is1*(is1+1)/2)
         c2 = anorm(is2*(is2+1)/2)
         c3 = anorm(is3*(is3+1)/2)
c        skip weak sites to avoid contamination
         if (c1+c2+c3.gt.cut_lvl) goto 130
         read (line,*) vme,vce,vmn,vcn,vcrr
         vme = vme*1.0d-3
         vmn = vmn*1.0d-3
         i1 = i1+1
         indx(i1) = is1
         i1 = i1+1
         indx(i1) = is2
         i1 = i1+1
         indx(i1) = is3
         isit = isit+1
         iesit(isit) = ib
         dx = x(ib)-xs
         dy = y(ib)-ys
         dz = z(ib)-zs
         call sph_ca(a1,a2,dx,dy,dz,de,dn,du,2)
         adjx = bnorm(is1)+vx(ib)
         adjy = bnorm(is2)+vy(ib)
         adjz = bnorm(is3)+vz(ib)
         soxy = soxy+dn**2+de**2
         call sph_ca(slat(ib),slon(ib),adjx,adjy,adjz,dve,dvn,dvu,2)
         call sph_ca(slat(ib),slon(ib),dx,dy,dz,de,dn,du,2)
         bs = dsin(slat(ib))
         bc = dcos(slat(ib))
         bso = dsin(slon(ib))
         bco = dcos(slon(ib))
         tfm(1) = tfm(1)+adjx+vme*bso+vmn*bs*bco
         tfm(2) = tfm(2)+adjy-vme*bco+vmn*bs*bso
         tfm(3) = tfm(3)+adjz-vmn*bc
         tfm(4) = tfm(4)-dn*(dve-vme)+de*(dvn-vmn)
 130  continue
c
 50   fac = 1.0d0/dsqrt(dble(isnum))
      tfm(1) = tfm(1)*fac
      tfm(2) = tfm(2)*fac
      tfm(3) = tfm(3)*fac
      tfm(4) = tfm(4)/dsqrt(soxy)
      print*,' tfm: ',(tfm(j),j=1,4)
      goto 200
c
c     junk list
 70   print*,' unknown site list !!! '
      return
c
 200  continue
c
      return
      end
c------------------------------------------------------------------
      subroutine get_iom_site(ifile,isnum,indx,xs,ys,zs,job)
c
c     get site for inner/outer/model coordinate solution
c
c     job = 1:   inner coordinate site list
c     job = 2:   outer coordinate site list
c     job = 3:   model coordinate site list
c
      implicit real*8(a-h,o-z)
      include 'solvem.fti'
c
      character*128 line
      integer job,indx,i1,i0,ic,is1,ifile,iread,isit,j
      integer ib,isnum,is2,match_name
      integer is3,count_arg,lift_arg
      character*8 sitnam
      dimension indx(maxsit*18)

c     discard junk job
      if (job.le.0.or.job.ge.4) goto 200
c
      cut_lvl = 0.3d-2

      i1 = 0
      xs = 0.0d0
      ys = 0.0d0
      zs = 0.0d0
c     read site name
      isnum = 0
      isit = 0
c     different format
      if (job.eq.3) goto 110
      do 100 iread = 1,maxsit
         read (ifile,'(a)',end=50) line
         ic = count_arg(line)
         if (ic.le.0) goto 100
         do 10 j = 1,ic
            i0 = lift_arg(line,sitnam,j)
            if (i0.le.0) goto 10
            ib = match_name(nsit,i0,sname,sitnam)
            if (ib.le.0) then
               print*,' GET_IOM_SITE name mismatch:',sitnam
               goto 10
            endif
            is1 = map((ib-1)*6+4)
            is2 = map((ib-1)*6+5)
            is3 = map((ib-1)*6+6)
            if (is1.le.0.or.is1.gt.nlive) goto 10
            if (is2.le.0.or.is2.gt.nlive) goto 10
            if (is3.le.0.or.is3.gt.nlive) goto 10
            c1 = anorm(is1*(is1+1)/2)
            c2 = anorm(is2*(is2+1)/2)
            c3 = anorm(is3*(is3+1)/2)
c           skip weak sites to avoid contamination
            if (c1+c2+c3.gt.cut_lvl) goto 10
            i1 = i1+1
            indx(i1) = is1
            i1 = i1+1
            indx(i1) = is2
            i1 = i1+1
            indx(i1) = is3
            xs = xs+x(ib)
            ys = ys+y(ib)
            zs = zs+z(ib)
            isit = isit+1
            iesit(isit) = ib
 10      continue
 100  continue
      goto 50
c
c     model coordinate solution
 110  do 130 iread = 1,maxsit
         read (ifile,'(a)',end=50) line
         ic = count_arg(line)
         if (ic.le.1) goto 50
         i0 = lift_arg(line,sitnam,6)
         if (i0.le.0) goto 130
         ib = match_name(nsit,i0,sname,sitnam)
         if (ib.le.0) then
            print*,' GET_IOM_SITE name mismatch:',sitnam
            goto 130
         endif
         is1 = map((ib-1)*6+4)
         is2 = map((ib-1)*6+5)
         is3 = map((ib-1)*6+6)
         if (is1.le.0.or.is1.gt.nlive) goto 130
         if (is2.le.0.or.is2.gt.nlive) goto 130
         if (is3.le.0.or.is3.gt.nlive) goto 130
         c1 = anorm(is1*(is1+1)/2)
         c2 = anorm(is2*(is2+1)/2)
         c3 = anorm(is3*(is3+1)/2)
c        skip weak sites to avoid contamination
         if (c1+c2+c3.gt.cut_lvl) goto 130
         xs = xs+x(ib)
         ys = ys+y(ib)
         zs = zs+z(ib)
         isit = isit+1
 130  continue
c
 50   isnum = isit
      if (isnum.le.0) goto 200
      xs = xs/dble(isit)
      ys = ys/dble(isit)
      zs = zs/dble(isit)
c
 200  continue
      return
      end
