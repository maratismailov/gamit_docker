      subroutine put_link(wcmd,lcmd,chi,job)
c
c     perform link constraint
c
c     job = 1:   link multiple sites
c     job = 2:   link baseline direction
c     job = 3:   link baseline length
c     job = 4:   link velocity gradient
c     job = 5:   link episodic parameters
c
c     type = 2:  link adjusted values
c     type = 3:  link total values (only for velocity)
c     type = 4:  
c
      implicit real*8(a-h,o-z)
      include 'solvem.fti'
c
      character*(*) wcmd
      character*10 code
      character*16 title
      integer job,indx,i1,i0,ic,is1,icol,mc,type,mdm
      integer ib,isnum,is2,i,match_name,i2,is3,imatch
      integer lcmd,count_arg,lift_arg,lmode(6),sit_id(20)
      dimension indx(12),coef(12),cjaco(9),tmp(9),tmp1(9)
c      dimension cobs(9)
      logical linkx,linkv

c  
c     decomposite command line
      ic = count_arg(wcmd)
c     not enough arguments
      if (ic.le.1) goto 200
c     pointer to coordinate and velocity
      ib = lift_arg(wcmd,code,1)
      if (ib.le.0) goto 200
      do i = 1,6
         lmode(i) = 0
      enddo
      do i = 1,ib
         if (code(i:i).eq.'x') lmode(1) = 1
         if (code(i:i).eq.'y') lmode(2) = 1
         if (code(i:i).eq.'z') lmode(3) = 1
         if (code(i:i).eq.'u') lmode(4) = 1
         if (code(i:i).eq.'v') lmode(5) = 1
         if (code(i:i).eq.'w') lmode(6) = 1
      enddo
      linkx = .false.
      linkv = .false.
      if (lmode(1)+lmode(2)+lmode(3).gt.0) linkx = .true.
      if (lmode(4)+lmode(5)+lmode(6).gt.0) linkv = .true.
c
c     set the "hardness" of the link
      ib = lift_arg(wcmd,code,2)
      read (code,*,err=999) sigm
999   continue
c     print*, 'code ',code
c     print*, 'sigm= ',sigm
      sigm= sigm/1.0d+3
c     print*, 'sigm= ',sigm
      sigm = sigm**2
c
c     identify link type
      ib = lift_arg(wcmd,code,3)
      type = 2
      if (code(1:3).eq.'tot') type = 3
      if (job.eq.5) then
         type = 3
         read (code,*) e_time
      endif
      if (job.eq.4) then
         type = 4
         mdm = 1
         if (code.eq.'dvxdx') mdm = 1
         if (code.eq.'dvxdy') mdm = 2
         if (code.eq.'dvydx') mdm = 3
         if (code.eq.'dvydy') mdm = 4
         ib = lift_arg(wcmd,code,4)
         read (code,*) cita
         cita = cita*dtor
      endif
c
c     read site name and match their index
      i1 = 0
      do 10 i = type+1,ic
         ib = lift_arg(wcmd,code,i)
         if (ib.le.0) goto 10
         i1 = i1+1
         i0 = match_name(nsit,ib,sname,code)
         if (i0.le.0) then
            print*,' PUT_LINK name mismatch: ',code
            i1 = i1-1
            goto 10
         endif
         sit_id(i1) = i0
 10   continue
      if (i1.le.1) goto 200
      isnum = i1
c
      goto (20,50,100,150,250,200) job
c
c     force coordinate or velocity to be the same
 20   continue
      i0 = sit_id(1)
c     ---- not regorous, modify it later
      is1 = map(i0*6-5)-1
      if (is1.lt.0) print*,' PUT_LINK: missing site ',sname(i0)
      if (is1.lt.0) goto 200
      coef(1) = 1.0d0
      coef(2) = -1.0d0
      do 40 ib = 2,isnum
         i1 = sit_id(ib)
         is2 = map(i1*6-5)-1
         if (is2.lt.0) goto 40
c         print*,' --- force site & index:',jtoi(i0),jtoi(i1),is1,is2
         chio = chi
         if (id_frame.eq.1) then
            do 30 i = 1,6
               if (lmode(i).ge.1) then
                  indx(1) = is1+i
                  indx(2) = is2+i
                  obsc = 0.0d0
                  if (type.ge.2.and.i.gt.3) then
                     if (i.eq.4) obsc = vx(i1)-vx(i0)
                     if (i.eq.5) obsc = vy(i1)-vy(i0)
                     if (i.eq.6) obsc = vz(i1)-vz(i0)
                  endif
                  call exert_link1(nlive,2,indx,coef,obsc,chi,sigm)
               endif
 30         continue
         else
            sl = dsin(slon(i0))
            cl = dcos(slon(i0))
            sf = dsin(slat(i0))
            cf = dcos(slat(i0))
            do i = 1,3
               indx(i) = is1+i
               indx(i+3) = is2+i
            enddo
            do 35 i = 1,6
               if (lmode(i).lt.1) goto 35
               if (i.gt.3) then
                  do i2 = 1,3
                     indx(i2) = is1+i2+3
                     indx(i2+3) = is2+i2+3
                  enddo
               endif
               if (i.eq.1.or.i.eq.4) then
                  coef(1) = -sl
                  coef(2) =  cl
                  coef(3) = 0.0d0
                  coef(4) = -coef(1)
                  coef(5) = -coef(2)
                  coef(6) = -coef(3)
               endif
               if (i.eq.2.or.i.eq.5) then
                  coef(1) = -cl*sf
                  coef(2) = -sl*sf
                  coef(3) = cf
                  coef(4) = -coef(1)
                  coef(5) = -coef(2)
                  coef(6) = -coef(3)
               endif
               if (i.eq.3.or.i.eq.6) then
                  coef(1) = cl*cf
                  coef(2) = sl*cf
                  coef(3) = sf
                  coef(4) = -coef(1)
                  coef(5) = -coef(2)
                  coef(6) = -coef(3)
               endif
               obsc = 0.0d0
               if (type.ge.2.and.i.gt.3) then
                  if (i.eq.4) obsc = ve(i1)-ve(i0)
                  if (i.eq.5) obsc = vn(i1)-vn(i0)
                  if (i.eq.6) obsc = vu(i1)-vu(i0)
               endif
               call exert_link1(nlive,6,indx,coef,obsc,chi,sigm)
 35         continue
         endif
         write(6,'(a,a8,a2,a8,f18.4,f12.4)') 
     .      ' LINK (site) chi2, dchi2: ',
     .      sname(i0),'--',sname(i1),chi,chi-chio
         if (iomode(6).gt.0) then
            title = ' link site  ,  :'
            if (linkx) title(12:12) = 'x'
            if (linkv) title(14:14) = 'v'
            write(20,'(a16,2x,a8,2x,a8,10x,f18.4,f12.4)') 
     .      title,
     .      sname(i0),sname(i1),chi,chi-chio
         endif
 40   continue
      goto 200
c
c     constrain baseline direction
 50   continue
      i0 = sit_id(1)
      is1 = map(i0*6-5)-1
      i1 = sit_id(2)
      is2 = map(i1*6-5)-1
      if (is1.lt.0.and.is2.lt.0) 
     .   print*,' PUT_LINK bad index: ',i0,i1,is1,is2
      if (is1.lt.0.and.is2.lt.0) goto 200
      sl = dsin(slon(i0))
      cl = dcos(slon(i0))
      sf = dsin(slat(i0))
      cf = dcos(slat(i0))
      fac1 = -(x(i1)-x(i0))*cl*sf-(y(i1)-y(i0))*sl*sf
     .       +(z(i1)-z(i0))*cf
      fac1 = fac1/(-(x(i1)-x(i0))*sl+(y(i1)-y(i0))*cl)
      icol = 6
      mc = 3
      if (is1.lt.0.or.is2.lt.0) icol = 3
      coef(1) = sf*cl-fac1*sl
      coef(2) = sf*sl+fac1*cl
      coef(3) = -cf
c
      if (lmode(1).ge.1.and.lmode(2).ge.1) then
c        constrain dn - k1*de = 0
         obsc = 0.0d0
         chio = chi
         do i = 1,mc
            indx(i) = is1+i
            indx(i+mc) = is2+i
            coef(i+mc) = -coef(i)
         enddo
         if (is1.lt.0) then
            do i = 1,mc
               indx(i) = indx(i+mc)
               coef(i) = coef(i+mc)
            enddo
         endif
         call exert_link1(nlive,icol,indx,coef,obsc,chi,sigm)
         write(6,'(a,a8,a2,a8,f19.4,f12.4)') 
     .      ' LINK (x dir) chi2, dchi2: ',
     .      sname(i0),'--',sname(i1),chi,chi-chio
         if (iomode(6).gt.0)
     .      write(20,'(a16,2x,a8,2x,a8,10x,f18.4,f12.4)') 
     .      ' link site cdir:',
     .      sname(i0),sname(i1),chi,chi-chio
      endif
c
      if (lmode(4).ge.1.and.lmode(5).ge.1) then
c        constrain dVn - k1*dVe = 0
         obsc = 0.0d0
         if (type.ge.2) obsc = vn(i0)-vn(i1)-fac1*(ve(i0)-ve(i1))
         chio = chi
         do i = 1,mc
            indx(i) = is1+i+3
            indx(i+mc) = is2+i+3
            coef(i+mc) = -coef(i)
         enddo
         if (is1.lt.0) then
            do i = 1,mc
               indx(i) = indx(i+mc)
            enddo
         endif
         call exert_link1(nlive,icol,indx,coef,obsc,chi,sigm)
         write(6,'(a,a8,a2,a8,f19.4,f12.4)') 
     .      ' LINK (v dir) chi2, dchi2: ',
     .      sname(i0),'--',sname(i1),chi,chi-chio
         if (iomode(6).gt.0)
     .      write(20,'(a16,2x,a8,2x,a8,10x,f18.4,f12.4)') 
     .      ' link site vdir:',
     .      sname(i0),sname(i1),chi,chi-chio
      endif
c
c ----- currently only suitable for geocentric frame
      if (lmode(1).ge.1.and.lmode(3).ge.1) then
c        constrain dz-k2*dx = 0
         fac2 = (z(is2)-z(is1))/(x(is2)-x(is1))
         obsc = 0.0d0
         coef(1) = fac2
         coef(2) = -1.0d0
         coef(3) = -fac2
         coef(4) = 1.0d0
         chio = chi
         indx(1) = is1+1
         indx(2) = is1+3
         indx(3) = is2+1
         indx(4) = is2+3
         call exert_link1(nlive,4,indx,coef,obsc,chi,sigm)
         write(6,'(a,a8,a2,a8,f19.4,f12.4)') 
     .      ' LINK (3dx dir) chi2, dchi2: ',
     .      sname(i0),'--',sname(i1),chi,chi-chio
         if (iomode(6).gt.0)
     .      write(20,'(a16,2x,a8,2x,a8,10x,f18.4,f12.4)') 
     .      ' link site udir:',
     .      sname(i0),sname(i1),chi,chi-chio
      endif
      if (lmode(4).ge.1.and.lmode(6).ge.1) then
c        constrain dVz-k2*dVx = 0
         fac2 = (z(is2)-z(is1))/(x(is2)-x(is1))
         obsc = 0.0d0
         if (type.ge.2) obsc = vz(i0)-vz(i1)-fac2*(ve(i0)-ve(i1))
         coef(1) = fac2
         coef(2) = -1.0d0
         coef(3) = -fac2
         coef(4) = 1.0d0
         chio = chi
         indx(1) = is1+4
         indx(2) = is1+6
         indx(3) = is2+4
         indx(4) = is2+6
         call exert_link1(nlive,4,indx,coef,obsc,chi,sigm)
         write(6,'(a,a8,a2,a8,f19.4,f12.4)') 
     .      ' LINK (3dv dir) chi2, dchi2: ',
     .      sname(i0),'--',sname(i1),chi,chi-chio
         if (iomode(6).gt.0)
     .      write(20,'(a16,2x,a8,2x,a8,10x,f18.4,f12.4)') 
     .      ' link site wdir:',
     .      sname(i0),sname(i1),chi,chi-chio
      endif
      goto 200
c
c     constrain baseline length (currently 3-d only, modify it later)
 100  continue
      i0 = sit_id(1)
      i1 = sit_id(2)
      bl = baslen(i0,i1,0.0d0,3,2)
      is1 = map(i0*6-5)-1
      is2 = map(i1*6-5)-1
      if (is1.lt.0.and.is2.lt.0) print*,' check index: ',i0,i1,is1,is2
      if (is1.lt.0.and.is2.lt.0) goto 200
      coef(1) = (x(i0)-x(i1))/bl
      coef(2) = (y(i0)-y(i1))/bl
      coef(3) = (z(i0)-z(i1))/bl
      icol = 6
      if (is1.lt.0.or.is2.lt.0) icol = 3
      if (lmode(1).ge.1) then
         obsc = 0.0d0
         chio = chi
         do i = 1,3
            indx(i) = is1+i
            indx(i+3) = is2+i
            coef(i+3) = -coef(i)
         enddo
         if (is1.lt.0) then
            do i = 1,3
               indx(i) = indx(i+3)
               coef(i) = coef(i+3)
            enddo
         endif
         call exert_link1(nlive,icol,indx,coef,obsc,chi,sigm)
         write(6,'(a,a8,a2,a8,f19.4,f12.4)') 
     .      ' LINK (len) chi2, dchi2: ',
     .      sname(i0),'--',sname(i1),chi,chi-chio
         if (iomode(6).gt.0)
     .      write(20,'(a16,2x,a8,2x,a8,10x,f18.4,f12.4)') 
     .      ' link base lenx:',
     .      sname(i0),sname(i1),chi,chi-chio
      endif
      if (lmode(4).ge.1) then
         obsc = 0.0d0
         if (type.ge.2) obsc = -coef(1)*(vx(i0)-vx(i1))-
     .      coef(2)*(vy(i0)-vy(i1))-coef(3)*(vz(i0)-vz(i1))
         chio = chi
         do i = 1,3
            indx(i) = is1+i+3
            indx(i+3) = is2+i+3
            coef(i+3) = -coef(i)
         enddo
         if (is1.lt.0) then
            do i = 1,3
               indx(i) = indx(i+3)
            enddo
         endif
         call exert_link1(nlive,icol,indx,coef,obsc,chi,sigm)
         write(6,'(a,a8,a2,a8,f19.4,f12.4)') 
     .      ' LINK (vel) chi2, dchi2: ',
     .      sname(i0),'--',sname(i1),chi,chi-chio
         if (iomode(6).gt.0)
     .      write(20,'(a16,2x,a8,2x,a8,10x,f18.4,f12.4)') 
     .      ' link base lenv:',
     .      sname(i0),sname(i1),chi,chi-chio
      endif
      goto 200
c
c     constraint velocity gradient (geodetic frame currently)
 150  continue
      i0 = sit_id(1)
      i1 = sit_id(2)
      i2 = sit_id(3)
      is1 = map(i0*6-3)+1
      is2 = map(i1*6-3)+1
      is3 = map(i2*6-3)+1
      if (is1.lt.0.and.is2.lt.0) 
     .   print*,' PUT_LINK bad index: ',i0,i1,is1,is2,
     .       +(z(i1)-z(i0))*cf
      if (is1.lt.0.and.is2.lt.0) goto 200
      dx = x(i0)-x(i1)
      dy = y(i0)-y(i1)
      dz = z(i0)-z(i1)
      alo = slon(i1)
      ala = slat(i1)
      call sph_ca(ala,alo,dx,dy,dz,de1,dn1,du1,2)
      dx = x(i2)-x(i1)
      dy = y(i2)-y(i1)
      dz = z(i2)-z(i1)
      call sph_ca(ala,alo,dx,dy,dz,de2,dn2,du2,2)
      s1 = dsin(cita)
      c1 = dcos(cita)
      dec1 = de1*c1+dn1*s1
      dnc1 = dn1*c1-de1*s1
      dec2 = de2*c1+dn2*s1
      dnc2 = dn2*c1-de2*s1
c     Ve gradient
      do i = 1,3
         indx(i) = is1+i-1
         indx(i+3) = is2+i-1
         indx(i+6) = is3+i-1
      enddo
      a1 = slon(i0)
      a2 = slat(i0)
      a3 = radius
      call getjac(a1,a2,a3,cjaco,5)
c     mdm = 1,2
      if (mdm.eq.1.or.mdm.eq.2) then
         if (mdm.eq.1) then
            da1 = dec1
            da2 = dec2
         endif
         if (mdm.eq.2) then
            da1 = dnc1
            da2 = dnc2
         endif
      tmp(1) = da2*c1
      tmp(2) = da2*s1
      call axb(1,2,3,tmp,cjaco,tmp1,1,0)
      do i = 1,3
         coef(i) = tmp1(i)
      enddo
      tmp(1) = (da1-da2)*c1
      tmp(2) = (da1-da2)*s1
      call axb(1,2,3,tmp,cjaco,tmp1,1,0)
      do i = 1,3
         coef(i+3) = tmp1(i)
      enddo
      tmp(1) = -da1*c1
      tmp(2) = -da1*s1
      call axb(1,2,3,tmp,cjaco,tmp1,1,0)
      do i = 1,3
         coef(i+6) = tmp1(i)
      enddo
      obsc = (ve(i0)*da2+ve(i1)*(da1-da2)-ve(i2)*da1)*c1
      obsc = obsc+(vn(i0)*da2+vn(i1)*(da1-da2)-vn(i2)*da1)*s1
      chio = chi
      call exert_link1(nlive,9,indx,coef,obsc,chi,sigm)
      write(6,'(a,a8,2(a2,a8),f19.4,f12.4)') 
     .      ' LINK (Ve gradi) chi2, dchi2: ',
     .      sname(i0),'--',sname(i1),'--',sname(i2),chi,chi-chio
         if (iomode(6).gt.0)
     .      write(20,'(a16,3(2x,a8),f18.4,f12.4)') 
     .      ' link Ve gradie:',
     .      sname(i0),sname(i1),sname(i2),chi,chi-chio
      endif

c     Vn gradient
c     mdm = 3,4
      if (mdm.eq.3.or.mdm.eq.4) then
         if (mdm.eq.3) then
            da1 = dec1
            da2 = dec2
         endif
         if (mdm.eq.4) then
            da1 = dnc1
            da2 = dnc2
         endif
      tmp(1) = -da2*s1
      tmp(2) = da2*c1
      call axb(1,2,3,tmp,cjaco,tmp1,1,0)
      do i = 1,3
         coef(i) = tmp1(i)
      enddo
      tmp(1) = -(da1-da2)*s1
      tmp(2) = (da1-da2)*c1
      call axb(1,2,3,tmp,cjaco,tmp1,1,0)
      do i = 1,3
         coef(i+3) = tmp1(i)
      enddo
      tmp(1) = da1*s1
      tmp(2) = -da1*c1
      call axb(1,2,3,tmp,cjaco,tmp1,1,0)
      do i = 1,3
         coef(i+6) = tmp1(i)
      enddo
      obsc = -(ve(i0)*da2+ve(i1)*(da1-da2)-ve(i2)*da1)*s1
      obsc = obsc+(vn(i0)*da2+vn(i1)*(da1-da2)-vn(i2)*da1)*c1
      chio = chi
      call exert_link1(nlive,9,indx,coef,obsc,chi,sigm)
      write(6,'(a,a8,2(a2,a8),f19.4,f12.4)') 
     .      ' LINK (Vn gradi) chi2, dchi2: ',
     .      sname(i0),'--',sname(i1),'--',sname(i2),chi,chi-chio
         if (iomode(6).gt.0)
     .      write(20,'(a16,3(2x,a8),f18.4,f12.4)') 
     .      ' link Vn gradie:',
     .      sname(i0),sname(i1),sname(i2),chi,chi-chio
      endif
      goto 200
c
c     
 250  continue
      i0 = sit_id(1)
      call chk_epi_list(i0,e_time,is1,imatch)
      if (imatch.le.0) print*,' PUT_LINK: missing epi_site ',sname(i0)
      if (imatch.le.0) goto 200
      coef(1) = 1.0d0
      coef(2) = -1.0d0
      do 210 ib = 2,isnum
         i1 = sit_id(ib)
         call chk_epi_list(i1,e_time,is2,imatch)
         if (imatch.lt.0) goto 210
c         print*,' --- force site & index:',jtoi(i0),jtoi(i1),is1,is2
         chio = chi
         do 220 i = 1,6
            if (lmode(i).ge.1) then
               indx(1) = is1+i-1
               indx(2) = is2+i-1
               obsc = 0.0d0
               call exert_link1(nlive,2,indx,coef,obsc,chi,sigm)
            endif
 220     continue
         write(6,'(a,a8,a2,a8,f18.4,f12.4)') 
     .      ' LINK (epi_site) chi2, dchi2: ',
     .      sname(i0),'--',sname(i1),chi,chi-chio
         if (iomode(6).gt.0)
     .      write(20,'(a16,2x,a8,2x,a8,10x,f18.4,f12.4)') 
     .      ' link epi_site :',
     .      sname(i0),sname(i1),chi,chi-chio
 210  continue
      goto 200
c
 200  continue
c
      return
      end


