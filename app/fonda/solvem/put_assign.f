      subroutine put_assign(wcmd,lcmd,chi,job)
c
c     perform assign constraint
c
c     job = 1:   assign network center constraint
c     job = 2:   assign network rotation constraint
c     job = 3:   assign site direction constraint
c     job = 4:   assign inner coordinate constraint
c     job = 5:   assign outer coordinate constraint
c     job = 6:   assign model coordinate constraint
c     job = 7:   assign episodic site displacement
c
      implicit real*8(a-h,o-z)
      include 'solvem.fti'
c
      character*(*) wcmd
      character*80 code
      integer job,indx,i1,i0,ic,is1,i2,i3,is3,lm
      integer ib,isnum,is2,i,match_name,imatch
      integer lcmd,count_arg,lift_arg,lmode(6),sit_id(20)
      dimension indx(maxsit*18),coef(maxsit*18)

c     discard junk job
      if (job.lt.1.or.job.gt.8) goto 70
c  
c     decomposite command line
      ic = count_arg(wcmd)
c     not enough arguments
      if (ic.le.1) goto 70
c     pointer to coordinate and velocity
      ib = lift_arg(wcmd,code,1)
      if (ib.le.0) goto 70
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
c
c     read site name and match their index
      ib = lift_arg(wcmd,code,2)
      if (ib.le.0) goto 70
      if (job.eq.4) then
         open (34,file=code,status='old')
         call get_iom_list(34,isnum,indx,xs,ys,zs,slxy,1)
         close (34)
         goto 90
      endif
      if (job.eq.5) then
         open (35,file=code,status='old')
         ib = lift_arg(wcmd,code,3)
         read(code,*) obst
         azio = obst*dtor
         call get_iom_list(35,isnum,indx,xs,ys,zs,soxy,2)
         close (35)
         goto 90
      endif
      if (job.eq.6) then
         open (36,file=code,status='old')
         call get_iom_list(36,isnum,indx,xs,ys,zs,slxy,3)
         close (36)
         goto 90
      endif
      if (job.eq.7) then
         i0 = match_name(nsit,ib,sname,code)
         if (i0.le.0.or.ic.le.4) then
            print*,' PUT_ASSIGN name mismatch: ',code,' ',ic
            goto 70
         endif
         sit_id(1) = i0
         ib = lift_arg(wcmd,code,3)
         read(code,*) e_time
         goto 90
      endif
         
      if (code(1:3).eq.'net'.or.code(1:3).eq.'all') then
         isnum = nsit
      else
         i0 = match_name(nsit,ib,sname,code)
         if (i0.le.0) then
            print*,' PUT_ASSIGN name mismatch:',code
            goto 200
         endif
         sit_id(1) = i0
         isnum = 1
      endif
c
c     get constraint value
      if (ic.le.2) then
         obst = 0.0d0
      else 
         ib = lift_arg(wcmd,code,3)
         read(code,*) obst
      endif
c
c     get constraint sigma
      if (ic.le.3) then
         sigc = 0.0d0
      else 
         ib = lift_arg(wcmd,code,4)
         read(code,*) sigc
      endif
 90   goto (20,50,100,150,250,300,350,200) job
c
c     network center movement constraint
 20   continue
      do 30 i = 1,6
         if (lmode(i).le.0) goto 30
         i1 = 0
         chio = chi
         do 40 ib = 1,isnum
            is2 = map((ib-1)*6+i)
            if (is2.lt.0.or.is2.gt.nlive) goto 40
            c1 = anorm(is2*(is2+1)/2)
            if (c1.le.0.0d0) goto 40
            i1 = i1+1
            indx(i1) = is2
            coef(i1) =  1.0d0/c1
 40      continue
         if (i1.le.0) goto 30
         call exert_link1(nlive,i1,indx,coef,obst,chi,0.0d0)
         print*,' ASSIGN: (center) chi2, dchi2: ',chi,chi-chio
         if (iomode(6).gt.0)
     .      write(20,'(a16,30x,f18.4,f12.4)')
     .      ' fix net centr:',chi,chi-chio
 30   continue
      goto 200
c
c     network rotation constraint
 50   continue
      do 60 i = 1,6
         if (lmode(i).le.0) goto 60
         i1 = 0
         i2 = i+1
         i3 = i+2
         if (i.eq.2) i3 = 1
         if (i.eq.3) i2 = 1
         if (i.eq.3) i3 = 2
         if (i.eq.5) i3 = 4
         if (i.eq.6) i2 = 4
         if (i.eq.6) i3 = 5
         call geocnt(xs,ys,zs,slxy,slxz,slyz,so2y,1)
         chio = chi
         do 80 ib = 1,isnum
            is1 = map((ib-1)*6+i2)
            is2 = map((ib-1)*6+i3)
            if (is1.lt.0.or.is1.gt.nlive) goto 80
            if (is2.lt.0.or.is2.gt.nlive) goto 80
            c1 = anorm(is1*(is1+1)/2)
            c2 = anorm(is2*(is2+1)/2)
            if (c1.le.0.0d0) goto 80
            if (c2.le.0.0d0) goto 80
            c1 = 1.0d0
            c2 = 1.0d0
            i1 = i1+1
            indx(i1) = is1
            if (i.eq.1.or.i.eq.4) fac = y(ib)-ys
            if (i.eq.2.or.i.eq.5) fac = z(ib)-zs
            if (i.eq.3.or.i.eq.6) fac = x(ib)-xs
            coef(i1) =  fac/c1
            i1 = i1+1
            indx(i1) = is2
            if (i.eq.1.or.i.eq.4) fac = z(ib)-zs
            if (i.eq.2.or.i.eq.5) fac = x(ib)-xs
            if (i.eq.3.or.i.eq.6) fac = y(ib)-ys
            coef(i1) =  -fac/c2
 80      continue
         if (i1.le.0) goto 60
         call exert_link1(nlive,i1,indx,coef,obst,chi,0.0d0)
         print*,' ASSIGN: (rotation) chi2, dchi2: ',chi,chi-chio
         if (iomode(6).gt.0)
     .      write(20,'(a16,30x,f18.4,f12.4)')
     .      ' fix net rotat:',chi,chi-chio
 60   continue
      goto 200
c
c     site direction constraint
 100  continue
      i0 = sit_id(1)
      is1 = map((i0-1)*6+4)
      if (is1.lt.0.or.is1.gt.nlive) goto 200
      is2 = map((i0-1)*6+5)
      if (is2.lt.0.or.is2.gt.nlive) goto 200
      is3 = map((i0-1)*6+6)
      if (is3.lt.0.or.is3.gt.nlive) goto 200
      fac = dtan(obst*dtor)
      c1 = dcos(slon(i0))
      s1 = dsin(slon(i0))
      c2 = dcos(slat(i0))
      s2 = dsin(slat(i0))
      indx(1) = is1
      coef(1) = -s1+c1*s2*fac 
      indx(2) = is2
      coef(2) = c1+s1*s2*fac 
      indx(3) = is3
      coef(3) = -c2*fac
      ob = fac*vn(i0)-ve(i0)
      chio = chi
      sigc2 = sigc**2
      call exert_link1(nlive,3,indx,coef,ob,chi,sigc2)
      print*,' ASSIGN: (site dir) chi2, dchi2: ',
     .   sname(i0),' ',chi,chi-chio
      if (iomode(6).gt.0)
     .   write(20,'(a16,2x,a8,2x,f7.1,11x,f18.4,f12.4)')
     .   ' assign sit dir:',sname(i0),obst,chi,chi-chio
      goto 200
c
c     inner coordinate constraint
c     we always perform 3-d translation and 2-d rotation
 150  continue
      lm = isnum*3
      call inner_sln(4,lm,indx,xs,ys,zs,slxy)
      print*,' ASSIGN: inner coordinate solution',lm
      if (iomode(6).gt.0)
     .   write(20,'(a16,2x,a9,19x,f18.4,f12.4)')
     .   ' inner coor sln:','inner sln',chi,chi-chio
      goto 200
c
c     outer coordinate constraint
c     we always perform 2-d outer coordinate solution in 
c     local topocentric frame
 250  continue
c     velocity     
      lm = isnum*3
      call outer_sln(4,lm,indx,xs,ys,zs,soxy)
      print*,' ASSIGN: outer coordinate solution',lm,azio
      if (iomode(6).gt.0)
     .   write(20,'(a16,2x,a9,19x,f18.4,f12.4)')
     .   ' outer coor sln:','outer sln',chi,chi-chio
      goto 200
c
c     model coordinate constraint
c     we always perform 3-d translation and 2-d rotation
 300  continue
      lm = isnum*3
      call model_sln(4,lm,indx,xs,ys,zs,slxy)
      print*,' ASSIGN: model coordinate solution',lm
      if (iomode(6).gt.0)
     .   write(20,'(a16,2x,a9,19x,f18.4,f12.4)')
     .   ' model coor sln:','model sln',chi,chi-chio
      goto 200
c
c     
 350  continue
      i0 = sit_id(1)
      call chk_epi_list(i0,e_time,is1,imatch)
      if (imatch.le.0) print*,' PUT_ASSIGN: missing epi_site ',sname(i0)
      if (imatch.le.0) goto 70
      sl = dsin(slon(i0))
      cl = dcos(slon(i0))
      sf = dsin(slat(i0))
      cf = dcos(slat(i0))
      i1 = 0
      indx(1) = is1
      indx(2) = is1+1
      indx(3) = is1+2
      do 370 i = 1,3
         if (lmode(i).lt.1) goto 370
         i1 = i1+1
         ib = lift_arg(wcmd,code,3+i1)
         if (ib.le.0) goto 70
         read(code,*) obst
         i1 = i1+1
         ib = lift_arg(wcmd,code,3+i1)
         if (ib.le.0) goto 70
         read(code,*) sigc
         sigc2 = sigc**2
         if (i.eq.1) then
            coef(1) = -sl
            coef(2) = cl
            i2 = 2
         endif
         if (i.eq.2) then
            coef(1) = -sf*cl
            coef(2) = -sf*sl
            coef(3) = cf
            i2 = 3
         endif
         if (i.eq.3) then
            coef(1) = cf*cl
            coef(2) = cf*sl
            coef(3) = sf
            i2 = 3
         endif
         chio = chi
         call exert_link1(nlive,i2,indx,coef,obst,chi,sigc2)
         write(6,'(a,a8,10x,f18.4,f12.4)') 
     .      ' ASSI (epi_site) chi2, dchi2: ',
     .      sname(i0),chi,chi-chio
         if (iomode(6).gt.0)
     .      write(20,'(a16,2x,a8,10x,10x,f18.4,f12.4)') 
     .      ' assi epi_site :',
     .      sname(i0),chi,chi-chio
 370  continue
      goto 200
c
c     junk assign
 70   print*,' unknown ASSIGN !!! ',wcmd(1:lcmd)
      return
c
 200  continue
c
      return
      end
