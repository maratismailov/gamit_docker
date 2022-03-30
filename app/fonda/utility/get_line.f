      program get_line
c
c     get single baseline observation information
c
c     options:
c        1 = get baseline length observations
c        2 = get residuals
c        3 = (anything else?)
c        4 = (anything else?)
c
c     unit:
c         x, y, z : m
c        vx,vy,vz : mm/year
c        slat,slon: radian
c        srad     : m
c        ve,vn,vu : mm/year
c        time     : year
c
      implicit real*8(a-h,o-z)
c
      character*8 sname1,sname2,tmpn
      character*64 infile,outfile,listfile
      character*120 command
      integer i,i1,i2,is,l1,l2,l4
      integer ios,lcmd,ioptn,len
      real*8 obs(200),time(200),sig(200),calc(200)
      character*120 line
      integer nblen,count_arg,lift_arg
      integer*4 iclarg
c
c     check if the runstring is complete.
      len = iclarg(1,infile)
      l2 = iclarg(2,outfile)
      l1 = iclarg(3,listfile)
      if (l1.le.0.or.l2.le.0.or.len.le.0) then
         print*,' Get length data from solvem.res file'
         print*,' Runstring: get_line <infile> <outfile>',
     .      ' <linelist> <option>'
         print*,'  option: 1= obs. data      2= residual'
         print*,'          3= (else?)        4 = (else?)'
         stop
      endif
      i = iclarg(4,tmpn)
      if (i.le.0) then
         ioptn = 1
      else
         read (tmpn,*) ioptn
         if (ioptn.le.0.or.ioptn.ge.3) ioptn = 1
      endif
c
c     open file
      open (12,file=outfile,status='unknown',err=1000)
      open (13,file=listfile,status='old',err=1000)
      print*,' Begin to subtract single baseline data .....'
c
c     get cleaned data file
      call blank(command)
      command(1:13) = 'grep -v "\*" '
      command(14:13+len) = infile(1:len)
      command(14+len:22+len) = ' >! tmp.1'
      l1 = system(command)
c
c     set up BIG DO loop to read until end of file...
      rewind (13)
      do 20 i = 1,1000
         call blank(line)
         ios = 0
         read (13,'(a)',iostat=ios,end=100,err=20) line
         lcmd = count_arg(line)
c        need 2 sites 
         if (lcmd.le.1) goto 20
         if (line(1:1).eq.'*') goto 20
         i1 = lift_arg(line,sname1,1)
         i2 = lift_arg(line,sname2,2)
         print*,' get baseline: ',sname1,' -- ',sname2
         call blank(command)
         command(1:8) = 'grep -i '
         command(9:9+i1) = sname1(1:i1)
         command(10+i1:26+i1) = ' tmp.1 | grep -i '
         command(27+i1:26+i1+i2) = sname2(1:i2)
         command(27+i1+i2:35+i1+i2) = ' >! tmp.2'
         l1 = system(command)
         open (11,file='tmp.2',status='old',err=20)
         l4 = 0
         t0 = 0.0d0
         t1 = 0.0d0
         wx = 0.0d0
         wv = 0.0d0
         w0 = 0.0d0
         do 30 is = 1,1000
            call blank(line)
            ios = 0
            read (11,'(a)',iostat=ios,end=25,err=20) line
            lcmd = count_arg(line)
            lcmd = nblen(line)
            if (lcmd.le.50) goto 30
            l4 = l4+1
            if (ioptn.eq.1) then
               read (line(1:lcmd),14,iostat=ios,err=30) 
     .         time(l4),obs(l4),calc(l4),sig(l4)
            else
               read (line(1:lcmd),15,iostat=ios,err=30) 
     .         time(l4),sig(l4),obs(l4)
            endif
 14         format (23x,f8.3,2x,f14.5,f17.5,2x,f11.5)
 15         format (23x,f8.3,35x,f11.5,2x,f13.5)
            if (l4.eq.1) then
               tsave = time(l4)
               osave = obs(l4)
               if (ioptn.eq.2) osave = 0.0d0
            endif
            time1 = time(l4)-tsave
            obs1 = obs(l4)-osave
            wt = 1.0d0/(sig(l4)**2)
            w0 = w0+wt
            t0 = t0+time1*wt
            wx = wx+obs1*wt
            t1 = t1+time1*time1*wt
            wv = wv+time1*obs1*wt
 30      continue
 25      close (11)
         if (l4.le.1) goto 20
         t0 = t0/w0
         wx = wx/w0
         wv = (wv-wx*t0*w0)/(t1-t0*t0*w0)
         fac = 1.0d3
         if (ioptn.eq.2) fac = 1.0d0
c
c        statistics
         sigma = 0.0d0
         do 50 is = 1,l4
            time1 = time(is)-tsave-t0
            obs1 = obs(is)-osave
            wt = 1.0d0/(sig(is)**2)
            dev = (obs1-wx-wv*time1)*fac
            sigma = sigma+dev*dev*wt
 50      continue
         scat = dsqrt(scat/dble(l4))
         if (l4.le.2) then
            rms = 0.0d0
            scat = 0.0d0
         else
            scat = dsqrt(sigma*dble(l4)/dble(l4-2)/w0)
            rms = dsqrt(sigma/dble(l4-2))
         endif
         if (ioptn.eq.1) then
            write (12,'(''* Baseline : '',a8,'' -- '',a8,3x,
     .      ''t(y)='',f8.3,2x,''l(m)='',f13.4,2x,
     .      ''v(mm/y)='',f8.2,2x,''s(mm)='',f8.2,2x,''nrms='',f6.2)')
     .      sname1,sname2,t0+tsave,wx+osave,wv*fac,scat,rms
         else
            write (12,'(''* Residual : '',a8,'' -- '',a8,4x,
     .      ''t(y)='',f8.3,2x,''r(mm)='',f11.3,2x,
     .      ''v(mm/y)='',f8.2,2x,''s(mm)='',f8.2,2x,''nrms='',f6.2)')
     .      sname1,sname2,t0+tsave,wx+osave,wv*fac,scat,rms
         endif
         do is = 1,l4
            write (12,'(3x,f8.3,3x,f14.5,3x,f12.5,3x,f14.5)')
     .         time(is),(obs(is)-wx-osave)*fac,sig(is),
     .         (calc(is)-wx-osave)*fac
         enddo
         write (12,40) t0+tsave,wx+osave,wv*fac,scat
         cri = 3.0d0
         l1 = 0
         do is = 1,l4
            time1 = time(is)-tsave-t0
            obs1 = obs(is)-osave
            dev = (obs1-wx-wv*time1)*fac
            if (dabs(dev).gt.cri*scat) then
               write (12,60) is,dev
               l1 = l1+1
            endif
         enddo
         if (l1.le.0) write (12,'(''* no outlier(3 sigma).'')')

 20   continue  

      goto 100

 1000 print *,' Can not open the file :'
      stop 'in NET_UPDATE. '

 100  continue
 40   format ('* solution : t=',f10.4,'  x(m)=',f14.5,'  v(mm/y)=',
     .        f10.3,'  sigma(mm)=',f10.3)
 60   format ('* outlier(2 sigma) : obs. ID =',i3,'  dev(mm)=',f10.3)
      close (12)
      close (13)
c
      stop
      end
