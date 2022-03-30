      program get_scater
c
c     get baseline scatters
c
c     options:
c        1 = (undecided)
c        2 = (undecided)
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
      integer i,i1,i2,is,l1,l2,l4,isa
      integer ios,lcmd,ioptn,len
      integer*4 iclarg
      real*8 obs(200),time(200),sig(200),calc(200)
      real*8 scx(200),scy(200),scs(200)
      character*120 line
      integer nblen,count_arg,lift_arg
      real   upb,lob
      real*8 up_rms,lo_rms
c
c     check if the runstring is complete.
      len = iclarg(1,infile)
      l2 = iclarg(2,outfile)
      l1 = iclarg(3,listfile)
      if (l1.le.0.or.l2.le.0.or.len.le.0) then
         print*,' Calculate length scatter from solvem.res file'
         print*,' Runstring: get_scater <infile> <outfile>',
     .      ' <linelist> <option>'
         print*,'  option: 1= raw data      2= postfit    '
         print*,'  option: 3= best of 1,2   4= (undecided)'
         stop
      endif
      i = iclarg(4,tmpn)
      if (i.le.0) then
         ioptn = 1
      else
         read (tmpn,*) ioptn
         if (ioptn.le.0.or.ioptn.ge.4) ioptn = 1
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
      isa = 0
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
         print*,' scatter for baseline: ',sname1,' -- ',sname2
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
         as = 0.0d0
         sigma1 = 0.0d0
         fac = 1.0d3
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
               read (line(1:lcmd),14,iostat=ios,err=30) 
     .         time(l4),obs(l4),calc(l4),sig(l4)
            endif
 14         format (23x,f8.3,2x,f14.5,f17.5,2x,f11.5)
 15         format (23x,f8.3,35x,f11.5,2x,f13.5)
            if (l4.eq.1) then
               tsave = time(l4)
               osave = obs(l4)
            endif
            time1 = time(l4)-tsave
            obs1 = obs(l4)-osave
            wt = 1.0d0/(sig(l4)**2)
            as = as+sig(l4)
            w0 = w0+wt
            t0 = t0+time1*wt
            wx = wx+obs1*wt
            t1 = t1+time1*time1*wt
            wv = wv+time1*obs1*wt
            if (ioptn.eq.2.or.ioptn.eq.3) then
               dev = (obs(l4)-calc(l4))*fac
               sigma1 = sigma1+dev*dev*wt
            endif 
 30      continue
 25      close (11)
         if (l4.le.2) goto 20
         as = as/dble(l4)
         t0 = t0/w0
         wx = wx/w0
         wv = (wv-wx*t0*w0)/(t1-t0*t0*w0)
         if (ioptn.eq.2) then
            sigma = sigma1
            goto 55
         endif
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
         if (ioptn.eq.3.and.sigma.gt.sigma1) sigma = sigma1
 55      rms = dsqrt(sigma*dble(l4)/(w0*dble(l4-2)))
c        68.24% confidence level is equivalent to 1 sigma
         call chi2_dis(l4-2,0.6824,upb,lob)
         aw = as
         up_rms = aw*dabs(dsqrt(dble(upb)/dble(l4-2))-1.0d0)
         lo_rms = aw*dabs(1.0d0-dsqrt(dble(lob)/dble(l4-2)))
         sigma = 0.5d0*(up_rms+lo_rms)
         write (12,40) as,rms,sigma
         sigma = (0.5d0*(dble(upb)-dble(lob))/dble(l4-2))*(aw**2)
         isa = isa+1
         scx(isa) = obs(1)**2
         scy(isa) = rms**2
         scs(isa) = sigma
 20   continue
      goto 100

 1000 print *,' Can not open the file :'
      stop 'in GET_SCATER. '

 100  continue
c
c     estimate error model sigma**2 = a**2 + (b*L)**2
      if (isa.ge.2) then 
         wx = 0.0d0
         wv = 0.0d0
         w0 = 0.0d0
         w1 = 0.0d0
         ws = 0.0d0
         do 130 i = 1,isa
            s1 = 1.0d0/(scs(i)**2)
            ws = ws+s1
            wx = wx+scx(i)*s1
            wv = wv+scy(i)*s1
            w0 = w0+scx(i)**2*s1
            w1 = w1+scx(i)*scy(i)*s1
 130     continue
         dlt = w0*ws-wx**2
         If (dabs(dlt).lt.1.0d-10) then
            aterm = 0.0d0
            bterm = 0.0d0
            a1 = 0.0d0
            b1 = 0.0d0
         else
            aterm = (wv*w0-w1*wx)/dlt
            bterm = (ws*w1-wx*wv)/dlt
            fac = 1.0d0
            if (aterm.lt.0.0d0) fac = -1.0d0
            aterm = fac*dsqrt(dabs(aterm))
            fac = 1.0d3
            if (bterm.lt.0.0d0) fac = -1.0d3
            bterm = fac*dsqrt(dabs(bterm))
            a1 = dsqrt(w0/dlt)
            b1 = dsqrt(ws/dlt)*1.0d6
         endif
      else
         aterm = 0.0d0
         bterm = 0.0d0
         a1 = 0.0d0
         b1 = 0.0d0
      endif

c     write (12,'(''*  n ='',i4,''  a ='',f6.2,''+-'',f6.2,
c    .   '' mm   b ='',f6.2,''+-'',f6.2,'' ppm'')') 
c    .   isa,aterm,a1,bterm,b1
      write (12,'(''*  n ='',i4,''  a ='',f6.2,
     .   '' mm   b ='',f6.2,'' ppm'')') 
     .   isa,aterm,bterm
      l1 = system('RM tmp.1 tmp.2')
 40   format (3f10.3)
      close (12)
      close (13)
c
      stop
      end
