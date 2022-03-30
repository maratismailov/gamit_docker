      subroutine fonda_out_compress(ifil1,ifil2,ifil3)
c
c
      implicit real*8(a-h,o-z)
      include 'maked.fti'

      character*5 netid(maxsit)
      character*8 name1,name2,sitnam(maxsit)
      integer  ist1(maxobs),ist2(maxobs)
      integer nblen, ie1,it,match_name,ibase1,ibase2,i2
      integer iobs1,i1
      character*30 note
      character*22 notice(2)
      dimension otime(maxobs)
      character*128 line
      integer i,j,icount,isdx,nlen,isit1,isit2,llen
      integer iline,irate,itest,ifil1,ifil2,ifil3,ks1,ks2
      integer iadd
      equivalence (ist1,idobs(1,1))
      equivalence (ist2,idobs(1,2))
      data notice/' Abnormal velocity:   ',
     .            ' Abnormal observation:'/
c
c     default format
      if (infmt(1:1).eq.'*')
     . infmt(1:31) = '(1x,f10.4,f14.5,f14.8,2(2x,a8))'
c
      iadd = 0
c     
c     First, get site names in the network
      read (ifil1,*) nsit
      note = '{network site number          '
      write (ifil2,30) nsit,note
      do i = 1,nsit
         read (ifil1,'(a8,1x,a5)') sitnam(i),netid(i)
c        read (ifil1,'(a5,1x,a8)') netid(i),sitnam(i)
         write (ifil2,'(a8,1x,a5)') sitnam(i),netid(i)
      enddo
      print*,' Begin to compress data of --- ',infil
c
c     get experiment number 
      read (ifil1,*) iexp
      note = '{experiment number            '
      write (ifil2,30) iexp*2,note
      close (ifil2)
      open (24,file='poor.obs',status='unknown')
      write (24,'(a)') ' ----- Table of poor observations -----'
c
c     experiment loop
      llen = nblen(outfil)
      note = '{exp. index, obs. type, number'
      do 90 i = 1,iexp
c
c        open auxiluary files
         open (21,file='rate.obs',status='unknown')
         open (22,file='length.obs',status='unknown')
         open (23,file='label.obs',status='unknown')
c        read experiment index, obs. type and number
         read (ifil1,*) ie1,it,iobs1
c
         if (it.gt.10) then
            iadd = iadd+1
            write (23,140) ie1,it,iobs1,note
            do j = 1,iobs1
               call blank(line)
               read (ifil1,'(a)') line
               nlen = nblen(line)
               write (23,'(a)') line(1:nlen)
            enddo
            close (23)
            call blank(line)
            line(1:17+llen) = 'cat label.obs >> ' // outfil(1:llen)
            j = system(line)
            goto 90
         endif   
c        read observation data
         do 50 j = 1,iobs1
            read(ifil1,infmt,err=300) 
     .         otime(j),data(j),erd(j),name1,name2 
            nlen = nblen(name1)
            isit1 = match_name(nsit,nlen,sitnam,name1)
            if (isit1 .le. 0) goto 100
            nlen = nblen(name2)
            isit2 = match_name(nsit,nlen,sitnam,name2)
            if (isit2 .le. 0) goto 100
            ist1(j) = isit1
            ist2(j) = isit2
            igobs(j) = j
 50      continue
c
c        rearrange the data baseline by baseline
         call isort2(iobs1,ist1,ist2,igobs)
c
c        estimate baseline length and baseline length rate
         ibase1 = ist1(1)
         ibase2 = ist2(1)
         icount = 0
         t0 = 0.0d0
         t1 = 0.0d0
         w0 = 0.0d0
         wx = 0.0d0
         wv = 0.0d0
         iline = 0
         irate = 0
         do 40 j = 1,iobs1
            isit1 = ist1(j) 
            isit2 = ist2(j) 
c           same baseline
            if (isit1.eq.ibase1.and.isit2.eq.ibase2) then
               icount = icount+1
               isdx = igobs(j)
               if (icount.eq.1) then
                  tsave = otime(isdx)
                  osave = data(isdx)
               endif
               time1 = otime(isdx)-tsave
               obs1 = data(isdx)-osave
               wt = 1.0d0/erd(isdx)**2
               t0 = t0+time1*wt
               w0 = w0+wt
               wx = wx+obs1*wt
               t1 = t1+time1*time1*wt
               wv = wv+time1*obs1*wt
               if (j.lt.iobs1) goto 40
            endif
c           change baseline
            if (isit1.ne.ibase1.or.isit2.ne.ibase2.or.j.eq.iobs1) then
               iline = iline+1
               t0 = t0/w0
               wx = wx/w0
               cerr = dsqrt(1.0d0/w0)
c              put the squeezed data into temporary file
               write (22,150)
     .            t0+tsave,wx+osave,cerr,sitnam(ibase1),
     .            sitnam(ibase2),icount
c              one observation can not get rate estimate
               if (icount.le.1) goto 20
c              output mapping file with GMT format
               if (iomode(5).gt.0) then
                  ks1 = match_name(maxnet,8,sname,sitnam(ibase1))
                  if (ks1.le.0) 
     .               print*,' Mismatch site: ',sitnam(ibase1)
                  ks2 = match_name(maxnet,8,sname,sitnam(ibase2))
                  if (ks2.le.0) 
     .               print*,' Mismatch site: ',sitnam(ibase2)
c                 use western longitude to meet GMT default
                  wlo1 = slo(ks1)*rtod
                  if (wlo1.gt.180.0d0) wlo1 = wlo1-360.0d0
                  wlo2 = slo(ks2)*rtod
                  if (wlo2.gt.180.0d0) wlo2 = wlo2-360.0d0
                  if (ks1.gt.0.and.ks2.gt.0) 
     .               write (ifil3,160) sitnam(ibase1),sitnam(ibase2),
     .               wlo1,sla(ks1)*rtod,
     .               wlo2,sla(ks2)*rtod
               endif
               wv = (wv-wx*t0*w0)/(t1-t0*t0*w0)
               cverr = dsqrt(1.0d0/(t1-t0*t0*w0))
               irate = irate+1
c              store abnormal observations into poor.obs file
c              check sigma and outliers
               sigma = 0.0d0
               if (icount.gt.2) then
                  i2 = j-1
                  i1 = j-icount
                  if (j.eq.iobs1) i2 = j
                  if (j.eq.iobs1) i1 = j-icount+1
                  do itest = i1,i2
                     isdx = igobs(itest)
                     time1 = otime(isdx)-tsave-t0
                     obs1 = data(isdx)-osave
                     dev = (obs1-wx-wv*time1)*1.0d3
                     sigma = sigma+dev*dev
                  enddo
                  sigma = dsqrt(sigma/dble(icount))
                  cri = 3.0d0
                  if (sigma.gt.1.0d2) cri = 2.0d0
                  if (sigma.gt.1.0d3) cri = 1.0d0
                  if (sigma.lt.1.0d0) cri = 5.0d0
                  if (sigma.lt.1.0d-1) cri = 10.0d0
                  do itest = i1,i2
                     isdx = igobs(itest)
                     time1 = otime(isdx)-tsave-t0
                     obs1 = data(isdx)-osave
                     dev = (obs1-wx-wv*time1)*1.0d3
                     if (dabs(dev).gt.cri*sigma) then
                        write(24,'(a22,2(5x,a4,f11.4))') 
     .                     notice(2),'sig=',sigma,'dev=',dev
                        write(24,infmt) 
     .                     otime(isdx),data(isdx),erd(isdx),
     .                     sitnam(ibase1),sitnam(ibase2)
                     endif
                  enddo
               endif      
c              using mm/year as the velocity unit
               write (21,150)
     .            t0+tsave,wv*1.0d3,cverr,sitnam(ibase1),
     .            sitnam(ibase2),icount,sigma
               write (6,150)
     .            t0+tsave,wv*1.0d3,cverr,sitnam(ibase1),
     .            sitnam(ibase2),icount,sigma
               if (dabs(wv*1.0d3).gt.0.5d2) then 
                  write (24,'(a)') notice(1)
                  write (24,150)
     .            t0+tsave,wv*1.0d3,cverr,sitnam(ibase1),
     .            sitnam(ibase2),icount,sigma
               endif
 20            icount = 1
               isdx = igobs(j)
               tsave = otime(isdx)
               osave = data(isdx)
               wt = 1.0d0/erd(isdx)**2
               t0 = 0.0d0
               t1 = 0.0d0
               w0 = wt
               wx = 0.0d0
               wv = 0.0d0
               ibase1 = isit1
               ibase2 = isit2
            endif
 40      continue
         write (23,140) i*2-iadd-1,it+10,irate,note
         write (23,140) i*2-iadd,it,iline,note
         close (21)
         close (22)
         close (23)
         call blank(line)
         line(1:21+llen) = 'head -1 label.obs >> ' // outfil(1:llen)
         isdx = system(line)
         call blank(line)
         line(1:16+llen) = 'cat rate.obs >> ' // outfil(1:llen)
         isdx = system(line)
         call blank(line)
         line(1:21+llen) = 'tail -1 label.obs >> ' // outfil(1:llen)
         isdx = system(line)
         call blank(line)
         line(1:18+llen) = 'cat length.obs >> ' // outfil(1:llen)
         isdx = system(line)
 90   continue
      isdx = system("\\rm -f label.obs rate.obs length.obs")
      close (24)
      goto 120
c 
c     error
 100  print*,' Site name mismatch at FONDA_OUT_COMPRESS: ',
     .   name1,' ',name2
      goto 120
 300  print*,' Read file error at FONDA_OUT_COMPRESS'
c
 30   format (i5,25x,a30)
 140  format (3i5,15x,a30)
 150  format (1x,f10.4,f14.5,f14.8,2(2x,a8),2x,i4,2x,f12.2)
 160  format (">",a8," to ",a8,/,2(1x,f12.6),/,2(1x,f12.6))
c
 120  continue
      if (iadd.gt.0) print*,' There are uncompressible data. ',
     .  'Total # of experiments should add ',iadd

      return
      end
