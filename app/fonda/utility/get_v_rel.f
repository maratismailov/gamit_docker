      program get_v_rel
c
c     get relative velocity, the variance keeps the same.
c     reference site should be in the first effective line.
c
c     options:
c        1 = FONDA mapping file format
c        2 = FONDA apriori file format
c        3 = GLOBK output velocity file format
c        4 = GENRELREF output file format
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
      character*64 infile,outfile
      character*8 tmpn
      integer i,l1,l2
      integer ios,lcmd,ioptn,len,iyr,imon,iday,ihr,imn,isec
      character*128 line
      integer nblen
      integer*4 iclarg
c
c     check if the runstring is complete.
      len = iclarg(1,infile)
      l2 = iclarg(2,outfile)
      if (l2.le.0.or.len.le.0) then
         print*,' Runstring: get_v_rel <infile> <outfile>',
     .      ' <option>'
         print*,'  option: 1=FONDA map file  2=FONDA apriori file'
         print*,'          3=GLOBK vel file  4=GENRELREF file'
         stop
      endif
      i = iclarg(3,tmpn)
      if (i.le.0) then
         ioptn = 1
      else
         read (tmpn,*) ioptn
         if (ioptn.le.0.or.ioptn.gt.4) ioptn = 1
      endif
c
c     open file
      open (11,file=infile,status='old',err=1000)
      open (12,file=outfile,status='unknown',err=1000)
      print*,' Begin to compare the site velocities ...'
c
c     set up BIG DO loop to read until end of file...
      rewind (11)
      l1 = 0
      do 20 i = 1,1000
c        read site coordinte and velosity
         ios = 0
         read (11,'(a)',iostat=ios,end=100,err=20) line
         lcmd = nblen(line)
         if (lcmd.le.0) goto 20
c        copy comment lines
         if (line(1:1).ne.' ') then
            if (i.eq.1) then
               call getdat(iyr,imon,iday)
               call gettim(ihr,imn,isec,ios)
               write (12,'(a18,a,a8,i4,2(a1,i2),3x,i2,2(a1,i2))')
     .            '* GET_V_REL from: ',infile(1:len),'  time: ',
     .            iyr,'/',imon,'/',iday,ihr,':',imn,':',isec
            endif
            write (12,'(a)') line(1:lcmd)
            goto 20
         endif
         if (ioptn.eq.1) 
     .      read (line,*,err=20) alon,alat,ve,vn
         if (ioptn.eq.2) 
     .      read (line(1:lcmd),'(68x,3f8.4)',err=20) ve,vn,vu
         if (ioptn.eq.3) 
     .      read (line(1:lcmd),'(16x,2f8.2,40x,f8.2)',err=20) ve,vn,vu
         if (ioptn.eq.4) 
     .      read (line(1:lcmd),'(19x,2f8.2,40x,f8.2)',err=20) ve,vn,vu
         l1 = l1+1
         if (l1.eq.1) then
            v1 = ve
            v2 = vn
            v3 = vu
         endif
         ve = ve-v1
         vn = vn-v2
         vu = vu-v3
         if (ioptn.eq.1) write (line(24:41),'(2f9.3)') ve,vn
         if (ioptn.eq.2) 
     .      write (line(69:92),'(3f8.4)') ve,vn,vu
         if (ioptn.eq.3) write (line(17:32),'(2f8.2)') ve,vn
         if (ioptn.eq.3) write (line(73:80),'(f8.2)') vu
         if (ioptn.eq.4) write (line(20:35),'(2f8.2)') ve,vn
         if (ioptn.eq.4) write (line(76:83),'(f8.2)') vu
         write (12,'(a)') line(1:lcmd)
c
 20   continue  

      goto 100

 1000 print *,' Can not open the file :'
      stop 'in NET_UPDATE. '

 100  continue
      close (11)
      close (12)
c
      stop
      end
