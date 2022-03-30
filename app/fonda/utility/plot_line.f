      program plot_line
c
c     plot multiple figures for baseline observations
c
c     options:
c        1 = 1 column plot format
c        2 = 2 columns plot format
c        3 = 3 columns plot format
c        4 = (anything else?)
c        5 = (anything else?)
c        6 = (anything else?)
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
      character*8 tmpn,sit1,sit2,ref1(400),ref2(400)
      character*64 infile,outfile,listfile,gmtfile
      character*120 command
      integer i,type,l1,l2,n,len,lrow,lplt,l3,j,j1
      integer ios,lcmd,ioptn,list,lpage,j2
      character*120 line
      character*8 tmp_d, tmp_t
      integer nblen,lift_arg
      integer*4 iclarg
      logical pass
c
c     check if the runstring is complete.
      len = iclarg(1,infile)
      l2 = iclarg(2,outfile)
      l3 = iclarg(4,gmtfile)
      l1 = iclarg(3,listfile)
      if (l1.le.0.or.l2.le.0.or.len.le.0.or.l3.le.0) then
         print*,' Runstring: plot_line <infile> <psfile>',
     .      ' <listfile> <gmtfile> <option>'
         print*,'  <listfile> = <all>  means default (all lines)'
         print*,'  option: 1= 1 column   2= 2 columns'
         print*,'  option: 3= 3 columns  4= (undecided)'
         stop
      endif
      i = iclarg(5,tmpn)
      if (i.le.0) then
         ioptn = 1
      else
         read (tmpn,*) ioptn
         if (ioptn.le.0.or.ioptn.ge.4) ioptn = 1
      endif
c
c     open file
      open (11,file=infile,status='old',err=1000)
      list = 1
      if (listfile(1:3).eq.'all') list = 0
      if (list.gt.0) then
         open (13,file=listfile,status='old',err=1000)
         j1 = 0
         do 60 j = 1,1000
            call blank(command)
            read (13,'(a)',end=40) command
            if (command(1:1).eq.'*') goto 60
            j2 = lift_arg(command,sit1,1)
            if (j2.le.0) goto 60
            j2 = lift_arg(command,sit2,2)
            if (j2.le.0) goto 60
            j1 = j1+1
            ref1(j1) = sit1
            ref2(j1) = sit2
 60      continue
 40      list = j1
         close (13)
      endif   
      print*,' Begin to plot multiple baseline data .....'
c
c     set plot number per page
      lrow = 4
      lplt = ioptn
      if (ioptn.eq.3) lrow = 5
      lpage = lrow*lplt
c
c     set up BIG DO loop to read until end of file...
      n = 0
      tmp_d(1:6) = 'tmp_d.'
      tmp_t(1:6) = 'tmp_t.'
      pass = .true.
      do 20 i = 1,100000
         call blank(line)
         call blank(command)
         ios = 0
         read (11,'(a)',iostat=ios,end=50,err=20) line
         lcmd = nblen(line)
         if (lcmd.le.1) goto 20
         if (line(3:10).eq.'Baseline'.or.line(3:10).eq.'Residual') then
            if (list.gt.0) then
               read(line(14:33),'(a8,4x,a8)') sit1,sit2
               do 30 j = 1,list
                  if (sit1.ne.ref1(j).or.sit2.ne.ref2(j)) goto 30
                  pass = .true.
                  goto 70
 30            continue
               pass = .false.
               goto 20
            endif
 70         n = n+1
            type = 1
            if (line(3:10).eq.'Residual') type = 2
            if (n.le.9) then
               write(tmp_d(7:8),'(i1,'' '')') n
               write(tmp_t(7:8),'(i1,'' '')') n
            else
               write(tmp_d(7:8),'(i2)') n
               write(tmp_t(7:8),'(i2)') n
            endif
            open (14,file=tmp_d,status='unknown',err=1000)
            open (15,file=tmp_t,status='unknown',err=1000)
            if (ioptn.le.2) then
               write (15,'(a,a)') '0.9 1. 7 0 0 11 ',line(13:69)
               write (15,'(a,a)') '0.9 0.95 7 0 0 11 ',line(70:lcmd)
            else
               write (15,'(a,a)') '0.9 1. 5 0 0 11 ',line(13:69)
               write (15,'(a,a)') '0.9 0.95 5 0 0 11 ',line(70:lcmd)
            endif
            close (15)
            goto 20
         endif
         
         if (.not. pass) goto 20
         if (line(1:1).eq.' ') then
            write (14,'(a)') line
         else
            if (line(3:10).ne.'solution') goto 20
            close (14)
            if (n-n/lpage*lpage.eq.0) then
               call blank(command)
               command = gmtfile(1:l3)// ' ' //outfile(1:l2)// ' '
               write (command(3+l3+l2:8+l3+l2),'(i1,i2,i3)') 
     .            lrow,lplt,lpage
               command(9+l3+l2:12+l3+l2) = ' dat'
               if (type.eq.2) command(9+l3+l2:12+l3+l2) = ' res'
               l1 = system(command)
               n = 0 
            endif
            goto 20
         endif

 20   continue  

 50   if (n.gt.0) then
         i = n/lplt
         if (n-i*lplt .gt. 0) i = i+1
         call blank(command)
         command = gmtfile(1:l3)// ' ' //outfile(1:l2)// ' '
         write (command(3+l3+l2:8+l3+l2),'(i1,i2,i3)') i,lplt,n
         command(9+l3+l2:12+l3+l2) = ' dat'
         if (type.eq.2) command(9+l3+l2:12+l3+l2) = ' res'
         l1 = system(command)
      endif

      goto 100

 1000 print *,' Can not open the file :'
      stop 'in PLOT_LINE. '

 100  continue
      close (11)
      if (list.gt.0) close (13)
      l1 = system('\\rm -f tmp_d.* tmp_t.*')
c
      stop
      end
