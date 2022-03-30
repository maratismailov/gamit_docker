      program rem_model
c
c     remove model velocity solution from regular velocity solution
c
c     options:
c        1 = solvem mapping file
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
      character*1 sym1
      character*8 sit_nam, tmpn
      character*64 infile,outfile,modfile
      integer i,is,l1,l2,i1,lift_arg
      integer ios,lcmd,ioptn,len,iyr,imon,iday,ihr,imn,isec
      character*128 line
      integer nblen
      integer*4 iclarg
c
c     check if the runstring is complete.
      len = iclarg(1,infile)
      l2 = iclarg(2,outfile)
      l1 = iclarg(3,modfile)
      if (l1.le.0.or.l2.le.0.or.len.le.0) then
         print*,' Runstring: rem_model <infile> <outfile>',
     .      ' <modfile> <option>'
         print*,'  option: 1=remove model vel. from map file'
         print*,'          2=velocity res. + covariance diff.'
         print*,'          3=(undecided)'
         stop
      endif
      i = iclarg(4,tmpn)
      if (i.le.0) then
         ioptn = 1
      else
         read (tmpn,*) ioptn
         if (ioptn.le.0.or.ioptn.gt.2) ioptn = 1
      endif
c
c     open file
      open (11,file=infile,status='old',err=1000)
      open (12,file=outfile,status='unknown',err=1000)
      open (13,file=modfile,status='old',err=1000)
      print*,' Begin to remove model velocity solution ...'
      f1 = 2.45d0
c
c     set up BIG DO loop to read until end of file...
      rewind (11)
      do 20 i = 1,1000
         rewind (13)
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
               write (12,'(a16,a,a5,a,a8,i4,2(a1,i2),3x,i2,2(a1,i2))')
     .            '* dif. between: ',infile(1:len),' and ',
     .            modfile(1:l1),'  time: ',
     .            iyr,'/',imon,'/',iday,ihr,':',imn,':',isec
            endif
            write (12,'(a)') line(1:lcmd)
            goto 20
         endif
         read (line(1:lcmd),*,err=20) 
     .      slon,slat,ve,vn,sige,sign,corr
         is = lift_arg(line,sit_nam,8)
 14         format (a1,4f14.8,2f10.3,f9.4,3x,a8)
c        check if the name in the model velocity file
         do 30 is = 1,1000
            ios = 0
c           read site coordinte and velosity
            read (13,'(a)',iostat=ios,end=40,err=30) line
            lcmd = nblen(line)
c           skip comment lines
            if (line(1:1).ne.' '.or.lcmd.le.1) goto 30

            read (line(1:lcmd),*,err=30) 
     .      slo2,sla2,ve2,vn2,sige2,sign2,corr2
            i1 = lift_arg(line,tmpn,8)

         if (tmpn .eq. sit_nam) then
            sym1 = ' '
            ve = ve-ve2
            vn = vn-vn2
            if (ioptn.eq.2) then
               if (sige.lt.sige2.or.sign.lt.sign2) sym1 = '?'
               fe = sige**2-sige2**2
               fn = sign**2-sign2**2
               if (fe.gt.0.0d0) fe = dsqrt(fe)
               if (fn.gt.0.0d0) fn = dsqrt(fn)
               if (fe.lt.1.0d-2) fe = 1.0d-2
               if (fn.lt.1.0d-2) fn = 1.0d-2
               if (dabs(ve).gt.fe*f1.or.dabs(vn).gt.fn*f1)
     .            sym1 = '!'
               temp = corr*sige*sign-corr2*sige2*sign2
               if (fe.le.1.0d-2.or.fn.le.1.0d-2) then
                  corr = corr-corr2
               else
                  corr = temp/fe/fn
               endif
               sige = fe
               sign = fn
            endif
c
            write (12,14)
     .      sym1,slon,slat,ve,vn,sige,sign,corr,sit_nam
            goto 20
         endif

 30      continue
 40      if (ioptn.eq.1) sym1 = '*'
         if (ioptn.eq.2) sym1 = '$'
         write (12,14)
     .      sym1,slon,slat,ve,vn,sige,sign,corr,sit_nam

 20   continue  

      goto 100

 1000 print *,' Can not open the file :'
      stop 'in REM_MODEL. '

 100  continue
      close (11)
      close (12)
      close (13)
c
      stop
      end
