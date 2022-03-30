      subroutine read_fondah_dat(ifil1,ifil2,mode)
c
c     1. read multiple FONDA-output hfiles and create FONDA-readable hfiles
c     2. create a summary list to output file
c
c     ifil1 :  input data file 
c     ifil2 :  output observation data file
c     ifil3 :  site list file
c
      implicit real*8(a-h,o-z)
      include 'maked.fti'

      character*64 subfile,line,subhfl
      character*6  fmt1,fmt2
      integer mode,len,ifile,k,j,i,lift_arg,j1,i1,jcut
      integer ifil1,ifil2,ifil3,ilst,iff,lcmd,k1,n1,n2
      integer itp,jobs,iyear,month,iday,ihour,imin 
      integer junk_sit,junk_lst(200),j2,jk,j3,jl,isit
      integer nblen,match_name
      integer*4 julday
      character*8  st1,st_list(200)
      character*8  sitnam(100),stnm1,glb_st(200)
      character*30 note
      character*80 wcmd
      dimension    temp(600)
      logical      old

      print*,' Begining to extract site adjusts from h-files ....'
      mode = 33
      if (in_opt(1:3).eq.'uvw') mode = 34
c
c     get all explicitely designed site names
      ilst = 0
      if (site_list(1:1).ne.'*') then
         ifil3 = 10
         open (10,file=site_list,status='old',err=300)
         do 10 i = 1,1000
            read (ifil3,'(a)',end=71) stnm1
            j = lift_arg(stnm1,st_list(i),1)
            ilst = ilst+1
 10      continue
         close (10)
      endif
c     
 71   continue

 15   nsit = 0
      rewind (ifil1)
c     get input file name list
      do 20 iff = 1,200
         read (ifil1,'(a)',err=300,end=200) subfile
         len = nblen(subfile)
c        screen display to avoid sleep
         print*,'   Processing  <<<  ',subfile(1:len)
         ifile = ifile+1
         open (21,file=subfile(1:len),status='old',err=300)
c        get site names
         call getcmd(21,'track',wcmd,lcmd,2)
         isit = 0
         do 35 j = 1,200
            read(21,'(a64)',end=20) line
            read (line,'(i4,12x,a8)',err=40) i1,st1
c           check site name list
            if (ilst.gt.0) then
               k1 = match_name(ilst,8,st_list,st1)
               if (k1.le.0) goto 35
            endif
            isit = isit+1
            sitnam(isit) = st1
 35      continue
c        add site name to global site list
 40      if (ifile.eq.1) then
            do j = 1,isit
               glb_st(j) = sitnam(j)
            enddo
            nsit = nsit+isit
         else
            do j = 1,isit
               old = .false.
               do k = 1,nsit
                  if (sitnam(j).eq.glb_st(k)) old = .true.
               enddo
               if (.not.old) then
                  nsit = nsit+1
                  glb_st(nsit) = sitnam(j)
               endif
            enddo
         endif
c        output site names to sub_hfile
         call blank(fmt1)
         call blank(fmt2)
         call blank(subhfl)
         write(fmt2,'(i3)') ifile
         n2 = lift_arg(fmt2,fmt1,1)
         n1 = nblen(outfil)
         subhfl(1:n1+n2+1) = outfil(1:n1) // '_' // fmt1(1:n2)
         open (22,file=subhfl,status='unknown',err=300)
         note = '{network site number          '
         write (22,50) isit,note
         write (22,'(a8)') (sitnam(j),j=1,isit)
c        experiment number (always 1)
         iexp = 1
         itp = mode
         jobs = isit
         note = '{experiment number            '
         write (22,50) iexp,note
         note = '{exp. index, obs. type, number'
         write (22,140) iexp,itp,jobs,note
c        use start time as the observation time
         call getcmd(21,'Start',wcmd,lcmd,1)
         if (lcmd.gt.0) then
            read(wcmd,*) iyear,month,iday,ihour,imin,sec
            time1 = 1900.0d0+julday(month,iday,iyear,1)/365.2422d0
            time1 = time1 + 
     .              (dble(ihour)+dble(imin)/60.0)/24.0d0/365.2422d0
         else
            print*,' Can not get the observation time.'
         endif
         call getcmd(21,'Prefit',wcmd,lcmd,2)
c        skip next two lines
         read (21,'(2x,/,2x)')
c        transfer coordinate and velocity adjustments
         junk_sit = 0
         do 60 j = 1,200
c           skip coordinate lines in 'uvw' mode
            if (mode.eq.34) then
               read(21,'(2x/2x/2x)')
               goto 52
            endif
            read (21,250,err=90) i1,symbl,st1,note(1:16),stx,sadj
            stx = stx+sadj
            read (21,250,err=90) i1,symbl,st1,note(1:16),sty,sadj
            sty = sty+sadj
            read (21,250,err=90) i1,symbl,st1,note(1:16),stz,sadj
            stz = stz+sadj
c           skip velocity lines in 'xyz' mode
            if (mode.eq.33) then
               read(21,'(2x/2x/2x)')
               goto 54
            endif
 52         read (21,250,err=90) i1,symbl,st1,note(1:16),stu,sadj
            stu = stu+sadj
            read (21,250,err=90) i1,symbl,st1,note(1:16),stv,sadj
            stv = stv+sadj
            read (21,250,err=90) i1,symbl,st1,note(1:16),stw,sadj
            stw = stw+sadj
c           check site name list
 54         k1 = j
            if (ilst.gt.0) then
               k1 = match_name(isit,8,sitnam,st1)
               if (k1.le.0) then
                  junk_sit = junk_sit+1
                  junk_lst(junk_sit) = j
                  goto 60
               endif
            endif
            if (mode.eq.33) write (22,160) time1,stx,sty,stz,sitnam(k1)
            if (mode.eq.34) write (22,160) time1,stu,stv,stw,sitnam(k1)
 60      continue
c        transfer covariance submatrix
 90      call getcmd(21,'Covariance',wcmd,lcmd,2)
         do 70 j = 1,600
            j1 = j - (j-1)/6*6
            j2 = j
            read (21,110,end=20,err=20) i1,(temp(k),k=1,j)
c           remove coordinate part or velocity part
            if (mode.eq.33.and.j1.ge.4) goto 70
            if (mode.eq.34.and.j1.lt.4) goto 70
            j3 = 0
            do 130 jk = 1,j
               jcut = jk-(jk-1)/6*6
               if (mode.eq.33.and.jcut.lt.4) goto 130
               if (mode.eq.34.and.jcut.ge.4) goto 130
               do jl = jk-j3,j-1-j3
                  temp(jl) = temp(jl+1)
               enddo 
               j3 = j3+1
 130        continue
            j2 = j2-j3
            if (j2.le.0) goto 70
c           remove part of junk sites
            if (junk_sit.gt.0) then
               j3 = 0
               jcut = j2
               do 170 jk = 1,jcut
                  do 120 k = 1,junk_sit
                     i = junk_lst(k)*3
                     if (jcut.gt.i-3.and.jcut.le.i) goto 70
                     if (jk.gt.i-3.and.jk.le.i) then
                        do jl = jk-j3,jcut-1-j3
                           temp(jl) = temp(jl+1)
                        enddo 
                        j3 = j3+1
                     endif
 120              continue
 170           continue
               j2 = j2-j3
            endif
            write(22,110) j2,(temp(k),k=1,j2)
 70      continue
 20   continue
c
c     output global site name list and sub-hfile list
 200  note = '{network site number          '
      write (ifil2,50) nsit,note
      write (ifil2,'(a8)') (glb_st(j),j=1,nsit)
      note = '{sub-hfile number            '
      write (ifil2,50) ifile,note
      rewind (ifil1)
c     output file name list
      n1 = nblen(outfil)
      call blank(subhfl)
      subhfl(1:n1+1) = outfil(1:n1) // '_' 
      do 80 iff = 1,ifile
         call blank(fmt1)
         call blank(fmt2)
         write(fmt2,'(i3)') iff
         n2 = lift_arg(fmt2,fmt1,1)
         subhfl(n1+2:n1+n2+1) = fmt1(1:n2)
         write (ifil2,'(a)') subhfl
 80   continue

 50   format (i5,25x,a30)
 140  format (3i5,15x,a30)
 110  format (1x,i3,'. ',50(5(1x,d23.16),:,/,6x))
 150  format (f9.4,3(1x,f23.16,1x),1x,a8)
 160  format (f9.4,2(1x,d23.16,1x),1x,d23.16,2x,a8)
 250  format (1x,i4,a1,a8,a16,1x,d23.16,4x,d23.16)
      goto 1000

 300  print *,'READ_HFILE_DAT: error in reading file: ',subfile
      stop 'MAKED: READ_HFILE_DAT ...'
         
 1000 continue
      return
      end
