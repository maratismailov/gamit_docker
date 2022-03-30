      subroutine read_gipsy_dat(ifil1,ifil2,mode)
c
c     1. read multiple GIPSY sta.cov files and create FONDA files
c     2. create a summary list to output file
c
c     ifil1 :  input data file 
c     ifil2 :  output observation data file
c     ifil3 :  site list file
c
      implicit real*8(a-h,o-z)
      include 'maked.fti'

      character*64 subfile,subgfl
      character*6  fmt1,fmt2
      integer mode,len,ifile,k,j,i,lift_arg,j1,i1
      integer ifil1,ifil2,ifil3,loop,ilst,iff,ipp,iyear,imonth,iday
      integer ipar,k1,n2,n1,itp,jobs
      integer nblen,match_name
      integer junk_par,junk_lst(600),j2,jk,j3,isit
      integer*4 julday
      character*3  month
      character*4  tail,st1,st_list(200)
      character*8  sitnam(100),stnm1,glb_st(200)
      character*30 note
      character*80 line
      dimension    temp(600),sigm(600)
      logical      old

      print*,' Begining to extract site adjusts from gipsy files ....'
      mode = 33
      loop =1
      ifile = 0
c
c     get all explicitely designed site names
      ilst = 0
      tail = '_GPS'
      if (site_list(1:1).ne.'*'.and.site_list(1:1).ne.'') then
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
         if (len.le.1) goto 20
c        screen display to avoid sleep
         print*,'   Processing  <<<  ',subfile(1:len)
         ifile = ifile+1
         open (21,file=subfile(1:len),status='old',err=300)
c        get parameter number
         read(21,'(a80)',end=20) line
         len = nblen(line)
         read(line(1:len),'(i5)',err=20) ipp
c get date and convert to fonda format
         read (line,'(20x,i2,a3,i2)',err=700) iyear,month,iday
c modify yymondd to epoch(year.partofyear)
         if (month.eq.'JAN') then
            imonth=1
         else if (month.eq.'FEB') then
            imonth=2
         else if (month.eq.'MAR') then
            imonth=3
         else if (month.eq.'APR') then
            imonth=4
         else if (month.eq.'MAY') then
            imonth=5
         else if (month.eq.'JUN') then
            imonth=6
         else if (month.eq.'JUL') then
            imonth=7
         else if (month.eq.'AUG') then
            imonth=8
         else if (month.eq.'SEP') then
            imonth=9
         else if (month.eq.'OCT') then
            imonth=10
         else if (month.eq.'NOV') then
            imonth=11
         else
            imonth=12
         end if
c     calling function julday  ~/fonda/com/. to do the conversion
         stime = 1900.0d0+julday(imonth,iday,iyear,1)/365.2422d0
c        parameter loop
         isit = 0
         ipar = 0
         junk_par = 0
         do 40 j = 1,ipp
            read(21,'(a80)',end=20) line
            read (line,250,err=40) i1,st1,note(1:13),add1,add2
c           check site name list
            if (ilst.gt.0) then
               k1 = match_name(ilst,4,st_list,st1)
               if (k1.le.0) then
                  junk_par = junk_par+1
                  junk_lst(junk_par) = j
                  goto 40
               endif
            endif
            if (isit.gt.0) then
               k1 = match_name(isit,4,sitnam,st1)
               if (k1.gt.0) goto 35
            endif
            isit = isit+1
            sitnam(isit) = st1 // tail
 35      continue
c        store adjustment and sigma
         ipar = ipar+1
         temp(ipar) = add1
         sigm(ipar) = add2
c        add site name to global site list
         do j1 = 1,isit
            old = .false.
            do k = 1,nsit
               if (sitnam(j1).eq.glb_st(k)) old = .true.
            enddo
            if (.not.old) then
               nsit = nsit+1
               glb_st(nsit) = sitnam(j1)
            endif
         enddo
 40      continue
c        output site names to sub_gfile
         call blank(fmt1)
         call blank(fmt2)
         call blank(subgfl)
         write(fmt2,'(i3)') ifile
         n2 = lift_arg(fmt2,fmt1,1)
         n1 = nblen(outfil)
         subgfl(1:n1+n2+1) = outfil(1:n1) // '_' // fmt1(1:n2)
         open (22,file=subgfl,status='unknown',err=300)
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
         do 60 j1 = 1,ipar/3
            j = (j1-1)*3+1
            write (22,160) stime,temp(j),temp(j+1),temp(j+2),sitnam(j1)
 60      continue
c        transfer covariance submatrix
         do 70 j2 = 1,ipp
            do 90 j = 1,j2
               if (j.eq.j2) then
                  temp(j) = sigm(j)**2
               else
                  read (21,'(13x,e22.15)',err=70) temp(j)
                  temp(j) = temp(j)*sigm(j2)*sigm(j)
               endif
 90         continue
c           remove part of junk parameters
            if (junk_par.gt.0) then
               j3 = 0
               do 130 jk = 1,j2
                  do 120 k = 1,junk_par
                     i = junk_lst(k)
                     if (jk.eq.i) goto 130
 120              continue
                  j3 = j3+1
                  temp(j3) = temp(jk)
 130           continue
            else
               j3 = j2
            endif
            write(22,110) j3,(temp(k),k=1,j3)
 70      continue
          
 20   continue
c
c     output global site name list and sub-hfile list
 200  note = '{network site number          '
      write (ifil2,50) nsit,note
      write (ifil2,'(a8)') (glb_st(j),j=1,nsit)
      note = '{sub-gipsy_file number       '
      write (ifil2,50) ifile,note
      rewind (ifil1)
c     output file name list
      n1 = nblen(outfil)
      call blank(subgfl)
      subgfl(1:n1+1) = outfil(1:n1) // '_' 
      do 80 iff = 1,ifile
         call blank(fmt1)
         call blank(fmt2)
         write(fmt2,'(i3)') iff
         n2 = lift_arg(fmt2,fmt1,1)
         subgfl(n1+2:n1+n2+1) = fmt1(1:n2)
         write (ifil2,'(a)') subgfl(1:n1+n2+1)
 80   continue

 50   format (i5,25x,a30)
 140  format (3i5,15x,a30)
 110  format (1x,i3,'. ',50(5(1x,d23.16),:,/,6x))
 160  format (f9.4,3(1x,f22.13,1x),1x,a8)
 250  format (1x,i4,2x,a4,a13,1x,e22.15,6x,e21.15)
      goto 1000

 300  print *,'READ_GIPSY_DAT: error in reading file: ',subfile
      stop 'MAKED: READ_GIPSY_DAT ...'
 700  print *,'700 error in reading file: '
         
 1000 continue
      return
      end
