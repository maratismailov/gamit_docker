      subroutine read_hfile_dat(ifil1,ifil2,mode)
c
c     1. read multiple GAMIT hfiles and create FONDA hfiles
c     2. create a summary list to output file
c
c     ifil1 :  input data file 
c     ifil2 :  output observation data file
c     ifil3 :  site list file
c
      implicit real*8(a-h,o-z)
      include 'maked.fti'

      character*64 subfile,line,subhfl
      character*128 command,command2
      character*6  fmt1,fmt2,l1,l2
      integer mode,len,ifile,k,j,i,lift_arg,j1,i1
      integer ifil1,ifil2,ifil3,loop,ilst,ierr
      integer istart,icol,n1,n0,n2,k1,itp,jobs,lcmd
      integer iyear,month,iday,ihour,imin,iff
      integer julday,match_name,nblen,l
      integer junk_sit,junk_lst(200),j2,jk,j3,jl,isit
      character*4  tail,st1,st_list(200)
      character*8  sitnam(100),stnm1,glb_st(200)
      character*30 note
      character*80 wcmd
      dimension    temp(600)
      logical      old

      print*,' Begining to extract site adjusts from h-files ....'
      mode = 31
      loop =4
      if (in_opt(1:8).eq.'tight_fr') loop = 1
      if (in_opt(1:8).eq.'tight_fi') loop = 2
      if (in_opt(1:8).eq.'loose_fr') loop = 3

c
c     get all explicitely designed site names
      ilst = 0
      tail = '_GPS'
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
         call blank(command)

c        get the number counter of the line before the station list
         command(1:17) = 'grep -n -e track '
         command(18:18+len) = subfile(1:len)
         command(19+len:30+len) = ' > tmp_file'
         
         ierr = system(command)

c        get number of stations from hfile
         open (22,file=subfile(1:len),status='old',err=300)
         call getcmd(22,'stations',wcmd,lcmd,1) 
         read (wcmd,'(i4)') n0
         close (22)

         open (21,file='tmp_file',status='old',err=300)
         istart = 0


 30      call blank(command)
         call blank(fmt1)
         call blank(fmt2)
     
c        get start line
         read(21,'(a64)',end=30) line
         close (21)
         icol = index(line,':')
         read(line(1:icol-1),'(i4)') n1
         
c        calc end

         n2 = n1 + n0  
  
c        convert to ascii and remove spaces
         write(fmt1,'(i4)') n1
         write(fmt2,'(i4)') n2
         l = lift_arg(fmt1,l1,1)
         l = lift_arg(fmt2,l2,1)

c        command = 'head -'//l2//' '//subfile(1:len)//' | '
         command2(1:1+6+nblen(l1)+len+12) = 'tail +'//l1(1:nblen(l1))//
     .     ' '//subfile(1:len)//' > tmp_fil2'

         ierr = system(command2)

c        get site names
         open (21,file='tmp_fil2',status='old',err=300)
c        skip first line
         read (21,'(2x)')
         isit = 0
         do 35 j = 1,200
            read(21,'(a64)',end=20) line
            read (line,'(i4,2x,a4)',err=40) i1,st1
c           check site name list
            if (ilst.gt.0) then
               k1 = match_name(ilst,4,st_list,st1)
               if (k1.le.0) goto 35
            endif
            isit = isit+1
            sitnam(isit) = st1 // tail
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
         itp = 31
         jobs = isit
         note = '{experiment number            '
         write (22,50) iexp,note
         note = '{exp. index, obs. type, number'
         write (22,140) iexp,itp,jobs,note
c        use start time as the observation time
         call getcmd(21,'Start',wcmd,lcmd,1)
         if (lcmd.gt.0) then
            read(wcmd,*) iyear,month,iday,ihour,imin,sec
            iyear = iyear + 1900
            time1 = 1900.0d0+julday(month,iday,iyear,1)/365.2422d0
            time1 = time1 + 
     .              (dble(ihour)+dble(imin)/60.0)/24.0d0/365.2422d0
         else
            print*,' Can not get the observation time.'
         endif
         call getcmd(21,'Prefit',wcmd,lcmd,2)
c        skip next two lines
         read (21,'(2x,/,2x)')
c        transfer coordinate adjustments
c        change the unit of radius from km to meter
         junk_sit = 0
c        do 60 j = 1,200
         do 60 j = 1,n0
            read (21,250,err=90) i1,symbl,st1,note(1:16),stla,sadj
            stla = stla+sadj
            read (21,250,err=90) i1,symbl,st1,note(1:16),stlo,sadj
            stlo = stlo+sadj
            read (21,250,err=90) i1,symbl,st1,note(1:16),stra,sadj
c           check site name list
            k1 = j
            if (ilst.gt.0) then
               k1 = match_name(isit,4,sitnam,st1)
               if (k1.le.0) then
                  junk_sit = junk_sit+1
                  junk_lst(junk_sit) = j
                  goto 60
               endif
            endif
            stra = (stra+sadj)*1.0d3
            write (22,160) time1,stla,stlo,stra,sitnam(k1)
 60      continue
c        transfer covariance submatrix
c        rescaling the radius part to fit the unit of meter
 90      call getcmd(21,'Covariance',wcmd,lcmd,2)
c        do 70 j = 1,600
c        three vcv lines for each site in hfile
         do 70 j = 1,(n0*3)
            j1 = j - j/3*3
            read (21,110,end=180,err=180) i1,(temp(k),k=1,j)
c           scale vcv up as well
            do k = 3,j,3
               temp(k) = temp(k)*1.0d3
            enddo
c           scale radius vcv to mm as well
            if (j1.eq.0) then
               do k = 1,j
                  temp(k) = temp(k)*1.0d3
               enddo
            endif
c           remove part of junk sites
            j2 = j
            if (junk_sit.gt.0) then
               j3 = 0
               do 130 jk = 1,j
                  do 120 k = 1,junk_sit
                     i = junk_lst(k)*3
                     if (j.gt.i-3.and.j.le.i) goto 70
                     if (jk.gt.i-3.and.jk.le.i) then
                        do jl = jk-j3,j-1-j3
                           temp(jl) = temp(jl+1)
                        enddo 
                        j3 = j3+1
                     endif
 120              continue
 130           continue
               j2 = j2-j3
            endif
            write(22,110) j2,(temp(k),k=1,j2)
 70      continue
c        remove temporary files
 180     call blank(command)
         command(1:24) = '\rm -f tmp_file tmp_fil2'
         ierr = system(command)
          
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
 160  format (f9.4,2(1x,f23.16,1x),1x,f23.13,2x,a8)
 250  format (1x,i4,a1,a4,a16,1x,f23.16,4x,d23.16)
      goto 1000

 300  print *,'READ_HFILE_DAT: error in reading file: ',subfile
      stop 'MAKED: READ_HFILE_DAT ...'
         
 1000 continue
      return
      end
