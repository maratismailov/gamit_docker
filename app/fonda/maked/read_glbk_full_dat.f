      subroutine read_glbk_full_dat(ifil1,ifil2,mode)
c
c     read GLORG output full covariance file and create FONDA file
c
c     ifil1 :  input data file 
c     ifil2 :  output observation data file
c     ifil3 :  site list file
c
      implicit real*8(a-h,o-z)
      include 'maked.fti'

      character*6  fmt1
      integer mode,len,k,j,i,lift_arg,i1,idx
      integer iobs,ll,k1,j1,j2,isit
      integer ifil1,ifil2,lcmd,ilst,ifil3,itp,jobs 
      integer nblen,match_name
      character*8  st1,st_list(400)
      character*8  sitnam(400),stnm1
      character*30 note
      character*80 wcmd
      character*80 line
      dimension    temp(2400)
      logical      find

      print*,' Begining to extract site velocity adjusts ....'
c     mode = 34
c     get reference time as the observation time
      call getcmd(ifil1,'Solution refers',wcmd,lcmd,1)
      if (lcmd.gt.0) then
         j = nblen(wcmd)
c        GLOBK format has changed with the addition of [Seconds tag  59.080] 
c         or the like. Need to check length of line first
         if (j.gt.30) then
            note(1:9) = wcmd(j-31:j-23) 
            read (note,'(f9.4)') time1
         else
            note(1:9) = wcmd(j-8:j)
            read (note,'(f9.4)') time1
         endif
         print*,' get reference time:',time1
      else
         print*,' Can not get the observation time.'
      endif
c
c     get all explicitely designed site names
      ilst = 0
      if (site_list(1:1).ne.'*') then
         ifil3 = 10
         open (10,file=site_list,status='old',err=300)
         do 10 i = 1,1000
            read (ifil3,'(a)',end=15) stnm1
            j = lift_arg(stnm1,st_list(i),1)
            ilst = ilst+1
 10      continue
         close (10)
      endif
c     
 15   nsit = 0
      rewind (ifil1)
      call blank(note)
      note(1:12) = 'X coordinate'
      len = nblen(infil)
c     screen display to avoid sleep
      print*,'   Processing  <<<  ',infil(1:len)
c     find the line with coordinate adjustment
      find = .false.
      do i = 1,1000
         read (ifil1,'(a80)') wcmd
         if (wcmd(16:27).eq.note(1:12)) then
            find = .true.
            goto 17
         endif
      enddo
      if (.not.find) goto 400
 17   open (21,file='tmp_file',status='unknown',err=300)
      idx = 0
c     magic number 10000 is the max number of lines in the GLORG file
      do 20 isit = 1,100000
         read (wcmd,'(a6,a8)') fmt1,st1
         if (fmt1.eq.' CORRE') goto 40
         if (wcmd(16:27).ne.note(1:12)) goto 35
c        check site name list
         if (ilst.gt.0) then
            len = nblen(st1)
            k1 = match_name(ilst,len,st_list,st1)

c           no match with st_list, so try next line
            if (k1.le.0) then
               read (ifil1,'(2x,/,2x)')
c              read (ifil1,'(2x,/,2x,/,2x,/,2x,/,2x,/,2x,/,2x)')
c              read (ifil1,'(2x,/,2x,/,2x,/,2x,/,2x,/,2x,/,2x)')
c              read (ifil1,'(2x,/,2x,/,2x,/,2x,/,2x)')
               goto 35
            endif

         endif
         nsit = nsit+1
         if (len.lt.8) then
            do k = len+1,8
               st1(k:k) = '_'
            enddo
         endif
         sitnam(nsit) = st1 
c        skip next two lines to reach velocity line
         read (ifil1,'(2x,/,2x)')

cmk      need to consider case of no velocity for a site 
         call blank(line)
         read (ifil1,'(a80)') line
         if (line(1:4).ne.'Unc.') then
             read (line,60) igobs(idx+1),stvx,erx
             read (ifil1,60) igobs(idx+2),stvy,ery
             read (ifil1,60) igobs(idx+3),stvz,erz
c        skip remaining lines
             read (ifil1,'(2x,/,2x,/,2x,/,2x,/,2x,/,2x,/,2x)')
             read (ifil1,'(2x,/,2x,/,2x,/,2x,/,2x,/,2x,/,2x)')
c            read (ifil1,'(2x,/,2x,/,2x)')
         else
c            velocities are zero
             igobs(idx+1) = 0.0d0
             igobs(idx+2)=0.0d0
             igobs(idx+3)=0.0d0
             stvx=0.0d0
             stvy=0.0d0
             stvz=0.0d0
             erx=0.0d0
             ery=0.0d0
             erz=0.0d0
             
c            skip remaining lines
             read (ifil1,'(2x,/,2x,/,2x,/,2x,/,2x)')
         endif

c        scaling (see Feigl et al. 1993)
c         scl = 2.0
c        scaling (see Dong 1993)
         scl = 3.0
         erd(idx+1) = erx*scl
         erd(idx+2) = ery*scl
         erd(idx+3) = erz*scl
         idx = idx+3
         write (21,160) time1,stvx,stvy,stvz,st1
 35      read (ifil1,'(a80)') wcmd
 20   continue
c
      print*,' Can not find CORELATION but ...',wcmd
c     
c     output site names 
 40   call blank(note)
      note = '{network site number          '
      write (ifil2,50) nsit,note
      write (ifil2,'(a8)') (sitnam(j),j=1,nsit)
c     experiment number (always 1)
      iexp = 1
c     Cartesian velocity
      itp = 34
      jobs = nsit
      note = '{experiment number            '
      write (ifil2,50) iexp,note
      note = '{exp. index, obs. type, number'
      write (ifil2,140) iexp,itp,jobs,note
c     transfer coordinate adjustments
      rewind (21)
      do 90 i = 1,nsit
         read (21,'(a80)') wcmd
         write (ifil2,'(a80)') wcmd
 90   continue
      close (21)

c     transfer covariance submatrix
      iobs = igobs(idx)
      ll = 1
      j2 = igobs(ll)
      erx = erd(ll)
      do 70 j = 1,iobs
         if (j.lt.j2) then
            j1 = (j-1)/10+1
            do k = 1,j1
               read (ifil1,'(2x)') 
            enddo
            goto 70
         else
            read (ifil1,110,end=70,err=70) i1,(temp(k),k=1,j)
         endif
c        remove junk terms 
         do 130 j1 = 1,j
            do 120 k = 1,ll
               i = igobs(k)
               ery = erd(k)
               if (j1.ne.i) goto 120
               temp(k) = temp(j1)*ery*erx
 120        continue
 130     continue
         write(ifil2,250) ll,(temp(k),k=1,ll)
         ll = ll+1
         j2 = igobs(ll)
         erx = erd(ll)
 70   continue
c
 50   format (i5,25x,a30)
 60   format (i4,39x,f13.4,13x,f13.4)
 140  format (3i5,15x,a30)
 110  format(i5,'. ', 10(1x,f7.4),:/,200(7x,10(1x,f7.4),:/) )
 150  format (f9.4,3(1x,f23.16,1x),1x,a8)
 160  format (f9.4,3(1x,f13.5,1x),2x,a8)
 250  format(i5,'. ', 10(1x,d11.4),:/,200(7x,10(1x,d11.4),:/) )
      goto 1000

 300  print *,'READ_GLOBK_FULL_DAT: error in reading file: ',infil
      stop 'MAKED: READ_GLOBK_FULL_DAT ...'

 400  print *,'READ_GLOBK_FULL_DAT: missing site: ',
     .         note(1:12),wcmd(16:27)
      stop 'MAKED: READ_GLOBK_FULL_DAT ...'
         
 1000 continue
      return
      end
