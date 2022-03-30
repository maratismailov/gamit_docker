      subroutine read_ofile_dat(ifil1,ifil2,mode)
c
c     read GAMIT o-file and create FONDA format data file
c     
c     ifil1 :  input data file 
c     ifil2 :  output observation data file
c     ifil3 :  site list file
c
      implicit real*8(a-h,o-z)
      include 'maked.fti'

      character*4 tail,st1,stbas,st2,stmp,st_list(maxsit)
      character*8 sitnam(maxsit),stnm1,stnm2
      character*9 tbase
      character*30 note
      integer mode,iyr,iscal
      dimension sigma(4),rho(3),obs(4),iscal(maxsit)
      integer i,itp,jobs,issn,ifil1,ifil2,ifil3,ilst
      integer j,k1,k2
      integer lift_arg,match_name
      logical old1,old2
c
      small = 1.0d-5
c     
      print*,' Begining o-file data transfering ...'
      ifil3 = 0
      if (site_list(1:1).ne.'*') then
         ifil3 = 10  
         open (ifil3,file=site_list,status='old',err=102)
      endif
c
c     default format
      if (infmt(1:1).eq.'*')
c    . infmt(1:46) = '(a9,1x,f8.3,1x,4(3x,f13.4,4x,f8.4),32x,3f10.5)'
*    . infmt(1:46) = '(a9,1x,f8.3,1x,4(2x,f14.4,4x,f8.4),32x,3f10.5)'
     . infmt(1:46) = '(a9,1x,f8.3,1x,4(2x,f14.4,4x,f8.4),31x,3f10.5)'
c
c     get all designed site names
      ilst = 0
      if (ifil3.gt.0) then
         do 10 i = 1,3000
            read (ifil3,'(a)',end=70,err=70) stnm1
            j = lift_arg(stnm1,st_list(i),1)
            ilst = ilst+1
 10      continue
 70      close (ifil3)
      endif
c     
c     First, sort out all working sites and experiment number.
      nsit = 0
      tail = '_GPS'
      issn = 0
      jobs = 0
      time0 = 0.0d0
         print*,'now reading file'
      do 20 i = 1,30000
         read (ifil1,fmt=infmt,end=50,err=50) 
     .      tbase,time1,obs(1),sigma(1),obs(2),sigma(2),obs(3),sigma(3),
     .      obs(4),sigma(4),rho(1),rho(2),rho(3)
c        decomposite to get site name
         st1 = tbase(1:4)
         st2 = tbase(6:9)
c        check site name list
         if (ilst.gt.0) then
            k1 = match_name(ilst,4,st_list,st1)
            k2 = match_name(ilst,4,st_list,st2)
            if (k1.le.0.or.k2.le.0) goto 20
         endif
         if (issn.eq.0.or.time0.ne.time1) then
            stbas = st1
            time0 = time1
            issn = issn+1
            iscal(issn) = 1
         endif
         jobs = jobs+1
c        only use independent baseline vector
         if (st1.ne.stbas.and.time1.eq.time0) goto 20 
         iscal(issn) = iscal(issn)+1
         if (jobs.eq.1) then
            nsit = nsit+1
            sitnam(nsit) = st1 // tail
         endif
c        update effective site name list
         if (issn.eq.1) then
            nsit = nsit+1
            sitnam(nsit) = st2 // tail
         else
            old1 = .false.
            old2 = .false.
            do 40 j = 1,nsit
               stmp = sitnam(j)(1:4)
               if (st1.eq.stmp) old1 = .true.
               if (st2.eq.stmp) old2 = .true.
               if (old1.and.old2) goto 20
 40         continue
            if (.not.old1) then
               nsit = nsit+1
               sitnam(nsit) = st1 // tail
            endif
            if (.not.old2) then
               nsit = nsit+1
               sitnam(nsit) = st2 // tail
            endif
         endif
c        
 20   continue

 50   print*,' nsit,st_list,session=',nsit,ilst,issn
      if (nsit.lt.1) goto 120
c
c     site id
      note = '{network site number          '
      write (ifil2,30) nsit,note
      write (ifil2,'(a8)') (sitnam(i),i=1,nsit)
c
c     experiment number (always 1)
      i = 1
      itp = mode
      note = '{experiment number            '
      write (ifil2,30) i,note
      note = '{exp. index, obs. type, number'
      write (ifil2,140) i,itp,jobs,note

c     second, shift the format one by one
      rewind (ifil1)
      time0 = 0.0d0
      issn = 0
      do 60 i = 1,30000
         read (ifil1,fmt=infmt,end=120) 
     .      tbase,time1,obs(1),sigma(1),obs(2),sigma(2),obs(3),sigma(3),
     .      obs(4),sigma(4),rho(1),rho(2),rho(3)
c        decomposite to get site name
         st1 = tbase(1:4)
         st2 = tbase(6:9)
c        check site name list
         if (ilst.gt.0) then
            k1 = match_name(ilst,4,st_list,st1)
            k2 = match_name(ilst,4,st_list,st2)
            if (k1.le.0.or.k2.le.0) goto 60
         endif
         if (i.eq.1.or.time0.ne.time1) then
            time0 = time1
            issn = issn+1
            scalor = sqrt(dble(iscal(issn))/2.0d0)
         endif
c        Adopt Zhengkang's approach:
c           Using all baseline combinations, and enlarge the sigma
c           by a factor of sqrt(n/2).  
c           Where n is the effective site number
         sigma(1) = scalor*sigma(1)
         sigma(2) = scalor*sigma(2)
         sigma(3) = scalor*sigma(3)
         stnm1 = st1 // tail
         stnm2 = st2 // tail
         iyr = int(time1)
         doy = (time1-dble(iyr))/0.3652422d0
         time1 = dble(iyr)+doy
c        o-file output N,E,U. FONDA uses E,N,U.
         if (mode.eq.27) then
            temp1 = obs(1)
            temp2 = sigma(1)
            obs(1) = obs(2)
            obs(2) = temp1
            sigma(1) = sigma(2)
            sigma(2) = temp2
            temp1 = rho(2)
            rho(2) = rho(3)
            rho(3) = temp1
         endif
         write (ifil2,150) time1,obs(1),sigma(1),obs(2),sigma(2),
     .      obs(3),sigma(3),rho(1),rho(2),rho(3),stnm1,stnm2
c---------------------------------------------------
c        temporary change 
c---------------------------------------------------
c         write (ifil2,'(1x,f10.4,f14.5,f14.8,2(2x,a8))')
c     .      time1,obs(4),sigma(4)*1.0d3,stnm1,stnm2
 60   continue
      goto 120
c
 102  print*,' can not open the file : ',site_list
      stop ' stop at READ_OFILE_DAT ...'
 
 30   format (i5,25x,a30)
 140  format (3i5,15x,a30)
 150  format (f9.4,3(1x,f13.4,1x,f7.4),3f9.5,2(1x,a8))
c
 120  continue
      return
      end
