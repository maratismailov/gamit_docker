      subroutine read_usgs_dat(ifil1,ifil2,mode,type)
c
c     read USGS data file and create FONDA format data file
c     error model:  sigma**2 = a**2 + (b*L)**2
c     a = 3mm;   b = 0.2ppm
c     (Savage & Prescott 1973, JGR, 78, pp.6001-6008)
c     mode = 1: distinguesh network
c     mode = 2: merge same site name from different network
c     type = 1: mark to mark distance
c     type = 2: Clarke-arc distance
c
      implicit real*8(a-h,o-z)
      include 'maked.fti'

      character*5 netid(maxsit),netido,netid1,upperc,net_abt(50)
      character*8 name1,name2,sitnam(maxsit),name_ns1,name_ns2
      character*8 abort_name(50)
      integer  id_list(50)
      character*30 note
      character*128 line(50)
      integer iyear,iday,month,id1,id2,mode,type,k
      integer ifil1,ifil2,nsub,isub,jobs,itp,igp
      integer*4 julday
      logical old1,old2,dmg,hold
      character*5 sbntm(20)
      character*2 hp
      integer i,j,remedy_space,idup,unique_name,iho
      common/subnet/nsub,sbntm
c
      small = 1.0d-5
c
c     default format
      if (infmt(1:1).eq.'*')
     . infmt(1:54) = 
     . '(a5,1x,a8,1x,a8,1x,3i2,6x,f6.2,2(f11.4,2x),f6.4,2x,a2)'
c     
c     First, sort out all working sites and experiment number.
      nsit = 0
      iexp = 0
      netido = 'xxxxx'
      print*,' nsub,sbntm:',nsub,(' ',sbntm(i),i=1,nsub)
      do 20 i = 1,30000
         read (ifil1,fmt=infmt,end=50) 
     .   netid1,name1,name2,iyear,month,iday,
     .   azmth,obs1,obs2,flag,hp
         obs = obs2
         if (type.eq.1) obs = obs1
         if (dabs(obs).lt.small.and.obs2.gt.small) obs = obs2
         if (dabs(obs).lt.small) goto 20
c        only record selected subnetwork sites
         do 80 isub = 1,nsub
            if (upperc(netid1).eq.upperc(sbntm(isub))) goto 88
 80      continue
         goto 20
c        replace space inside name1 and name2 with character '_'
 88      j = remedy_space(name1,name_ns1,8,'_',1)
         j = remedy_space(name2,name_ns2,8,'_',1)
         if (nsit.eq.0.or.netid1.ne.netido) then
            nsit = nsit+1
            sitnam(nsit) = name_ns1
            netid(nsit) = netid1
            nsit = nsit+1
            sitnam(nsit) = name_ns2
            netid(nsit) = netid1
            netido = netid1
            iexp = iexp+1
            goto 20
         endif 
c        check if it is a new site
         old1 = .false.
         old2 = .false.
         do 40 j = 1,nsit
            if (mode.eq.1) then
               if (name_ns1.eq.sitnam(j).and.netid1.eq.netid(j))
     .            old1 = .true.
               if (name_ns2.eq.sitnam(j).and.netid1.eq.netid(j))
     .            old2 = .true.
            endif
            if (mode.eq.2) then
               if (name_ns1.eq.sitnam(j))
     .            old1 = .true.
               if (name_ns2.eq.sitnam(j))
     .            old2 = .true.
            endif
            if (old1.and.old2) goto 60
 40      continue
         if (.not.old1) then
            nsit = nsit+1
            sitnam(nsit) = name_ns1
            netid(nsit) = netid1
         endif
         if (.not.old2) then
            nsit = nsit+1
            sitnam(nsit) = name_ns2
            netid(nsit) = netid1
         endif
c        
 60      iexp = iexp+1

 20   continue
 50   print*,' nsit,iexp=',nsit,iexp
      if (nsit.lt.1) goto 120
      jobs = iexp
c
c     check and change duplicated site name
      idup = unique_name(nsit,netid,sitnam,abort_name,
     .          net_abt,id_list)
      if (idup.gt.0) print*,' Total ',idup,' site names'
     .   ,' changed.',
     .   '  Please modify net and map files manually.'
c
c     site id
      note = '{network site number          '
      write (ifil2,30) nsit,note
      write (ifil2,'(a8,1x,a5)') (sitnam(i),netid(i),i=1,nsit)
c
c     experiment number (always 1)
      i = 1
      itp = 4
      note = '{experiment number            '
      write (ifil2,30) i,note
      note = '{exp. index, obs. type, number'
      write (ifil2,140) i,itp,jobs,note

c     second, shift the format one by one
      rewind (ifil1)
      igp = 1
      iho = 0
      aterm = 3.0d0
      bterm = 0.2d-3
      do 90 i = 1,iexp
 85      read (ifil1,fmt=infmt,end=120) 
     .   netid1,name1,name2,iyear,month,iday,
     .   azmth,obs1,obs2,flag,hp
         obs = obs1
         if (type.eq.2) obs = obs2
         hold = .false.
         dmg = .false.
         if (dabs(obs).lt.small.and.obs2.gt.small) then
            obs = obs2
            dmg = .true.
         endif
         if (dabs(obs).lt.small) goto 85
c        only record selected subnetwork sites
         do isub = 1,nsub
            if (upperc(netid1).eq.upperc(sbntm(isub))) goto 95
         enddo
         goto 85
 95      btm = bterm
         if (hp.eq.'hp'.or.hp.eq.'HP') btm = 2.0d-3
         err = dsqrt(aterm**2+(btm*obs)**2)
         if (dmg) err = dsqrt(3.6d1+(0.7d-3*obs)**2)
c         if (flag.gt.9.0d0) err = err*1.4142d0
         time1 = 1900.0d0+julday(month,iday,iyear,1)/365.2422d0
         j = remedy_space(name1,name_ns1,8,'_',1)
         j = remedy_space(name2,name_ns2,8,'_',1)
c        check if this site belongs to duplicated site name
         if (idup.gt.0) then
            do j = 1,idup
               if (name_ns1.eq.abort_name(j).and.netid1.eq.
     .             net_abt(j)) name_ns1 = sitnam(id_list(j))
               if (name_ns2.eq.abort_name(j).and.netid1.eq.
     .             net_abt(j)) name_ns2 = sitnam(id_list(j))
            enddo
         endif
       
         old1 = .false.
         old2 = .false.
         do 70 j = 1,nsit
            if (mode.le.1) then
               if (name_ns1.eq.sitnam(j).and.netid1.eq.netid(j)) then
                  old1 = .true.
                  id1 = j
               endif
               if (name_ns2.eq.sitnam(j).and.netid1.eq.netid(j)) then
                  old2 = .true.
                  id2 = j
               endif
            endif
            if (mode.eq.2) then
               if (name_ns1.eq.sitnam(j)) then
                  old1 = .true.
                  id1 = j
               endif
               if (name_ns2.eq.sitnam(j)) then
                  old2 = .true.
                  id2 = j
               endif
            endif
            if (old1.and.old2) then
               if (dmg) then
                  iho = iho+1
                  write (line(iho),150) 
     .               time1,obs,err,name_ns1,name_ns2,azmth
                  goto 90
               endif
               if (.not.dmg.and.iho.gt.0) then
                  dif = obs1-obs2
                  do k = 1,iho
                     read (line(k),'(11x,f14.5)') obs2
                     write (line(k)(12:25),'(f14.5)') obs2+dif
                     write (ifil2,'(a)') line(k)(1:68)
                     call blank(line(k))
                  enddo
                  iho = 0
               endif
               write (ifil2,150) time1,obs,err,name_ns1,name_ns2,azmth
               goto 90
            endif
 70      continue
         print*,'  missing case:',i,name1,name2

 90   continue

c
c10   format (a4,2x,a8,1x,a8,1x,3i2,12x,f11.4)
 30   format (i5,25x,a30)
 140  format (3i5,15x,a30)
 150  format (1x,f10.4,f14.5,f14.8,2(2x,a8),4x,f5.1)
c
 120  continue
      return
      end
