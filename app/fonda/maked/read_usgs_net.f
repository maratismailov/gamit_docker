      subroutine read_usgs_net(ifil1,ifil2,ifil3,nnet,netfmt,
     .   net_opt,rtime)
c
c     read USGS site-table file and create FONDA format site-table file
c     ( SV6 format).
c
c     unit:
c         x, y, z : km
c        vx,vy,vz : mm/year
c        time     : year
c
      implicit real*8(a-h,o-z)

      character*5 netid,subnam
      character*8 name, name_ns
      character*80 netfmt,net_opt
      character*12 net_temp, net_sub
      character*1 id1,id2
      character*5 sbntm(10)
      integer i,isub,j,remedy_space
      integer ifil1,ifil2,ifil3,nnet,nsub,ld,lm,mm,md
      common/subnet/nsub,sbntm
c     
      write (ifil3,'(a)') ' Network distribution'
      ve = 0.0d0
      vn = 0.0d0
      vu = 0.0d0
      sigx = 0.0d0
      sigy = 0.0d0
      sigz = 0.0d0
c
c     decode to get subnetwork selection (moved to read_driv)
c      read (net_opt,*) nsub,(sbntm(i),i=1,nsub)
c
c     usgs default format
      if (netfmt(1:1).eq.'*') 
     .   netfmt(1:41) = '(a5,1x,a8,4x,2i3,f9.5,i4,i3,f9.5,4x,f8.2)'
c
c     write FONDA style network file header
      call net_head(ifil2)
c
c     subtract sites according to selected subnetworks
      do 80 isub = 1,nsub
         subnam = sbntm(isub)
         rewind (ifil1)
         do 20 i = 1,nnet
c           shift format line by line
            read (ifil1,fmt=netfmt) 
     .      netid,name,ld,lm,sla,md,mm,slo,ht
c     later add case insensive function
            if (netid.ne.subnam) goto 20
            s1 = dabs(dble(ld))+dble(lm)/60.0d0+sla/3600.0d0
c           westward longitude
            s2 = (dabs(dble(md))+dble(mm)/60.0d0+slo/3600.0d0)
            id1 = 'N'
            if (ld.lt.0) then
               id1 = 'S'
               ld = -ld
               s1 = -s1
            endif
            id2 = 'E'
            if (md.lt.0) then
               id2 = 'W'
               md = -md
               s2 = -s2
            endif
c           replace space inside name with character '_'
            j = remedy_space(name,name_ns,8,'_',1)
c           expand usgs 5-character subnet id code to 8-character
            net_temp(1:5) = netid
            net_temp(6:9) = '_net'
            j = remedy_space(net_temp,net_sub,9,'_',2)
            write (ifil2,30) name_ns,net_sub(1:9),
     .         id1,ld,lm,sla,id2,md,mm,slo,ht,ve,vn,vu,
     .         rtime,sigx,sigy,sigz
            write (ifil3,40) s2,s1,ve,vn,name_ns
c            if(i.lt.nnet) read (ifil1,'(2x)')
 20      continue
 80   continue
c
 30   format (1x,a8,1x,a9,4x,a1,2(i2,1x),f8.5,1x,
c    .        a1,i3,1x,i2,1x,f8.5,f13.4,3f8.4,f9.3,3f8.4)
     .   a1,i3,1x,i2,1x,f8.5,f13.4,2f8.2,f7.3,f9.3,3f8.4)
 40   format (1x,2f11.5,2f9.3,2x,a8)
c
      return
      end
