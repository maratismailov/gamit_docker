      subroutine read_usgs_new_net(ifil1,ifil2,ifil3,nnet,netfmt)
c
c     read USGS new data file and create FONDA format site-table file
c     ( SV6 format).
c
c     unit:
c         x, y, z : km
c        vx,vy,vz : mm/year
c        time     : year
c
      implicit real*8(a-h,o-z)

      character*4 netid
      character*50 netfmt
      character*8 name
      character*1 id1,id2
c     
      write (ifil3,'(a)') ' Network distribution'
      ve = 0.0d0
      vn = 0.0d0
      vu = 0.0d0
c
c     default format
      if (netfmt(1:1).eq.'*')
     . netfmt(1:33) = '(a4,2x,a8,1x,a8,1x,3i2,12x,f11.4)'
c
c     write reference epoch
      write (ifil2,'(3x,a16,2x,f9.4)') 'reference epoch:',rtime
c
      read (ifil1,50) nobs,nnet,name,name2,id1,id2
      read (ifil1,60) sla,i,j,ld,slo
 50   format (11x,i3,11x,i3,10x,a8,13x,a8,4x,a1,6x,a1)
 60   format (f8.1,3i2,f8.1)
      do 20 i = 1,nnet
c        shift format line by line
         read (ifil1,fmt=netfmt) 
     .      netid,name,id1,ld,lm,sla,md,mm,slo,id2,ht
         s1 = dabs(dble(ld))+dble(lm)/60.0d0+sla/3600.0d0
c        westward longitude
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
         write (ifil2,30) netid,name,
     .         id1,ld,lm,sla,id2,md,mm,slo,ht,ve,vn,vu
         write (ifil3,40) s2,s1,ve,vn,name
c         if(i.lt.nnet) read (ifil1,'(2x)')
 20   continue
c
c10   format (a4,2x,a8,3x,a1,2i3,f9.5,i4,i3,f9.5,3x,a1,f8.2)
 30   format (a4,1x,a8,5x,a1,i2,1x,i2,1x,f8.5,
     .        1x,a1,i3,1x,i2,1x,f8.5,f13.4,3f8.2)
 40   format (1x,2f11.5,2f9.3,2x,a8)
c
      return
      end
