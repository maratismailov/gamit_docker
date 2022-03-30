      subroutine read_glbk_full_net(ifil1,ifil2,imfil,
     .   frame)
c
c     transfer GLORG coordinate solution to FONDA net file
c     reads the 'Unc.' line - already grepped to a temp file 
c      in the calling routine
c
c     ifil1 :  input network file 
c     ifil2 :  output FONDA network file
c     ifil3 :  site list file
c     imfil :  output mapping file
c
      implicit real*8(a-h,o-z)
      include 'maked.fti'

      character*6  frame
      integer len,k,j,i,lift_arg,idx,imfil,match_name
      integer isit,k1,ifil1,ifil2,ifil3,ilst,nblen
      integer id1,im1,id2,im2
      character*1 sym1,sym2
      character*8  st1,st_list(400),stnm1
      real*8 fct

      print*,' Begining to transfer coordinates ....'
c     
      if (iomode(5).gt.0) 
     .   write (imfil,'(a)') ' Network distribution'

      

      if (iomode(8).eq.1) then
         fct=mmtom
         write (imfil,'(a)') ' Velocities in m/yr '
      else 
         fct=1.0d0
         write (imfil,'(a)') ' Velocities in mm/yr '
      endif

c
c     default format
      if (netfmt(1:1).eq.'*') then
       netfmt(1:31) = '(5x,a8,3f14.4,3f8.4,f9.3,3f8.4)'
c      netfmt(1:32) = '(5x,a8,3f14.4,3f10.4,f9.3,3f8.4)'
          print*,'Using default file format : ',netfmt
      endif
c
c     write FONDA style network file header
      call net_head(ifil2)
c
c     get all explicitely designed site names
      ilst = 0
      if (site_list(1:1).ne.'*') then
         ifil3 = 14
         open (ifil3,file=site_list,status='old',err=1000)
         do 10 i = 1,1000
            read (ifil3,'(a)',end=15) stnm1
            j = lift_arg(stnm1,st_list(i),1)
            ilst = ilst+1
 10      continue
         close (ifil3)
      endif
c     
 15   nsit = 0
      rewind (ifil1)
c     screen display to avoid sleep
      print*,'   Processing  <<<  ',infil
      idx = 0
      stnm1 = 'glbk_sln'
      call geotab(frame,1,semi,finv,tx,ty,tz)
      do 20 isit = 1,100000
         read (ifil1,fmt=netfmt,end=1000) 
     .      st1,bx,by,bz,vxt,vyt,vzt,time1,sx,sy,sz
c        check site name list
         len = nblen(st1)
         k1 = 1
         if (ilst.gt.0) then
            k1 = match_name(ilst,len,st_list,st1)
         endif
         nsit = nsit+1
         if (len.lt.8) then
            do k = len+1,8
               st1(k:k) = '_'
            enddo
         endif
c        from geocentric to geodetic frame
         call geoxyz(semi,finv,tx,ty,tz,alat,alon,ht,
     .      bx,by,bz,2,heit)
         call rtodms(alat,id1,im1,sela,1)
         sym1 = 'N'
         if (alat.lt.0.0d0) then
            sym1 = 'S'
            id1 = iabs(id1)
            im1 = iabs(im1)
            sela = dabs(sela)
         endif
         call rtodms(alon,id2,im2,selo,1)
c        use western longitude (follow the tradition)
         sym2 = 'E'
         if (alon.gt.pi) then
            sym2 = 'W'
            id2 = 359-id2
            im2 = 59-im2
            selo = 60.0d0-selo
            if (selo.ge.60.0d0) then
               selo = selo-60.0d0
               im2 = im2+1
            endif
            if (im2.ge.60) then
               im2 = im2-60
               id2 = id2+1
            endif
         endif

c        output mapping file
         if (iomode(5).gt.0.and.k1.gt.0) then
            if (alon.gt.pi) alon=alon-2*pi
            write (imfil,40) alon*rtod,alat*rtod,ve1*fct,vn1*fct,st1
         end if


         call sph_ca(alat,alon,vxt,vyt,vzt,v1,v2,v3,2)
c
c
         write (ifil2,30) st1,stnm1,
     .      sym1,id1,im1,sela,sym2,id2,im2,selo,ht,v1,v2,v3,
     .      time1,sx,sy,sz
 
 20   continue
c
 30   format (1x,a8,2x,a8,4x,a1,2(i2,1x),f8.5,1x,
c    .        a1,i3,1x,i2,1x,f8.5,f13.4,3f8.4,f9.3,3f8.4)
     .   a1,i3,1x,i2,1x,f8.5,f13.4,2f8.2,f7.3,f9.3,3f8.4)
 40   format (1x,2f11.5,2f9.3,2x,a8)
 
 1000 continue
      return
      end
