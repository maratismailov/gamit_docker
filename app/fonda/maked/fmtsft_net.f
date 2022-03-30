      subroutine fmtsft_net(frame)
c
c     coordinate system : (comode)
c         1. geocentric spherical coordinate (lat,lon,radius)
c         2. geodetic (ellipsoidal) coordinate (lat,lon,height)
c         3. geocentric Cartesian coordinate (x,y,z)
c         4. topocentric coordinate (x(east),y(north),z(height))
c
c     network file style : (net_style)
c         1. GAMIT l-file 
c         2. GLOBK 
c         3. FONDA
c         4. BLUE BOOK
c         5. USGS statab.lis
c         6  UCSD
c         7. IPGP 
c         8. GLORG full
c
      implicit real*8(a-h,o-z)
      include 'maked.fti'

      character*128 line
      integer l1,l2,nblen
c     character*5 sbntm(10)
      character*6 frame
c     common/subnet/nsub,sbntm

      print*,' Begin to transfer net-file ...'

      open (9,file=newnet,status='unknown')
      if (iomode(5).gt.0) open (10,file=mapfil,status='unknown')
      open (8,file=netfil,status='old',err=1000)

      if (net_style.eq.1) call read_lfile(8,9,10,nnet,netfmt,rtime)
      if (net_style.eq.2) call read_glbk_net(8,9,10,nnet,netfmt,rtime)
      if (net_style.eq.4) 
     .   call read_bbook_net(8,9,10,nnet,netfmt,iesit,sname,rtime)
      if (net_style.eq.5) 
     .   call read_usgs_net(8,9,10,nnet,netfmt,net_opt,rtime)
      if (net_style.eq.6) call read_ucsd_net(8,9,10,nnet,netfmt,rtime)
      if (net_style.eq.7) call read_ipgp_net(8,9,10)
      if (net_style.eq.8) then
         close (8)
         l1 = nblen(netfil)
         call blank(line)
         line = 'grep -e "Unc." ' // netfil(1:l1) // 
     .    ' | grep -v "[\*][NaN]" >! tmp_file'
         l2 = system(line)
c        print lines containing * to the screen   
         print*,'The following sites have problems ... '
         line = 'grep -e "Unc." ' // netfil(1:l1) //
     .    ' | grep "*" '
         
         l2 = system(line)
         open (8,file='tmp_file',status='old',err=1000)
         call read_glbk_full_net(8,9,10,frame)
      endif
      if (net_style.eq.9) 
     .   call read_newbb_net(8,9,10,nnet,netfmt,iesit,sname,rtime)
      if (net_style.eq.10) 
     .   call read_gipsy_net(8,9,10,nnet,netfmt,iesit,sname,rtime)
      if (net_style.eq.12) call read_free_net(8,9,frame)
      close (8)
      close (9)
      if (iomode(5).gt.0) close (10)
      goto 100

 1000 continue
      print*,'error in openning ',netfil
      stop

 100  continue
      return
      end
