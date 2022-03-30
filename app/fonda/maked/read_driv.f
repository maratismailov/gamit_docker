      subroutine read_driv(frame)
c
c     all driving informations
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
c         8. GLORG 
c         9. new BLUE BOOK
c        10. GIPSY 
c        12. GLOBK free format (xyz)
c        13. (waiting list)
c
c     input data file style : (in_style)
c         1. GAMIT o-file 
c         2. GLOBK glorg output
c         3. FONDA maked.out file
c         4. BLUE BOOK 
c         5. USGS data
c         6. UCSD data (little different from USGS)
c         7. IPGP data
c         8. GAMIT h-file 
c         9. FONDA output h-file
c        10. GLOBK glorg output with full covariance matrix
c        11. new BLUE BOOK
c        12. GIPSY sta.cov file
c        13. SINEX
c
      implicit real*8(a-h,o-z)
      include 'maked.fti'
c
      character*6 frame
      character*80 wcmd
      character*30 task(5)
      integer lcmd, nsub,i,iyear,iday,month
      character*5 sbntm(10),vunit
      integer*4 julday
      common/subnet/nsub,sbntm
c
      data task/' simulate data                ',
     .          ' transfer data to FONDA format',
     .          ' extract specified data only  ',
     .          ' compress repeated observation',
     .          ' merge several maked.out files'/

c     task identification:
      call getcmd(8,'task',wcmd,lcmd,1)
      if (lcmd.gt.0) then
         if (wcmd(1:1).eq.'1'.or.wcmd(1:1).eq.'s') iomode(1) = 1
         if (wcmd(1:1).eq.'2'.or.wcmd(1:1).eq.'t') iomode(1) = 2
         if (wcmd(1:1).eq.'3'.or.wcmd(1:1).eq.'e') iomode(1) = 3
         if (wcmd(1:1).eq.'4'.or.wcmd(1:1).eq.'c') iomode(1) = 4
         if (wcmd(1:1).eq.'5'.or.wcmd(1:1).eq.'m') iomode(1) = 5
         print*, ' task of this running: ',task(iomode(1))
      endif
      if (iomode(1).eq.3) goto 10
      if (iomode(1).le.0.or.iomode(1).gt.5) then
         wcmd = ' No task has been identified. What do you want?'
         goto 1000
      endif

c     network coordinate file name and its format
      call getcmd(8,'netfil',wcmd,lcmd,1)
      if (lcmd.gt.0) then
         netfil = wcmd(1:lcmd)
         iomode(2) = 1
         print*, ' network coordinate file name: ',netfil
         call getcmd(8,'netfmt',wcmd,lcmd,1)
         if (lcmd.gt.0) netfmt = wcmd(1:lcmd)
         call getcmd(8,'net_style',wcmd,lcmd,1)
         if (lcmd.gt.0) then
            if (wcmd(1:5).eq.'lfile'.or.wcmd(1:5).eq.'LFILE') 
     .         net_style = 1
            if (wcmd(1:5).eq.'globk'.or.wcmd(1:5).eq.'GLOBK') 
     .         net_style = 2
            if (wcmd(1:5).eq.'fonda'.or.wcmd(1:5).eq.'FONDA') 
     .         net_style = 3
            if (wcmd(1:5).eq.'bbook'.or.wcmd(1:5).eq.'BBOOK') 
     .         net_style = 4
            if (wcmd(1:4).eq.'usgs'.or.wcmd(1:4).eq.'USGS') 
     .         net_style = 5
            if (wcmd(1:4).eq.'ucsd'.or.wcmd(1:4).eq.'UCSD') 
     .         net_style = 6
            if (wcmd(1:4).eq.'ipgp'.or.wcmd(1:4).eq.'IPGP') 
     .         net_style = 7
            if (wcmd(1:5).eq.'glorg'.or.wcmd(1:5).eq.'GLORG') 
     .         net_style = 8
            if (wcmd(1:5).eq.'newbb'.or.wcmd(1:5).eq.'NEWBB') 
     .         net_style = 9
            if (wcmd(1:5).eq.'gipsy'.or.wcmd(1:5).eq.'GIPSY') 
     .         net_style = 10
            if (wcmd(1:4).eq.'free'.or.wcmd(1:4).eq.'FREE') 
     .         net_style = 12
         endif
         print*, ' network file style: ',wcmd(1:lcmd)
      else
c        if no netfil, then check netput file list
         call getcmd(8,'net_list',wcmd,lcmd,1)
         if (lcmd.gt.0) then
            net_list = wcmd(1:lcmd)
            print*, ' network file list: ',net_list
            iomode(2) = 1
         else
            iomode(2) = 0
         endif
      endif
c
c     subnetwork option and output net file
      call getcmd(8,'net_opt',wcmd,lcmd,1)
      if (lcmd.gt.0) then
         net_opt = wcmd(1:lcmd)
         read (net_opt,*) nsub,(sbntm(i),i=1,nsub)
      endif
c
      if (iomode(2).gt.0) then
         call getcmd(8,'comode',wcmd,lcmd,1)
         if (lcmd.gt.0) then
            if (wcmd(1:5).eq.'spher'.or.wcmd(1:5).eq.'SPHER') 
     .         comode = 1
            if (wcmd(1:5).eq.'geode'.or.wcmd(1:5).eq.'GEODE') 
     .         comode = 2
            if (wcmd(1:5).eq.'carte'.or.wcmd(1:5).eq.'Carte') 
     .         comode = 3
            if (wcmd(1:5).eq.'topoc'.or.wcmd(1:5).eq.'local') 
     .         comode = 4
         endif
c        output new network file
         call getcmd(8,'newnet',wcmd,lcmd,1)
         if (lcmd.gt.0) then
            newnet = wcmd(1:lcmd)
            print*, ' new network file name: ',newnet
            iomode(3) = 1
         else
            iomode(3) = 0
         endif
      endif

c     reference network file for alignment
      call getcmd(8,'ref_net',wcmd,lcmd,1)
      if (lcmd.gt.0) then
         ref_net = wcmd(1:lcmd)
         iomode(6) = 1
         print*, ' reference network file name: ',ref_net
      else
         iomode(6) = 0
      endif

c     rename file for correct site name and adjustment
      call getcmd(8,'rename',wcmd,lcmd,1)
      if (lcmd.gt.0) then
         renam_fil = wcmd(1:lcmd)
         iomode(7) = 1
         print*, ' rename file name: ',renam_fil
      else
         iomode(7) = 0
      endif
c
c     input data file
 10   call getcmd(8,'infil',wcmd,lcmd,1)
      if (lcmd.gt.0) then
         infil = wcmd(1:lcmd)
         print*, ' input data file name: ',infil
         iomode(4) = 1
      else
         call getcmd(8,'in_list',wcmd,lcmd,1)
         if (lcmd.gt.0) then
            in_list = wcmd(1:lcmd)
            print*, ' input file list: ',in_list
            call getcmd(8,'in_opt',wcmd,lcmd,1)
            if (lcmd.gt.0) in_opt = wcmd(1:lcmd)
            iomode(4) = 2
         else
            if (iomode(2).gt.0) goto 30
            wcmd = ' No input file. No net file. What do you want?'
            goto 1000
         endif
      endif
c
      call getcmd(8,'infmt',wcmd,lcmd,1)
      if (lcmd.gt.0) infmt = wcmd(1:lcmd)
      call getcmd(8,'in_style',wcmd,lcmd,1)
      if (lcmd.gt.0) then
         if (wcmd(1:5).eq.'ofile'.or.wcmd(1:5).eq.'OFILE') 
     .      in_style = 1
         if (wcmd(1:5).eq.'globk'.or.wcmd(1:5).eq.'GLOBK') 
     .      in_style = 2
         if (wcmd(1:5).eq.'fonda'.or.wcmd(1:5).eq.'FONDA') 
     .      in_style = 3
         if (wcmd(1:5).eq.'bbook'.or.wcmd(1:5).eq.'BBOOK') 
     .      in_style = 4
         if (wcmd(1:4).eq.'usgs'.or.wcmd(1:4).eq.'USGS') 
     .      in_style = 5
         if (wcmd(1:4).eq.'ucsd'.or.wcmd(1:4).eq.'UCSD') 
     .      in_style = 6
         if (wcmd(1:4).eq.'ipgp'.or.wcmd(1:4).eq.'IPGP') 
     .      in_style = 7
         if (wcmd(1:5).eq.'hfile'.or.wcmd(1:5).eq.'HFILE') 
     .      in_style = 8
         if (wcmd(1:6).eq.'fondah'.or.wcmd(1:6).eq.'FONDAH') 
     .      in_style = 9
         if (wcmd(1:5).eq.'gkful'.or.wcmd(1:5).eq.'GKFUL') 
     .      in_style = 10
         if (wcmd(1:5).eq.'newbb'.or.wcmd(1:5).eq.'NEWBB') 
     .      in_style = 11
         if (wcmd(1:5).eq.'gipsy'.or.wcmd(1:5).eq.'GIPSY') 
     .      in_style = 12
         if (wcmd(1:5).eq.'sinex'.or.wcmd(1:5).eq.'SINEX')
     .      in_style = 13
         if (wcmd(1:3).eq.'glf'.or.wcmd(1:3).eq.'GLF') 
     .      in_style = 15
         if (wcmd(1:3).eq.'v_s'.or.wcmd(1:3).eq.'V_S') 
     .      in_style = 16
      endif
      print*, ' input data file style: ',wcmd(1:lcmd)
c
c     Although not all in_style files use site_list yet,
c     it is harmless to read this command
c     Currently ofile, gkful, gipsy & sinex
c      if (in_style.eq.1.or.in_style.eq.10.or.in_style.eq.12
c    .     .or.in_style.eq.13) then
 30      call getcmd(8,'site_list',wcmd,lcmd,1)
         if (lcmd.gt.0) then
            site_list = wcmd(1:lcmd)
            print*, ' site selection file: ',site_list
         else
            site_list = '*'
         endif
c      endif
c
c     xyu and/or uvw
      if (in_style.eq.15) then
         call getcmd(8,'in_opt',wcmd,lcmd,1)
         if (lcmd.gt.0) in_opt = wcmd(1:lcmd)
      endif
c
c     output file name
      call getcmd(8,'outfil',wcmd,lcmd,1)
      if (lcmd.gt.0) then
         outfil = wcmd(1:lcmd)
         print*, ' output file name: ',outfil
      else
         wcmd = ' No output file.'
         goto 1000
      endif

      if (iomode(1).eq.3) goto 200
c
c     mapping file name
      call getcmd(8,'mapping',wcmd,lcmd,1)
      if (lcmd.gt.0) then
         mapfil = wcmd(1:lcmd)
         print*, ' mapping file name: ',mapfil
         iomode(5) = 1
      else
         iomode(5) = 0
      endif
c
c     reference frame
      call getcmd(8,'frame',wcmd,lcmd,1)
      if (lcmd.gt.0) frame = wcmd(1:6)
      print*,' reference frame: ',frame
c
c     reference time and network site number
      call getcmd(8,'rtime',wcmd,lcmd,1)
      if (lcmd.gt.0) then
         read(wcmd,*) iyear,iday,nnet
         rtime = 1900.0d0+julday(month,iday,iyear,2)/365.2422d0
      endif
c
c     dimension of deformation field
      call getcmd(8,'dim',wcmd,lcmd,1)
      if (lcmd.gt.0) then
         read(wcmd,*) dim
      endif

c     
c     velocity units m/yr or mm/yr (default)
      iomode(8)=0 
      call getcmd(8,'Velo',wcmd,lcmd,1)
      if (lcmd.gt.0) then
         read(wcmd,'(a)') vunit
         if (vunit.eq.'m/yr') then
             iomode(8)=1
         else
             iomode(8)=0
         endif
      endif 
c 
      goto 200

c     troubled stop
 1000 print*,'suicide due to .. ',wcmd
      stop ' stop at READ_DRIV.'
c
 200  continue
      return
      end
