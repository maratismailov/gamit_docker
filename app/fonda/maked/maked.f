      program maked
c
c     this program can perform 4 tasks (for the time being):
c         1 = simulation
c         2 = transfer data to FONDA format
c         3 = extract specified data only
c         4 = compress data as obs + obs rate
c         5 = merge multiple data files into single data file
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
c         6. UCSD 
c         7. IPGP
c         8. (waiting list)
c
c     input data file style : (in_style)
c         1. GAMIT o-file 
c         2. GLOBK glorg output
c         3. FONDA maked.out file
c         4. BLUE BOOK 
c         5. USGS data
c         6. UCSD data (little different from USGS)
c         7. IPGP data
c         8. (waiting list)
c
      implicit real*8(a-h,o-z)
      include 'maked.fti'
c
      character*6 frame
      character*16 version
      character*32 drvfil
      integer system,iclarg,i1,i2,get_path,ierr,mode,noisl,nobs
c 
c     If the driving file name is missing, the program will give
c     the description only.
      i1 = iclarg(1,drvfil)
      if (i1.le.0) then
         i2 = get_path('FONDA_H',netfmt)
         infmt(1:24+i2) = 
     .   'head -13 ' // netfmt(1:i2) // 'help/maked.help'
         ierr = system(infmt)
         stop
      endif
c
c     get version ID
      call museum(version)

      print*, ' --------------------------------------------------'
      print *,'              MAKED version: ',version
      print*, ' --------------------------------------------------'

c     read driving file
      open (8,file=drvfil,status='old',err=1000)
      call read_driv(frame)
      close (8)
      print*,' comode: ',comode
c
c     get reference site coordinates
      if (iomode(6).gt.0) then
         open (17,file=ref_net,status='old',err=1000)
         call read_ref_net(17)
         close (17)
      endif
c
c     transfer other formats to standard FONDA format
      if (iomode(1).eq.2) then
         if (iomode(2).gt.0) call fmtsft_net(frame)
         if (iomode(4).gt.0) call fmtsft_dat()
         goto 100
      endif
      if (iomode(4).le.0) goto 100
c
c     merge several maked.out files to create a list file
      if (iomode(1).eq.5) then
c         require in_list (iomode(4).eq.2) for this option so
         if (iomode(4).le.1) goto 100
         open (8,file=in_list,status='old',err=1000)
         open (9,file=outfil,status='unknown')
         call merge_list(8,9,mode)
         goto 100
      endif
c
c     get geodetic parameters
      call geotab(frame,1,radius,finv,tx,ty,tz)
c
c     get network site coordinates
      if (iomode(2).gt.0) then
         open (8,file=netfil,status='old',err=1000)
         call getnet(8)
         close (8)
      endif
c
c     compress repeated observations
      if (iomode(1).eq.4) then
         call compress_dat()
         goto 100
      endif
c
c     read input informations
      open (10,file=infil,status='old',err=1000)
      call redinf(10,noisl)
      close (10)
c
c     construct the pure simulated observables
      call frmobs(nobs)
c
c     add observation noise (option)
c      if (noisl.gt.0) then
c      call addnoi(sdata,nmode,...)
c      endif
c
c     output results by standard format
      open (10,file=outfil,status='unknown',err=1000)
      if (iomode(5).gt.0)
     .   open (11,file=mapfil,status='unknown',err=1000)
      call outdat(10,11)
      close (10)
      if (iomode(5).gt.0) close (11)
      goto 100
c
c     troubled stop
 1000 print*,'suicide due to ..'
      stop
c
c     normal stop
 100  print*,' Normal stop. Congratulations!'
      print*, ' --------------------------------------------------'
      stop
      end
