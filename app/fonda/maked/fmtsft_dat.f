      subroutine fmtsft_dat()
c
c     coordinate system : (comode)
c         1. geocentric spherical coordinate (lat,lon,radius)
c         2. geodetic (ellipsoidal) coordinate (lat,lon,height)
c         3. geocentric Cartesian coordinate (x,y,z)
c         4. topocentric coordinate (x(east),y(north),z(height))
c
c     input data file style : (in_style)
c         1. GAMIT o-file 
c         2. GLOBK glorg output
c         3. FONDA maked.out file
c         4. BLUE BOOK 
c         5. USGS data
c         6. UCSD data (slightly different from USGS)
c         7. IPGP data
c         8. GAMIT h-file 
c         9. FONDA output h-file
c        10. GLOBK glorg output with full covariance matrix
c        11. NEW BLUE BOOK format
c        12. GIPSY data
c        13. SINEX file
c

      implicit real*8(a-h,o-z)
      include 'maked.fti'

      integer mode,ier,np
      character*256 line
c     character*5 sbntm(10)
c     common/subnet/nsub,sbntm
      
      if (iomode(4).eq.1) open (8,file=infil,status='old',err=101)
      if (iomode(4).eq.2) open (8,file=in_list,status='old',err=101)
      open (9,file=outfil,status='unknown')

      if (in_style.eq.1) then
         call grep_data(8,in_opt,mode)
         close (8)
         open (8,file='tmpofile',status='old',err=103)
         call read_ofile_dat(8,9,mode)
      else if (in_style.eq.2) then
          call read_glbk_dat(8,9,1)
      else if (in_style.eq.4) then
          call read_bbook_dat(8,9,1,1)
      else if (in_style.eq.5) then
          call read_usgs_dat(8,9,1,1)
      else if (in_style.eq.6) then
          call read_ucsd_dat(8,9)
      else if (in_style.eq.7) then
          call read_ipgp_dat(8,9)
      else if (in_style.eq.8) then
         call read_hfile_dat(8,9,mode)
      else if (in_style.eq.9) then
         call read_fondah_dat(8,9,mode)
      else if (in_style.eq.10) then
          call read_glbk_full_dat(8,9,1)
      else if (in_style.eq.11) then
          call read_newbb_dat(8,9,1,1)
      else if (in_style.eq.12) then
          call read_gipsy_dat(8,9,1)
      else if (in_style.eq.13) then
          call read_sinex_dat(8,9,line,np)
      else
          print *,'FMTSFT_DAT: unknown in_style ',in_style
          stop 'FMTSFT_DAT: unknown in_style '
      endif

      close (8)
      close (9)
      if (in_style.eq.1) ier = system('\rm -f tmpofile')
      goto 1000
 101  print*,' can not open the file 1: ',infil
      stop ' stop at FMTSFT_DAT ...'
 103  print*,' can not open the file 2: ','tmpofile'
      stop ' stop at FMTSFT_DAT ...'

 1000 continue
      return
      end
