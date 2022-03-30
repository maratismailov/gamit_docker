      subroutine compress_dat()
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
c         9. (waiting list)
c

      implicit real*8(a-h,o-z)
      include 'maked.fti'

      if (iomode(4).eq.1) open (8,file=infil,status='old',err=101)
      if (iomode(4).eq.2) open (8,file=in_list,status='old',err=101)
      open (9,file=outfil,status='unknown')
      if (iomode(5).gt.0) open (10,file=mapfil,status='unknown')

      if (in_style.eq.1) then
          continue
      else if (in_style.eq.2) then
          continue
      else if (in_style.eq.3) then
          call fonda_out_compress(8,9,10)
      else if (in_style.eq.5) then
          continue
      else if (in_style.eq.6) then
          continue
      else if (in_style.eq.7) then
          continue
      else if (in_style.eq.8) then
          continue
      else
          print *,'FMTSFT_DAT: unknown in_style ',in_style
          stop 'FMTSFT_DAT: unknown in_style '
      endif

      close (8)
      if (iomode(5).gt.0) close (10)

      goto 1000
 101  print*,' can not open the file 1: ',infil
      stop ' stop at COPMRESS_DAT ...'
 103  print*,' can not open the file 2: ','tmpofile'
      stop ' stop at COPMRESS_DAT ...'

 1000 continue
      return
      end
