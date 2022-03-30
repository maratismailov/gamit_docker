      program make_res
 
      implicit none 

*     Program to make the weekly sinex combinated residual file
*
*         max_site      - Maximum number of sites
*           max_center  - Maximum number of centers
 
      integer*4 max_site, max_center
 
      parameter ( max_site = 500 )
      parameter ( max_center = 10 )
 
*   len_run     - Length of string
*   ierr, jerr      - IOSTAT errors
*   i,j,k,l     - Loop indices
*   indx            - position in string
*   num_site        - Number of sites
*   num_center      - Number of centers
 
*   list(max_site)  - List of sites in ascending order after sort
*   iel, jel        - Site and center index numbers
 
      integer*4 len_run, ierr, jerr, i,j,k,l, indx, num_site,
     .    num_center, list(max_site), iel, jel, rcpar, trimlen
 
*   dpos(3)     - NEU adjuists
*   spos(3)     - Sigmas for adjustments
 
*   pdiff(3,max_site, max_center)   - position adjustments
*   sdiff(3,max_site, max_center)   - Sigma for adjustments
 
*   site_stats(4,3,max_site)        - Number, Mean, RMS, Chi by
*                   - coordinate and site
*   center_stats(4,3,max_center)    - Number, Mean, RMS, Chi by
*                   - coordinate and center
 
      real*8 dpos(3), spos(3), pdiff(3,max_site, max_center),
     .    sdiff(3,max_site, max_center), site_stats(4,3,max_site),
     .    center_stats(4,3,max_center)
 
*   infile          - Input file name
*   outfile     - Output file name
*   headfile        - File with description information
*   line            - Line read/written from/to files
 
      character*256 infile, outfile, line, headfile

*   date_start      - Character String with start date
*   date_end        - Character string wiht end date

      character*8 date_start, date_end

*   GPS_week        - Week for which this file refers

      integer*4 GPS_week
 
*   site_names(max_site)    - Site names
*   center_names(max_center)    - Center names
 
      character*4 site_names(max_site), center_names(max_center),
     .            name, center
 
*   type(3)     - Position type (NEU)
 
      character*1 type(3)
 
      data type  / 'N', 'E', 'U' /
 
***** Set up the center names
      num_center = 7
      center_names(1) = 'COD ' 
      center_names(2) = 'JPL ' 
      center_names(3) = 'SIO ' 
      center_names(4) = 'EMR ' 
      center_names(5) = 'ESA ' 
      center_names(6) = 'GFZ ' 
      center_names(7) = 'NGS ' 
 
****  Get the runstring
      len_run = rcpar(1,infile)
      if( len_run.eq.0 ) then
          call proper_runstring('make_res.hlp', 'make_res', 1)
      end if
 
      open(100,file=infile, status='old', iostat=ierr)
      call report_error('IOSTAT',ierr,'open',infile,1,'make_res')
 
*     Get the output file name
      len_run = rcpar(2,outfile)
      if( len_run.eq.0 ) then
          call proper_runstring('make_res.hlp', 'make_res', 1)
      end if
 
      open(200,file=outfile, status='unknown', iostat=ierr)
      call report_error('IOSTAT',ierr,'creat',outfile,1,'make_res')

*     Get GPS week
      len_run = rcpar(3,line)
      if( len_run.gt.0 ) then
          read(line,*) GPS_week
      else
          GPS_week = 900
      end if

*     Get the start and end dates
      call gpsw_date(gps_week, 0, date_start)
      call gpsw_date(gps_week, 6, date_end)

*     Write header to output file
      write(200,100,iostat=ierr) gps_week, date_start, date_end
100   format(/,80('*'),/,
     .       'MIT RESIDUAL REPORT GPS-WEEK ',i4.4,20x,
     .       'Dates: ',a8,' to ',a8,/,80('*'))

***** Now see if we can read the header file
      len_run = rcpar(4,headfile)
      if( len_run.gt.0 ) then
          open(101,file=headfile, iostat=jerr, status='old')
          if( jerr.eq.0 ) then
              do while (jerr.eq.0 )
                 read(101,'(a)', iostat=jerr) line
                 if( jerr.eq.0 ) then
                     write(200,'(a)') line(1:max(1,trimlen(line)))
                 end if
              end do
          end if
      end if
 
 
****  Initialize and start reading
      num_site = 0
 
      do i = 1,max_center
          do j = 1, max_site
              do k = 1,3
                  pdiff(k,j,i) = 999.9
                  sdiff(k,j,i) = 999.9
              end do
          end do
      end do
 
      do while ( ierr.eq.0 )
          read(100,'(a)',iostat=ierr) line
          if( ierr.eq.0 .and. index(line,'wrt COMB').gt.0 ) then
 
*             Found a valid entry; read the name and the site
*             adjustments
              read(line,120,iostat=jerr) name, (dpos(i),spos(i),i=1,3),
     .                center
 120          format(1x,a4,2x,6F7.1,1x,a4)
 
*             See if this is a site we have
              indx = 1
              call get_cmd(name, site_names, num_site, iel, indx)
              if( iel.le.0 ) then
                  num_site = num_site + 1
                  iel = num_site
                  site_names(iel) = name
              end if
 
*             Check for center number
              indx = 1
              call get_cmd(center, center_names, num_center,
     .                jel, indx)
              if( jel.le.0 ) then
                  num_center = num_center + 1
                  jel = num_center
                  center_names(jel) = center
              end if
 
****          Now add entry to list
              do i = 1,3
                  pdiff(i,iel, jel) = dpos(i)
                  sdiff(i,iel, jel) = spos(i)
              end do
          end if
      end do
 
****  Now arrange the values and compute the statistics
      do i = 1, num_site
          list(i) = i
      end do
      call nsort(num_site, site_names, list )
 
***** Now do the stats on each site and for each center
      do i = 1, num_site
          do j = 1,3
              call get_stats( site_stats(1,j,i), pdiff, sdiff,
     .                max_site, max_center, num_site, num_center,
     .                i, j, 'SITE')
          end do
      end do
 
      do i = 1, num_center
          do j = 1,3
              call get_stats( center_stats(1,j,i), pdiff, sdiff,
     .                max_site, max_center, num_site, num_center,
     .                i, j, 'CENTER')
          end do
      end do
 
*     Now write our the results
      write(200,300)
 300  format('SITE AND CENTER RESIDUALS TO COMBINED SOLUTION')
      write(200,320) (center_names(i), i = 1, num_center)
 320  format(' SITE T #   RMS   Chi ',10(5x,a4,6x))
      write(200,340) (' dPOS    +- ',i = 1, num_center)
      write(200,340) (' (mm)   (mm)',i = 1, num_center)
 340  format('                     ',10(a14,1x))
 
      do i = 1, num_site
          j = list(i)
          do k = 1,3
 
              write(line,400) site_names(i), type(k),
     .            nint(site_stats(1,k,j)), site_stats(3,k,j),
     .            site_stats(4,k,j), (pdiff(k,j,l), sdiff(k,j,l),
     .            l = 1, num_center)
400           format(1x,a4,1x,a1,i2,2F6.1,10(1x,F6.1,1x,F6.1,1x))
              call sub_char(line,'999.9','  -- ')
              write(200,'(a)') line(1:trimlen(line))
          end do
          write(200,'(1x)')
      end do

*     Write center statistics
      write(200,500) 
 500  format(/'CENTER SUMMARY STATISTICS')
      write(200,520) 
 520  format(' CENTER   #       North           East           Up',/,
     .       '               RMS    Chi     RMS    Chi    RMS     Chi',/,
     .       '              (mm)           (mm)          (mm)')
      do i = 1, num_center
          write(200,540) center_names(i), nint(center_stats(1,1,i)),
     .        (center_stats(3,k,i), center_stats(4,k,i), k = 1,3)
 540      format(1x,a4,3x,i3,10(1x,F6.1,1x,F6.1,1x))
      end do
 
****  Thats all
      end
 
CTITLE GET_STATS
 
      subroutine get_stats( stats, pdiff, sdiff,
     .                max_site, max_center, num_site, num_center,
     .                ent, cmp, type)
 
      implicit none 

*     Rouitne to compute statistics
 
*         max_site      - Maximum number of sites
*           max_center  - Maximum number of centers
 
      integer*4 max_site, max_center
 
*   i,j,k,l     - Loop indices
*   num_site        - Number of sites
*   num_center      - Number of centers
*   ent             - Entry number (either site or center)
*   cmp             - Component (1=N, 2=E, 3=U)
 
*   num             - number of things to loop over
 
      integer*4 i, num_site, num_center, ent, cmp, num
 
*   pdiff(3,max_site, max_center)   - position adjustments
*   sdiff(3,max_site, max_center)   - Sigma for adjustments
 
*   stats(4)        - Number, Mean, RMS, Chi
*   dpos, spos      - Position and sigma estimate to use
*   summ,sumv,sumw  - Summation variables for stats/
 
      real*8 pdiff(3,max_site, max_center),
     .    sdiff(3,max_site, max_center), stats(4), dpos, spos,
     .    summ,sumv,sumw
 
 
      character*(*) type
 
****  Loop over data doing the statistics
      do i = 1,4
          stats(i) = 0
      end do
      summ = 0.d0
      sumv = 0.d0
      sumw = 0.d0
 
      if( type(1:1).eq.'S') then
          num = num_center
      else
          num = num_site
      end if
 
      do i = 1, num
          if( type(1:1).eq.'S') then
              dpos = pdiff(cmp,ent,i)
              spos = sdiff(cmp,ent,i)
          else
              dpos = pdiff(cmp,i,ent)
              spos = sdiff(cmp,i,ent)
          end if
 
          if( dpos.ne.999.9d0 ) then
              stats(1) = stats(1) + 1
              summ = summ + dpos/spos**2
              sumv = sumv + dpos**2/spos**2
              sumw = sumw + 1.d0/spos**2
          end if
      end do
 
****  See if have enough values
      if( stats(1).gt.1 ) then
          stats(2) = summ/sumw
          stats(4) = (sumv - stats(2)**2*sumw)/(stats(1)-1)
          stats(3) = sqrt((1.d0/sumw)*stats(1)*stats(4))
          stats(4) = sqrt(stats(4))
      end if
 
****  Thats all
      return
      end
 
 
 
CTITLE NSORT
 
      subroutine nsort( num, names, list)
 
      implicit none 

*     This routine uses an exchande sort algormithm to sort
*     the list names into ascending order.  There are num values
*     in names.  list starts as a sequence 1-n and returns with the
*     in the order of the sort.
 
*   num     - Number of values to be sorted
*   names(num)  - List to be sorted in to ascending order.
 
      integer*4 num, list(num)
 
 
      character*4 names(num)
 
* LOCAL VARIABLES
 
*   i,j     - Loop counters
*   smallest_one    - Smallest integer in current pass.
*   iswap   - Value used to swap integers
 
 
      integer*4 i,j, smallest_one, iswap
 
 
      character*4 swap
 
****  Start loop using exchange sort
 
      do i = 1, num
          smallest_one = i
          do j = i+1, num
              if( names(j).lt. names(smallest_one) ) then
                  smallest_one = j
              end if
          end do
 
*****     See if we should swap
          if( smallest_one.gt. i ) then
              swap = names(smallest_one)
              names(smallest_one) = names(i)
              names(i) = swap
 
              iswap = list(smallest_one)
              list(smallest_one) = list(i)
              list(i) = iswap
          end if
      end do
 
***** Thats all.  Now sorted in ascending order
      return
      end
 
CTITLE GPSW_DATE

      subroutine gpsw_date( gps_week, gps_dow, date_str )

      implicit none 

*     Convert gps_week and seconds to string with date

      integer*4 gps_week, gps_dow
      character*(*) date_str

*   gps_start_mjd - MJD of the start of the GPS week count (80/1/5)
      integer*4 gps_start_mjd, date(5)
      real*8 mjd, sectag

      character*3 months(12)


      data months / 'Jan','Feb','Mar','Apr','May','Jun',
     .              'Jul','Aug','Sep','Oct','Nov','Dec'  /
      data gps_start_mjd / 44243 /


      mjd = gps_start_mjd + gps_week*7 + gps_dow + 1
      call mjd_to_ymdhms( mjd, date, sectag)
      if( date(1).gt.500 ) date(1) = date(1) - 1900
      write(date_str,200) date(1), months(date(2)), date(3)
 200  format(i2,a3,i2.2)

****  Thats all
      return
      end


