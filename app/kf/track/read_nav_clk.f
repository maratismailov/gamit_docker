CTITLE READ_nav
 
      subroutine read_nav_clk

      implicit none
 
*     Thuis routine will read Nav file and get the clock information
*     from this file.  The polynomials are read and then computed
*     at the epochs of the sp3 file epochs.
 
      include '../includes/const_param.h'
      include 'track_com.h'
 
* LOCAL VARIABLES
 
*   ierr        - IOSTAT error
*   i           - Loop counter
*   rp          - Read Prn number
*   trimlen     - length of string
*   ll          - Length of line string
*   date(5)     - ymdhm
*   nn          - Entry for nav clock number
 
      integer*4 ierr, i, j, trimlen,  date(5), pn, nn
 
*   sectag      - Seconds tag in date
*   af(3)       - Coeffocients read from file
*   tsp3        - Time of SP3 entries
*   dt          - Time difference from nav clk to tsp3 (seconds)
*   clk         - Clk estimate from Nav file
*   dclk        - Difference in clock time (usec)
*   mjd          - Julian date of measurement
*   tsp3        - Time from SP3 entries
 
      real*8 sectag, af(3), dt, clk, dclk, mjd, tsp3
 
*   still_header    - Indicates that we are sill in the header
*                   - part of file.
 
      logical still_header
 
*   line            - line read from file
 
      character*256 line
 
*   cr - Carriage return (for handling dos files)
 
      character*1 cr
 
      cr = char(13)
 
****  Open the sp3_file
 
      open(100, file=sv_clk_file, status='old', iostat=ierr)
      call report_error('IOSTAT',ierr,'open',nav_file,1,
     .            'maked/read_nav')
 
 
      still_header = .true.
 
      do while ( still_header )
          read(100,'(a)', iostat=ierr) line
          if( index(line,'END OF HEADER').gt.0 ) still_header = .false.
          if( trimlen(line).eq.0 ) still_header = .false.
          if( ierr.ne.0 ) then
              write(*,120)sv_clk_file(1:trimlen(sv_clk_file))
 120          format(' Error reading ',a,' before end of ',
     .               'header found')
              stop 'read_nav_clk: End of file before end of header'
           end if
      end do
 
      do i = 1, num_sat
         num_nav(i) = 0
      end do
 
*     Now start reading the entries
      do while ( ierr.eq.0 )
 
*         Read down the elements of the satellites
          read(100,'(a)',iostat=ierr) line
          if( line(1:3).ne.'   ' .and. ierr.eq.0 ) then
 
*             This must be a PRN line, so get clock info
              read(line, 220) pn, date, sectag, af
 220          format(i2,5i3,f5.1,3d19.8)
 
*             Get the MJD of the clock reference
              call ymdhms_to_mjd( date, sectag, mjd)
 
*             Save the value for this satellite
              num_nav(pn) = num_nav(pn) + 1
c              nav_clk(1,pn,num_nav(pn)) = mjd
              nav_clk(1,num_nav(pn),pn) = mjd
c	write(*,220) date, sectag, pn,  num_nav(pn)
	      do i = 1,3
c                  nav_clk(i+1,num_nav(pn),pn) = af(i)
                  nav_clk(i+1,num_nav(pn),pn) = af(i)
              end do
          end if
      end do

****  Now sort each the nav entries into time order
      do i = 1, num_sat
c	 write(*,*)  i, num_nav(i), nav_clk(1,i,1)
	          call sort_nav(num_nav(i), nav_clk(1,1,i))
c	 write(*,*)  num_nav(i), nav_clk(1,i,1)
*        Now loop over the sp3 file entries interpolating the clocks
         if( num_nav(i).gt.0 ) then
             do j = 1, num_sp3
                tsp3 = sp3_time(j)
 
*               Loop over entries until we find one after the desired time.
*               At which time we use the one before it.
                nn = 1
                if( tsp3.gt.nav_clk(1,nn,i) ) then ! Search
                    do while ( tsp3.lt.nav_clk(1,nn,i) .and.
     .                         nn.lt.num_nav(i) )
                       nn = nn + 1
                    end do
                end if
 
                dt = (tsp3-nav_clk(1,nn,i))*86400.d0
                clk = nav_clk(2,nn,i) +
     .                nav_clk(3,nn,i)*dt +
     .                nav_clk(4,nn,i)*dt*dt/2.d0
c               write(*,300)j, i, tsp3, clk,dt/86400.0
c     .                       ,(nav_clk(k,nn,i),k=1,4)
 300                    format(i3,' Prn ',i2,' Time ',f10.3,
     .                      1x,   f20.10,' sec ',4f20.6)

*               See how big difference is
                if( sp3_clk(i,j).lt. 99999.999d-6 ) then
c                     dclk = (sp3_clk(i,j)-clk)
                    dclk = (sp3_clk(i,j)-clk)*1.d6
                    if( abs(dclk).gt.1.d0 ) then
                        call jd_to_ymdhms(tsp3,date, sectag)
                        write(*,400) i, date, sectag, dclk
 400                    format('CLKDIFF: Prn ',i2,' Date ',i5,4i3,
     .                         f8.3,1x,f10.2,' usec')
                    end if
                end if
                sp3_clk(i,j) = clk
              end do
           end if
      end do
 
****  Thats all
      close(100)
      return
      end
 
CTITLE SORT_NAV
 
      subroutine sort_nav( num, nav_clk)
 
      implicit none

*     Routine to time sort the nav clk entries.
 
*
* num    - Number of entries
 
 
      integer*4 num
 
* nav_clk(4,num)  - Values of the clocks
 
 
      real*8 nav_clk(4,num)
 
*LOCAL
      integer*4 i, j, smallest_one
 
 
      real*8 swap
 
*     Sort the values
      do i = 1, num
          smallest_one = i
          do j = i+1, num
              if( nav_clk(1,j).lt. nav_clk(1,smallest_one) ) then
                  smallest_one = j
              end if
          end do
 
*****     See if we should swap
          if( smallest_one.gt. i ) then
              do j = 1,4
                  swap = nav_clk(j,smallest_one)
                  nav_clk(j,smallest_one) = nav_clk(j,i)
                  nav_clk(j,i) = swap
              end do
          end if
      end do
 
***** Thats all
      return
      end
 
