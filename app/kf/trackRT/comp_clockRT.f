
CTITLE EPH_TO_XYZRT
 
      subroutine eph_to_xyzRT(d, t, i, sys, ep)

      implicit none
 
*     Routine to compute XYZ coordinates of satellite at time t for
*     sattellite number i.  SYS is 'I' or 'E' for inertial or Earth
*     fixed.
* MOD TAH 110504: Added check to exclude PRN if orbit is too low an
*     altitude
 
      include 'trackRT.h'

*  order - Order of the interpolating polynomial.  Must be even 
*          value

      integer*4 order

      parameter ( order = 14 )

*   d       - MJD of the day 
*   t       - Time for computation (Seconds from start of sp3 file).
  
      real*8 d, t
 
*   i       - Satellite number
*   ep      - Epoch numbers
      integer*4 i, ep
 
*   sys     - System for results.
 
 
      character*(*) sys
 
* LOCAL VARIABLES
 
*   j           - Loop counter in Kepler's equation
 
 
      integer*4 j, k,  bin, ein
 
* int_xyz(3) - Inertial coordinates of satellite.  Copied back to
*              return variable

      real*8 err, rot_mat(3,3), xa(100), ya(100), za(100), 
     .       int_xyz(3)
      real*8 bado_err(3)   ! Error in orbit maked as bad

      real*8 rad    ! Radius of SVS 

      logical small_rad  ! Set true if any radius values are small
      logical inlist     ! Set true if SVS in exclude list
      logical bad_orbit  ! Set true if orbits is bad (in exclude) or
                         ! has interpolation errors.
      integer*4 pn       ! PRN number


* BE_err(32) -- Set true once an interpolation error has been 
*     reported.
      logical BE_err(32)
  
      logical found

      data BE_err / 32*.false. /
      save BE_err

****  Find out which two points are closest to the epoch
      pn = i   ! Assign PRN number
      j = 0
      found = .false.
      do while ( .not.found .and. j.le.num_sp3)
          j = j + 1
          if( sp3_sec(j).gt. t  ) then
              found = .true.
          end if
      end do
 
*     Get the start index for interpolation
      bin = j - order/2
      if( bin.le.0 ) bin = 1

*     Check we have a valid entry
      do while ( sp3_xyz(1,i,bin).eq.0 .and. bin.lt.num_sp3 )
          bin = bin + 1
          if( .not. BE_err(i) ) then
             print *,'**WARNING** Adjusting Start Bin for PRN ',i,
     .               ' Bin ',bin, ' Seconds past SP3 start ', t
             print *,'**WARNING** Message only printed once'
             BE_err(i) = .true.
          endif 

      end do

      ein = bin + order - 1
      if( ein.gt.num_sp3 ) then
          ein = num_sp3
          bin = ein - (order-1)
      end if

****  Check we have a valid end of sp3 records
      do while ( sp3_xyz(1,i,ein).eq.0 .and. ein.gt.1 )
          ein = ein - 1
          bin = bin - 1
          if( .not. BE_err(i) ) then
             print *,'**WARNING** Adjusting End Bin for PRN ',i,
     .               ' Bin ',bin, ' Seconds past SP3 start ', t
             print *,'**WARNING** Message only printed once'
             BE_err(i) = .true.
          endif 
      end do
      
 
***** Now do interporlation
      small_rad = .false.
      bad_orbit = BE_err(pn)

      do j = 0, order-1
          xa(j+1) = sp3_xyz(1,i,bin+j)
          ya(j+1) = sp3_xyz(2,i,bin+j)
          za(j+1) = sp3_xyz(3,i,bin+j)
          rad = sqrt(xa(j+1)**2+ya(j+1)**2+za(j+1)**2)
          if( rad < 6.4d6 ) small_rad = .true.
      end do

****  See if we should exclude this satellite or re-include it if
*     is now OK
*     See if have a bad orbit due to interpolation that may be
*     good again
      if( BE_err(pn) ) then

          call polint(sp3_sec(bin), xa, order, t,
     .               svs_xyz(1,i), bado_err(1) )
          call polint(sp3_sec(bin), ya, order, t,
     .               svs_xyz(2,i), bado_err(2) ) 
          call polint(sp3_sec(bin), za, order, t,
     .               svs_xyz(3,i), bado_err(3))
*         See if orbit looks OK
          if( bado_err(1)**2+bado_err(2)**2+bado_err(3).lt.0.01 ) then
              inlist = .false.
              do j = 1, num_exclude
                  if( -pn.eq.ss_exclude(j) ) inlist = .true.
              end do

!             See if already excluded: Good now so add back
              if ( inlist ) then
                 do k = j+1, num_exclude
                    ss_exclude(k-1) = ss_exclude(k)
                 end do 
                 num_exclude = num_exclude - 1
                 call report_stat('warning','trackrt','eph_to_xyzRT',
     .             ' ','PRN included again after bad orbit PRN ',pn)
                 BE_err(pn) = .false.
              endif
          end if
       end if

****  See if other probems with orbit
      if( small_rad .or. BE_err(pn) ) then
!         See if in list
          inlist = .false.
          do j = 1, num_exclude
              if( pn.eq.abs(ss_exclude(j)) ) inlist = .true.
          end do
!         See if already excluded
          if ( .not. inlist ) then
              num_exclude = num_exclude + 1
              ss_exclude(num_exclude) = -pn  ! Use -ve to show automatic
              call report_stat('warning','trackrt','eph_to_xyzRT',
     .             ' ','PRN excluded due to bad orbit PRN ',pn)
          endif
          bad_orbit = .true.
      else
!         Normal case, see of excluded and should be put back
          inlist = .false.
          do j = 1, num_exclude
              if( -pn.eq.ss_exclude(j) ) inlist = .true.
          end do

!         See if already excluded: Good now so add back
          if ( inlist ) then
              do k = j+1, num_exclude
                 ss_exclude(k-1) = ss_exclude(k)
              end do 
              num_exclude = num_exclude - 1
              call report_stat('warning','trackrt','eph_to_xyzRT',
     .             ' ','PRN included again after bad orbit PRN ',pn)
          endif
      end if

! MOD TAH 110506: See if PRN excluded
      do j = 1, num_exclude
          if( pn.eq.abs(ss_exclude(j)) ) bad_orbit = .true.
      end do
      
      if( ep.ge.debug(1) .and. ep.le.debug(2) ) then
          write(*,320) ep, d,t,i, bin, sp3_sec(1), bad_orbit, small_rad
 320      format('MJD, Sec, PRN ',i6,1x,F15.6,1x,F16.3,1x,I2,1x,i4,
     .                           1x,F16.3,' Status ',2L3)
          write(*,340) xa(1),ya(1),za(1), xa(order),ya(order),za(order)
 340      format('SP3 Pos', 6F16.3)
          if( num_exclude.gt.0 )
     .    write(*,360) num_exclude,(ss_exclude(k),k=1, num_exclude)
 360      format('Exluded Satellites # ',I3,' List ',32(I4.2))
      end if

****  Only process if orbit OK             
      if( .not. small_rad .and. .not. bad_orbit ) then 

         call polint(sp3_sec(bin), xa, order, t,
     .               svs_xyz(1,i), err)
         if( abs(err).gt.1.0 ) write(*,500) i,ep, d,t, 'X', err
         call polint(sp3_sec(bin), ya, order, t,
     .               svs_xyz(2,i), err)
         if( abs(err).gt.1.0 ) write(*,500) i,ep,d,t, 'Y', err
 
         call polint(sp3_sec(bin), za, order, t,
     .               svs_xyz(3,i), err)
         if( abs(err).gt.1.0 ) write(*,500) i,ep,d,t, 'Z', err
  500    format(' Interpolation error PRN ',i2,' Epoch ',i6,
     .       ' MJD ', f12.4,' Sec ',e16.6,
     .       ' Component ',a1,' Magnitude ', e16.6,' m')

         if( abs(err).gt.10.0 ) then
             write(*,*) 'Interpolation errors too large. Check SP3 ',
     .        'file or used EXCLUDE_SVS command to remove satellite'
             do j = 0, order-1
                write(*,450) j, bin+j, xa(j+1), ya(j+1), za(j+1)
 450            format('Values: ',2i4,' XYZ ',3F15.3)
             end do 
! MOD TAH 110506: Removed stop and added ss_exclude
!             stop 'TRACK: SP3 interpolation errors'
              num_exclude = num_exclude + 1
              ss_exclude(num_exclude) = -pn
              call report_stat('warning','trackrt','eph_to_xyzRT',
     .             ' ','PRN exclude to Interpolation error PRN ',pn)

         end if
     
*****    Check to see if Inertial coordinates wanted.
         if( sys(1:1).eq.'I' ) then
             call earth_to_inert(d, svs_xyz(1,i), int_xyz, 'E','I')
*            Copy the interial coordinates to the return position
*            variable
             do j = 1, 3
                svs_xyz(j,i) = int_xyz(j)
             end do
         end if   

*        Now convert to latitude and longitude
         call xyz_to_geod( rot_mat, svs_xyz(1,i), svs_loc(1,i))
      else
         svs_xyz(1,i) = 26000.0d3
         svs_xyz(2,i) =     0.0d0
         svs_xyz(3,i) =     0.0d0
      endif
 
****  Thats all
      return
      end
 
CTITLE POLINT
 
      subroutine polint(xa,ya,n,x,y,dy)

      implicit none
 
*     Routine to do an 11 point interpolation of the SP3
*     file entries.  (From Numerical Recipes)
 
      integer*4 nmax, n, i,m, ns
 
      parameter (nmax=21)
 
      real*8 xa(n), ya(n), c(nmax),d(nmax), dif, dift,
     .    den, w, hp, ho, x,y, dy
 
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.) then
             print *,'*ERROR* Interpolation error in polint: Den is 0'
             den = 1.d-10
          endif
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      end
 
      
CTITLE comp_svs_clkRT

      subroutine comp_svs_clkRT( epoch )

      implicit none
      
*     Rouitne to compute satellite clock coorection at time epoch

      include 'trackRT.h'
 
* epoch - MJD at which clocks are needed

      real*8 epoch
      
* LOCAL VARIABLES
*   i,j - Loop counters
*   nn  - Tabular point immediately after data_epoch
*       - in the SP3 file.
 
      integer*4 i,nn
 
*   dt  - Time between the data_epoch and SP3 entry
*       - immediately before it (days)
*   dt_sp3  - Spacing between SP3 entries (days)
*   dc_sp3  - Difference in the SP3 clocks (sec)

      real*8 dt, dt_sp3, dc_sp3
      
***** First evaluate the satellie clock corrections at this time
*     Find the SP3 entries that surrond the current data_epoch
      nn = 2
      do while ( sp3_time(nn).lt. epoch .and.
     .        nn.lt.num_sp3 )
          nn = nn + 1
      end do
 
****  The two entries needed are nn-1 and nn.  Warn user if these
*     are actually outside the range of the data
      if( sp3_time(nn-1).gt.epoch .or.
     .    sp3_time(nn).lt. epoch ) then
          write(*,120) epoch, nn, num_sp3
 120      format('WARNING: Data epoch ',f12.4,' before or after ',
     .        'SP3 file range. Using point ',i4,' of ',i4)
      end if
 
****  Now get the interpolated satellite time
      dt = epoch - sp3_time(nn-1)
      dt_sp3 = sp3_time(nn) - sp3_time(nn-1)
      do i = 1, num_sat
          dc_sp3 = sp3_clk(i,nn) - sp3_clk(i,nn-1)
          svs_clk(i) = sp3_clk(i,nn-1) + (dc_sp3/dt_sp3)*dt
      end do
      
***  Thats all
      return
      end 
