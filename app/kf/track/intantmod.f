CTITLE IntAntMod

      real*8 function IntAntMod(num_dph, dzn,
     .                   dphs, zd, az, max_zen, max_az, ob)

      implicit none 

*     Routine to linearly interpolate phase center values from
*     tables

* PASSED VARIABLES
      
      integer*4 num_dph(2)  ! Number of entries in ZD and Azimuth
      integer*4 max_zen, max_az  ! Dimenisioning in ZD and AZ directions
      real*8 dzn(3,2)       ! Limits in Zenith distance and Azimuth
*                           ! Z1, Z2 and dZD and 0, 360, dAZ (if dAz
*                           ! is zero then no az dependence)
      real*8 zd, az         ! Zenith distance or Nadir angle to be
                            ! interpolated to.
      real*4 dphs(max_zen, max_az) ! Table to be interpolated (m)

      logical ob            ! Set true if value out of bounds


* LOCAL VARIABLES
      integer*4 inzd, inaz  ! bins that point to tabular point 
                            ! at lower end of current value for ZD and AZ
      real*8 Z1, Z2, ZA     ! Interpolation in ZD at two azimuths and 
                            ! with azimuth



****  See if we have tables
      ob = .false.
      if( num_dph(1) .eq. 0 ) then
        IntAntMod = 0.d0
        RETURN
      end if

****  Compute the ZD index to interpolate on:
      if( zd.lt.dzn(2,1) ) then
          inzd = int((zd-dzn(1,1))/dzn(3,1))+1
      else
C          write(*,120) zd, dzn
 120      format('**WARNING** ZD ',F8.2,' Exceeds table ',6F8.2)
          inzd = int((dzn(2,1)-dzn(1,1))/dzn(3,1))
*         Only report problem is more than 0.5 degrees off
          if( zd-dzn(2,1).gt.0.5 ) ob = .true.
      endif

****  See if we need to interpolate in ZD and AZ or just ZD
      if( dzn(3,2).gt.0 ) then
         inaz = int((az-dzn(1,2))/dzn(3,2))+1
      else
         inaz = 0
      endif

****  OK: Do the intepolation
      if( inaz.eq. 0 ) then 
         Z1 = dphs(inzd,1)+(dphs(inzd+1,1)-dphs(inzd,1))/dzn(3,1)*
     .                                     (zd-(inzd-1)*dzn(3,1))
         IntAntMod = Z1
      else
         Z1 = dphs(inzd,inaz)+
     .              (dphs(inzd+1,inaz)-dphs(inzd,inaz))/dzn(3,1)*
     .                                     (zd-(inzd-1)*dzn(3,1))
 
         Z2 = dphs(inzd,inaz+1)+
     .          (dphs(inzd+1,inaz+1)-dphs(inzd,inaz+1))/dzn(3,1)*
     .                                     (zd-(inzd-1)*dzn(3,1))
         ZA = Z1 + (Z2-Z1)/dzn(3,2)*(az-(inaz-1)*dzn(3,2))
         IntAntMod = ZA
      endif

****  Thats all
      return
      end

