Ctitle Save_head
      subroutine  save_head(i)

      implicit none

*     Subroutine to Save head information for each data filef
*     GANG CHEN   of MIT, 1997      

      include '../includes/xfile_def.h'
      include 'track_com.h'
      
***   i --- the index of file
      integer*4 i
***   Local valuables
      integer*4 j, k, trimlen

* int_uset -- Logical that is true if user set interval
      logical int_uset

      data  int_uset / .true. /

      save int_uset

*     When not found in bat file , use site position at input file as default
      if(.not.read_site_apr(i)) then
           do j=1, 3
               site_apr(j,i) = xf_pos(j)
               site_int(j,i) = xf_pos(j)
               site_vel(j,i) = 0.0d0
           end do 
           site_ep(i) =  51543.d0
      endif

****  Rinex file: use original data search time internal
*     If the user has not passed an interval to use, use the one from the
*     rinex file                                                               
      if( usr_interval.eq.0 ) then
          usr_interval = xf_ircint
          int_uset = .false. 
      else
* MOD TAH 070108: Check to see if this interval is more ot less
*         currently set.  Largest interval will be ser by default
          if( xf_ircint .gt. usr_interval .and. .not. int_uset ) then
              write(*,100) usr_interval,  xf_ircint, 
     .                     obs_file(i)(1:trimlen(obs_file(i)))
 100          format('**WARNIMG** Replacing ',f6.2,' sampling rate',
     .               ' with ',F6.2,' sec from file ',a,/,
     .               '** Use INTERVAL command to override')
              usr_interval = xf_ircint
          end if
      end if

*     Check for RX files without usr_interval
      if( xf_ircint.le.0 ) then
          write(*,120) usr_interval
 120      format('**WARNING** Interval from RX files is zero, using ',
     .           f6.2,' s user interval')
          xf_ircint = usr_interval
          if( usr_interval.eq.0 ) then
              write(*,140) 
 140          format('**DISASTER** user inverval is zero,',
     .               ' use INTERVAL command to set')
              stop 'TRACK: No data spacing interval'
          end if
      end if
      
      sv_ndat(i) = xf_ndat
      sv_mdat(i) = xf_mdat

* MOD TAH 180322: Save for each GNSS (xf_mdat now limited to 6)
      do j=1, xf_mdat
         do k = 1, 7
            sv_dattyp(j,k,i)= xf_dattyp(j,k)
         end do
      enddo
      
      do j=1, 3
         sv_offarp(j,i)=xf_offarp(j)
      enddo

****  Thats all
      return
      end
      
c______________________________________________________________________

Ctitle get_headinfor
      subroutine  get_headinfor(i)

      implicit none

*     Subroutine to re-gain head information from sv_*
*     GANG CHEN   of MIT, 1997
      include '../includes/xfile_def.h'
      include 'track_com.h'
      
***   i --- the index of file
      integer*4 i
***   Local valuables
      integer*4 j, k


      xf_ndat=sv_ndat(i) 
      xf_mdat=sv_mdat(i)
      
* MOD TAH 180322: Restore the header data types for the RINEX file
*     being read
      do j=1, xf_mdat
         do k = 1, 7
            xf_dattyp(j,k) = sv_dattyp(j,k,i) 
         end do
      enddo

C     do j=1, xf_ndat
C       xf_dattyp(j)= sv_dattyp(j,i) 
C     enddo
      
      return
      end
     
