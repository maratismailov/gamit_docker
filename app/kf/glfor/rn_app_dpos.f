CTITLE RN_APP_DPOS
 
      subroutine rn_app_dpos(sol_obs)

      implicit none 
 
*     Routine to apply the change in position associated woth
*     site renames.  This feature is to many handle the instrument
*     height changes that have not been correctly put into the binary
*     h-files.
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
      include '../includes/glb_hdr_def.h'
 
* PASSED Variables
 
*      sol_obs(cnum_parn)       - Solution vector read from
*                       - binary hfile
 
      real*8 sol_obs(cnum_parn)
 
* LOCAL variables.
 
*   i,j     - Loop counters
*   ns, ls  - Gloval and local site numbers
*   indx    - Postion in string and site numer
*   type    - Type of parameter
*   np      - parameter number to be changed
 
      integer*4 i,j, ns, ls, indx, type, np
 
*   loc_coord(3)        - Local coordinates of the site
*   rot_matrix(3,3)     - Rotation matrix from rn_types system
*                       - to XYZ
*   dxyz(3)             - XYZ coordinate shift (m)
 
      real*8 loc_coord(3), rot_matrix(3,3), dxyz(3)
 
*       done            - Set true when we have applied the
*                       - Z shift to the coordinates (Error
*                       - printed if Z was not estimated)
*       kbit            - Checks bit is set
 
 
      logical done, kbit
 
****  First see if the times match
      do i = 1, num_renames
        if( cepoch_start.ge. rn_times(1,i) .and.
     .      cepoch_end.le.rn_times(2,i)     ) then
 
*         This experiment falls in the correct time range.
*         See if we need to change position
          if( rn_dpos(1,i).ne.0.0d0 .or. rn_dpos(2,i).ne.0.d0 .or.
     .        rn_dpos(3,i).ne.0.d0 ) then
 
*             A position change is needed, so see if station used.
              indx = 1
              call get_cmd(rn_codes(2,i),gsite_names, gnum_sites,
     .                ns, indx)

              if( ns.eq.-2 ) then
                  write(*,100) rn_codes(2,i), (gsite_names(j),
     .                         j = 1, gnum_sites)
 100              format('WARNING: Duplicate name for ',a8,
     .                   ' Site names:',1000(/,10(a8,1x)))
              end if
 
*             See if the site is used this time
              ls = -1
              do j = 1, cnum_sites
                  if( ltog_sites(j).eq.ns ) then
                      ls = j
                  end if
              end do

*             Only change the position if the name of this site
*             has been changed
* MOD TAH 001004: Added additional check to make sure that the rename
*             has been used before applying the correction.  This
*             treats the case when two different sites are renamed
*             and moved to the same site name in overlapping intervals
*             of time.
              if( .not.kbit(rn_name_changed,ns) .or.
     .            .not.kbit(rn_used,i)              ) ls = 0
 
*             If ls is greater than zero then the site is used so
*             apply correction
              if( ls .gt.0 ) then
 
*                 Convert the postion correction to XYZ if it is in
*                 NEU (we save this to now so that we will use the most
*                 accurate coordinates)
                  call rotate_geod(rn_dpos(1,i), dxyz, rn_types(i),
     .                    'XYZ', apr_val_site(1,1,ns), loc_coord,
     .                     rot_matrix)
 
*                 OK Now find the parameter and apply the correction
                  done = .false.
                  np = 0
                  do while ( .not.done .and.np.le.cnum_parn)
                      np = np + 1
                      call decode_code(gpar_codes(np), type, indx )
                      if( type.eq.7 .and. indx.eq.ls ) then
                          sol_obs(np) = sol_obs(np) + dxyz(1)
                      end if
                      if( type.eq.8 .and. indx.eq.ls ) then
                          sol_obs(np) = sol_obs(np) + dxyz(2)
                      end if
                      if( type.eq.9 .and. indx.eq.ls ) then
                          sol_obs(np) = sol_obs(np) + dxyz(3)
                          write(*,9000) rn_codes(1,i), gsite_names(ns),
     .                                  dxyz
9000                      format('Updating from ',a8,' to ',a8,' by ',
     .                           3F10.4,' m')
                          done = .true.
                      end if
                  end do
 
*                 Print error if we had a problem i.e., if done is
*                 still false.
                  if( .not.done ) then
                      write(*,120) i, rn_codes(2,i), ns, ls
 120              format('**WARNING** Problem applying rename ',
     .                    i4,1x,A8,' Global site ',i4,' Local site ',
     .                    i4,/,
     .                    12x,' Did not find station coordinates in',
     .                        ' parameter list')
                  end if
              end if
          end if
        end if
      end do
 
****  Thats all
      return
      end
 
 
