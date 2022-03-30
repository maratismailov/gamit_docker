CTITLE INC_RENAME    
 
      subroutine inc_rename( ns, os )

      implicit none
 
*     This routine will check sites that have position changes
*     associated with their renames and add new entries to the
*     rename list so that after an earthquake occurrs the position
*     change will continue to be applied.  (This basically involves
*     adding the earthquake named site to the rn_code(2,*) entry
*     can copying the remainder of the rn_name entries.
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
      include '../includes/glb_hdr_def.h'

* PASSED VARIABLES

*   ns      - New site number in the global list of sites after
*             and earthquake occurs
*   os      - Old site number for the site name that actually 
*             appears in the binary hfiles.

      integer*4 ns, os
 
* LOCAL VARIABLES
 
*   i,j     - Loop counters
 
*   trimlen - Length of string
*   len_rnhf - Length of hfile name string.
*   org_num_renames - Initial number of renames.  We may change
*     number in this routine.
*   indx     - Position in string
 
      integer*4 i,j, k, trimlen, len_rnhf, org_num_renames, indx
 
*   use_rnh     - Set true when rn_hfiles set and the name matches.
*   update_needed - Set true to indicate that we need to add an
*                 entry to the rename list.

      logical use_rnh, update_needed

*     Run down the list of renames and see if any have position
*     changes at a site effected by an earthquake.  Run loop on
*     org_num_renames since num_renames may change
      org_num_renames = num_renames
     
      do i = 1, org_num_renames

****      See if there is a position change with this earthquake
*         We only need to update in this case so that the position
*         change will be later applied in rn_app_dpos

          if( abs(rn_dpos(1,i)) +  abs(rn_dpos(2,i)) +  
     .        abs(rn_dpos(3,i)) .gt. 0.d0 ) then

*             OK, position change for this one, now see if sites
*             involved are for our stations and this named h-file.
*             Check to see if hfile name match needed.
              use_rnh = .true.
              if( trimlen(rn_hfiles(i)).gt.0 ) then
                  len_rnhf = trimlen(rn_hfiles(i))
                  indx = index(glb_inp_file,rn_hfiles(i)(1:len_rnhf)) 
                  if( indx.eq.0 ) then
                     use_rnh = .false.
                  end if
              endif

*             Now see if epoch range and site names match (plus that
*             we are using this hfile)
* MOD TAH 130506: Changed overlap test so that mid-point name changes
*              will be assigned to either the rename before after so that
*             it is not left in an un-renamed state.
C             if( cepoch_start.ge.rn_times(1,i)   .and.
C    .            cepoch_end.le.rn_times(2,i)     .and. 
C    .            gsite_names(os).eq.rn_codes(1,i) .and. use_rnh ) then
              if( cepoch_start-rn_times(1,i).gt. -0.5d0  .and.
     .            cepoch_end-rn_times(2,i).le. +0.5d0  .and.
     .            gsite_names(os).eq.rn_codes(1,i) .and. use_rnh ) then

*                 OK, site effected, now see if we already have the
*                 entry we need (either automatically entered or 
*                 because the user put the entry in).
                  update_needed = .true.
                  do j = 1, num_renames

*                    Check all parts of the rename entry.  If they all
*                    match then set the update_needed to false (use 
*                    new site name). 
                     if( gsite_names(os).eq. rn_codes(1,j) .and.
     .                   gsite_names(ns).eq. rn_codes(2,j) .and.
     .                   rn_times(1,i) .eq. rn_times(1,j) .and.
     .                   rn_times(2,i) .eq. rn_times(2,j) .and.
     .                   rn_hfiles(i)  .eq. rn_hfiles(j)  .and.
     .                   rn_dpos(1,i)  .eq. rn_dpos(1,j)  .and.
     .                   rn_dpos(2,i)  .eq. rn_dpos(2,j)  .and.
     .                   rn_dpos(3,i)  .eq. rn_dpos(3,j)  .and.
     .                   rn_types(i)   .eq. rn_types(j)        ) then
                         update_needed = .false.
                     end if
                  end do
 
*                 If an update needed then add new entry
                  if( update_needed ) then
                     if( num_renames+1.gt.max_rn ) then
                        write(*,200) max_rn, gsite_names(ns)
 200                    format('**WARNING** Too many renames (',I4,
     .                         ') need to be added for ',a,
     .                         ' post-earthquake.',
     .                         ' Additional entry ignored')
                        if( log_unit.ne.6 ) 
     .                  write(*,200) max_rn, gsite_names(ns)
                     else

                        num_renames = num_renames + 1
                        k = num_renames
                        rn_codes(1,k) = gsite_names(os)
                        rn_codes(2,k) = gsite_names(ns)
                        if( rn_times(1,k).lt. cepoch_start ) then
                            rn_times(1,k) = cepoch_start  
                        else
                            rn_times(1,k) = rn_times(1,i)
                        end if
                        if( rn_times(2,k).gt. cepoch_end ) then
                            rn_times(2,k) = cepoch_end   
                        else
                            rn_times(2,k) = rn_times(2,i)
                        end if
                        rn_hfiles(k)  = rn_hfiles(i) 
                        rn_dpos(1,k)  = rn_dpos(1,i) 
                        rn_dpos(2,k)  = rn_dpos(2,i)
                        rn_dpos(3,k)  = rn_dpos(3,i)
                        rn_types(k)   = rn_types(i) 
                     end if
                  end if
              end if
          end if
      end do

****  Thats all 
      return
      end

