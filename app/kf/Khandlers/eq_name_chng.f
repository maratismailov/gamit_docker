CTITLE EQ_NAME_CHANGE
 
      subroutine eq_name_change(opt)

      implicit none
 
*     This routine will check the earthquake list and rename
*     lists to see if we should rename the sites in the current
*     list of sites

* MOD TAH 950613: Changed to keep track of renames when an earthquake
*     also occurrs.  If a site is renamed and then affected by an
*     earthquake, the rn_name_changed will be propagated forwarded
*     to the new earthquake renamed site name.
* MOD TAH 080910: Added check to make sure that multi-day h-files
*     that have had earthquake renames applied do not have the
*     renamed sites converted back to origonal name.  In original
*     code, this happens when an earthquake occurrs during the 
*     interval of the h-file.
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
      include '../includes/glb_hdr_def.h'
      
* PASSED VARIABLES
*   opt     - Option set to update coordinates in this routine.
*             For the first set of calls in glinit and glist
*             the coordinates should be set to update, after that
*             the NO option should be used so that coordinates are
*             not changed

       character*(*) opt     
 
* LOCAL VARIABLES
 
*   i,j     - Loop counters
*   ns, os  - New and old site numbers when we rename
*   indx    - Position in string
 
*   last_eq     - Earth number for the last earthquake to
*           - to effect the site being considered.
*   trimlen - Length of string
*   len_rnhf - Length of hfile name string.
*   len_rns  - Length of rename site (used to allow short names)
 
      integer*4 i,j, k, ns, os, indx, last_eq, trimlen,
     .          len_rnhf, len_rns
 
*   last_eqep   - Epoch of the last earthquake
*   dist        - Distance from epicenter to site
 
      real*8 last_eqep, dist

*   kbit        - Checks status of bit
*   duplicate   - Set true if the rename of site will generate a
*                 site already in this experiment.
*   use_rnh     - Set true when rn_hfiles set and the name matches.
*   ud          - Logical set true to update coordinates (based on
*                 opt passed to routine)

      logical kbit, duplicate, use_rnh, ud

*   new_name    - New name of site after an earthquake.
 
      character*8 new_name

*   norename(max_eq) -- List of earthquakes codes not to
*     have the rename applied for because they fall inside
*     hfile.
*   num_norename -- Number of no-renames
*   noeqrn  -- Logical set true if no EQ rename needed

      integer*4 num_norename
      character*4 norename(max_eq)
      logical noeqrn 

* MOD TAH 000902: See what option passed by user and set the
*     update logical
      ud = .true.
      if( opt(1:2).eq.'NO' .or. opt(1:2).eq.'no' ) ud = .false.
      
****  First clear the array that says the name has changed
      do i = 1, max_glb_site_wrds
         rn_name_changed(i) = 0
      end do

* MOD TAH 080910: See how many earthqukes fall inside this hfile
      num_norename = 0
      do i = 1, num_eq
         if ( eq_epoch(i).ge.cepoch_start .and.
     .        eq_epoch(i).lt.cepoch_end        ) then
*          EQ is inside hfile.  Do no rename sites with this name
           num_norename = num_norename + 1
           norename(num_norename) = eq_codes(i)
         endif
      enddo

 
****  Now start the specific rename operations.  We do
*     this first so that the site can be further renamed
*     due to earthquake processes.
 
      do i = 1, num_renames
 
*         See if this experiment is in the correct time
*         range of the rename
* MOD TAH 971112: Check to see if hfile name match needed.
          use_rnh = .true.
          if( trimlen(rn_hfiles(i)).gt.0 ) then
              len_rnhf = trimlen(rn_hfiles(i))
              indx = index(glb_inp_file,rn_hfiles(i)(1:len_rnhf)) 
              if( indx.gt.0 ) then
                 use_rnh = .true.
              else
                 use_rnh = .false.
              end if
          endif

* MOD TAH 081028: See if short first name
          len_rns = trimlen(rn_codes(1,i))
          if( rn_codes(1,i)(len_rns:len_rns).eq.'@' .or.
     .        rn_codes(1,i)(len_rns:len_rns).eq.'*' ) then
              len_rns = len_rns -1 
              if( len_rns.eq.0 ) then
                  call report_stat('WARNING','GLOBK','eq_name_chng',
     .                 'rename code','Only wild-card @ name',0)
                  len_rns = 8
              endif
          endif

* MOD TAH 130506: Changed overlap test so that mid-point name changes
*         will be assigned to either the rename before after so that
*         it is not left in an un-renamed state.
C         if( cepoch_start.ge.rn_times(1,i) .and.
C    .        cepoch_end.le.rn_times(2,i)   .and. use_rnh ) then
          if( cepoch_start-rn_times(1,i).gt. -0.5d0  .and.
     .        cepoch_end-rn_times(2,i).le. +0.5d0  .and. use_rnh ) then
 
*             OK time range good.  See if we have this
*             station.
              do j = 1, cnum_sites
* MOD TAH 131205: Only allow the rename if the gsite_name is not an
*                 _X.. site.  Covers cases when XPS rename starts 
*                 before later renames (rare).  RESET will be needed
*                 allow saved _X.. sites to be renamed.
                  if( gsite_names(ltog_sites(j))(1:len_rns).eq.
     .                rn_codes(1,i)(1:len_rns) .and.
     .                gsite_names(ltog_sites(j))(6:6).ne.'X' ) then
 
*                     Yes, local site j needs to be renamed.
*                     See if we already have the new name
*                     in the list of stations.  If we do just
*                     update the name, other wise add the new
*                     name to the list
* MOD TAH 980514: Set bit to show rename used.
                      call sbit(rn_used,i,1)

                      indx = 1
                      call get_cmd(rn_codes(2,i), gsite_names,
     .                    gnum_sites, ns, indx)

*                     Now scan the list of sites in this h-file and
*                     make sure that the new name does not match
*                     a site already in this file.  (This can happen
*                     in combined globals).
                      duplicate = .false.
                      do k = 1, cnum_sites 
                          if( gsite_names(ltog_sites(k)).eq.
     .                        rn_codes(2,i) .and. j.ne.k ) 
     .                                              duplicate = .true.
                      end do
                      if ( duplicate ) then
                          write(*,105) rn_codes(1,i), rn_codes(2,i)
                          if( log_unit.ne.6 .and. log_unit.ne.0 ) 
     .                    write(log_unit,105) rn_codes(1,i),
     .                                        rn_codes(2,i)
 105                      format('**DUPLICATE** site names in',
     .                           ' rename from ', a,' to ',a)
                      end if

*                     Only continue if not a duplicate
                      if( .not.duplicate ) then
                        if( ns.gt.0 ) then
                            os = ltog_sites(j)
                            if( os.le.0 ) then
                                write(*,110) ns, os
  110                           format('ERROR in times_used, ns, os ',
     .                                  2i6)
                            end if
                            ltog_sites(j) = ns
                            times_used(os) = times_used(os) - 1
                            times_used(ns) = times_used(ns) + 1
* MOD TAH 000902: If we have not used the original site, set its
*                           apriori coordinates back to zero.
* MOD TAH 101228: Only reset the coordinates if the rename is not an
*                           edit.  _XCL _XPS edits on a site that 
*                           already has eq-renames in the h-file can
*                           reduce number to 0 and reset the coordinates                           
                            if( times_used(os).eq.0 .and. ud .and.
     .                          gsite_names(ns)(5:6).ne.'_X' ) then
                                do k = 1,3
                                   apr_val_site(k,1,os) = 0.d0
                                end do
                             end if  
                             call sbit(rn_name_changed,ns,1)
                        else if ( ns.lt.-1 ) then
 
*****                       Major problem.  We have a non unique name:
                            write(*,120) rn_codes(2,i), gnum_sites,
     .                          (gsite_names(k),k=1, gnum_sites)
                            stop 'GLINIT/EQ_NAME_CHN: Non unique names'
 
                        else
 
*****                       This the first time we have seen site so added
*                           new site name:
                            gnum_sites = gnum_sites + 1
                            ns = gnum_sites
                            gsite_names(ns) = rn_codes(2,i)
                            os = ltog_sites(j)
                                if( os.le.0 ) then
                                write(*,130) ns, os
  130                           format('ERROR: times_used, new ns, os ',
     .                                 2i6)
                            end if
                            ltog_sites(j) = ns
                            times_used(os) = times_used(os) - 1
                            times_used(ns) = times_used(ns) + 1
                            call sbit(rn_name_changed,ns,1)

*                           Copy station coordinate
* MOD TAH 000902: Added extra logic to saving apriori coordinates
*                           Changes are:
*                           (a) If apriori already exists, do not
*                               save a new one.
*                           (b) If this is the first time a site has
*                               been used, reset the original aprioris
*                               to zero so that they will be updated
                            do k = 1,3
                                if( apr_val_site(k,1,ns).eq.0 
     .                             .and. ud ) then
                                   apr_val_site(k,1,ns) =
     .                                     apr_val_site(k,1,os)
                                   apr_val_site(k,2,ns) = 
     .                                     apr_val_site(k,2,os)
                                end if
* MOD TAH 101228: Only reset the coordinates if the rename is not an
*                               edit.  _XCL _XPS edits on a site that 
*                               already has eq-renames in the h-file can
*                               reduce number to 0 and reset the coordinates                           
                                if( times_used(os).eq.0 .and.ud .and.
     .                              gsite_names(ns)(5:6).ne.'_X' ) then
                                   apr_val_site(k,1,os) = 0.d0
                                end if  
                            end do
* MOD TAH 070902: Copy zenith delay as well
                            apr_val_atm(ns) = apr_val_atm(os)
                        end if
*                     end of not a duplicate
                      end if  
                  end if
              end do
          end if
      end do
 
 
****  Start with earthquakes:  Loop over all the stations getting
*     the last earthquake to affect the station name
*     in list and see if we need doing anything
      do i = 1, cnum_sites
 
*         Now check the effect of each earthquake
          last_eq = 0
          last_eqep = 0.d0

* MOD TAH 080910: Make sure not one of the sites not to be renamed
          noeqrn = .false.
          do j = 1, num_norename
             if ( gsite_names(ltog_sites(i))(7:8).eq.
     .            norename(j)(1:2) ) noeqrn = .true.
          end do
 
* MOD TAH 090910: Only check EQ if we know needed
          if ( .not. noeqrn ) then 
            do j = 1, num_eq
                if( eq_rename(j) ) then
 
*                   We are renaming sites for this earthquake.  See
*                   if epoch falls in correct range
                    if( cepoch_start.ge.eq_epoch(j) ) then
 
*                       We are after the earthquake.  See if any of the
*                       sites in the current set need to have their
*                       names changed.
 
                        call eval_dist( eq_pos(1,j),
     .                      apr_val_site(1,1,ltog_sites(i)), dist)
 
***                     Now see if close enough.
* MOD TAH 011126: Check that the name is not _XCL before renaming
* MOD TAH 110122: Added XPS to exclude
* MOD TAH 100602: Only add name on first pass when ud is true.
                        if( dist.le.eq_rad(j) .and.
     .                      gsite_names(ltog_sites(i))(5:8).ne.
     .                                                '_XCL' .and.                         
     .                      gsite_names(ltog_sites(i))(5:8).ne.
     .                                                '_XPS' ) then                          
**    .                      '_XCL' .and. ud ) then 
 
*                           Save the new site code and eq. Number
*                           if this is later than any others.
*                           (We need to check other earthquakes in
*                           there is a latter on which will effect
*                           the station.
                            if( eq_epoch(j).gt.last_eqep ) then
                                last_eq = j
                                last_eqep = eq_epoch(j)
                            end if
                        end if
                    end if
                end if
            end do
          end if
 
****      If the last_eq value is non-zero then we have an earthquake
*         that we need to rename for.
          if( last_eq.gt.0 ) then

* MOD TAH 980514: Set bit to show eq has been used.
              call sbit(eq_used, last_eq, 1)
 
*             Generate the new station name
              new_name = gsite_names(ltog_sites(i))
              new_name(7:8) = eq_codes(last_eq)(1:2)
 
*             Make sure there are no blanks
              call sub_char(new_name,' ','_')
 
*             see if have already used this new name before
              indx = 1
              call get_cmd(new_name, gsite_names, gnum_sites,
     .                ns, indx)
 
*             If new site number is not zero, then we have used
*             before so save under new new name.
              if( ns.gt.0 ) then

*                 Decrement count on old site name
                  os = ltog_sites(i)
                  times_used(os) = times_used(os) - 1
*                 All we need do is change ltog_sites.
                  ltog_sites(i) = ns
                  times_used(ns) = times_used(ns) + 1

*                 Copy the rename forward to the post-earthquake
*                 site name
                  if( kbit(rn_name_changed,i) ) then
                      call sbit( rn_name_changed, ns,1)
                  end if

* MOD TAH 980329: Check to see if we should add entry to rename
*                 list to move site after an earthquake.
                  call inc_rename( ns, os )
              else if ( ns.lt.-1 ) then
 
*****             Major problem.  We have a non unique name:
                  write(*,120) new_name, gnum_sites,
     .                (gsite_names(j),j=1, gnum_sites)
 120              format('**DISASTER** Non-unqiue name found when ',
     .                    a,' added to ',i4,' existing sites',/,
     .                    'CURRENT SITE LIST:',200(/,8(a8,2x)))
                  stop 'GLINIT/EQ_NAME_CHANGE: Non unique names'
 
              else
 
*****             This the first time we have seen site so added
*                 new site name:
* MOD TAH 100616: Only add the new site name if this is the update
*                 called (ud == .true )
                  if( ud ) then 
                     gnum_sites = gnum_sites + 1
                     ns = gnum_sites
                     gsite_names(ns) = new_name
                     os = ltog_sites(i)
                     ltog_sites(i) = ns
                     times_used(os) = times_used(os) - 1
                     times_used(ns) = times_used(ns) + 1

*                    Copy the rename forward to the post-earthquake
*                    site name.
                     if( kbit(rn_name_changed,os) ) then
                         call sbit( rn_name_changed, ns,1)
                     end if

*                    Copy station coordinate
                     do j = 1,3
                         apr_val_site(j,1,ns) =
     .                           apr_val_site(j,1,os)
                         apr_val_site(j,2,ns) = 
     .                           apr_val_site(j,2,os)
                     end do
* MOD TAH 070902:    Copy zenith delay as well
                     apr_val_atm(ns) = apr_val_atm(os)

* MOD TAH 980329:    Check to see if we should add entry to rename
*                    list to move site after an earthquake.
                     call inc_rename( ns, os )
                  else
                     write(*,500) new_name, gnum_sites
 500                 format('**WARNING** Unexpected need for name of ',
     .                      a,' There are ',i6,' sites')
                  end if
 
              end if
          end if
      end do
C     print *,'Num Site ',gnum_sites
C     write(*,999) gnum_sites,
C    .       (j,gsite_names(j),times_used(j),j=1, gnum_sites)
C999  format('Next ',i4,/,(i4,2x,a8,2x,i4))
C     write(*,998) (i,ltog_sites(i),i=1,cnum_sites)
C998  format(5(i4,':'i5,' '))

****  Thats all
      return
      end 
 
