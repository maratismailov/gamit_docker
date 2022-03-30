CTITLE REMOVE_PARAMS
 
      subroutine remove_params( cov_obs, sol_obs, part_pnt, a_part)

      implicit none  
 
*     Routine to remove the parameters from COV_OBS etc. which are
*     not being estimated in the global solution.  If this is not
*     done, then these parameters are, in effect, set to zero.
*     (And option for the future if we wish to do such things)
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
      include '../includes/glb_hdr_def.h'
 
*   des         - Row number for destination position
*   i,j,k       - Loop counters
*   ir          - pointer to row to be compressed
*   part_pnt(2,max_glb_deriv,cnum_parn) - Pointers to which
*               - partials are being used.
*   zeros       - Number of zeros in the indx_pnt (i.e., number
*               - of not being estimated so far)
 
      integer*4 des, i,j, ir, part_pnt(2,max_glb_deriv,cnum_parn),
     .    zeros
 
*   a_part(max_glb_deriv,cnum_parn) - Compressed partial derivatives
*   cov_obs(cnum_parn,cnum_parn)    - Covariance matrix of input
*               - parameter estimates
*   sol_obs(cnum_parn)  - Solution vector
 
      real*8 a_part(max_glb_deriv,cnum_parn),
     .    cov_obs(cnum_parn,cnum_parn), sol_obs(cnum_parn)

* LOCAL Variables to implement adding noise
* gs -- Global site number
* len_ssh -- Length of hfile pattern string
* indx    -- Position in string
* ns      -- Return from trying to match global site names to
*            list of sites with wild cards (1 mean match found)
* lwc     -- Length of wild card string
* pwc     -- Position of wild card pattern.  Must be at end of
*            string.
* trimlen -- Length of string
      
      integer*4 gs, len_ssh, indx, ns, lwc, pwc, trimlen

* dmn     -- One minute in days.  Used to add a little tolerance
*            to time checks

      real*8 dmn

* use_ssh -- Set true of the hfile name matches the pattern
      logical use_ssh 

* wcard   -- Wild card pattern from ss_codes that start wiht @
      character*8 wcard

      data dmn / 0.0007d0 / 

* MOD TAH 020603: Added feature to weight stations by the number of times
*     they are used.
      uni_wght = .false. 
      if( uni_wght .and. gepoch_end-gepoch_start.le.1.d0 ) then
          call uni_weigh_sites(cov_obs, sol_obs)
      endif

* MOD TAH 000903: Added code to add noise to a site component based
*     on list given by user.
      
      do i = 1, num_ss

*        See if times and hf code match.
         use_ssh = .true.
         if( trimlen(ss_hfiles(i)).gt.0 ) then
             len_ssh = trimlen(ss_hfiles(i))
             indx = index(glb_inp_file,ss_hfiles(i)(1:len_ssh))
             if( indx.gt.0 ) then
                use_ssh = .true.
             else
                use_ssh = .false.
             end if
         endif
         if( cepoch_expt.ge.ss_times(1,i)-dmn .and.
     .       cepoch_expt.le.ss_times(2,i)+dmn .and. use_ssh ) then

*            OK: Times and hfile name match.  Now see if have any
*            sites.
             do gs = 1, gnum_sites

*               Does this sites name match the one for this sig_site.
*               (Use of get_cmd will allow short names for sites
*               and get multiple matchs)
                indx = 0
                call get_cmd(ss_codes(i),gsite_names(gs),1,ns,indx)
*               If we did not match, see if special code used for
*               this ss command (i.e,@_XXX)
                if( ss_codes(i)(1:1).eq.'@' ) then
                    wcard = ss_codes(i)(2:)
                    lwc = trimlen(wcard)
*                   See if the end of the string matches the wild card
                    if( lwc.gt.0 ) then
                        pwc = index(gsite_names(gs),wcard(1:lwc))
                        if( pwc.eq. 9-lwc) ns = 1
                    else
*                       Automatic match since only @ given
                        ns = 1
                    end if
                end if
*               See if we should continue
                if( ns.eq.1 ) then
*                   OK: We match on everything so report to user and
*                   add the noise
                    call add_ssnoise(gs, i, cov_obs)
                end if
             end do
         end if
      end do   
       
***** Look for all parameters for which the number of partials is
*     zero i.e., they are not used in the solution
 
*                             ! Set to all used
      cnum_used = cnum_parn
*                             ! initialize the row number to be compressed
      ir        = 0
*                             ! initialize destination row number
      des       = 0
*                             ! initialize the number of deleted parameters
      zeros     = 0
 
      do i = 1, cnum_parn
          ir = ir + 1
*                                         ! Parameter not used
          if( indx_pnt(ir).eq.0 ) then
              zeros = zeros + 1
              ! print *,'NON-USED PARAM ',i,cnum_parn, zeros
*                                         ! parameter found compress
          else
*                                         ! any earlier parameters
*                                         ! indcrement destiantin
              des = des + 1
*                                         ! we need to compress
              if( des.lt.ir ) then
 
*                 Compress the indx_pnt and parameter codes
                  do j = ir, cnum_used
                      indx_pnt  (des+j-ir) = indx_pnt  (j)
                      gpar_codes(des+j-ir) = gpar_codes(j)
                      sol_obs   (des+j-ir) = sol_obs   (j)
 
*                     Now used VIS to compress matrices (move V1 to V2)
                      call DWMOV(a_part(1,       j),1,
     .                           a_part(1,des+j-ir),1, max_glb_deriv)
 
*                     Since part_pnt is in pairs move as a R*8 vector
*                     (pairs of I*4 values)
                      call DWMOV(part_pnt(1,1, j),1,
     .                     part_pnt(1,1,des+j-ir),1, max_glb_deriv)
 
*                     Move column of cov_obs, followed by row
                      call DWMOV(cov_obs(1,       j),1,
     .                           cov_obs(1,des+j-ir),1,cnum_parn)
                      call DWMOV(cov_obs(j ,1),cnum_parn,
     .                     cov_obs(des+j-ir,1),cnum_parn,cnum_parn)
*                             ! Compressing
                  end do
 
*                 Get new size for compressed matrix
                  cnum_used = cnum_parn - zeros
*                             ! new destination
                  ir = des
*                             ! we needed to compress
              end if
*                             ! Current parameter used
          end if
*                             ! Looping over all parameters
      end do
 
***** Thats all, make sure we have the new number of parameters
      cnum_used = cnum_parn - zeros
      return
      end

CTITLE ADD_SSNOISE

      subroutine add_ssnoise(gs, i, cov_obs)

      implicit none 

*     Routine to add additional noise onto the data covariance
*     matrix based on user input using the sig_neu command

      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
      include '../includes/glb_hdr_def.h'

* PASSED VARIABLES
* gs -- Global site number of site to be changed
* i  -- Number of site sig entry to add

      integer*4 gs, i

* cov_obs(cnum_parn,cnum_parn) -- data covariance matrix

      real*8 cov_obs(cnum_parn,cnum_parn)

* LOCAL VARIABLES
* j,k,l -- Loop countes
* ls    -- Local site site number
* nl(3) -- pointers to XYZ positions in cov_obs
* lh    -- Length of hfile name for output
* trimlen -- Length of string
* type, indx -- Parameter type and index number (station
*               number for types 7,8,9 (XYZ position)
      integer*4 j,k,l, ls, nl(3), lh, trimlen, type, indx

* rot_mat(3,3) -- Rotation matrix from NEU to XYZ
* loc_coord(3) -- Local coorinates
* cov_neu(3,3), cov_xyz(3,3), temp_cov(3,3) -- Covariance
*     matrices for NEU, XYZ and temp version for comps.

      real*8 rot_mat(3,3), loc_coord(3), 
     .       cov_neu(3,3), cov_xyz(3,3), temp_cov(3,3) 

****  Find the local site (basically, when called we know
*     that if the site is observed, it should be down weighted)
      ls = -1
      do j = 1, cnum_sites
         if( ltog_sites(j).eq.gs ) then
             ls = j
         end if
      end do

****  If we found the site locally, continue
      if( ls.gt.0 ) then

****      Tell user what we are doing (at least for the moment)
          lh = trimlen(glb_inp_file)
          write(*,110) gsite_names(gs), ss_codes(i), 
     .       (ss_sig(j,i)*1000,j=1,3),  glb_inp_file(1:lh)
          if( log_unit.ne.6 ) 
     .    write(log_unit,110) gsite_names(gs), ss_codes(i), 
     .        (ss_sig(j,i)*1000,j=1,3),  glb_inp_file(1:lh)
 110      format('SS_SIG Applied to ',a8,' Code ',a8,' Sig NEU ',
     .           3(1x,F10.1),' mm. Hfile ',a)

****      OK: Based on coordinates of site get the NEU to XYZ
*         transformation
* MOD TAH 030116: Changed call to NEU from GEOD 
          call XYZ_to_NEU( rot_mat, apr_val_site(1,1,gs),
     .                      loc_coord)

*         Transpose the matrix
          do k = 1,2
              call dvswp(rot_mat(k,k+1),3,
     .                   rot_mat(k+1,k),1, 3-k)
          end do

*         Set up the diagonal NEU matrix
          do k = 1,3
             do l = 1,3
                 cov_neu(k,l) = 0.d0
             end do
             cov_neu(k,k) = ss_sig(k,i)**2
          end do

          call var_comp(rot_mat, cov_neu, cov_xyz,
     .                  temp_cov, 3,3,1 )

*****     OK, Now we need to add to the correct places in
*         cov_obs
          do j = 1,3
             nl(j) = 0
          end do
          do j = 1, cnum_parn
             call decode_code(gpar_codes(j),type,indx)
             if( indx.eq.ls .and. type.ge.7 .and.
     .           type.le.9 ) nl(type-6) = j
          end do

****      OK: Now add extra noise
          do j = 1, 3
             do k = 1,3
                if( nl(j).gt.0 .and. nl(k).gt.0 ) then
                    cov_obs(nl(j),nl(k)) = 
     .                   cov_obs(nl(j),nl(k)) + cov_xyz(j,k)
                end if
             end do
          end do
      end if

***** Thats all 
      return
      end

CTITLE UNI_WEIGH_SITES
 
      subroutine uni_weigh_sites( cov_obs, sol_obs )

      implicit none  
 
*     Routine to re-scale the site coordinates assoicated with sites
*     used multiple times in a network analysis
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
      include '../includes/glb_hdr_def.h'
 

* PASSED VARIABLES
* cov_obs(cnum_parn,cnum_parn) -- Data covariance matrix

      real*8 cov_obs(cnum_parn,cnum_parn)
     .,      sol_obs(cnum_parn) 

* LOCAL VARIABLES 
      real*8 wght   ! Weight for station (multiplier on variance)
      real*8 scale(max_glb_parn)  ! Scale for inversion

      integer*4 cs  ! Loop over local stations
     .,         gs  ! Corresponding global site number
     .,         i,j,k  ! Loop counter
     .,         type, indx  ! Type and index of local parameter types
     .,         ipivot(max_glb_parn)   ! Pivot elements


****  Loop over all the local sites, checking how many times 
*     they have been used and re-scaling site covariances 
*     appropriately 
                    

***   CODE DOES NOT WORK:
      uni_wght = .false.
c      RETURN

c      call invert_vis(cov_obs,sol_obs, scale, ipivot,
c     .                cnum_parn,cnum_parn,1)
      
      do cs = 1, cnum_sites

*        See if used more than once
         gs = ltog_sites(cs)
         if( times_used(gs)-1.gt.1 ) then

*            Find all parameters associated with site and
*            rescale the variances
             wght = 1.d0/dble(times_used(gs))
             do i = 1, cnum_parn
                call decode_code(gpar_codes(i),type,indx)
                if( indx.eq.cs .and. type.eq.7  ) then
                    do j = 0, 2
                       do k = 0, 2
                          cov_obs(i+j,i+k) = cov_obs(i+j,i+k)/wght
                       end do
                    end do
                end if
             end do
         end if
      end do

c      call invert_vis(cov_obs,sol_obs, scale, ipivot,
c    .                cnum_parn,cnum_parn,1)

****  Thats all
      return
      end
 

           
