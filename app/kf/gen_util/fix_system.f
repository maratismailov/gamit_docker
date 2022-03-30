CTITLE FIX_SYSTEM
 
      subroutine fix_system( site_pos_main, site_pos_ema, dims, type,
     .    translation, sol_vec, parn, dimp, use_site, num_sites,
     .    work, ltog, pmu_parts )
 

      implicit none 
 
*     Routine to compute the "optiminal" site coordinate translation
*     and rotation to be used so that the horizonal adjustment to the site
*     positions is minimized.  The sites to be used are given by
*     the bits in USE_SITE.  The routine may be passed site poisitions
*     in either main memory or ema.  The type passed is given by the
*     TYPE variable. If PMU_parts is passed with zeros then only the
*     translation will be determined.
 
*   dimp        - Second dimension of the PARN array which gives
*               - the parameter numbers of the site position
*               - parameters.  The first dimension is assumed to be
*               - three (for XYZ).  Dimp indicates if rates are given
*   dims        - Similar to dimp, but for the site positions and
*               - rates which are passed.
 
*   i,j,k,l     - Loop counters
*   is          - Correct site number to use for the coordinates of
*               - the site. (Depends on LTOG values)
*   iuse(6)     - Pivot rows from SMNV8
*   JEL         - Function to return index in lower triangular
*               - matrix.
*   ltog(1)     - Conversion from local site number to global site
*               - number.  If first entry is -1 then local and global
*               - site numbers are the same.
*   ltog_map    - integer*4 function which converts local to global
*               - site number.
 
*   num_iterations  - Number of iterations, limited to 200 to get
*               - convergence.
 
*   num_sites   - Total number of sites in this solution
*   num_terms   - Either 3 if only translations are done, or
*               - 6 if translations and rotations are done.
*   num_used    - Number of sites actually used
*   parn(3,dimp,1)  - Array giving the parameter numbers of the site
*               - positions or rates to have their origin fixed.
*   use_site(1) - Bit mapped words with bits set if this site is to
*               - be used.
 
      integer*4 dimp, dims, i,j,k,l, is, iuse(6), JEL, ltog(1),
     .    ltog_map, num_iterations, num_sites, num_terms, num_used,
     .    parn(3,dimp,1), use_site(1)
 
*   a(21),b(6)  - Normal equations and solution vector for solving
*               - quadatic form. (Normal equations are lower diagonal)
*   beta        - intermediate value used in computing the new
*               - origin values.
*   C0          - Chi**2 value at current origin
*   chi(-1:1,6) - Chi**2 of the horizontal displacements for choices
*               - of the origin and rotation
*   corr(6)     - Corrections to current origin and translation.
*               - Used to check convergence.
*   corr_scale  - Scaling factor for corrections. Used too speed
*               - up convergence.
*   ds1, ds2    - Changein in Chi**2 bewteen origin trials
*   horizontal  - Horizonal displacement at a site
*   org(6)      - Working values of the translations and rotations.
*   pos(3)      - Working values for site adjustment after origin
*               - removed
*   pmu_parts(3,3,1)    - Polar motion/UT1 partials for each site
*               - NOTE: This array is changed if not all sites are
*               - used.
*   prev_corr(6)    - Previous correction.  Used to speed up
*               - convergence by scaling correction
*   pvrow(6), pvcol(6), pvrwb   - Scratch arrays used by SMNV8
*   radius      - Radius for current site
*   radial      - Radial displacement of the current site
*   site_pos_ema(3,dims,1)  - Site positions stored in ema
*   site_pos_main(3,dims,1) - Site positions stored in main memory
*   step        - Step to be used in testing trial values
*   total           - Total motion of the site
*   translation(6)  - Values of the three translations and 3 rotations
*                   - which should be applied to all site position
*                   - site position adjustments
*   sol_vec(1)      - solution vector for all parameters estmated
*   work(3,2,1)     - Partials of site position with respect to
*                   - radial displacement and the adjustment to the
*                   - site positions (second index) (STORED IN EMA)
      real*8 a(21),b(6), C0, chi(-1:1,6), corr(6), corr_scale,
     .    org(6), pmu_parts(3,3,1), pvrow(6), pvcol(6), pvrwb, radius,
     .    site_pos_ema(3,dims,1), site_pos_main(3,dims,1), step,
     .    translation(6), sol_vec(1), work(3,2,1)
 
*   all_est         - Indicates that all position components have
*                   - been estimated. If they have not been estimated
*                   - then routine returns with zeros for the
*                   - translations
*   converged       - Indicates values have converged (site origin
*                   - Changes of < 0.0001 m or m/yr)
*   kbit            - Function for checking if a bit is set
 
      logical all_est, converged, kbit
 
*   type            - Indicates type of storage for site positions
*                   - Either 'EMA' or 'MAIN' memory storage.
 
      character*4 type
 
 
****  Start, check that all selected sites have been estimated.
 
      step = 0.1
      num_iterations = 0
*                         ! Set for doing only translation. If pmu_parts is
      num_terms = 3
*                         ! Non-zero then value will be changed to 6
*                                                ! This way a single 0 can
      if( pmu_parts(1,1,1).ne.0 ) num_terms = 6
*                         ! by passed
 
*                 ! Counter for number of site used
      k = 0
 
      do i = 1, num_sites
 
*                                         ! We are to use
          if( kbit( use_site,i ) ) then
 
*             Check that all components of site position have been
*             estimated
              all_est = .true.
*                                          ! Loop over XYZ
              do j = 1,3
                  if( parn(j,1,i).eq.0 ) then
                      all_est = .false.
                  end if
              end do
 
              k = k + 1
 
*                                         ! Cant fix origin since all sites
              if( .not.all_est ) RETURN
*                                         ! not estimated.
 
*             Now get radial partials and save
*                                              ! Ema site position
              IF( type(1:3).eq.'EMA' ) THEN
                  is = ltog_map( ltog, i)
                  radius = sqrt( site_pos_ema(1,1,is)**2 +
     .                           site_pos_ema(2,1,is)**2 +
     .                           site_pos_ema(3,1,is)**2 )
                  do j = 1, 3
                      work(j,1,k) = site_pos_ema(j,1,is)/radius
                      work(j,2,k) = sol_vec( parn(j,1,i) )
                  end do
*                                             ! Main memory
              ELSE
                  is = ltog_map( ltog, i)
                  radius = sqrt( site_pos_main(1,1,is)**2 +
     .                           site_pos_main(2,1,is)**2 +
     .                           site_pos_main(3,1,is)**2 )
                  do j = 1, 3
                      work(j,1,k) = site_pos_main(j,1,is)/radius
                      work(j,2,k) = sol_vec( parn(j,1,i) )
                  end do
              END IF
 
*             Now move the PMU partials
              if( k.ne.i .and. num_terms.eq.6 ) then
                  do j = 1, 3
                      do l = 1,3
                          pmu_parts(j,l,k) = pmu_parts(j,l,i)
                      end do
                  end do
              end if
 
*                         ! Site used
          end if
 
*                         ! Looping over the sites
      end do
 
*                         ! Save number of sites
      num_used = k
 
*     Now start iterating to get the site position translations
 
      converged = .false.
*                         ! Use half step at first to stop growing
      corr_scale = 0.5
*                         ! oscillations
*     Clear the translation and rotations
      do i = 1, num_terms
          translation(i) = 0.0d0
      end do
 
C     do while ( .not. converged )
C
*         Loop over XYZ and PMU generating trial values
C         do i = 1,num_terms
C
*             Copy current translation into working origin
C             do j = 1,num_terms
C                 org(j) = translation(j)
C             end do
C
*             Now compute trial chi**2 at +-1 step about current origin
C             do j = -1, 1
C                 org(i) = translation(i) + step*j
C                 chi(j,i) = 0.d0
C
*                 Now sum up horizontal displacements for this choice
*                 of origin
C                 do k = 1, num_used
C
C                     pos(1) = work(1,2,k) - org(1)
C                     pos(2) = work(2,2,k) - org(2)
C                     pos(3) = work(3,2,k) - org(3)
C
*                     Now do the rotation
C                     do l = 1,3      ! Loop over x,y and UT1 values
C                         pos(1) = pos(1) - pmu_parts(l,1,k)*org(3+l)
C                         pos(2) = pos(2) - pmu_parts(l,2,k)*org(3+l)
C                         pos(3) = pos(3) - pmu_parts(l,3,k)*org(3+l)
C                     end do
C
C                     radial = work(1,1,k)*pos(1) + work(2,1,k)*pos(2) +
C    .                         work(3,1,k)*pos(3)
C
C                     total =  pos(1)**2 + pos(2)**2 + pos(3)**2
C                     horizontal = abs(total - radial**2)
C                     chi(j,i) = chi(j,i) + horizontal
C                 end do
C             end do          ! Looping over -1,0,+1 steps
C         end do              ! Looping over XYZ origin.
C
*         Now get the new origin estimate
C         converged = .true.
C
*         See if we should try to speed convergence
C         corr_scale = 0.5
C         if( num_iterations.gt.25 ) then
C             corr_scale = 0
C             do i = 1,num_terms
C                 corr_scale = corr_scale +
C    .                         abs(corr(i)/(prev_corr(i)-corr(i)))
C             end do
C             corr_scale = corr_scale/num_terms
C             if( corr_scale.gt. 50 ) then
C                 corr_scale =  50
C             end if
C             if( corr_scale.lt.1 ) corr_scale = 0.5 ! Make sure we dont
C                                                    ! oscillate
C         end if
C
*         Now compute corrections
C
C         do i = 1,num_terms  ! Loop over XYZ and xy and UT1
C             ds1 = chi(-1,i) - chi(0,i)
C             ds2 = chi( 1,i) - chi(0,i)
C             beta = ds1/ds2
C
C             prev_corr(i) = corr(i)
C             corr(i) = corr_scale * (beta-1)/(2*(beta+1)) * step
C
*             Check convergence
C             if( abs(corr(i)).gt.1.d-5 ) converged = .false.
C             translation(i) = translation(i) + corr(i)
C         end do
C
* DEBUG
cd        write(1,999,iostat=l) num_iterations, translation, corr
cd  999     format(' Iter ',i3,' Values ',3(f12.7,1x),/,
cd     .                    t18,          3(f12.7,1x),/,
cd     .                    t9,' Corrs  ',3(f12.7,1x),/,
cd     .                    t18,          3(f12.7,1x)    )

cd         write(1,998) sqrt(chi(0,1)/(num_used-3)), corr_scale
cd  998     format(' Horizontal RMS ',f12.7,' corr_scale ',f12.7)

*         Check number of iterations to get convergence
C         num_iterations = num_iterations + 1
C         if( num_iterations.gt. 10 ) then
C             write(1,100) num_iterations
C100          format(' ** WARNING ** Fix_origin failed to converge',
C    .               ' after ',i3,' iterations')
C             converged = .true.
C         end if
C
C     end do                  ! Looping until we are convered
 
***** Now go to solution for quadtratic form and see what happens
 
      converged = .false.
      do while ( .not. converged )
 
*         Get chi for current origin
          call comp_hchi( C0, work, translation, pmu_parts, num_used)
 
*         Now do chi's for all of the unit vectors and the cross terms
*                             ! Loop over XYZ and xy UT1 values
          do i = 1, num_terms
 
*             Do main chi unit vector
 
              call tran_to_org( translation, org, i, -step)
              call comp_hchi( chi(-1,i), work, org, pmu_parts, num_used)
 
              call tran_to_org( translation, org, i, +step)
              call comp_hchi( a(jel(i,i)), work, org, pmu_parts,
     .                        num_used)
 
*                                 ! Do the cross vector
              do j = 1, i-1
 
                  call tran_to_org( translation, org, i, +step)
                  call tran_to_org(         org, org, j, +step)
 
                  call comp_hchi( a(jel(i,j)), work, org, pmu_parts,
     .                            num_used )
              end do
          end do
 
****      Now construct the normal equations and b vector for solution
          do i = 1, num_terms
              b(i)        = ( a(jel(i,i)) - chi(-1,i) )/2
              a(jel(i,i)) = ( a(jel(i,i)) + chi(-1,i) )/2 - C0
 
              do j = 1, i-1
                  a(jel(i,j)) = ( a(jel(i,j)) - a(jel(i,i)) -
     .                            a(jel(j,j)) - b(i) - b(j) - C0 )/2
              end do
          end do
 
*****     Now get the solution
          call smnv8(a,b, num_terms, 1, pvrow, pvcol, iuse, pvrwb)
 
*         Now get corrections
          converged = .true.
          num_iterations = num_iterations + 1
          do i = 1, num_terms
              corr(i) = -b(i)*step*corr_scale
              translation(i) = translation(i) + corr(i)
              if( abs(corr(i)).gt.1.d-5 ) converged = .false.
          end do
 
          if( num_iterations.gt. 200 ) converged = .true.
 
*                 ! Looping until convereged
      end do
 
C     write(1,999, iostat=l) num_iterations, translation, corr
C     write(1,998) sqrt(C0/num_used), corr_scale
 
      work(1,1,1) = sqrt(C0/num_used)
 
***** Thats all
      return
      end
 
