CTITLE FIX_ORIGIN
 
      subroutine fix_origin( site_pos_main, site_pos_ema, dims, type,
     .    translation, sol_vec, parn, dimp, use_site, num_sites,
     .    work, ltog)
 

      implicit none 
 
*     Routine to compute the "optiminal" site coordinate translation
*     to be used so that the horizonal adjustment to the site
*     positions is minimized.  The sites to be used are given by
*     the bits in USE_SITE.  The routine may be passed site poisitions
*     in either main memory or ema.  The type passed is given by the
*     TYPE variable
 
 
*   dimp        - Second dimension of the PARN array which gives
*               - the parameter numbers of the site position
*               - parameters.  The first dimension is assumed to be
*               - three (for XYZ).  Dimp indicates if rates are given
*   dims        - Similar to dimp, but for the site positions and
*               - rates which are passed.
 
*   i,j,k       - Loop counters
*   is          - Correct site number to use for the coordinates of
*               - the site. (Depends on LTOG values)
*   ltog(1)     - Conversion from local site number to global site
*               - number.  If first entry is -1 then local and global
*               - site numbers are the same.
*   ltog_map    - integer*4 function which converts local to global
*               - site number.
 
*   num_iterations  - Number of iterations, limited to 200 to get
*               - convergence.
 
*   num_sites   - Total number of sites in this solution
*   num_used    - Number of sites actually used
*   parn(3,dimp,1)  - Array giving the parameter numbers of the site
*               - positions or rates to have their origin fixed.
*   use_site(1) - Bit mapped words with bits set if this site is to
*               - be used.
 
      integer*4 dimp, dims, i,j,k, is, ltog(1), ltog_map,
     .    num_iterations, num_sites, num_used, parn(3,dimp,1),
     .    use_site(1)
 
*   beta        - intermediate value used in computing the new
*               - origin values.
*   chi(-1:1,3) - Chi**2 of the horizontal displacements for choices
*               - of the origin
*   corr(3)     - Corrections to current origin.  Used to check
*               - convergence
*   ds1, ds2    - Changein in Chi**2 bewteen origin trials
*   horizontal  - Horizonal displacement at a site
*   org(3)      - Working values of the translations
*   pos(3)      - Working values for site adjustment after origin
*               - removed
*   radius      - Radius for current site
*   radial      - Radial displacement of the current site
*   site_pos_ema(3,dims,1)  - Site positions stored in ema
*   site_pos_main(3,dims,1) - Site positions stored in main memory
*   step        - Step to be used in testing trial values
*   total           - Total motion of the site
*   translation(3)  - Values of the three translations which should
*                   -be applied to all site position adjustments
*   sol_vec(1)      - solution vector for all parameters estmated
*   work(3,2,1)     - Partials of site position with respect to
*                   - radial displacement and the adjustment to the
*                   - site positions (second index) (STORED IN EMA)
      real*8 beta, chi(-1:1,3), corr(3), ds1, ds2, horizontal, org(3),
     .    pos(3), radius, radial, site_pos_ema(3,dims,1),
     .    site_pos_main(3,dims,1), step, total, translation(3),
     .    sol_vec(1), work(3,2,1)
 
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
 
      do i = 1,3
          translation(i) = 0.d0
      end do
 
      step = 0.1d0
      num_iterations = 0
 
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
*                         ! Site used
          end if
 
*                         ! Looping over the sites
      end do
 
*                         ! Save number of sites
      num_used = k
 
*     Now start iterating to get the site position translations
 
      converged = .false.
      do while ( .not. converged )
 
*         Loop over XYZ generating trial values
          do i = 1,3
 
*             Copy current translation into working origin
              do j = 1,3
                  org(j) = translation(j)
              end do
 
*             Now compute trial chi**2 at +-1 step about current origin
              do j = -1, 1
                  org(i) = translation(i) + step*j
                  chi(j,i) = 0.d0
 
*                 Now sum up horizontal displacements for this choice
*                 of origin
                  do k = 1, num_used
 
                      pos(1) = work(1,2,k) - org(1)
                      pos(2) = work(2,2,k) - org(2)
                      pos(3) = work(3,2,k) - org(3)
 
                      radial = work(1,1,k)*pos(1) + work(2,1,k)*pos(2) +
     .                         work(3,1,k)*pos(3)
 
                      total =  pos(1)**2 + pos(2)**2 + pos(3)**2
                      horizontal = abs(total - radial**2)
                      chi(j,i) = chi(j,i) + horizontal
                  end do
*                             ! Looping over -1,0,+1 steps
              end do
*                             ! Looping over XYZ origin.
          end do
 
*         Now get the new origin estimate
          converged = .true.
*                             ! Loop over XYZ
          do i = 1,3
              ds1 = chi(-1,i) - chi(0,i)
              ds2 = chi( 1,i) - chi(0,i)
              beta = ds1/ds2
              corr(i) = (beta-1)/(2*(beta+1)) * step
 
*             Check convergence
              if( abs(corr(i)).gt.1.d-5 ) converged = .false.
              translation(i) = translation(i) + corr(i)
          end do
 
*         Check number of iterations to get convergence
          num_iterations = num_iterations + 1
          if( num_iterations.gt.200 ) then
              write(*,100) num_iterations
 100          format(' ** WARNING ** Fix_origin failed to converge',
     .               ' after ',i3,' iterations')
              converged = .true.
          end if
 
*                             ! Looping until we are convered
      end do
 
***** Thats all
      return
      end
 
