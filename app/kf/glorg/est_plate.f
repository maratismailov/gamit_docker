CTITLE EST_PLATE      

      subroutine est_plate(iout, cov_parm, sol_parm, cov_obs,
     .                     sol_obs, plate_parts, atv) 

      implicit none 

*     Routine to estimate plate rotation vector based on velocities
*     of sites selected with the plate command.

      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include '../includes/globk_common.h'
      include '../includes/glorg_common.h'

*   iout                - Output unit
 
      integer*4 iout
 
*   cov_parm(num_glb_parn,num_glb_parn) - Covariance matrix.
*   sol_parm(num_glb_parn,numc) - Solution vector
*   cov_obs(3*num_plate_sites, 3*num_plate_sites)  - Covariances
*             for the sites used in the plate vector solution
*   sol_obs(3*num_plate_sites) - Solution vector for velocities
*             of sites.
*   plate_parts(3*num_plate_sites,3*num_plates)  - Partials of
*             site velocties with respect to the plates.
*   atv(3*num_plate_sites,3*num_plates) - ATV-1 matrix (to speed
*             computation
 
      real*8  cov_parm(num_glb_parn,num_glb_parn), 
     .        sol_parm(num_glb_parn),
     .        cov_obs(3*num_plate_sites, 3*num_plate_sites) ,
     .        sol_obs(3*num_plate_sites),
* MOD TAH 030513: Added three parameters for translation
     .        plate_parts(3*num_plate_sites,3*num_plates+3),
     .        atv(3*num_plate_sites,3*num_plates+3)
 
* LOCAL VARIABLES
 
*   i,j,k,l         - Loop counters
*   np, mp          - Parameter numbers for estimated sites
*   nu, mu          - Incremental counters for the sites used
*                     in the plate solution
*   pn, pm          - Plate numbers for sites
*   ns, ps, ms      - Parameter numbers in the plate solution
*                     for sites, plate and second site
*   ipivot(3*max_plates) - Pivot elements for inversion

*   ntran            - Number of translation parameters
 
      integer*4 i,j,k,l,ipivot(3*max_glb_sites), np, mp, nu, mu,
     .          pn, pm, ns, ps, ms, ntran
 
*   scale(3*max_plates)            - Used by invert_vis
*   norm_eq(3*max_plates,3*max_plates) - Normal eqautions
*   bvec(3*max_plates)             - Solution vectors
*   res                            - result from dot product
*   assign_part(3,3)  -- Partials of XYZ with respect to the
*     pole of rotaton.
 
      real*8 scale(3*max_glb_sites), 
     .       norm_eq(3*max_plates,3*max_plates),
     .       bvec(3*max_plates), res,
     .       assign_part(3,3)

***** See if we are estimating translation with the plate poles
      if( PlateTrans ) then
          ntran = 3
      else
          ntran = 0
      end if

***** First loop over the sites and find which plate they
*     belong to (if any) and form the partials and covariance
*     matrix
      nu = 0
      do i = 1, gnum_sites
         pn = plate_number(i)
         np = parn_site(1,2,i)
         if( pn.gt.0 .and. np.gt.0 ) then 

****         OK found a site with a plate.  Clear
*            partial array and form partials for this site.
*            Increment site number used for plates.
             nu = nu + 1
             ns = 3*(nu-1)
             do j = 1, num_plates+1   ! Add one for translation
                 do k = 1,3
                    do l = 1,3
                       plate_parts(ns+l,3*(j-1)+k) = 0.d0
                    end do
                 end do
             end do

*****        Now form the partials
             ps = 3*(pn-1)
             plate_parts(ns+1,ps+2) =  apr_val_site(3,1,i)*1.d-6
             plate_parts(ns+1,ps+3) = -apr_val_site(2,1,i)*1.d-6

             plate_parts(ns+2,ps+1) = -apr_val_site(3,1,i)*1.d-6
             plate_parts(ns+2,ps+3) =  apr_val_site(1,1,i)*1.d-6

             plate_parts(ns+3,ps+1) =  apr_val_site(2,1,i)*1.d-6
             plate_parts(ns+3,ps+2) = -apr_val_site(1,1,i)*1.d-6

* MOD TAH 030513: Add partials for translation
             if( PlateTrans ) then
                plate_parts(ns+1,3*num_plates+1) = 1.d0
                plate_parts(ns+2,3*num_plates+2) = 1.d0
                plate_parts(ns+3,3*num_plates+3) = 1.d0
             end if
 

*****        Put the velcoity estimates into the solution vector
             do j = 1,3
                sol_obs(ns+j) = apr_val_site(j,2,i) +
     .                          sol_parm(np+j-1)
             end do

*****        Now move the covariance elements into cov_obs
             mu = 0
             do j = 1, gnum_sites
                pm = plate_number(j)
                mp = parn_site(1,2,j)
                if( pm.gt.0 .and. mp.gt.0 ) then 
                    mu = mu + 1
                    ms = 3*(mu-1)
                    call mov_cov(parn_site(1,2,i), parn_site(1,2,j),3,
     .                     cov_obs(ns+1,ms+1), 3*num_plate_sites,
     .                     cov_parm, num_glb_parn)
                end if
             end do
          end if
      end do

***** We have partials, covariance matrix and solution vector.  Now
*     invert cov_obs and form LSQ estimator

      call invert_vis(cov_obs, sol_obs,scale,ipivot,3*num_plate_sites,
     .                3*num_plate_sites, 0 )

*               T -1       T -1
****  Now form A V  A and A V  f
*
*           T -1
*     Form A V   first
      do i = 1, 3*num_plates+ntran
         do j = 1, 3*num_plate_sites
            call dwdot(res, plate_parts(1,i),1, cov_obs(1,j),1,
     .                3*num_plate_sites)
            atv(j,i) = res 
         end do
      end do

****  Now form the normal equations and the solution vector
      do i = 1, 3*num_plates+3
         call dwdot(res, atv(1,i),1, sol_obs(1),1, 3*num_plate_sites)
         bvec(i) = res
         do j = 1, 3*num_plates+ntran
            call dwdot(res, atv(1,i),1, plate_parts(1,j),1, 
     .                 3*num_plate_sites)
            norm_eq(i,j) = res
         end do
      end do

****  Now get the solution
      call invert_vis(norm_eq, bvec,scale, ipivot, 3*num_plates+ntran,
     .                3*max_plates, 1 )

****  Now write out the solution

      call write_plate(iout, norm_eq, bvec )

****  Now compute the residuals and save in sol_parm.  (We also
*     update the apriori values for the rates)
      nu = 0
      do i = 1, gnum_sites
         pn = plate_number(i)
         np = parn_site(1,2,i)
         if( pn.gt.0 .and. np.gt.0 ) then 
             ps = 3*(pn-1)
             nu = nu + 1
             ns = 3*(nu-1)
             do j = 1,3
                call dwdot(res, 
     .                     plate_parts(ns+j,ps+1),3*num_plate_sites,
     .                     bvec(ps+1),1, 3)
                sol_parm(np+j-1) = apr_val_site(j,2,i) + 
     .                             sol_parm(np+j-1) - res
                apr_val_site(j,2,i) = res
             end do
         else
* MOD TAH 980611:
*            See if a site has been assigned to a plate
             pn = assign_number(i)
             np = parn_site(1,2,i)
             if( pn.gt.0 .and. np.gt.0 ) then 
                 ps = 3*(pn-1)
                 assign_part(1,1) = 0.d0
                 assign_part(1,2) =  apr_val_site(3,1,i)*1.d-6
                 assign_part(1,3) = -apr_val_site(2,1,i)*1.d-6
      
                 assign_part(2,1) = -apr_val_site(3,1,i)*1.d-6
                 assign_part(2,2) = 0.d0
                 assign_part(2,3) =  apr_val_site(1,1,i)*1.d-6
       
                 assign_part(3,1) =  apr_val_site(2,1,i)*1.d-6
                 assign_part(3,2) = -apr_val_site(1,1,i)*1.d-6
                 assign_part(3,3) = 0.d0

                 do j = 1,3
                    call dwdot(res, 
     .                         assign_part(j,1),3,
     .                         bvec(ps+1),1, 3)
                    sol_parm(np+j-1) = apr_val_site(j,2,i) + 
     .                                 sol_parm(np+j-1) - res
                    apr_val_site(j,2,i) = res
                 end do
              end if
         end if
      end do

* MOD TAH 030513: Remove the translation from all of the sites
      if( PlateTrans ) then
         do i = 1, gnum_sites
            do j = 1,3 
               np = parn_site(j,2,i)
               if(np.gt.0 ) then 
                   sol_parm(np) = sol_parm(np) - bvec(3*num_plates+j)
               end if
            end do
         end do
*        Add the translation to any estimate of the translation
         do j = 1,3
            np = parn_tran(j,2)
            if( np.gt.0 ) then
                sol_parm(np) = sol_parm(np) + bvec(3*num_plates+j)
            end if
         end do
      end if

         

****  Thats all
      return
      end

CTITLE WRITE_PLATE

      subroutine write_plate( iout, norm_eq, bvec)

      implicit none 

*     Routine to write out the estimated plate rotation vectors.

      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include '../includes/globk_common.h'
      include '../includes/glorg_common.h'

*   iout                - Output unit
 
      integer*4 iout

*   norm_eq(3*max_plates,3*max_plates) - Normal eqautions
*   bvec(3*max_plates)             - Solution vectors
 
      real*8 norm_eq(3*max_plates,3*max_plates),
     .       bvec(3*max_plates)

* LOCAL VARIABLES

* i,j   - Loop counters
* pn    - Plate number for site i
* nu    - Counter for number of sites for plate i modulo
*         number of sites names written to line
* is,js - Start elements for plate i and j
* ip(3), jp(3) - Parameter numbers to be used by mov_cov

      integer*4 i,j, k, pn, nu, is,js, ip(3), jp(3) 

* rot_mat(3,3)  - Rotation matrix from XYZ to NEU
* trans(3,6)    - Transformation from two Euler poles to
*                 difference (XYZ) and NEU of difference
* loc_coord(3)  - Colat, Long and length of vector
* xyz_adj(3)    - Dummy XYZ adjust (for call compatability)
* neu_adj(3)    - Dummy NEU adjust (for call compatability)

* cov_pole(6,6) - Covariance matrix for pair of Euler poles
* cov_sing(3,3) - covariance matrix for single Euler pole.
* temp_cov(6,6) - Working storage.

* est(3)        - Estimate for output
* sig(3)        - Sigma's for output
* rho(3)        - Correlations between Lat, long and rate
* cov_out(3,3)  - Derived covariance matrix

* temp          - Introduced to HP compiler warning 


      real*8 rot_mat(3,3), trans(3,6), loc_coord(3), xyz_adj(3),
     .       neu_adj(3) , cov_pole(6,6), 
     .       cov_sing(3,3), temp_cov(6,6), est(3), sig(3), 
     .       rho(3), cov_out(3,3)

* first  - True if first line of sites not written

      logical  first

* list   - Line with site names written

      character*128 list

      j = 0

      write(iout,100)
 100  format(/,' PLATE ROTATION VECTOR RESULTS',/,
     .         ' -----------------------------')

      write(iout,120)
 120  format('SITES USED TO ESTIMATE POLE POSITIONS',/,
     .       '-------------------------------------',/,
     .       ' PLATE      SITES')

      do i = 1, num_plates

****     Build up list of sites
         first = .true.
         nu = 0
         do j = 1, gnum_sites
            pn = plate_number(j)
            if( pn.eq.i ) then
                nu = nu + 1
                list(9*(nu-1)+1:) = gsite_names(j)
                if( nu.eq.6 ) then
                    if( first ) then
                        write(iout,150) plate_names(i), list(1:60)
 150                    format(1x,a8,' : ',a)
                        first = .false.
                    else
                        write(iout,'(12x,a)') list(1:60)
                    end if
                    nu = 0
                end if
             end if
          end do

*****     See if anything left in list:
          if( nu.gt.0 ) then
              if( first ) then
                  write(iout,150) plate_names(i), list(1:9*nu)
              else
                  write(iout,'(12x,a)') list(1:9*nu)
              end if
          end if
      end do

*  MOD TAH 980612: List the assigned sites to plates.
      write(iout,160)
 160  format('SITES ASSIGNED TO PLATE BUT NOT USED IN ESTIMATE',/,
     .       '------------------------------------------------',/,
     .       ' PLATE      SITES')

      do i = 1, num_plates

****     Build up list of sites
         first = .true.
         nu = 0
         do j = 1, gnum_sites
            pn = assign_number(j)
            if( pn.eq.i ) then
                nu = nu + 1
                list(9*(nu-1)+1:) = gsite_names(j)
                if( nu.eq.6 ) then
                    if( first ) then
                        write(iout,150) plate_names(i), list
                        first = .false.
                    else
                        write(iout,'(12x,a)') list
                    end if
                    nu = 0
                end if
             end if
          end do

*****     See if anything left in list:
          if( nu.gt.0 ) then
              if( first ) then
                  write(iout,150) plate_names(i), list(1:9*nu)
              else
                  write(iout,'(12x,a)') list(1:9*nu)
              end if
          end if
      end do

****  MOD TAH 030513: Output the translation estimate
      if ( PlateTrans ) then
         j = 3*num_plates
         do i = 1,3
            sig(i) = sqrt(norm_eq(j+i,j+i))
         end do
         rho(1) = norm_eq(j+1,j+2)/(sig(1)*sig(2))
         rho(2) = norm_eq(j+1,j+3)/(sig(1)*sig(3))
         rho(3) = norm_eq(j+2,j+3)/(sig(2)*sig(3))
         write(iout, 180)
 180     format(/,' TRANSLATION  X (mm/yr)    +-     Y  (mm/yr)   +-',
     .         '        Z  (mm/yr)   +-       RhoXY   RhoXZ   RhoYZ')
         write(iout, 190) (bvec(j+i)*1000,sig(i)*1000, i=1,3),
     .                 (rho(i), i=1,3)
 190     format(10x,3(F10.2,1x,F10.2,1x),1x,3(F6.3,1x),' TRAN')
      endif


****  Now write the results
      write(iout, 200)
 200  format(/,' PLATE       Wx (deg/My)   +-    Wy (deg/My)   +-',
     .         '       Wz (deg/My)   +-       RhoXY   RhoXZ   RhoYZ')
      do i = 1, num_plates

         do j = 1,3
            do k = 1,3
               cov_out(j,k) = norm_eq(3*(i-1)+j,3*(i-1)+k)
            end do
         end do
         do j = 1, 3
            sig(j) = sqrt(cov_out(j,j))
         end do
         do j = 1,2
           rho(j) = cov_out(1,j+1)/(sig(1)*sig(j+1))
         end do
         rho(3) = cov_out(2,3)/(sig(2)*sig(3))

         write(iout, 220) plate_names(i), 
     .           (bvec(3*(i-1)+j)*180.d0/pi,
     .           sig(j)*180.d0/pi, j=1,3), (rho(j),j=1,3)
 220     format(1x,a8,1x,3(F10.6,1x,F10.6,1x),1x,3(F6.3,1x),' XYZ')
      end do

****  Now compute the lat, long, and length of each vector
      do i = 1,3
         xyz_adj(j) = 0.d0
      end do
      write(iout, 300)
 300  format(/,' PLATE       Lat. (deg)    +-    Long (deg)    +-',
     .         '       Mag (deg/My)   +-   RhoLtLg RhoLtMg RhoLgMg')

      do i = 1, num_plates
         is = 3*(i-1)+1
         call rotat_crd(xyz_adj, neu_adj, 'XYZ', 'NEU',
     .                   bvec(3*(i-1)+1), loc_coord, rot_mat)

*        Get the covariance matrix elements
         do k = 0,2
            ip(k+1) = is+k
         end do
         call mov_cov(ip, ip, 3, cov_sing,3, norm_eq, 3*max_plates)

*        Now compute var_covariance matrix 
         call var_comp(rot_mat, cov_sing, cov_out, temp_cov, 3,3,1)

*        Now convert the results:
         est(1) = 90.d0 - loc_coord(1)*180.d0/pi
         est(2) = loc_coord(2)*180.d0/pi
         est(3) = loc_coord(3)*180.d0/pi

*        Now get sigmas
         sig(1) = sqrt(cov_out(1,1))/loc_coord(3)*180.d0/pi
         sig(2) = sqrt(cov_out(2,2)/sin(loc_coord(1))**2)/
     .                 loc_coord(3)*180.d0/pi
         sig(3) = sqrt(cov_out(3,3))*180.d0/pi

*        Now get the correlations. NOTE: Correlations do not depend
*        units provided I use the sqrt(cov_out) diagonal elements
*        as sigmas.
         rho(1) = cov_out(1,2)/sqrt(cov_out(1,1)*cov_out(2,2))
         rho(2) = cov_out(1,3)/sqrt(cov_out(1,1)*cov_out(3,3))
         rho(3) = cov_out(2,3)/sqrt(cov_out(2,2)*cov_out(3,3))

         write(iout, 320) plate_names(i), (est(j), sig(j), j=1,3),
     .                    (rho(j),j=1,3)
 320     format(1x,a8,1x,2(F10.3,1x,F10.3,1x),(F10.6,1x,F10.6),
     .          2x,3(f6.3,1x),' LLM')
      end do


****  Now do the "differences between plates"
      if( num_plates.gt.1 ) write(iout, 400)
 400  format(/,' PLATE -  PLATE        Wx (deg/my)   +-    ',
     .         'Wy (deg/my)   +-       Wz (deg/My)   +-    ',
     .         'RhoXY   RhoXZ   RhoYZ')

*     Form the transformation
      do i = 1,3
         do j = 1,3
            trans(i,j)   = 0.d0
            trans(i,j+3) = 0.d0
         end do
         trans(i,i)   =  1.d0
         trans(i,i+3) = -1.d0

      end do

      do i = 1,num_plates - 1
         is = 3*(i-1)+1
         do j = 0,2
            ip(j+1) = is + j
         end do
         call mov_cov(ip,ip, 3, cov_pole(1,1), 6, norm_eq,3*max_plates)
         do j = i+1, num_plates

*           Move the covariance matrix elements into plate
            js = 3*(j-1)+1
            do k = 0,2
               jp(k+1) = js + k
            end do
            call mov_cov(ip,jp,3,cov_pole(1,4),6,norm_eq,3*max_plates)
            call mov_cov(jp,ip,3,cov_pole(4,1),6,norm_eq,3*max_plates)
            call mov_cov(jp,jp,3,cov_pole(4,4),6,norm_eq,3*max_plates)

*           Now form the solution
            do k = 1,3
               est(k) = bvec(is+k-1) - bvec(js+k-1)
            end do

*           Now form the covariance matrix
            call var_comp(trans, cov_pole, cov_out, temp_cov, 3,6,1)
            do k = 1,3
               sig(k) = sqrt(cov_out(k,k))
            end do
            do k = 1,2
               rho(k) = cov_out(1,k+1)/(sig(1)*sig(k+1))
            end do
            rho(3) = cov_out(2,3)/(sig(2)*sig(3))

            write(iout,420) plate_names(i), plate_names(j), 
     .                     (est(k)*180.d0/pi, sig(k)*180.d0/pi, k=1,3), 
     .                     (rho(k), k=1,3)
 420        format(1x,a8,'-',a8,1x,3(F10.6,1x,F10.6,1x),
     .             1x,3(F6.3,1x),' XYZ')
         end do
      end do

****  Now do the "differences between plates".  Write the results as
*     lat/long/and magnitude.  Here we re-form the differences and
*     then convert to lat, long and magnitude
      if( num_plates.gt.1 ) write(iout, 500)
 500  format(/,' PLATE -  PLATE        Lat (deg)     +-    ',
     .         'Long (deg)    +-       Mag (deg/My)   +-   ',
     .         'RhoLaLg RhoLaMa RhoLgMa')

      do i = 1,num_plates - 1
         is = 3*(i-1)+1
         do j = 0,2
            ip(j+1) = is + j
         end do
         call mov_cov(ip,ip, 3, cov_pole(1,1), 6, norm_eq,3*max_plates)

         do j = i+1, num_plates

*           Move the covariance matrix elements into plate
            js = 3*(j-1)+1
            do k = 0,2
               jp(k+1) = js + k
            end do
            call mov_cov(ip,jp,3,cov_pole(1,4),6,norm_eq,3*max_plates)
            call mov_cov(jp,ip,3,cov_pole(4,1),6,norm_eq,3*max_plates)
            call mov_cov(jp,jp,3,cov_pole(4,4),6,norm_eq,3*max_plates)

*           Now form the solution
            do k = 1,3
               est(k) = bvec(is+k-1) - bvec(js+k-1)
            end do

*           Now form the covariance matrix
            call var_comp(trans, cov_pole, cov_sing, temp_cov, 3,6,1)

****        Rotate system into NEU components.
            call rotat_crd(xyz_adj, neu_adj, 'XYZ', 'NEU',
     .                     est, loc_coord, rot_mat)

*           Now compute var_covariance matrix 
            call var_comp(rot_mat, cov_sing, cov_out, temp_cov, 3,3,1)

*           Now convert the results:
            est(1) = 90.d0 - loc_coord(1)*180.d0/pi
            est(2) = loc_coord(2)*180.d0/pi
            est(3) = loc_coord(3)*180.d0/pi

*           Now get sigmas
            sig(1) = sqrt(cov_out(1,1))/loc_coord(3)*180.d0/pi
            sig(2) = sqrt(cov_out(2,2)/sin(loc_coord(1))**2)/
     .                 loc_coord(3)*180.d0/pi
            sig(3) = sqrt(cov_out(3,3))*180.d0/pi

*           Now get the correlations. NOTE: Correlations do not depend
*           units provided I use the sqrt(cov_out) diagonal elements
*           as sigmas.
            rho(1) = cov_out(1,2)/sqrt(cov_out(1,1)*cov_out(2,2))
            rho(2) = cov_out(1,3)/sqrt(cov_out(1,1)*cov_out(3,3))
            rho(3) = cov_out(2,3)/sqrt(cov_out(2,2)*cov_out(3,3))
            write(iout,520) plate_names(i), plate_names(j), 
     .                     (est(k), sig(k), k=1,3), (rho(k),k=1,3)
 520        format(1x,a8,'-',a8,1x,2(F10.3,1x,F10.3,1x),
     .               (F10.6,1x,F10.6),2x,3(F6.3,1x),' LLM')
         end do
      end do

****  Thats all
      return
      end 

         
