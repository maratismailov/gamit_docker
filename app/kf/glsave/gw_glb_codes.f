CTITLE GW_GLB_CODES

      subroutine gw_glb_codes

      implicit none 
 
 
*     This routine will generate the parameter codes for both the
*     apriori site and source positions and the estimated parameters
*     It will also determine the number of the first parameter to be
*     saved.
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
      include '../includes/obs_header.h'
 
      include '../includes/glsave_comm.h' 

      integer*4 na, np, ne, i,j,k
 
 
***** Countup the list of apriori site/source positions which were
*     not estimated and make list.
 
      num_apr_codes = 0
 
      call glb_apr_codes(parn_site,   3, 2, -gnum_sites,   7, 
     .                   guse_site, 1 )
      call glb_apr_codes(parn_atm,    1, 1, -gnum_sites,   61, 
     .                   guse_site, 0 )
      call glb_apr_codes(parn_axo,    1, 2, -gnum_sites,   10, 
     .                   guse_site, 0 )
      call glb_apr_codes(parn_source, 2, 2, gnum_sources, 11,0, 1)
 
      call glb_apr_codes(parn_site(1,2,1),   3,2,-gnum_sites,    42, 
     .                   guse_site, 1)
*     Do not save axis offset rate of change.
      call glb_apr_codes(parn_source(1,2,1),2,2, gnum_sources, 45,0,1)

*     Save the apriori codes for the satellite orbits
      call glb_svs_codes(parn_svs(1,1), max_svs_elem,1, gnum_svs, 51,1)
 
****  Now generate codes for estimated parameters
 
      call glb_glb_codes(parn_site,    3,2, -gnum_sites,    7 , 
     .                   guse_site)
      call glb_glb_codes(parn_atm,     1,1, -gnum_sites,    61, 
     .                   guse_site)
      call glb_glb_codes(parn_axo,     1,2, -gnum_sites,    10, 
     .                   guse_site)
      call glb_glb_codes(parn_source,  2,2, gnum_sources,  11,0)
      call glb_glb_codes(parn_wob ,    1,1,            8,  13,0)
      call glb_glb_codes(parn_ut1,     1,1,            6,  14,0)
      call glb_glb_codes(parn_nut_ang, 1,1,            8,  15,0)

* MOD TAH 981104: Now add up the number of multi-day parameter epochs
*     needed.
      na = num_apr_codes
      ne  = 0 
      do i = 1, 3
         do j = 1,2
            do k = 1, num_mul_pmu
               if( parn_mul_pmu(j,i,k).ne.0 ) then
                  ne = ne + 1
                  na = na + 1
                  apr_codes(na) = 55+i + j*256 + ne*65536
*                 This is slightly out of place: Normally apriori values
*                 are normally saved in gw_aprioris 
                  apr_list(na) = apr_val_mul_pmu(j,i,k)
                  np = parn_mul_pmu(j,i,k)
                  glb_codes(np) = 55+i + j*256 + ne*65536
                  mul_par_ep(ne) = gmul_pmu_ep(k)
               end if
            end do
         end do
      end do 

*     Now save values
      num_apr_codes = na
      ent_par_ep = ne
 
      if( glb_glb_tides ) then
          call glb_glb_codes(parn_tid ,3,1,            0,  16,0)
      else
          call glb_glb_codes(parn_tid ,3,1,   -gnum_sites,  16, 
     .                   guse_site)
      end if
 
      call glb_glb_codes(parn_gamma,   1,1,             1, 41,0)
 
      call glb_glb_codes(parn_site(1,2,1),  3,2, -gnum_sites,   42, 
     .                   guse_site)

* MOD TAH 030615: Generate codes for log estimates
      call glb_glb_codes(parn_log(1,1), 3,1, -gnum_sites,   62, 
     .                   guse_site)


*     Do not save axis offset parameter rates of change
c     call glb_glb_codes(parn_axo(2,1),     1,2, gnum_sites,    0)
      call glb_glb_codes(parn_source(1,2,1),2,2, gnum_sources, 45,0)
 
      call glb_coeff_codes(parn_nut_coeff, 2, max_nut_coeff,    47,0)

*     Generate the codes the satellite orbit parameters
      call glb_svs_codes(parn_svs(1,1), max_svs_elem,1, gnum_svs, 51,0)

*     Generate the translation and translation rate of change codes 
      call glb_glb_codes(parn_tran(1,1), 1,1,            3,  52,0)
      call glb_glb_codes(parn_tran(1,2), 1,1,            3,  53,0) 

*     Generate the rotation and rotation rate of change codes 
      call glb_glb_codes(parn_rot (1,1), 1,1,            3,  59,0)
      call glb_glb_codes(parn_rot (1,2), 1,1,            3,  60,0)   

*     Generate codes for scale and scale rate of change
      call glb_glb_codes(parn_scale(1) , 1,1,            1,  54,0)
      call glb_glb_codes(parn_scale(2) , 1,1,            1,  55,0)

***** Thats all
      return
      end
 
