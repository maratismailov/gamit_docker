CTITLE site_rate_part
 
      subroutine site_rate_part( comp, pn, part_pnt, a_part, site)

      implicit none  
 
*     Routine to compute the site rate partial derivatives
*     These depend on:
*     i)  the global site rates   
*     ii) translation parameters (not yet)
*     iii)the PM/UT1 rate parameters in these were not estimated during the
*         the solution.
*
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
 
      include '../includes/glb_hdr_def.h'
 
 
*   comp        - Component of position (1=X,2=Y,3=Z)
*   i           - Loop counter
*   part_pnt(2,1)   - Partial pointers
*   pn          - local parameter number
*   site        - globak site number
 
      integer*4 comp, i, part_pnt(2,1), pn, site
 
*   a_part(1)   - Partial derivatives (non_zero values only)
*   deriv(3)    - Three derivatives for the PM/UT1 partials
 
      real*8 a_part(1), deriv(3)
 
*   kbit        - Bit checking function
 
      logical kbit
 
 
 
***** If this site is not being used then just return
 
      if( .not.kbit(guse_site, site) ) RETURN
 
****  See if global position turned on
 
      call add_partial(indx_pnt(pn), part_pnt, parn_site(comp,2,site),
     .                  a_part, 1.d0 )
 
****  See if translation is turned on , and not in input date
      if( .not.kbit(tran_est,comp+3) )
     .call add_partial(indx_pnt(pn), part_pnt, parn_tran(comp,2),
     .                 a_part, 1.d0 )

****  Compute the rotation partials in case needed.
      call pmu_part(site, deriv, apr_val_site(1,1,site), comp,
     .              cut1_apr) 

****  Check if rotation is turned on, and not turned on in the input data
      do i = 1, 3
         if( .not.kbit(rot_est,i+3) ) then
            call add_partial(indx_pnt(pn), part_pnt, parn_rot(i,2),
     .                 a_part, deriv(i) ) 
         end if
      end do

****  See if scale rate of change estimated.
      if( .not.kbit(scale_est,2) )
     .call add_partial(indx_pnt(pn), part_pnt, parn_scale(2),
     .                 a_part, apr_val_site(comp,1,site)*1.d-9 )

 
***** Now see if we need PM/UT1 parameters.  (Only needed if pmu were
*     not estimated along with the site coordinates)
 
*                             ! X,Y wobble
      do i = 1,2
	   if( .not.kbit(pmu_est,i+3) ) 
     .       call add_partial(indx_pnt(pn) ,part_pnt, parn_wob(i+2),
     .                     a_part, deriv(i) )
      end do
      if( .not.kbit(pmu_est,6) ) 
     .    call add_partial(indx_pnt(pn), part_pnt, parn_ut1(2),
     .                 a_part, deriv(3) )
 
***** Thats all
      return
      end
 
