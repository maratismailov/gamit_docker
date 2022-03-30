CTITLE SITE_PARTIAL
 
      subroutine site_partial( comp, pn, part_pnt, a_part, site)

      implicit none  
 
*     Routine to compute the site position partial derivatives
*     These depend on:
*     i)  the global site position
*     ii) translation parameters
*     iii)the PM/UT1 parameters in these were not estimated during the
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
 
      integer*4 comp, i, part_pnt(2,*), pn, site
 
*   a_part(1)   - Partial derivatives (non_zero values only)
*   deriv(3)    - Three derivatives for the PM/UT1 partials
 
      real*8 a_part(*), deriv(3)
 
*   kbit        - Bit checking function
 
      logical kbit
 
 
***** If this site is not being used then just return
 
      if( .not.kbit(guse_site, site) ) RETURN

 
****  See if global position turned on
 
      call add_partial(indx_pnt(pn), part_pnt, parn_site(comp,1,site),
     .                  a_part, 1.d0 )

* MOD TAH 131226: Only consider partials if the the coordinate of this
*     is being estimated.  If it is not estimated  or if it is an
*     excluded site (mostly would be used except _XPS on single day).
*     (Test after site position partial since XPS coordinates are estimated).
      if( parn_site(comp,1,site).eq.0 .or.
     .    gsite_names(site)(5:6).eq.'_X' ) RETURN

****  Add paritials for log after earthquake.  Only add partial if logs are
*     not included in the global file itself
      if( .not.kbit(tran_est,16) )
     .     call add_log(comp, pn, part_pnt, a_part, site)
 
****  See if translation is turned on, and not turned on in input data
      if( .not.kbit(tran_est,comp) )
     .call add_partial(indx_pnt(pn), part_pnt, parn_tran(comp,1),
     .                 a_part, 1.d0 )

****  Compute the rotation partials in case needed.
      call pmu_part(site, deriv, apr_val_site(1,1,site), comp,
     .              cut1_apr) 

****  Check if rotation is turned on, and not turned on in the input data
!      print *,'DEBUG: ROT_EST ',rot_est
      do i = 1, 3
         if( .not.kbit(rot_est,i) ) then
            call add_partial(indx_pnt(pn), part_pnt, parn_rot(i,1),
     .                 a_part, deriv(i) ) 
!            print *,'DEBUG: ROT  ',i,pn, indx_pnt(pn), parn_rot(i,1)
!            print *,'DEBUG: PARN ',part_pnt(:,1:2), deriv(i)
         end if
      end do

***** See if scale is turned on, and not estimated in input data
      if( .not.kbit(scale_est, 1) )
     .call add_partial(indx_pnt(pn), part_pnt, parn_scale(1),
     .                 a_part, apr_val_site(comp,1,site)*1.d-9 )
 
***** Now see if we need PM/UT1 parameters.  (Only needed if pmu were
*     not estimated along with the site coordinates).  If polar motion
*     or UT1 not estimated then compute partials. 
 
*                             ! X,Y wobble
      do i = 1,2
         if( .not.kbit(pmu_est,i) )
     .     call add_partial(indx_pnt(pn) ,part_pnt, parn_wob(i),
     .                     a_part, deriv(i) )
      end do
      if( .not.kbit(pmu_est,3) )  
     .     call add_partial(indx_pnt(pn), part_pnt, parn_ut1(1),
     .                     a_part, deriv(3) )
 
***** Thats all
      return
      end
 
CTITLE ADD_LOG

      subroutine add_log( comp, pn, part_pnt, a_part, site)

      implicit none 

*     Routine to add earthquake log partials

 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
 
      include '../includes/glb_hdr_def.h'
 
 
*   comp        - Component of position (1=X,2=Y,3=Z)
*   part_pnt(2,*)   - Partial pointers
*   pn          - local parameter number
*   site        - globk site number
 
      integer*4 comp, part_pnt(2,*), pn, site
 
*   a_part(1)   - Partial derivatives (non_zero values only)
 
      real*8 a_part(*)

* LOCAL VARIABLES
*   loc_coord(3) - Local lat,long and height
*   rot_mat(3,3) - rotation matrix from NEU to XYZ
*   dtl          - Log argument (1-dt/tau)
* MOD TAH 090811: Added accounting for a site being effected by 
*                 multiple earthquakes.
*   dpdt(max_eq)  - partial of position with respect to time (allow
*                 upto all earthquakes affecting site
*   eq_plst(3,max_eq) -- parameter numbers of log terms that affect
*                 this site
*   num_eqaff    - number of earthquakes affecting this site

      real*8 loc_coord(3), rot_mat(3,3), dtl, dpdt(max_eq)
      integer*4 eq_plst(3,max_eq), num_eqaff

*     i, ne -- looop counters
      integer*4 i, j, ne

      save rot_mat, dpdt, eq_plst, num_eqaff


***   For the first component we need to find which earthquake
*     MODTAH 090811: Removed because this site might be affected
*     by an earilier earthquake
C     if( parn_log(1,site) + parn_log(2,site)
C    .   +parn_log(3,site).eq.0               ) RETURN

      if( comp.eq.1 ) then
*        MODTAH 090811:See how many eq's affect this site
         num_eqaff = 0
         do ne = 1, num_eq
C            if( gsite_names(site)(7:8).eq.eq_codes(ne)(1:2) ) then
*
*               OK, found the earthquake.  Generate the time partial
C                dtl = (gepoch_expt-eq_epoch(ne))/eq_log_tau(ne)            
C                dpdt = log(1+dtl)
C
C            end if
*            MOD TAH 090811: Loop over possible sites to get parameters
             do i = 1, site  ! Could use a shorter list but is OK
*               See if we match site code and eq code
                if( gsite_names(i)(1:4).eq. gsite_names(site)(1:4) .and.
     .              gsite_names(i)(7:8).eq. eq_codes(ne)(1:2) ) then
*                   This earthquake affects this site, so add into list
*                   if there is a log parameter
                    if( parn_log(1,i).ne.0 ) then
                        num_eqaff = num_eqaff + 1
                        do j = 1, 3
                           eq_plst(j,num_eqaff) = parn_log(j,i)
                        end do
*                       Now get the time partial                        
                        dtl = (gepoch_expt-eq_epoch(ne))/eq_log_tau(ne)            
                        dpdt(num_eqaff) = log(1+dtl)
                    endif
                endif
             enddo   
*            END MODTAH 090811 
         end do
*        Now get the NEU to XYZ transformation
         if( num_eqaff.gt.0 ) then
            call xyz_to_neu(rot_mat, apr_val_site(1,1,site), loc_coord)

*           Now transpose rot_matrix ! See Page 9-15 VIS section of
*           relocable libaries
            do i = 1,2
               call dvswp(rot_mat(i,i+1),3,rot_mat(i+1,i),1,3-i)
            end do
         endif
c        write(*,999) site, gsite_names(site),num_eqaff,
c    .               (dpdt(i),eq_plst(1,i),i=1, num_eqaff)
c999     format('EQ Logs: ',i4,1x,a8,1x,'EQs' ,i3,20(F10.3,1x,I6))
      endif

*     MOD TAH 090811: If no earthquakes found, return
      if( num_eqaff.eq.0 ) RETURN

****  OK, Now add the partials
*     MOD TAH 090811: Loop over all affecting earthqukes
      do j = 1, num_eqaff
         do i = 1,3
C           call add_partial(indx_pnt(pn), part_pnt, parn_log(i,site),
C    .                    a_part, rot_mat(comp,i)*dpdt )
            call add_partial(indx_pnt(pn), part_pnt, eq_plst(i,j),
     .                    a_part, rot_mat(comp,i)*dpdt(j) )
         enddo
      end do

****  Thats all for a single earthquake effected site
      return
      end

 
