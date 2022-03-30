CTITLE CHANGE_LOAD

      subroutine change_load ( sol_obs )

      implicit none

*     Routine to change the loading applied to the solution (either 
*     add or remove as requested by user).

      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
      include '../includes/glb_hdr_def.h'

* PASSED 
*   sol_obs(cnum_parn)  - Solution vector.  

      real*8 sol_obs(cnum_parn) 

* LOCAL VARIABLES
* i,is  -- Loop counters, site number (global)
* dXYZ_atm(3), dXYZ_hyd(3) dNEU(3) -- load corrections in XYZ (m) for
*      atm and hydrology. Values computed from NEU in mm
* type, indx  -- Type and index for parameters in solution

      integer*4 i, type, indx, is 
      real*8 dXYZ_atm(3), dXYZ_hyd(3), dNEU(3)
      real*8 sign ! sign of correction : To apply the load, we remove it from
                  ! from the observed location (solobs) and sign is -1; to
                  ! remove the load coorection we add it back to the observed
                  ! values (sign +1)
      real*8 loc_coord(3), rot_matrix(3,3)  ! Need to transform NEU to XYZ
      logical kbit

****  Loop over the solvectr updating the values as needed.

      do i = 1, cnum_parn
         call decode_code(gpar_codes(i), type, indx )
         if( type.eq.7 ) then   ! This is a X-coordinate so
*             compute the changes needed.
*             See if atm load to be changed
              is = ltog_sites(indx)
              dXYZ_atm = 0.d0
              dXYZ_hyd = 0.d0
*             Check atmospheric load
              if( kbit(appload_mod,2) .neqv. kbit(cload_mod, 9) ) then
*                 Get the dXYZ (m) of load
                  sign = 1.d0  
*                 If bit set, we want to apply load correction.
                  if( kbit(appload_mod,2) ) sign = -1.d0
                  dNEU = sign*gatmload(:,is)/1000.d0   ! Convert to m and
*                      set sign so that added to sol_obs
                  call rotate_geod(dNEU, dXYZ_atm, 'NEU','XYZ',
     .                    apr_val_site(1,1,is), loc_coord, rot_matrix)
              endif
*             Check hydrologic load
              if( kbit(appload_mod,3) .neqv. kbit(cload_mod,25) ) then
*                 Get the dXYZ (m) of load
                  sign = +1.d0
                  if( kbit(appload_mod,2) ) sign = -1.d0
                  dNEU = sign*ghydload(:,is)/1000.d0   ! Convert to m and
*                      set sign so that added to sol_obs
                  call rotate_geod(dNEU, dXYZ_hyd, 'NEU','XYZ',
     .                    apr_val_site(1,1,is), loc_coord, rot_matrix)
              endif

              sol_obs(i) = sol_obs(i) + dXYZ_atm(1) + dXYZ_hyd(1)
              if( 1.eq. 2 )
     .        write(*,150) gsite_names(ltog_sites(indx)), 
     .                     dXYZ_atm(:)*1000, dXYZ_hyd(:)*1000
 150          format('Load Update Site ',a,' dXYZ_atm ',3f7.2,
     .               ' mm, dXYZ_hyd ',3F7.2,' mm')
         else if( type.eq.8 ) then
              sol_obs(i) = sol_obs(i) + dXYZ_atm(2) + dXYZ_hyd(2)
         else if( type.eq.9 ) then
              sol_obs(i) = sol_obs(i) + dXYZ_atm(3) + dXYZ_hyd(3)
         end if
      end do

****  Thats all
      return
      end

