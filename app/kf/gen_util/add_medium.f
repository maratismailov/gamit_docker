CTITLE ADD_MEDIUM
 
      subroutine add_medium(in)

      implicit none 
 
*     Routine to add the atmospheric propagation medium effects to
*     the theoretical delays.  When the user tries to apply a
*     particular calibration, the available of the information needed
*     is checked.  If all of the required information is not available
*     then the next 'best' calibration is used.  The cont_medium
*     values are updated to reflect this change.
*
*     Currently, the wet mapping function is used for the atmospheric
*     delay partial.
*
*     See MEDIUM_CONT_TYPES in OBS_VALUES.FTNI for list of options in the
*     calibrations.
 
      include '../includes/kalman_param.h'
      include '../includes/obs_values.h'
      include '../includes/obs_data.h'
 
*   in          - Indicates is site 1 or 2
 
      integer*4 in
 
*   kbit        - Bit checking function
 
      logical kbit
 
*   dry_map(2)  - Dry mapping function for delay and rate
*   dry_zen(2)  - Dry zenith delay term for delay and rate
*   wet_map(2)  - Wet mapping function for delay and rate
*   wet_zen(2)  - Wet zenith delay term for delay and rate.
*   sign        - Sign of contribution (plus for in=2, minus
*               - for in=1)
 
      real*8 dry_map(2), dry_zen(2), wet_map(2), wet_zen(2), sign
 
***** Start, work down all the possible choices and see what is to
*     be done.
 
***** DRY ZENITH TERM
*                                               ! Dry Sasstamonian delay
      if( kbit(cont_medium(site(in)),1) ) then
          call dry_saas( in, dry_zen )
      end if
 
*                                               ! Constant term given by user
      if( kbit(cont_medium(site(in)),2) ) then
          call dry_const( in, dry_zen )
      end if
 
*     if( kbit(cont_medium(site(in)),3) ) then  ! User function
*         call dry_user( in, dry_zen )
*     end if
 
*     if( kbit(cont_medium(site(in)),4) ) then  ! User dry term in data file
*         call dry_user_sup( in, dry_zen )
*     end if
 
*                                               ! Constant term based on
      if( kbit(cont_medium(site(in)),6) ) then
*                                               ! height
          call dry_hgt( in, dry_zen )
      end if
 
***** WET ZENITH TERM
 
*                                               ! Wet Sasstamonian delay
      if( kbit(cont_medium(site(in)),7) ) then
          call Wet_saas( in, Wet_zen )
      end if
 
C     if( kbit(cont_medium(site(in)),8) ) then  ! Wet delay from WVR antenna
*                                               ! values
C         call Wet_WVR_ant( in, Wet_zen )
C     end if
 
*                                               ! Wet delay from DataBase WVR
      if( kbit(cont_medium(site(in)), 9) ) then
*                                               ! values
          call Wet_WVR_DB( in, Wet_zen )
      end if
 
*     if( kbit(cont_medium(site(in)),10) ) then  ! User function
*         call Wet_user( in, Wet_zen )
*     end if
 
*     if( kbit(cont_medium(site(in)),11) ) then  ! User wet term in data file
*         call Wet_user_sup( in, Wet_zen )
*     end if
 
*                                                ! Constant term given by user
      if( kbit(cont_medium(site(in)),13) ) then
          call Wet_const( in, Wet_zen )
      end if
 
***** DRY MAPPING FUNCTION
 
*                                                ! Dry Marini mapping
      if( kbit(cont_medium(site(in)),14) ) then
          call Dry_mari( in, Dry_map , dry_zen, wet_zen )
      end if
 
*                                                ! CfA-2.2 mapping
      if( kbit(cont_medium(site(in)),15) ) then
          call Dry_CfA( in, Dry_map )
      end if
 
*                                                ! MIT mapping
      if( kbit(cont_medium(site(in)),16) ) then
          call Dry_MIT( in, Dry_map )
      end if

      if( kbit(cont_medium(site(in)),17) ) then
          call dry_mtt( in, Dry_map )
      end if
 
*     if( kbit(cont_medium(site(in)),17) ) then  ! User mapping function
*         call Dry_user_map( in, Dry_map )
*     end if
 
*                                                ! CHAO mapping
      if( kbit(cont_medium(site(in)),19) ) then
          call Dry_chao( in, Dry_map )
      end if
 
 
***** WET MAPPING FUNCTION
 
*                                                ! Wet Marini mapping
      if( kbit(cont_medium(site(in)),20) ) then
*                                                          ! Use Dry formula
          call Wet_mari( in, Wet_map , dry_zen, wet_zen )
      end if
 
*                                                ! CfA-2.2 mapping
      if( kbit(cont_medium(site(in)),21) ) then
*                                      ! Use Dry at the moment
          call Wet_CfA( in, Wet_map )
      end if
 
*                                                ! Cosecant mapping function
      if( kbit(cont_medium(site(in)),22) ) then
          call Wet_csc( in, Wet_map )
      end if
 
c     if( kbit(cont_medium(site(in)),23) ) then  ! Ifadis wet mapping function
c         call Wet_ifad( in, Wet_map )
c     end if
 
*                                              ! MIT Wet mapping function
      if (kbit(cont_medium(site(in)),24)) then
          call wet_mit(in,wet_map)
      end if
*                                              ! MTT Wet mapping function 
      if (kbit(cont_medium(site(in)),25)) then
          call wet_mtt(in,wet_map)
      end if
 
*                                                ! CHAO mapping
      if( kbit(cont_medium(site(in)),26) ) then
          call Wet_chao( in, Wet_map )
      end if
 
***** Now we have everything we need, update theoretical.
 
      sign = 2*in - 3
 
      theoretical(1) = theoretical(1) +
     .    sign*( dry_zen(1)*dry_map(1) + wet_zen(1)*wet_map(1) )
      theoretical(2) = theoretical(2) +
     .    sign*( dry_zen(1)*dry_map(1) + wet_zen(1)*wet_map(1) )
      theoretical(3) = theoretical(3) +
     .    sign*( dry_zen(1)*dry_map(1) + wet_zen(1)*wet_map(1) )
 
      theoretical(4) = theoretical(4) +
     .    sign*( dry_zen(2)*dry_map(1) + dry_zen(1)*dry_map(2) +
     .           wet_zen(2)*wet_map(1) + wet_zen(1)*wet_map(2))
 
***** Save the wet mapping function
 
      atm_part(in,1) = sign*wet_map(1)
      atm_part(in,2) = sign*wet_map(2)
 
***** Thats all
      return
      end
 
 
