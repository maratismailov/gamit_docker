CTITLE WET_WVR_DB
 
      subroutine Wet_WVR_db( in, wet_zen )

      implicit none 
 
 
*     Routine to return the values for the wet zenith delay as measured
*     by the WVR delays and rate stored in the database.
 
* JLD 870706 "Mapped" type 6 WVR codes (see OBS_DATA.FTNI) to type
*     2 WVR codes
 
      include '../includes/kalman_param.h'
      include '../includes/obs_values.h'
      include '../includes/obs_data.h'
 
*   in      - Index to site being processed (1 or 2)
*   ibitt   - Temporary bit index
*   temp_code   - Working version of WVR code
 
      integer*4 in, ibitt, temp_code
 
*   wet_zen(2)  - Wet zenith delay (ps and fs/s)
 
      real*8 wet_zen(2)
 
*   kbit    - Utility to see if bit turned on
 
      logical kbit
 
***** Get working version of wvr code
      temp_code = mod(wvr_code(in),100)
 
***** Map type 6 WVR codes to type 2
      if (temp_code .eq. 6) temp_code = 2
 
***** See if database WVR values are availble, if they are not turn
*     this calibration off.
 
*                                              ! Database WVR available
      IF( kbit(avail_met(site(in)),5) ) THEN
 
*         Now see if good for this observation
*                                                      ! WVR good
          if( .not.kbit(medium_flag(in),5) .and.
     .        (temp_code.ge.0 .and.
*                                                      ! Passes code
     .         temp_code.le.wvr_code_limit) ) then
*                                                      ! quality. Ignore
*                                                      ! number of samples in
*                                                      ! upper part of code
              wet_zen(1) = wvr_dbase(in,1)
              wet_zen(2) = wvr_dbase(in,2)
 
*                                              ! WVR bad for this obs.
          else
*                                              ! Flag data
              ibitt = 4
*                                              ! bad WVR delay
              call sbit( data_flag,   4 ,1)
 
          end if
 
*                                 ! WVR not available, resort to user
      ELSE
*                                 ! constant value
*                                                  ! turn off WVR dbase
          call sbit( cont_medium(site(in)), 9,0 )
*                                                  ! turn on UCON
          call sbit( cont_medium(site(in)),13,1 )
 
      END IF
 
***** Thats all
      return
      end
 
