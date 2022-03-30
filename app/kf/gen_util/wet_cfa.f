CTITLE WET_CFA
 
      subroutine wet_cfa( in, wet_map )
 
      implicit none 

 
*     Routine to compute the CfA-2.2 Dry mapping function.  If the pressure
*     temperature and humidity are not available, the mapping is set back
*     to CHAO.
 
      include '../includes/kalman_param.h'
      include '../includes/obs_values.h'
      include '../includes/obs_data.h'
 
*   in      - Index to site in baseline (1 or 2)
 
      integer*4 in
 
*   wet_map(2)  - Delay and rate mapping function
 
      real*8 wet_map(2)
 
*   kbit        - Ckeck to see if bit set
 
      logical kbit
 
***** Start, see if we pressure, temperature and rel.hum.
 
      if( kbit(avail_met(site(in)),1) .and.
     .    kbit(avail_met(site(in)),2) .and.
*                                                  ! all mets available
     .    kbit(avail_met(site(in)),3)       ) Then
 
*         For the moment, use the DRY_CfA mapping function (Since dry_CfA
*         is only called if all weather data is available, we do not
*         need to worry about the dry mapping function being changed
*         when this routine is called.
 
          call dry_CfA ( in, wet_map )
 
*                 ! Not enough data, go to CHAO
      ELSE
 
*                                                  ! Turn off CfA-2.2
          call sbit( cont_medium(site(in)),21,0 )
*                                                  ! Turn on  CHAO
          call sbit( cont_medium(site(in)),26,1 )
 
      END IF
 
***** Thats all
      return
      end
 
