CTITLE J2MOV
 
      subroutine j2mov( xp, yp )
 
*     Moves the pen to user coordinates xp, yp.  Uses the SPPS
*     routine FRSTPT.
 
* PASSED VARIABLES
 
*   xp, yp      - User coorinates to move to.
 
 
      real*4 xp, yp
 
*     Directly equivalent call.
 
      call frstpt( xp, yp )
 
*     Thats all
      return
      end
 
