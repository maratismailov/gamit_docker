CTITLE JWLOC
 
      subroutine jwloc( unit, num, character, xv, yv )
 
*     Supposed to return the position of the cursor.  Can't be
*     implemented at the 0A level of GKS.  Here we ask user to
*     type in values.
 
* PASSED VARIABLES
 
*   unit        - Graphics work station (ignored)
*   num         - number of characters (ignored)
 
 
      integer*4 unit, num
 
*   xv, yv      - Virtual coordinates of cursor
 
 
      real*4 xv, yv
 
*   character   - Charcter hit by user.
 
 
      character*(*) character
 
      write(*,100)
 100  format(' We cannot sense pen position at 0A level in GKS.',/,
     .       ' Enter the desired virtual coordinates ',
     .       ' (XV, YV each 0-1) ',$)
      read(*,*) xv, yv
 
 
 
      character = ' '
 
***** Thats all
      return
      end
 
