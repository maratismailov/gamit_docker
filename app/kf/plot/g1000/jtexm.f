CTITLE jtexM
 
      subroutine jtexm( ilen, text )
 
*     Emulates writing medium qulaity text using PWRITY.  We first need
*     to retreive the current pen position and then we use the JCENTER
*     and JOR_deg values saved from calls to jcori and jjust.
*
      include 'g1000.h'
 
* PASSED VARIABLES
 
*   ilen        - Length of the string to be written
 
 
      integer*4 ilen
 
*   text        - String to be written.
 
 
      character*(*) text
 
* LOCAL VARIABLES
 
*   px, py      - Pen poistion in user coordinates
 
 
      real*4 px, py
 
***** Get position of the text, accounting for height justification.
 
      call jsett( px, py )
 
*     Now use PWRITY
 
      call pwrity(px, py, text, ilen, JSZE, JOR_deg, JCENTER)
 
*     Thats all
      return
      end
 
