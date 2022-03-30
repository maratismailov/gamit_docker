CTITLE jtexH
 
      subroutine jtexh( ilen, text )
 
*     Emulates writing high qulaity text using PWRITX.  We first need
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
*   ksze        - Size rescaled for cartographic.
 
      integer*4 ksze
 
      real*4 px, py
 
***** Get position of the text, accounting for height justification.
 
      call jsett( px, py )
 
*     Now use PWRITX
*     Allow for rescaling of small cartographics font so that
*     it will be about the same as Principle size.
 
      ksze = jsze
      if( jfnt.eq.1 ) ksze = 2*jsze
 
*     see if we should set the font
      if( .not. jfntset ) then
          if( jfnt.eq.1 ) then
              call pwritx(px,py, '''K''',3, ksze, 0,0)
          else
              Call pwritx(0.,0.,'''P''', 3, ksze, 0,0)
          end if
          jfntset = .true.
      end if
 
      call pwritx(px, py, text, ilen, KSZE, JOR_deg, JCENTER)
 
*     Thats all
      return
      end
 
