CTITLE JCSIZ
 
      subroutine jcsiz ( width, height, gap )
 
*     Sets the character size.  Originally these were in world
*     coordinates (in G1000).  Here we need size in plotter units
*     therefore size will be set this way.  Width is used as the value
*     height and gap are ignored.
 
      include 'g1000.h'
 
* PASSED VARIABLES
 
*   width   - Character size to use . Value assumed in mm and is
*             converted to plotter units assuming 200 mm plot size.
*            (if value less than 0.8 mm then standard 0-3 values for
*             size will produce large letters)
*   height, gap - These were height and gap and are ignored
*   dwdh    - Ratio of world to mm on page
*   hwd     - Height of characters in world coods to use.

 
      real*4 width, height, gap, dwdh, hwd
 
      JSZE = nint(width/0.195)
      dwdh = (gwt-gwb)/(gvt-gvb)/180.0
c     hwd  = dwdh*width

*     Although documentation says this is world coordinates
*     it looks like virtual coordinates on the screen (i.e.1e-3
*     is a good size).
      hwd  = width/300.0
 
*     Set marker size as well
      call setusv('MS', JSZE*60 )
      call gschh( hwd )

*     Thats all
      return
      end
 
