CTITLE jsett
 
      subroutine jsett( px, py )
 
*     This rotuine will set the postion of text so as to emulate
*     height justification in strings.
 
      include 'g1000.h'
 
* PASSED VARIABLES
 
*   px, py      - Positions needded for the text to be justified
*               - correctly in height.
 
 
      real*4 px, py
 
* LOCAL VARIABLES
 
*   ix, iy    - Current plotter coordinates of pen
*   jx, jy    - Poistioned need to get height justifcation
 
      integer*4 ix, iy, jx, jy
 
*   ang         - Orientation in rads (used to get the correct
*               - direction to move.
*   cpux, cpuy  - Conversion from plotter to user coordinates
*   dxdh, dydh  - Change in x and y for change in height
 
      real*4 ang, cpux, cpuy, dxdh, dydh
 
*     Get current position of cursor
 
      call mxmy ( ix, iy )
 
*     Now based on current orientation get transformattion from heigt
*     to x and y
 
      ang = JOR_deg / 57.29577951
 
      dxdh = -sin(ang)
      dydh = -cos(ang)
 
*     Get the new position in plotter coords.
      jx = ix + dxdh*gheight*jsze
      jy = iy + dydh*gheight*jsze
 
*     Now convert to user coordinates
      px = cpux(jx)
      py = cpuy(jy)
 
*     Thats all
      return
      end
 
 
