CTITLE JCOLR
 
      subroutine jcolr( pen )
 
*     Routine to set the color.  Uses the SETUSV routine.  Form of the
*     pen is  Intensity|red|green|Blue where each value may take on
*     value between 0-9.
*
 
* PASSED VARIABLES
 
*   pen     - Pen to use (in sense that it defines colors)
 
 
      integer*4 pen
 
* LOCAL VARIABLES
 
*   ir, ig, ib  - Red, green, blue intensity
*   in          - Intensity to use scaled to max of 10000.
*   intensity   - Intensity of color.
*   nt          - Next entry in color table
 
 
      integer*4 ir, ig, ib, in, intensity, nt
      real*4 r,g,b
 
      data nt / 0 /
 
*     Split up the values
 
      intensity = pen/1000
      ir = (pen - intensity*1000)/100
      ig = (pen - intensity*1000 - ir*100)/10
      ib = mod(pen,10)
      in = min(intensity*1000,10000)

****  Set colors
      r = ir/9.01
      g = ig/9.01
      b = ib/9.01
 
*     Flush first
      call sflush
 
*     Now use the set routine
 
c     call setusv('IR', ir)
c     call setusv('IG', ig)
c     call setusv('IB', ib)
c
c     call setusv('IN', in )
 
*     Superceded by the new COLOR commands
      write(*,200) nt, r,g,b
 200  format(' Pen ',i3,' Colors ',3f5.2)
      call gscr( 1, nt, r, g, b )
      nt = nt + 1
 
*     Now set color for everything
      call gsplci(nt)
      call gspmci(nt)
      call gstxci(nt)
*     The following could be added, but we are not areas.
*     call gsfaci(1)
 
***** Thats all
      return
      end
 
