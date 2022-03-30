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
 
 
      integer*4 ir, ig, ib, in, intensity, nt, i

      real*4 rgbv(3,20)

      logical first_call

      data rgbv/  0.00,  0.00,  0.00 ,  !  1
     +            0.10 , 0.00 , 0.90 ,  !  2
     +            0.00 , 0.00 , 0.96 ,  !  3
     +            0.00 , 0.25 , 1.00 ,  !  4
     +            0.00 , 0.75 , 1.00 ,  !  5
     +            0.00 , 0.98 , 1.00 ,  !  6
     +            0.00 , 1.00 , 0.04 ,  !  7
     +            0.65 , 1.00 , 0.00 ,  !  8
     +            0.99 , 1.00 , 0.10 ,  !  9
     +            1.00 , 0.75 , 0.00 ,  ! 10
     +            0.79 , 0.71 , 0.00 ,  ! 11
     +            0.77 , 0.59 , 0.00 ,  ! 12
     +            0.99 , 0.00 , 0.00 ,  ! 13
     +            1.00 , 0.65 , 0.50 ,  ! 14
     +            0.99 , 0.00 , 0.99 ,  ! 15
     +            0.65 , 0.00 , 0.99 ,  ! 16
     +            0.65 , 0.30 , 0.75 ,  ! 17
     +            0.30 , 0.65 , 0.75 ,  ! 18
     +            0.75 , 0.65 , 0.30 ,  ! 19
     +            1.00 , 1.00 , 1.00 /  ! 20

      data first_call / .true. /
 
      call sflush

***** If this is the first call set up the color table
      if( first_call ) then
          call gscr(1,0,rgbv(1,pen),rgbv(2,pen),rgbv(3,pen))
c
          do i=1,20
             call gscr(1,i,rgbv(1,i),rgbv(2,i),rgbv(3,i))
          end do
          first_call = .false.
      end if

*     Now use the set routine
 
*     Superceded by the new COLOR commands
*     nt = nt + 1
*     call gscr( 1, nt, ir, ig, ib )
 
*     Now set color for everything
      nt = pen
      if( nt.gt.20 ) nt = 20
      call gsplci(nt)
      call gspmci(nt)
      call gstxci(nt)
*     The following could be added, but we are not areas.
*     call gsfaci(1)
 
***** Thats all
      return
      end
 
