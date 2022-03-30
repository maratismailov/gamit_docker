CTITLE    ...................................................................
 
      subroutine label_axis(ch_x,ch_y,  orient, label)
c
c     Routine to put the axis despcrition on the plot.
c     The cursor is initially positioned on ch_x and ch_y.  The
c     cursor is then diplaced by the amount need to account for
c     the labels on the axis.  (This depends on width and height
c     of the labels and the orientation of the labels which
c     have been output)
c
c Include files
c -------------
*                         ! the  parameter file
      include 'plot_param.h'
c
*                         ! the common block
      include 'plot_com.h'
c
c Variables
c ---------
c ch_x, ch_y -- the center of decription on the axis
c orient -- the orienation of the label
c     orient = 0 -- parrarell to x axis (top justification)
C            = 2 -- parrarell to x axis (bottom justifiction)
c     orient = 1 -- antiparrarell to y axis (right justification)
c              3 -- antiparrarell to y axis (left justification)
c label -- the description to be output
c
      real*4 ch_x, ch_y
 
c
      integer*4 orient
 
c
      character*(*) label
 
c
c Local variables
c ---------------
c dx, dy -- the character displacement converted world coordinates -- not used
c
*   blj     - Bottom left justification
*   trj     - Top/right justication value
      real*4 blj, trj
c     real*4 dx, dy
 
c
c.... Move to correct position
      call s2mov(sign_x*ch_x, sign_y*ch_y)
c
c.... Now move relative to acount for labels and save the orientation
 
*                                               ! X axes
      if( orient.eq.0 .or. orient.eq.1 ) then
*                     ! Top/right justification
          trj = 1.0
*                     ! Bottom/left justification
          blj = 0.0
*                                               ! Y axes
      else
          trj = 0.0
          blj = 1.0
      end if
 
*                                               ! Parrallel to x axis
      if( orient.eq.0 .or. orient.eq. 2 ) then
          if( sign_x*sign_y.gt.0 ) then
              call scj(sign_x*1.0,0.0,  0.5, trj)
*                                                   ! Shift for height
              call sr2mv(wh(2)*lf(1), wh(2)*lf(2))
          else
              call scj(sign_x*1.0,0.0,  0.5, blj)
*                                                   ! Shift for width
              call sr2mv(wh(1)*lf(1), wh(1)*lf(2))
          end if
 
      end if
 
*                                              ! Parrallel to y axis
      if( orient.eq.1 .or. orient.eq.3 ) then
          if( sign_x*sign_y.gt.0 ) then
              call scj(0.0, sign_y*1.0,  0.5, blj)
*                                                   ! Shift for width
              call sr2mv(wh(1)*lf(1), wh(1)*lf(2))
          else
              call scj(0.0, sign_y*1.0,  0.5, trj)
*                                                   ! Shift for height
              call sr2mv(wh(2)*lf(1), wh(2)*lf(2))
          end if
 
      end if
 
C     call conv_mmtow(charsz_x*sign_x*nc_x, 0.0, x_size_mm, dx)
C     call conv_mmtow(0.0, charsz_y*sign_y*nc_y, y_size_mm, dy)
C
C     call jr2mv( dx, dy)
C     call sr2mv( dx, dy)
C
C.... set up orientation and justification of descritption
C     if( orient.eq.0 ) then    ! Parralell to x axis
C        call jcori(sign_x*1.0,0.0,0.0, 0.0,sign_x*1.0,0.0)
C        call jjust( 0.5, 1.0)
C        call scj(sign_x*1.0, 0., 0.5, 1.0)
C     end if
C
C     if( orient.eq.1 ) then    ! parrarell to y axis
C        call jcori(0.0, sign_y*1.0, 0.0, -sign_y*1.0, 0.0,0.0)
C        call jjust( 0.5, 0.0)
C        call scj(0.0, sign_y*1.0, 0.5, 0.0)
C     end if
C
c.... Now we can output the label
 
      call save_label(label)
c
      return
      end
 
