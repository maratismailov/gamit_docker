CTITLE  ................................................................
 
      subroutine plts3

      implicit none 
c
c     This program does all of the high quality labelling and
c     axis drawing work.  Any routine used JTEXH should be called
c     from with in this segment
c
c Include files
c -------------
*                          ! the parameter file
      include 'plot_param.h'
c
*                          ! the common block
      include 'plot_com.h'
c
      character*256 copy_buf
c
c.... Go to the routine specified by pel
      if( pel.eq.10) then
         call dxaxis_mn
      end if
c
      if( pel.eq.11 ) then
         call dxaxis_mx
      end if
c
      if( pel.eq.12 ) then
         call dyaxis_mn
      end if
c
      if( pel.eq.13 ) then
         call dyaxis_mx
      end if
c
      if( pel.eq.17 ) then
         call put_label
      end if
c
*                           ! draw all axes
      if( pel.eq.28 ) then

* MOD TAH 131021: See if option not label some axes
         copy_buf = buffer
         call casefold(copy_buf)
c
c....    Set the buffer for defaults tics and spacing
         if( index(copy_buf,'-X').gt.0 ) then
             buffer(9: ) = ' -1 0 ' 
         else
             buffer(9: ) = ' -1 1 ' // xaxis_label
         endif 
         call dxaxis_mn
*        Y MIN AXIS
         if( index(copy_buf,'-Y').gt.0 ) then
             buffer(9: ) = ' -1 0 ' 
         else
             buffer(9: ) = ' -1 1 ' // yaxis_label
         end if
         call dyaxis_mn
c
c....    Now do rest of axes
         buffer(9: ) = ' -1 0 '
         call dxaxis_mx
         buffer(9: ) = ' -1 0 '
         call dyaxis_mx
c
      end if
c
c.... Thats all
*                   ! flush the plot buffers
      call jmcur
*                   ! return to main segment
      return
c
      end
 
