CTITLE DRAW_MAP
 
      subroutine draw_map
 
*     This routine draws a map based on the parameters set in the
*     map_set, map_proj commands.
 
      include 'plot_param.h'
      include 'plot_com.h'
 
***** CAll the EZMAP routines needed for making map.
      call maproj( map_proj, ppos(1), ppos(2), rota)
      call mapset( map_limit, plim(1,1), plim(1,2), plim(1,3),
     .             plim(1,4))
      call mappos( view_size(1), view_size(2), view_size(3),
     .             view_size(4) )
      call mapstc( 'OU', map_ou )
      call mapstr( 'GR', map_grid_space)
      if( map_grid_space.gt.0 ) then
          call mapstr( 'GD', map_grid_space/10.0 )
      end if
      call mapstl( 'LA', map_la)
      if( map_sa.gt.0 ) then
          call mapstr( 'SA', map_sa )
      end if
 
*     Initialze the map
      call mapint
 
*     Draw map
      call mapdrw
 
****  Thats all
      return
      end
 
 
 
