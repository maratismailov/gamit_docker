 
CTITLE GEN_ERR_ELLIPSE
 
      subroutine gen_err_ellipse ( sig_x, sig_y, rho, confid, xp, yp,
     .                        nump, errb_scale )
 
*     This routine generates an error ellipse of condfindence (confid)
*     sigmas in the x and y directions of sig_x and sig_y, correlation
*     rho.  The nump points in the ellipse are retured in xp and yp.
*
      include '../includes/const_param.h' 
 
* PASSED VARABLES
 
*   nump        - Number of points in the ellipse (supplied
*               - by user)
 
      integer*4 nump
 
*   sig_x, sig_y    - Sigmas in the x and y dirrections.
*   rho         - Correlation coefficient between x and y
*   confid      - Confidence interval wanted (0-1)
*   xp(nump), yp(nump)  - Points in the ellipse.  These values
*               - are relative to the end of the velocity
*               - vector. (When mapping we need an additional
*               - transformation, so we don't do that yet)
*   errb_scale(2) - Scaling for the error bars in x and y directions
 
      real*4 sig_x, sig_y, rho, confid, xp(nump), yp(nump),
     .       errb_scale(2)
 
* LOCAL VARIABLES
 
*   cov(2,2)        - Covariance matrix of data
*   a,b,c       - Constants used in getting eigenvalues
*   eigen1, eigen2  - The two eigenvalues
*   ang         - Orientation of ellipse relative to X
*               - axis (rads)
*   rot(2,2)        - Rotation matrix from XY to eignenvectors
*   conrad      - Radius for the confidence interval
*   xi,yi       - Point on the ellipse in the local frame
*               - orientated along eigen vector
 
      real*8 cov(2,2), a,b,c, eigen1, eigen2, ang, rot(2,2), conrad,
     .    xi,yi
 
*   i           - Counter in the xp and yp arrays.
 
      integer*4 i
 
****  Get the radius we need for correct confidence interval. Make
*     value is of confid is correct.
 
      if( confid.gt.1 .or. confid.le.0 ) then
          write(*,'(a,f6.2,a)') 
     .     ' INVALID Confidence interval in GEN_ERR_ELLIPSE (',confid
     .    ,'--value should be 0-1)'
          RETURN
      end if
 
      conrad = sqrt(-2.0*log(1.0-confid))
 
****  Generate covarinace matrix
      cov(1,1) = (sig_x*errb_scale(1))**2
      cov(2,2) = (sig_y*errb_scale(2))**2
      cov(1,2) = rho*sig_x*sig_y*errb_scale(1)*errb_scale(2)
      cov(2,1) = cov(1,2)
 
*     Now compute eigenvalues
      a = 1
      b = -(cov(1,1)+cov(2,2))
      c = (cov(1,1)*cov(2,2) - cov(1,2)**2)
 
*     Include confidence interval scaling here.  Treat the case
*     of zero correlation specially.
      if( abs(cov(1,2)).gt.1.d-5 ) then
          eigen1 = sqrt((-b - sqrt(b**2-4.d0*a*c))/(2.d0*a))* conrad
          eigen2 = sqrt((-b + sqrt(b**2-4.d0*a*c))/(2.d0*a))* conrad
 
*         Now get angle
          ang = atan(cov(1,2)/(cov(2,2)-(eigen1/conrad)**2))
      else
          eigen1 = sig_x*conrad
          eigen2 = sig_y*conrad
          ang    = 0.d0
      end if
 
*     Get the rotation matrix from XY to eigenvectors
      rot(1,1) = cos(ang)
      rot(1,2) = sin(ang)
      rot(2,1) = -rot(1,2)
      rot(2,2) =  rot(1,1)
 
*     Now loop forming the ellipse
 
      do i = 1, nump
          ang = 2.d0*pi * (i-1) / (nump-1)
 
          xi = eigen1*cos(ang)
          yi = sqrt(1.d0-(xi/eigen1)**2)*eigen2
          if( ang .lt. pi ) then
              xp(i) = rot(1,1)*xi + rot(1,2)*yi
              yp(i) = rot(2,1)*xi + rot(2,2)*yi
          else
              xp(i) = rot(1,1)*xi - rot(1,2)*yi
              yp(i) = rot(2,1)*xi - rot(2,2)*yi
          end if
*                     ! Looping over all of the angles
      end do
 
***** Thats all
      return
      end
 
CTITLE MAKE_ARROW
 
      subroutine make_arrow( dx, dy, arrow_len, dxp, dyp)
 
*     This routine makes and arrow returns the results in the
*     five vector dxp, dyp.  Arrow length in is same units as
*     dx,dy.  The first point in the vector is 0,0.
 
* PASSED VARIABLES
 
*   dx, dy  - Differential postion of arrow head
*   arrow_len   - Length of the arrow head (same units as
*           - dx,dy
*   dxp(5), dyp(5)  - Five points making at the arrow.
 
      real*4 dx, dy, arrow_len, dxp(5), dyp(5)
 
* LOCAL VARIABLES
 
*   i       - Loop counter
 
      integer*4 i
 
*   ang     - Angle made by dx, dy.  (used to get rotation
*           - matrix)
*   veclen  - Length of the vector
*   rot(4,4)    - Rotation matrix from xy to arrow dirrection
*   tx,ty   - Temporary storage of x and y position during
*           - rotation.
 
 
      real*8 ang, veclen, rot(4,4), tx,ty
 
***** Get the length and orientation of arrow
 
      veclen = sqrt(dx**2 + dy**2 )
*                                 ! Clear the arrays and get out
      if( veclen.eq.0.d0 ) then
          do i = 1,5
              dxp(i) = 0.
              dyp(i) = 0.
          end do
          RETURN
      end if
 
*     See if atan2 will overflow
      if( abs(dx/dy).lt.1.d-6 ) then
          dx = 1.d-6*dy
      end if
 
      ang = atan2(dy,dx)
 
*     Generate the rotation matrix
      rot(1,1) = cos(ang)
      rot(1,2) = -sin(ang)
      rot(2,1) = -rot(1,2)
      rot(2,2) = rot(1,1)
 
***** Generate the points in the arrow in a x directed frame
      dxp(1) = 0
      dyp(1) = 0
 
      dxp(2) = veclen
      dyp(2) = 0
 
      dxp(3) = dxp(2) - arrow_len
      dyp(3) = -arrow_len/4
 
      dxp(4) = dxp(3)
      dyp(4) = arrow_len/4
 
      dxp(5) = dxp(2)
      dyp(5) = dyp(2)
 
****  Now transform the vector into the correct frame
      do i = 1, 5
          tx = rot(1,1)*dxp(i) + rot(1,2)*dyp(i)
          ty = rot(2,1)*dxp(i) + rot(2,2)*dyp(i)
          dxp(i) = tx
          dyp(i) = ty
      end do
 
****  Thats all
      return
      end
 
CTITLE REPORT_MAP
 
      subroutine report_map
 
*     This routine reports the current values of the map parameters
*
 
      include 'plot_param.h'
      include 'plot_com.h'
 
* LOCAL VARIABLES
 
*   i           - Loop counter
*   proj_num        - Number of projection we are using
*   limit_num   - numner of limit we are using
 
      integer*4 i,j, proj_num, limit_num
 
*   proj_names(10)  - Full names of the projections
 
      character*32 proj_names(10)
 
*   limit_names(5)  - Full names for limit types

      character*8 limit_names(5)
 
*   proj_types(10)  - Codes for types of projections
*   limit_types(5)  - Codes for limit types
 
      character*2 proj_types(10), limit_types(5)
 
      data proj_types / 'LC', 'ST', 'OR', 'LE', 'GN', 'AE', 'SV',
     .                'CE', 'ME', 'MO' /
 
      data limit_types / 'MA', 'CO', 'PO', 'AN', 'LI' /
 
      data proj_names /
     .        'Lambert Conformal'
     .,       'STereographic'
     .,       'ORthographic'
     .,       'Lambert Equal area'
     .,       'GNomonic'
     .,       'Azimuthal Equidistant'
     .,       'Satellite-View'
     .,       'Cylindrical Equidistant'
     .,       'MErcator'
     .,       'MOllweide-type'            /
 
      data limit_names /
     .        'MAximum'
     .,       'COrners'
     .,       'POints'
     .,       'ANgles'
     .,       'LImits'  /
 
      write(*,100)
 100  format(' MAPPING PARAMETERS')
 
*     Find the projection we are using
      proj_num = -1
      do i = 1, 10
          if( proj_types(i).eq.map_proj ) then
              proj_num = i
          end if
      end do
 
*     See if we found it
*                                 ! Did not find it
      if( proj_num.eq.-1 ) then
          write(*,200) (proj_types(i), proj_names(i), i = 1,10)
  200     format('UNKNOWN MAP PROJECTION.  The choices are:',/,
     .           ' Code   Type ',
     .         10(/,2x,a2,4x,a))
      else
          write(*,220) Proj_names(proj_num)(1:23), ppos, rota
  220      format(' PROJ ',a,' LAT ',f6.2,' LONG ',f7.2,
     .           ' ROTA ',f7.2,' (deg)')
      end if
 
*     Write out the limits set.  Find the projection we are using
      limit_num = -1
      do i = 1, 5
          if( limit_types(i).eq.map_limit ) then
              limit_num = i
          end if
      end do
 
*     See if we found it
*                                     ! Did not find it
      if( limit_num.eq.-1 ) then
          write(*,300) (limit_types(i), limit_names(i), i = 1,10)
  300     format('UNKNOWN LIMIT TYPE.  The choices are:',/,
     .           ' Code   Type ',
     .         10(/,2x,a2,4x,a))
*                                 ! Yes we found.  Write limits in
*                                 ! appropriate format
      else
          if( map_limit.eq.'MA' ) then
 
              write(*,320) limit_names(limit_num)(1:8)
  320         format(' LIMIT ',a,' Full area shown')
 
          else if( map_limit.eq.'LI' ) then
 
              write(*,340) limit_names(limit_num)(1:8),
     .            ((plim(i,j),j=1,4), i=1,2)
  340         format(' LIMIT ',a,' LATS  (deg) ',4(f7.2,1x),/,
     .                15x, ' LONGS (deg)', 4(f7.2,1x))
 
*                     ! All the rest just have 4 arguments
          else
 
              write(*,360) limit_names(limit_num)(1:8),
     .            (plim(1,j),j=1,4)
  360         format(' LIMIT ',a,' Arguments (deg) ',4(f7.2,1x))
 
          end if
      end if
 
****  Write out the rest of the paramters
      write(*,400) map_ou, map_grid_space, map_sa
 400  format(' Outline',t20,'Grid Space',t40,'Satellite Alt.',/,
     .       4x,a2,t20,1x,f4.1,' deg',t40,f6.2,' Earth radii',$)
 
      if( map_la ) then
          write(*,420) 'Labels'
      else
          write(*,420) 'No Labels'
      end if
 420  format(10x,a)
 
****  Thats all return
      return
      end
 
CTITLE DECODE_MDRAW
 
      subroutine decode_mdraw( gbuffer, indx )
 
*     This routine decodes the various parameters which can be
*     set in the MAP_DRAW command.  These are the values whicj
*     are set with MAPSTx command.  We will decode them here
*     but they will set just befoe the map is drawn
 
      include 'plot_param.h'
      include 'plot_com.h'
 
* PASSED VARIABLES
 
*   indx    - Current pointer to position in buffer
 
      integer*4 indx
 
*   gbuffer - The buffer read from the command file
 
      character*(*) gbuffer
 
* LOCAL VARIABLES
 
*   jndx,kndx   - Working copy of indx
*   i           - dummy
*   ierr,jerr   - Error return from read_line for next
*               - word and  real*4 extraction
*               - (loop until ierr non zero)
 
      integer*4 jndx,kndx, i, ierr,jerr, trimlen
 
*   next_word   - Next word from gbuffer
*   cdum        - Dummy character string for read_line
 
      character*20 next_word, cdum
 
****  Initialize and start getting words until end
      ierr = 0
      jndx = indx
 
      do while ( ierr.eq.0 )
          call read_line(gbuffer,jndx,'CH',ierr, i, next_word)
 
*         See which parameter it is.
*                                     ! Decode option
          if( ierr.eq.0 ) then
              call casefold(next_word)
              kndx = jndx
 
*                                                     ! Outline next
              if( next_word(1:2).eq.'OU' ) then
                  call read_line(gbuffer,jndx,'CH',ierr, i, map_ou)
                  call casefold(map_ou)
 
*                                                     ! Grid spacing
              else if( next_word(1:2).eq.'GR' ) then
                  call read_line(gbuffer,jndx,'R4',jerr,
     .                 map_grid_space, cdum )
                  call report_error('IOSTAT',jerr,'decod',
     .                 gbuffer(kndx:),0,'decode_mdraw/Gridspace')
 
*                                                     ! Satellite Alt.
              else if( next_word(1:2).eq.'SA' ) then
                  call read_line(gbuffer,jndx,'R4',jerr,map_sa, cdum)
                  call report_error('IOSTAT',jerr,'decod',
     .                 gbuffer(kndx:),0,'decode_mdraw/Sat.Alt.')
 
*                                                     ! Labels
              else if( next_word(1:2).eq.'LA' ) then
                  call read_line(gbuffer,jndx,'CH',ierr,i, cdum)
                  call casefold(cdum)
                  if( cdum(1:1).eq.'Y' ) then
                      map_la = .true.
                  else
                      map_la = .false.
                  end if
              else
                  write(*,200) gbuffer(kndx:trimlen(gbuffer))
 200              format(' Could not decode > ',a,
     .                   ' in decode_mdraw')
              end if
*                     ! No error on getting next word
          end if
*                     ! Looping until end
      end do
 
****  Thats all
      return
      end
 
CTITLE SET_MAP
 
      subroutine set_map ( type )
 
*     Routine to tell PLOT and G1000 emulation which mode we are in
*     mapping (type <> 0) or regular (type = 0 )
 
      include 'plot_param.h'
      include 'plot_com.h'
 
*   type        - 0 for regular, <>1 for mapping
 
      integer*4 type
 
*     Test value of type and set plot mapping mode variable
      if( type.eq.0 ) then
          map_mode = .false.
      else
          map_mode = .true.
      end if
 
*     Now tell G1000 emulation (needed for drawing and labels)
      call jmapm( map_mode )
 
****  Thats all
      return
      end
 
CTITLE VDRAW
 
      subroutine vdraw( gbuffer, xp, yp, dxp, dyp, num_ellipse )
 
*     Routine to draw vectors on a plot.  Either a regular figure
*     or a map projection.  If a map projection is to be used, then
*     the map should be draw first to ensure that all of the needed
*     parameters are set.
*     The Error bar type is used to control the drawing of error
*     ellipses, and the point type is used to generate the point at
*     base of the vector.  The line type is used for the lines if
*     its value is greater than zero, otherwise a solid line is used.
*
      include 'plot_param.h'
      include 'plot_com.h'
 
* PASSED VALIABLES
 
*   num_ellipse - Number of line segments in the error
*               - ellipse
 
      integer*4 num_ellipse
 
*   xp(num_ellipse), yp(num_ellipse)    - The x and y
*               - cooridinates of the points in the error
*               - ellipse.  (UV values if in map mode)
*   dxp(num_ellipse), dyp(num_ellipse)  - The differential x and
*               - y cooridinates of the points in the error
*               - ellipse and error.  (UV values if in map mode)
 
      real*4 xp(num_ellipse), yp(num_ellipse),
     .    dxp(num_ellipse), dyp(num_ellipse)
 
*   gbuffer     - The command line.  (To get the scale and
*               - confidence interval (for error ellipses) for
*               - velocity vectors).
 
 
 
      character*(*) gbuffer
 
* LOCAL VARIABLES
 
*   i       - Loop index
*   symbol  - Point type to used.
 
      integer*4 i, symbol, trimlen
 
*   su, sv      - Postion of current point in u/x and v/y
*               - coordinates
*   en_to_mm(2,2)   - Transformation of e and n velcoties to
*               - mm on page (scaled by vscale)
*   arrow_length    - Length of the arrow head in east north
*               - coordinates
*   eu, ev      - End of the arrow
 
      real*4 su, sv, en_to_mm(2,2), arrow_length, eu, ev
 
****  Finish decoding the command line.  Do now so that values entered will
*     not be lost if there is no data.
      call decode_vdraw( gbuffer, vscale, vconfid )
 
****  Make sure we have some data
      if( num_velvec.eq.0 ) then
          write(*,100)
 100      format('VDRAW Error: No data read yet, use VREAD command')
          RETURN
      end if
 
*     If we are not in map mode set view and scales
      if( .not. map_mode ) then
          call set_view
          call set_scale
      end if
 
****  Now loop over the data
      do i = 1, num_velvec
 
*         Get transformation from e and n velocotity to mm on the
*         page scaled by the velocity scaling
 
          call get_map_tran( vpos(1,i), vpos(2,i),  su, sv,
     .                     en_to_mm, arrow_length )
 
*         Make the vector for the arrow
          call make_arrow( vvel(1,i), vvel(2,i), arrow_length,
     .                    dxp, dyp )
 
*         Convert the arrow into absolute cooridates and plot after
*         setting line style.
 
          call dxy_to_abs( su, sv, dxp, dyp, xp, yp, 5, en_to_mm )
 
*****     Now draw the line
          call j2ply( 5, xp, yp )

*         Save the end of the arrow
          eu = xp(2)
          ev = yp(2)
 
****      If point type is set plot the point type at start location
          if( point_type.gt.0 ) then
              call j2mrk( su, sv, -point_type )
*                                                 ! Plot character
          else if( point_type.eq.-1 ) then
              symbol = (vpnt(i)+64)
              call j2mrk( su, sv, symbol )
c             call points( su, sv, 1, symbol, 0 )
*                                                 ! Plotter symbol
          else if( point_type.eq.-2 ) then
              call j2mrk( su, sv, -vpnt(i) )
          end if
 
****      If the error bar type is set, then draw the error ellipse
          if( errb_type.gt.0 ) then
 
*             Generate the error ellipse.
              call gen_err_ellipse( vsig(1,i), vsig(2,i), vrho(i),
     .                vconfid, dxp, dyp, num_ellipse,errb_scale )
 
*             Now transform values.
              call dxy_to_abs( eu, ev, dxp, dyp, xp, yp, num_ellipse,
     .                        en_to_mm )
 
*             Now draw the ellipse
              call j2ply( num_ellipse, xp, yp )
          end if
 
*****     Finally see if we should label.  Tell the user where we are
*         placing it so that it can moved to a better spot later.
*                                 ! Label to right and centered in
*                                 ! height
          call jjust( 0.0, 0.5 )
          if( trimlen( vname(i) ).gt.0 ) then

*             Tell user where we are writing
c             write(*,300) vname(i), su, sv
c300          format(' Writing ',a,' X,Y ',2f10.4)
              call j2mov( su, sv )
              call write_label ( vname(i) )
          end if
 
*                     ! Looping over all of the vectors
      end do
 
***** Thats all
      return
      end
 
CTITLE READ_VFILE
 
      subroutine read_vfile( data_array )
 
*     This routine will loop over the data records read from the
*     velocity file and extract the information specified in the
*     v_field command.
 
      include 'plot_param.h'
      include 'plot_com.h'
 
* PASSED VARAIBLES
 
*   data_array(int_recl, num_epochs )   - This is an copy
*               - of the data records from the data file.
*               - (There are num_epochs of them)
 
      integer*4 data_array(int_recl, num_epochs )
 
* LOCAL VARIABLES
 
 
*   i,j         - Loop counters
*   idata(int_recl) - One line from data_array
*   largest_field   - Largest of the fields required.
*   id          - Dummy for readline
 
      integer*4 i,j, idata(int_recl), largest_field, id, indx,
     .    next_val_indx, ierr
 
*   values(50)  - Array to read data into before we set into
*               - correct arrays.
*   max_vec     - Length of the longest vector
*   vel_mag     - Length of vector for setting scales
*   dx, dy      - Length of x and y axes
 
 
      real*4 values(50), vel_mag, max_vec, dx, dy
 
*   cdata       - Data record as character string.
 
      character*(char_recl) cdata
 
*   next_word   - Word with the name of the site
 
      character*20 next_word
 
*   cd          - Dummy for readline
 
      character*4 cd

      equivalence ( cdata, idata )

****  Set to handle default scaling
      scale_set = .false.
      if( reset_scales ) then
         use_def_scale = .true.
      else
         use_def_scale = .false.
      end if

*     Set initial values
      do i = 1,3,2
         default_scale(i) = 1.d20
         default_scale(i+1) = -1.d20
      end do

*     See if file has been read
      if( .not.file_read ) then
          write(*,50) 
  50      format('*** ERROR IN VREAD *** No file given yet. Use the ',
     .           ' FILE command')
          do i = 1,3,2
             default_scale(i) = 0     
             default_scale(i+1) = 1
          end do
          return
      end if

 
****  Scan v_field and get the largest value of the reals
*     required
 
      largest_field = 0
*                     ! Skip 8 which is the name field.
      do i = 1, 7
          largest_field = max(largest_field, v_field(i))
          if( v_field(i).eq.0 ) v_field(i) = 50
      end do
 
      largest_field = max(largest_field, p_field(1))
 
      if( largest_field.gt.50 ) then
          write(*,100)
  100     format(' ** ERROR IN VREAD_DATA** Largest field in',
     .           ' V_FIELD too large.',/, 10x,
     .           ' Largest value allowed is 50.' )
          num_velvec = 0
          RETURN
      end if
 
****  Loop over all of the data records
      num_velvec = 0
      ref_valx = 0.0
      ref_valy = 0.0
      max_vec  = 0.0
      values(50) = 0.0
      
      do i = 1, num_epochs
 
*         Copy the data record to local storage
          do j = 1, int_recl
              idata(j) = data_array(j,i)
          end do

          next_word = cdata 
 
          num_velvec = num_velvec + 1
 
*         Now pull off the values.  Get all the reals upto the max
*         required and save values
 
*                             ! Character number in string
          indx = 1
*                             ! Will be reset if we read up to name
          next_val_indx = 1
 
*                                     ! Get reals up to this point
          if ( v_field(8).gt.0 ) then
              call multiread(cdata, indx, 'R4', ierr, values, cd,
     .                    v_field(8)-1)
              next_val_indx = v_field(8)
*             Get the site name
              call read_line(cdata, indx, 'CH', ierr, id, next_word)
 
*             Prepend a blank (to move name of point) and save
              vname(num_velvec) = ' ' // next_word
          end if
 
*         Now get the rest of reals
          call multiread(cdata, indx, 'R4', ierr,
     .        values( next_val_indx ), cd, largest_field-v_field(8))
 
*         Now save the reals
*                                     ! Save values, otherwise skip line
          if( ierr.eq.0 ) then
              vpos(1,num_velvec) = values(v_field(1))
              vpos(2,num_velvec) = values(v_field(2))
              vvel(1,num_velvec) = values(v_field(3))
              vvel(2,num_velvec) = values(v_field(4))
              vsig(1,num_velvec) = values(v_field(5))
              vsig(2,num_velvec) = values(v_field(6))
              vrho(num_velvec)   = values(v_field(7))
 
*****         Now get the point field (if set)
              if( p_field(1).gt.0 ) then
                  vpnt(num_velvec) = nint(values(p_field(1)))
              end if

*             check the scales
              default_scale(1) = min(vpos(1,i), default_scale(1))
              default_scale(2) = max(vpos(1,i), default_scale(2))
              default_scale(3) = min(vpos(2,i), default_scale(3))
              default_scale(4) = max(vpos(2,i), default_scale(4))

*             Get vector magnitude
              vel_mag = sqrt(vvel(1,i)**2 + vvel(2,i)**2)
              max_vec = max(vel_mag, max_vec)
          else
              call report_error('IOSTAT',ierr,'vread',cdata,0,
     .                          'read_vfile')
              num_velvec = num_velvec - 1
          end if
*                     ! Looping over data records
      end do

*     Now at end set the default scales.
      dx = default_scale(2) - default_scale(1)
      if( dx.ne.0 ) then
          default_scale(1) = default_scale(1) - 0.2*dx
          default_scale(2) = default_scale(2) + 0.2*dx
      else
          default_scale(1) = -max_vec * 1.2
          default_scale(2) =  max_vec * 1.2
      end if
      dy = default_scale(4) - default_scale(3)
      if( dy.ne.0 ) then
          default_scale(3) = default_scale(3) - 0.2*dy
          default_scale(4) = default_scale(4) + 0.2*dy
      else
          default_scale(3) = -max_vec * 1.2
          default_scale(4) =  max_vec * 1.2
      end if
 
***** Thats all
      return
      end
 
CTITLE DXY_TO_ABS
 
      subroutine dxy_to_abs( su, sv, dxp, dyp, xp, yp, num, en_to_mm)
 
*     Routine to add values computed relative to zero as origin to
*     to start point while applying the transformation to bring
*     e and n into uv plotting coordinate system,
 
* PASSED VARIABLES
 
*   num     - Number of points to be transformed
 
      integer*4 num
 
*   su, sv              - Start position
*   dxp(num), dyp(num)  - Differential positions
*   xp(num), yp(num)        - Absolute positions
*   en_to_mm(2,2)       - Transformation for e and n to uv
 
      real*4 su, sv, dxp(num), dyp(num), xp(num), yp(num),
     .    en_to_mm(2,2)
 
* LOCAL VARIABLES
 
*   i,j         - Loop counters
 
      integer*4 i
 
*   du,dv       - Change in u and v coordinates
 
      real*4 du,dv
 
****  Loop over points transforming dxp, dyp and adding to origin.
 
      do i = 1, num
 
*         Get the change in u and v from the dxp dyp
          du = en_to_mm(1,1)*dxp(i) + en_to_mm(1,2)*dyp(i)
          dv = en_to_mm(2,1)*dxp(i) + en_to_mm(2,2)*dyp(i)
 
*         Now add to origin
          xp(i) = su + du
          yp(i) = sv + dv
 
      end do
 
****  Thats all
      return
      end
 
CTITLE DECODE_VDRAW
 
      subroutine decode_vdraw( gbuffer, vscale, vconfid )
 
*     Routine to decode the remainder of the vdraw command line.
*     Extracts scale and confidience if they are given.
 
* PASSED VARIABLES
 
*   vscale  - Scale from units of vectors to mm on the page
*           - (for 200 mm wide page output)
*   vconfid  - Confidence interval (default is 0.95, if
*           - value passed in > 1 and <100 then converted to
*           - 0 to 1 value.
 
      real*4 vscale, vconfid
 
*   gbuffer - Command line read from user
 
      character*(*) gbuffer
 
* LOCAL VARIABLES
 
*   indx, jndx  - Pointers for reading through string.
*               - Arguments start at indx=9
*   id          - Dummy argument for read_line
 
      integer*4 indx, jndx, id, ierr
 
*   value       - Value read from string (used to stop
*               - over writing default)
 
      real*4 value
 
*   next_word   - Next word. Words pulled from buffer, to
*               - extract scale and confidence intervals.
 
      character*20 next_word
 
*   cd      - Dummy string for read_line
 
      character*4 cd
 
****  Get the scale argument from the string
 
      indx = 9
      call read_line(gbuffer, indx, 'CH', ierr, id, next_word )
*                                 ! Read the value from next_word
      if( ierr.eq.0 ) then
          jndx = 1
          call read_line(next_word, jndx, 'R4', ierr, value, cd)
          call report_error('IOSTAT',ierr,'decoding', gbuffer, 0,
     .        'DECODE_VDRAW/vector scale')
          if( ierr.eq.0 ) then
              vscale = value
          end if
      end if
 
****  Get the confidence interval argument from the string
 
      call read_line(gbuffer, indx, 'CH', ierr, id, next_word )
*                                 ! Read the value from next_word
      if( ierr.eq.0 ) then
          jndx = 1
          call read_line(next_word, jndx, 'R4', ierr, value, cd)
          call report_error('IOSTAT',ierr,'decoding', gbuffer, 0,
     .        'DECODE_VDRAW/confidence interval')
          if( ierr.eq.0 ) then
              vconfid = value
          end if
      end if
 
***** Thats all
      return
      end
 
CTITLE GET_MAP_TRAN
 
      subroutine get_map_tran( slong, slat, su, sv,
     .                       en_to_mm, arrow_length)
 
*     This routine determines the transformation between the east,north
*     cooorindate system and the uv system used by EZMAP at the point
*     being projected.
 
      include 'plot_param.h'
      include 'plot_com.h'
 
* PASSED VARIABLES
 
*   slong, slat - Start long and lat (or x and y) of the
*               - vector to be plotted.
*   su, sv      - U anf v cooridnates of start of vector (if no
*               - mapping than same as slong, slat)
*   en_to_mm(2,2)   - E and N to mm on the page transformation.
*   arrow_length    - Length of the arrow head in E coordinates
*               - units such that it will be half the character
*               - size on the page.
 
      real*4 slong, slat, su, sv, en_to_mm(2,2), arrow_length
 
* LOCAL VARIABLES
 
*   cufx, cufy  - Functions to convert user coordinates to the
*               - fractional system (0-1 accross the page)
*   dmmdu, dmmdv    - NUmber of millimeters on page for unit
*               - change in u and v coordinates
*   udlat, vdlat    -  u and v values when latitude is changed
*               - by  1 deg change
*   udlng, vdlng    -  u and v values when longitude is changed
*               - by  1 deg change
*   dudlat, dvdlat  - Changes in u and v for 1 deg change in
*               - latitude
*   dudlng, dvdlng  - Changes in u and v for 1 deg change in
*               - longitude
*   dl          - Length of change for making unit vectors.
 
 
      real*4 cufx, cufy, dmmdu, dmmdv, udlat, vdlat, udlng, vdlng,
     .    dudlat, dvdlat, dudlng, dvdlng, dl
 
****  First get the scaling of u and v to mm on the page.
 
      dmmdu = 1.0/( (cufx(1.0)-cufx(0.0))*x_size_mm )
      dmmdv = 1.0/( (cufy(1.0)-cufy(0.0))*x_size_mm )
 
*     Get the length of the arrow head now that we know E to mm on page.
*     We divide by vscale because the transformation en_to_mm will
*     scale the length back down thus giving true mm on the page.
 
      arrow_length = charsz_x * 0.5 / vscale
 
***** If we are in map_mode then get the transformation from east north
*     to uv using the MAPTRN routine (numerically)
 
      if( map_mode ) then
*                                             ! Start position
          call maptrn(slat, slong, su, sv )
          call maptrn(slat+1.0, slong, udlat, vdlat)
          call maptrn(slat, slong+1.0, udlng, vdlng)
 
*         Compute dudlat, dudlng, dvdlat, dvdlng
          dudlat = udlat - su
          dvdlat = vdlat - sv
          dudlng = udlng - su
          dvdlng = vdlng - sv
 
*         Make unit vectors for the long (e/x) and lat (n/y)
*         transformations
          dl = sqrt( dudlng**2 + dvdlng**2 )
          en_to_mm(1,1) = dmmdu*dudlng/dl * vscale
          en_to_mm(2,1) = dmmdv*dvdlng/dl * vscale
 
          dl = sqrt( dudlat**2 + dvdlat**2 )
          en_to_mm(1,2) = dmmdu*dudlat/dl * vscale
          en_to_mm(2,2) = dmmdv*dvdlat/dl * vscale
 
*                     ! We are working in a regular system
      else
 
          en_to_mm(1,1) = dmmdu * vscale
          en_to_mm(2,1) = 0.0
          en_to_mm(1,2) = 0.0
          en_to_mm(2,2) = dmmdv * vscale
 
*                         ! Save start coordinates of vector
          su = slong
          sv = slat
 
      end if
 
****  Thats all
      return
      end
 
 
