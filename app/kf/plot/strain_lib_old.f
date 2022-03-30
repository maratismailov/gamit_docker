CTITLE MAKE_CLOVER
 
      subroutine make_clover( eps1,eps2,theta,eps1sig,eps2sig,
     .                        sscale, dxp, dyp, nump)
 
c     Make a strain clover leaf
c     Kurt Feigl 910224
 
c     INPUT VARIABLES:
c        sscale               : scaling factor for size of cloverleaf
c        theta,thetasig       : azimuth of more compressive eigenvector (deg)
c        eps1,eps2            : eigenvalues of strain rate (1/yr)
c        eps1sig,eps2sig      : uncertainties of eigenvalues of strain rate (1/yr)
c     OUTPUT VARIABLES:
c       dxp, dyp                : arrays of X and Y 
c       nump                  : number of points
       real*4 sscale,theta,eps1sig,eps2sig,eps1,eps2
       real*4 dxp(*), dyp(*)
       integer nump,istep1,istep2,i

c     LOCAL VARIABLES:
c       x0, y0               : UV coordinates of origin
c       eps1snr,eps2snr:     : value/sigma of strain rate eigenvalues
c                            : eps1 is max extension 
c                            : eps2 is max compression
c       edot11,edot12,edot22 : components of strain rate (1/yr)
c                            : 1-axis points E, 2-axis points N
c                            : positive strain is extension!
      real*4 edot11,edot12,edot22,eps1snr,eps2snr, sscale
      real*4 degrad,th,eps,x0,y0,x1,x2,y1,y2,th0
      
c     degrees to radian conversion factor
      degrad = atan(1.0d0)/45.0d0

c     convert to radians
      th0 = theta * degrad

c     write (6,*) 'MAKE_CLOVER: sscale = ',sscale

c     change notation!
      edot11 = eps1*(cos(th0)**2) + eps2*(sin(th0)**2)
      edot22 = eps1*(sin(th0)**2) + eps2*(cos(th0)**2)
      edot12 = (eps2 - eps1)*sin(th0)*cos(th0)

c     set origin
      x0 = 0.
      y0 = 0.

c     number of lines depends on significance level of epsilon 1
      eps1snr = abs(eps1/eps1sig)
      if (eps1snr .gt. 3.0d0) then
         istep1 = 3
      else if (eps1snr .gt. 1.0d0) then
         istep1 = 1
      else
         istep1 = 0
      endif
c     number of lines depends on significance level of epsilon 2
      eps2snr = abs(eps2/eps2sig)
      if (eps2snr .gt. 3.0d0) then
         istep2 = 1
      else if (eps2snr .gt. 1.0d0) then
         istep2 = 3
      else
         istep2 = 0
      endif

      
      nump = 0
      do i = 0,360
         th = i * degrad
         eps = edot11 * (sin(th))**2
     .       + edot12 * sin(2.*th)
     .       + edot22 * (cos(th))**2
         x1 = x0 + eps * sin(th) * sscale
         y1 = y0 + eps * cos(th) * sscale

         nump = nump + 1
         dxp(nump) = x1
         dyp(nump) = y1

c        extensional quadrant
         if (eps .gt. 0d0 .and. istep1 .ne. 0) then
            if (mod(i,5*istep1) .eq. 0) then
c              go to zero and come back
               nump = nump + 1
               dxp(nump) = x0
               dyp(nump) = y0
               nump = nump + 1
               dxp(nump) = x1
               dyp(nump) = y1
            endif
         endif
c        compressional quadrant
         if (eps .le. 0d0 .and. istep2 .ne. 0) then
            if (mod(i,5*istep2) .eq. 0) then
               x2 = x0 + eps*sin(th0-(th-th0))*sscale
               y2 = y0 + eps*cos(th0-(th-th0))*sscale
c              go to opposite side and come back
               nump = nump + 1
               dxp(nump) = x2
               dyp(nump) = y2
               nump = nump + 1
               dxp(nump) = x1
               dyp(nump) = y1
            endif
         endif
      enddo

      return
      end

CTITLE make_strain_arrow
 
      subroutine make_strain_arrow( eps1,eps2,theta,eps1sig,eps2sig,
     .           arrow_len,sscale, dxp, dyp, nump)
 
c     Make a set of strain arrows
c     for the time being, neglect uncertainties.
c     Danan Dong 911203
 
c     INPUT VARIABLES:
c        sscale               : scaling factor for size of cloverleaf
c        theta,thetasig       : azimuth of more compressive eigenvector (deg)
c        eps1,eps2            : eigenvalues of strain rate (1/yr)
c        eps1sig,eps2sig      : uncertainties of eigenvalues of strain rate (1/yr)
c        arrow_len            : arrow length for the strain vector
c     OUTPUT VARIABLES:
c       dxp, dyp                : arrays of X and Y 
c       nump                  : number of points
c     LOCAL VARIABLES:
c       mode                  : 1 = extension;  2 = compression
c       dx, dy                : local coordinate of stain vector
c
       real*4 sscale,theta,eps1sig,eps2sig,eps1,eps2
       real*4 dxp(*), dyp(*),axp(5),ayp(5),arrow_len
       integer nump,istep1,istep2,i,j,mode

      real*4 edot11,edot12,edot22,eps1snr,eps2snr, sscale
      real*4 degrad,th,eps,x0,y0,x1,x2,y1,y2,th0,dx,dy
      
c     degrees to radian conversion factor
      degrad = atan(1.0d0)/45.0d0

c     convert to radians
      th0 = theta * degrad

c     strain vector has 4 arrows, every arrow has 5 points.
      nump = 20

c        extension component
         dx = eps1*cos(th0)*sscale
         dy = -eps1*sin(th0)*sscale
         mode = 1
         if (eps1.lt.0.0d0) mode = 2
         call make_arrow_new(dx,dy,arrow_len,dxp,dyp,mode)
         do j = 1,5
            dxp(j+5) = -dxp(j)
            dyp(j+5) = -dyp(j)
         enddo

c        compression component
         dx = eps2*sin(th0)*sscale
         dy = eps2*cos(th0)*sscale
         mode = 2
         if (eps2.gt.0.0d0) mode = 1
         call make_arrow_new(dx,dy,arrow_len,axp,ayp,mode)
         do j = 1,5
            dxp(j+10) = axp(j)
            dyp(j+10) = ayp(j)
            dxp(j+15) = -axp(j)
            dyp(j+15) = -ayp(j)
         enddo

      return
      end

CTITLE SDRAW
 
      subroutine sdraw(gbuffer)
 
*     Routine to strain cloverleaves on a plot.  Either a regular figure
*     or a map projection.  If a map projection is to be used, then
*     the map should be draw first to ensure that all of the needed
*     parameters are set.
*
      include 'plot_param.h'
      include 'plot_com.h'
 
*     PASSED VARIABLES
*     The command line to obtain scale factor for strain cloverleaves.
      character*(*) gbuffer
 
c     LOCAL VARIABLES
c     relative coordinates of cloverleaf 
      real*4 dxp(max_clover),dyp(max_clover)
c     absolute coordinates of cloverleaf 
      real*4 xp(max_clover),yp(max_clover)
c     scale factor
      real*4 sscale
c     number of points
      integer nump
c     loop index
      integer*4 i
c     type of plot (currently ignored)
      integer*4 istrain,ival

      integer ierr

      character*20 buff20
 
*   su, sv      - Postion of current point in u/x and v/y
*               - coordinates
*   en_to_mm(2,2)   - Transformation of e and n velcoties to
*               - mm on page (scaled by vscale)
*   arrow_length    - Length of the arrow head in east north
*               - coordinates
*   eu, ev      - End of the arrow
 
      real*4 su, sv, en_to_mm(2,2), arrow_length, eu, ev, value
 
****  Finish decoding the command line.  Do now so that values entered will
*     not be lost if there is no data.
c      call decode_sdraw( gbuffer, sscale, istrain)

      read (gbuffer,*,iostat=ierr) buff20,ival,value
      call report_error('IOSTAT',ierr,'decoding', gbuffer, 0,'SDRAW')
      if( ierr.eq.0 ) then
         sscale = value
         istrain = ival
      end if
         print *,'sscale, istrain ', sscale, istrain
 
 
****  Make sure we have some data
      if( num_strains.eq.0 ) then
          write(*,100)
 100      format('SDRAW Error: No data read yet, use SREAD command')
          RETURN
      end if
 
*     If we are not in map mode set view and scales
      if( .not. map_mode ) then
          call set_view
          call set_scale
      end if

c     do this for tom
      vscale = 1.d0

      print 17,'LAT','LON','EPS1','EPS1SIG','EPS2'
     .,'EPS2SIG','CW-SPIN','SPINSIG','THETA','THETASIG'
     .,'TRIANGLE'
      print 17,'(deg)','(deg)','(1/yr)','(1/yr)','(1/yr)'
     .,'(1/yr)','(1/yr)','(1/yr)','(deg)','(deg)','corners'
  17  format (12(a9,1x))


****  Now loop over the data
      do i = 1, num_strains
c        print out the strain rates:
         write (6,200) spos(2,i),spos(1,i),
     .      strain(1,i),strainsig(1,i),
     .      strain(2,i),strainsig(2,i),
     .      strain(3,i),strainsig(3,i),
     .      strain(4,i),strainsig(4,i)

200      format (2(1x,f9.4),6(1x,1pe9.2),2(1x,0pf9.4),1x,a16)

 
 
*        Get transformation from e and n velocotity to mm on the
*        page scaled by the velocity scaling
 
         call get_map_tran(spos(1,i), spos(2,i),  su, sv,
     .                     en_to_mm, arrow_length )

c        Make the cloverleaf in relative coordinates
         if (istrain.eq.1)
     .   call make_clover(strain(1,i),strain(2,i),strain(4,i),
     .         strainsig(1,i),strainsig(2,i),
     .         sscale, dxp, dyp, nump )

c        Make the cloverleaf in relative coordinates
         if (istrain.eq.2)
     .   call make_strain_arrow(strain(1,i),strain(2,i),strain(4,i),
     .         strainsig(1,i),strainsig(2,i),
     .         arrow_length,sscale, dxp, dyp, nump )

         call dxy_to_abs( su, sv, dxp, dyp, xp, yp, nump, en_to_mm )

c        Now draw the cloverleaf
         call j2ply( nump, xp, yp )

      end do
 
      return
      end
 
CTITLE READ_SFILE
 
      subroutine read_sfile( data_array )
 
*     This routine will loop over the data records read from the
*     strain rate file and extract the information specified in the
*     S_FIELD command.
 
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
*   strain_max     - Length of the longest vector
*   vel_mag     - Length of vector for setting scales
*   dx, dy      - Length of x and y axes
 
 
      real*4 values(50), vel_mag, strain_max, dx, dy
 
*   cdata       - Data record as character string.
 
      character*(char_recl) cdata
 
*   next_word   - Word with the name of the site
 
      character*20 next_word
 
*   cd          - Dummy for readline
 
      character*4 cd

      equivalence ( cdata, idata )

*     See if file has been read
      if( .not.file_read ) then
          write(*,50) 
  50      format('*** ERROR IN READ_SFILE *** ',
     .           'No file given yet. Use the ',
     .           ' FILE command')
          return
      end if

      largest_field = 0
      do i = 1, 10
          largest_field = max(largest_field, s_field(i))
      end do
 
      if( largest_field.gt.50 ) then
          write(*,100)
  100     format(' ** ERROR IN SREAD_DATA** Largest field in',
     .           ' S_FIELD too large.',/, 10x,
     .           ' Largest value allowed is 50.' )
          num_strains = 0
          RETURN
      end if
 
****  Loop over all of the data records
      num_strains = 0
      do i = 1, num_epochs
 
*         Copy the data record to local storage
          do j = 1, int_recl
              idata(j) = data_array(j,i)
          end do

          next_word = cdata 
 
          num_strains = num_strains + 1
 
*         Now pull off the values.  
          indx = 1
 
c         get the values
          call multiread(cdata,indx,'R4',ierr,values,cd,10)

c         save reals 
          if( ierr.eq.0 ) then
              spos(1,num_strains) = values(s_field(1))
              spos(2,num_strains) = values(s_field(2))
              strain(1,num_strains) = values(s_field(3))
              strain(2,num_strains) = values(s_field(5))
              strain(3,num_strains) = values(s_field(7))
              strain(4,num_strains) = values(s_field(9))
              strainsig(1,num_strains) = values(s_field(4))
              strainsig(2,num_strains) = values(s_field(6))
              strainsig(3,num_strains) = values(s_field(8))
              strainsig(4,num_strains) = values(s_field(10))
          else
              num_strains = num_strains - 1
          end if
c         Looping over data records
      end do

      print *,'Number of strains read: ',num_strains
      return
      end


CTITLE DECODE_SDRAW
 
      subroutine decode_sdraw( gbuffer, sscale, istrain)
 
*     Routine to decode the remainder of the sdraw command line.
*     Extracts scale if given.
*     extracts istrain if given
 
*     PASSED VARIABLES
 
*     Scale from microstrain to mm on the page (for 200 mm wide page output)
      real*4 sscale
c     type of strain plot 0 = cloverleaf
      integer*4 istrain
c     Command line read from user
      character*(*) gbuffer
 
* LOCAL VARIABLES
 
*   indx, jndx  - Pointers for reading through string.
*               - Arguments start at indx=9
*   id          - Dummy argument for read_line
 
      integer*4 indx, jndx, id, ierr, ival
 
*   value       - Value read from string (used to stop
*               - over writing default)
 
      real*4 value
 
*   next_word   - Next word. Words pulled from buffer, to
*               - extract scale
 
      character*20 next_word
 
      indx = 9
      call read_line(gbuffer, indx, 'CH', ierr, id, next_word )
*                                 ! Read the value from next_word
      if( ierr.eq.0 ) then
          read (next_word,*,iostat=ierr) ival,value
          call report_error('IOSTAT',ierr,'decoding', gbuffer, 0,
     .        'DECODE_SDRAW')
          if( ierr.eq.0 ) then
              sscale = value
              istrain = ival
          end if
      end if
 
      print *,'sscale, istrain ',sscale,istrain 
      return
      end
 
CTITLE MAKE_ARROW
 
      subroutine make_arrow_new(dx,dy,arrow_len,dxp,dyp,mode)
 
*     This routine makes and arrow returns the results in the
*     5 vector dxp, dyp.  Arrow length in is same units as
*     dx,dy.  The first point in the vector is 0,0.
 
* PASSED VARIABLES
 
*   dx, dy  - Differential postion of arrow head
*   arrow_len   - Length of the arrow head (same units as
*           - dx,dy
*   dxp(10), dyp(10)  - Five points making at the arrow.
 
      real*4 dx, dy, arrow_len, dxp(10), dyp(10)
 
* LOCAL VARIABLES
 
*   i       - Loop counter
 
      integer*4 i,mode
 
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
          do i = 1,10
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
    
      if (mode.eq.1) then
         dxp(2) = veclen
         dyp(2) = 0
    
         dxp(3) = dxp(2) - arrow_len
         dyp(3) = -arrow_len/4
    
         dxp(4) = dxp(3)
         dyp(4) = arrow_len/4
    
         dxp(5) = dxp(2)
         dyp(5) = dyp(2)
      endif
 
      if (mode.eq.2) then
         dxp(2) = arrow_len
         dyp(2) = -arrow_len/4
    
         dxp(3) = dxp(2)
         dyp(3) = arrow_len/4
    
         dxp(4) = 0
         dyp(4) = 0
    
         dxp(5) = veclen
         dyp(5) = 0
      endif
 
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
 
