CTITLE    .................................................................
 
      subroutine get_xy_scale(field, scale, window, ref_val)
c
c     Subroutine to read the xy scales in a format consistent with
c     the input of the data
c
c Include files
c -------------
*                         ! the parameter file
      include 'plot_param.h'
c
*                         ! the common block
      include 'plot_com.h'
c
c Variables
c ---------
c field  -- the field informatation which tells us how to decode the
c     the data field
c scale  -- the scale read from the buffer
c ref_val -- the reference values for the x or y data
c
      integer*4 field(3)
 
c
*   window(1)   - Window for fitting polynomials.  If reset scales
*               - is on (default) then window will also be reset.
      real*4 scale(2), window(2)
 
c
      real*8 ref_val
 
c
c Local variables
c ---------------
c ierr   -- error number
c values -- the values from the buffer (used to get scales)
c ivalues -- integer values read from buffer
c mjd_values(2) -- MJD of start and stop
c
      real*8 values(2), sectag(2)
 
c
      integer*4 ierr, ivalues(5,2), i,j
      integer*4 absf3
 
c
c
c.... Clear the values before we start
      do i = 1,2
         values(i) = 0.d0
         sectag(i) = 0.d0
c
         do j = 1,5
            ivalues(j,i) = 0
         end do
      end do
c
c.... See what sort of field this is (0 for time data, 1 for normal data)
*                                ! this is time data
* MOD TAH 200331: Added field 3 to be date type field.
      if( field(1).eq.0 .or. field(1).eq.3 ) then
c
c....    decode: read the number of values based on the fiels width
         if( field(3).gt.0 ) then
             absf3 = field(3)
         else
             absf3 = abs(field(3))+1
         end if
         if( absf3.lt.6 ) then  ! read upto YMD MM
             read(buffer(9:),*, iostat=ierr, err=1000, end=1000)
     .          ((ivalues(j,i),j=1,absf3),i=1,2)
         else
!            Read with floating point seconda tag
             read(buffer(9:),*, iostat=ierr, err=1000, end=1000)
     .          ((ivalues(j,i),j=1,5),sectag(i),i=1,2)
         endif
             
c
c....    Convert values to julian days
*                     ! lower and upper values
         do i = 1,2
            call ymdhms_to_mjd(ivalues(1,i),sectag(i),values(i))
         end do
c
*                                ! normal data, just read values
      else
c
         read(buffer(9:),*, iostat=ierr, err=1000, end=1000)
     .      (values(j),j=1,2)
c
      end if
      call report_error('IOSTAT',ierr,'decoding',buffer,0,
     .     'get_xy_scale')

      if( ierr.ne.0 ) RETURN
c
c.... Now remove the reference values
      scale(1) = values(1) - ref_val
      scale(2) = values(2) - ref_val
c
c.... Set the use_def_scale fale
      use_def_scale = .false.
*                           ! force rescaling of the plot
      scale_set = .false.
c
c.... Now set the scales
      call set_scale
c
c.... See if we should reset the window as well
 
      if( reset_scales ) then
          window(1) = scale(1)
          window(2) = scale(2)
      end if
c
c.... Thats all
      return
c
c.... Error address
 1000 continue
c
      call report_error('IOSTAT',ierr,'decod',buffer,0,'GET_XY_SCALE')
c
      return
      end
