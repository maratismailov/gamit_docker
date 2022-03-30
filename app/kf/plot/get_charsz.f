CTITLE    ................................................................
 
      subroutine get_charsz
c
c     Routine to get the character size on set Graphics 1000/II to
c     output this size charcaters
c
c Include files
c -------------
*                        ! the parameter file
      include 'plot_param.h'
c
*                        ! the common block
      include 'plot_com.h'
c
c Local variables
c ---------------
c ierr -- error number for readind label
c
      integer*4 ierr
 
c
c.... Get character size from buffer
      read(buffer(9:),*,iostat=ierr,err=1000, end=1000)
     .   charsz_x

      charsz_y = charsz_x
c
c.... The character size will be set when text is output by write_label
c
      return
c
c.... Error return
 1000 continue
      call report_error('IOSTAT',ierr,'decod',buffer,0,'GET_CHARSZ')
c
      return
      end
 
