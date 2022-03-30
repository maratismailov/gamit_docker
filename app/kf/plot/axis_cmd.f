CTITLE    .....................................................................
 
      subroutine axis_cmd(buffer, tic_space, label_mult, label)
c
c     Routine to read the axis information from the user buffer
c
c Variables
c ---------
c buffer -- the user input buffer
c tic_space -- the spacing between tics in the world coordinates
c label_mult -- the multiple of the tic spacing at which labels
c     will appear
c label -- the axis label to be used
c
      character*(*) buffer, label
 
c
      real*4 tic_space
 
c
      integer*4 label_mult, ierr
 
c
c.... get the spacing information
      read(buffer(9:),*,iostat=ierr, err=1000, end=1000 )
     .  tic_space, label_mult
c
c.... Now get the label
      call read_label(buffer,label)
c
      return
c
c.... Error return
 1000 continue
      call report_error('IOSTAT',ierr,'decod',buffer,0,'AXIS_CMD')
      return
c
      end
 
