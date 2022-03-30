CTITLE    ................................................................
 
      subroutine label_ax1(labform,tic, ref_val, idep)
c
c     Routine to output a type 1 data label
c
c Variables
c ---------
c labform -- the format for the label
c tic -- the value of the tic
c ref_val -- the reference value
c idep -- the depth of the label
c
      character*(*) labform
 
c
      real*4 tic
 
c
      real*8 ref_val
 
c
      integer*4 idep, ierr
 
c
c Local variables
c ---------------
c value -- actual value of tic mark
c label -- label to be written
c
      real*8 value
 
c
      character*20 label
 
c
c Functions
c ---------
c trimlen -- HP utility to return length of string
c
      integer*4 trimlen
 
c
c.... Compute tic mark value
c
      value = tic + ref_val
*      if( abs(value).lt.1.e-5 ) value = 0.0
c
      write(label,labform, iostat=ierr) value
      call report_error('IOSTAT',ierr,'us',labform,0,'LABEL_AX1')
c
c.... Remove any leading blanks
      call trim_lead(label, ierr)
c
      idep = max(idep,trimlen(label))
*                         ! Save the width and height of the label
      call swh( idep, 1)
c
c.... Output it
C     call write_label(label)
      call save_label(label)
c
      return
      end
 
