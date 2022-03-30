CTITLE    ..................................................................
 
      subroutine get_file
c
c     Routine to get the name of the data file
c
c Include files
c -------------
*                         ! the parameter file
      include 'plot_param.h'
c
*                         ! the common block
      include 'plot_com.h'
c
c Functions
c ---------
c trimlen -- HP utility
c
      integer*4 trimlen, iterm, ierr, jerr, NF
 
c
c.... Get the name of the file from the buffer
c     Trim blanks from beginning of file name. (The leading blanks
c     from the whole string were removed in get command.)
      call trim_lead(buffer(9:),ierr)
      if( ierr.ge.0 ) then
         read(buffer(9:),'(a)') input_file

* MOD TAH 131015: see if #N form
         if( input_file(1:1).eq.'#' ) then
             read(input_file(2:),*,iostat=jerr) NF
             if( NF.ge.1 .and. NF.le.num_files .and. jerr.eq.0 ) then
                 input_file = saved_files(NF)
             else
                 write(*,140) jerr, trim(input_file), NF, num_files
 140             format('** ERROR ',i4,' decoding file number from ',
     .                  a,' file number returns ',I2,' Max ',i2)
             endif
          endif
*                                            ! report error
      else
         iterm = 6                                                       
         write(iterm,'(a,i2,a,a)') ' TRIM_LEAD error ',ierr
     .        ,' occurred decoding ',buffer(1:trimlen(buffer))
c        below replaced by above to avoid splitting Hollerith  rwk 970920
c         write(iterm,'(" TRIM_LEAD error ",i2," ocurred",
c     .      " decoding ",a)') ierr, buffer(1:trimlen(buffer))
c
c....    Clear the file name so we will not try to open file
         input_file = ' '
c
      end if
c
      return
      end
 
