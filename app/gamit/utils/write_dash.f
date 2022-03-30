      subroutine write_dash(isnx)

c  put this into a routine 'cos I got sick of typing it all the time!!!!!
c  writes '*-----------------' etc out to column 80

      implicit none

      integer isnx
      character*80 dashes

      data dashes/'*----------------------------------------------------
     .---------------------------'/

      write(isnx,'(a80)')dashes
      return
      end
