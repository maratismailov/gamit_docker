ctitle
 
      subroutine set_use(name_of_file,use)
 
      implicit none 

c
c     routine to set 'use' false if the file name 'none' is given
c
c Variables
c ---------
c name_of_file -- file name to be checked
c use    -- logical variable which is set false if name_file equals
c        'none'.
c
      character*(*) name_of_file
 
c
      logical use
 
c
c.... see if name_of_file is none
      if( name_of_file(1:4).eq.'NONE' ) then
         use = .false.
      else
         use = .true.
      endif
c
      return
      end
 
c.........................................................................
