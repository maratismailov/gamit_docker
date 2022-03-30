CTITLE GET_INSTITUTE

      subroutine get_institute( creator, default )

      implicit none 

*     Get the users insitute from environment variable INSTITUTE
*
      character*(*)  creator, default

      integer*4 trimlen

      call getenv('INSTITUTE', creator )

      if ( trimlen(creator).eq.0 .or.
     .     ichar(creator(1:1)).eq.0   ) then
         write(*,120) default
 120     format('**WARNING** INSTITUTE enviroment variable ',
     .          'not set.  Default ',a,' used')
         creator = default
      end if

      return
      end

