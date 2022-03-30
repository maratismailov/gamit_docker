      program glout_main

      implicit none 

*     This is a "dummy" main program to allow the real glout
*     to be called as a subroutine.  MAIN is passed to the 
*     routine if it is being run as a main program

      call glout('MAIN','Dummy',0)

      end

 
*     Include the block data in the main program since some compilers/linkers (specifically OSX) will not correctly 
*     load the initialized block data routines from library archives.
      include '../globk/globk_cmd_bd.f' 
