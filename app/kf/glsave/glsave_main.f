      program glsave_main

      implicit none 

*     This is a "dummy" main program to allow the real glsave
*     to be called as a subroutine.  MAIN is passed to the 
*     routine if it is being run as a main program

      call glsave('MAIN','YES')

      end

