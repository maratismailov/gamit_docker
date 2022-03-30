      program glbak_main

      implicit none 

*     This is a "dummy" main program to allow the real glbak
*     to be called as a subroutine.  MAIN is passed to the 
*     routine if it is being run as a main program

      call glbak('MAIN')

      end

