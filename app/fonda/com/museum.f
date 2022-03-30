      subroutine museum(version)
c
c     Welcome to visit the museum of COM.
c     Museum stores all history of the evolution of COM.

      character*16 version

      version = '1.02            '

c
c------------------------------------------------------------------------------------
c   version   |  time   | designer   |       description 
c------------------------------------------------------------------------------------
c   1.00       04/10/93   D.N.Dong     finish the prototype of COM  
c    ---changes 04/10/93 - 10/05/95 not recorded

c   1.01       10/05/95   Ferhat     1. correct constants for Clarke ellipsoid
c                                       (GEOTAB )
c   1.02       15/02/01   VKotzev    1. numerous small fixes to allow compile on Linux

      return
      end
