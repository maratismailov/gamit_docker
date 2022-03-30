      SUBROUTINE PROPER(IUNIT)
c
c     Write copyright notice to logical unit iunit.

      integer*4 iunit

      WRITE(IUNIT,10)
   10 format(/,1x,'Copyright 1996 Massachusetts Institute of Technology'
     . ,/,' and The Regents of the University of California, San Diego.'
     . ,/,' All Rights Reserved',/ )
      RETURN
      END



