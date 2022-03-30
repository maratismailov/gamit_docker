      SUBROUTINE PUTCYL( NUMSIT,X,SITNAM,OFILE )
C
C       Write out a set or file of cylindrical coordinates.
C
      implicit none

      include '../includes/tform.h'

      INTEGER OFILE,numsit
      real*8 x(3,nsdim)         
      character*16 sitnam(nsdim)
                     
      write(iscrn,'(a)') 'PUTCYL not coded, returning without action'

c     dummy statement to avoid compiler warning of unused argument 
      if( numsit.eq.0 ) then
        print *,'PUTCYL numsit ',numsit,x,sitnam,ofile
      endif
 
      end
