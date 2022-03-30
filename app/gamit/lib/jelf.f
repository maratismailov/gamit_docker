      integer function jelf(i)

*     Function to compute i*(i-1)/2 and compare the
*     integer calculation with the floating value.
*     Written by T. Herring 980121 to cope with HP compiler
*     bug with this type of integer calculation.
* MOD TAH 020322: Converted back to integer calculation

      integer i
      
*     Local variables
c     integer icomp, ngood, nbad
      
c     save ngood, nbad

*     Method known to work with HPUX
c     jelf = nint((float(i)*float((i-1))/2))
      jelf = (i*(i-1))/2     
*     Standard
c     icomp = (i*(i-1))/2

*     See if results the same      
c     if( icomp.ne.jelf ) then
c         nbad = nbad + 1
c         write(*,100) i, jelf, icomp, ngood, nbad
c100     format('JELF Failure I = ',i5,' Results ',
c    .           2i12,' Ngood&bad ',2i8) 
c     else
c         ngood = ngood + 1    
c     end if
      
      end

