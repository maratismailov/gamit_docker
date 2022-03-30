CTITLE JBEGN
 
      subroutine jbegn( metafile )
 
*     Setup routine which opens the gkssystem.
 
      include 'g1000.h'
 
*   metafile        - NAme of metafile (saved in the
*                   - grahics common for later use)
 
      character*(*) metafile
 
      gmetaf = metafile
 
*     Start up system
 
      call opngks
 
*     Thats all
      return
      end
 
