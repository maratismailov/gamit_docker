      SUBROUTINE PUTCAR( NUMSIT,X,SITNAM,OFILE )

C       Write out a set or file of Cartesian coordinates.

c       --use apr-file convention of blank first column, 8-character site name


      implicit none

      include '../includes/tform.h'

         
      integer*4 ofile,numsit,i,j  
   
      real*8 x(3,nsdim)
                    
      character*8 siteid     
      character*12 sname
      character*16 sitnam(nsdim)
        
      do i=1,numsit
                 
        call getsname(siteid,sname,sitnam(i),2)  
        if( ofile.gt.0 ) then
          write(ofile,'(1x,a8,3f14.4)') siteid,(x(j,i),j=1,3) 
        else
         write(iscrn,'(1x,a8,3f14.4)') siteid,(x(j,i),j=1,3) 
        endif
     
      enddo

      return
      end     
