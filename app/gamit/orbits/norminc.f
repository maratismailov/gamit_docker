      Subroutine norminc( iepoch,ksat ) 

c     Form O-C for 3satellite coordinates at a single T-file epoch,
c     calculate the partials of coordinates wrt the estimated parameters
c     (ICs, radiation pressure, and rotations), and increment the normal
c     equations

      implicit none

c     Input:  iepoch         epoch number for indexing arrays   
c             ksat           index of SV on reference T-file

c    Output in commons /omcpart/ and /nrmcom/:
c             omc(3,maxsat,maxepc)  O-C's for x,y,z for each SV at each epoch
c             amat(mxoprm,mxoprm)  coefficient matrix for normal equations
c             bvec(mxoxprm)          right-hand-side of normal equations
         
      
      include '../includes/dimpar.h'
      include '../includes/orbits.h'
      include 'orbfit.h'

      integer*4 iepoch,ksat,i,j,k

      logical debug
 
      data debug/.false./      

c     Increment the left-hand side 
                         
       if( debug ) then
        if ( iepoch.le.22 ) then
          print *,'NORMINC iepoch ksat part '
     .      ,iepoch,ksat,(part(1,j,ksat,iepoch),j=1,nparam)
        endif  
       endif
      
      do j=1,nparam
        do k=1,nparam
          do i=1,3
            amat(j,k) = amat(j,k) + 
     .                     part(i,j,ksat,iepoch) * part(i,k,ksat,iepoch)
          enddo
        enddo
      enddo


c     Increment the right-hand side

      do j=1,nparam
        do i=1,3
          bvec(j) = bvec(j) + part(i,j,ksat,iepoch) * omc(i,ksat,iepoch)
        enddo
      enddo 
            
      if( debug ) then
        if ( iepoch.le.22 ) then
           print *,' bvec ',(bvec(i),i=1,nparam)  
        endif
      endif

      return
      end


