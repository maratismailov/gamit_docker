      Subroutine norm_solve

c     Solve the normal equations for ORBFIT.  

      implicit none
     
      include '../includes/dimpar.h'
      include '../includes/orbits.h'
      include 'orbfit.h'

      integer*4  job,info,i,j 

      real*8 scale(mxoprm),rcond0,z(mxoprm),dterm2(2)

      character*256 message

c     Invert the normal equations using linpack (blas 1) routines
c     assumng the that the matrix is positive definite (and thus
c     symmetric).  
           
      call report_stat('STATUS','ORBFIT','orbits/norm_solve',' '
     .                 ,'Solving the normal equations',0)
 
c     first scale the matrix
      do i=1,nparam
        scale(i) = 1.d0/dsqrt(amat(i,i))
      enddo
      do i=1,nparam
        do j=1,nparam
          amat(i,j) = scale(i)*scale(j)*amat(i,j)
        enddo 
      enddo
 
 
c     factor the matrix
      call dpoco( amat,mxoprm,nparam,rcond0,z,info )
      
      if ( info.eq.0 .and.rcond0.ge.1.d-16 ) then      

c       do the inversion (job = 1 means no determinant)
        job = 1                                   
c       dpodi takes a full matrix but returns the inverse in upper triangle
        call dpodi( amat,mxoprm,nparam,dterm2,job)
c       fill in the lower half and rescale 
        do i=1,nparam
          do j = 1, i
            amat(j,i) = scale(i)*scale(j)*amat(j,i)
c*              amat(j,i) = amat(j,i) / (scale(i)*scale(j))
            amat(i,j) = amat(j,i)
          enddo
        enddo

c       multiply inv(amat) * bvec to get the adjustment vector 
        call dgemv ('N', nparam, nparam, 1.d0, amat, mxoprm, bvec
     .             , 1, 0.d0, adjust, 1)

                                     
      else
                         
        write(message,'(a,d7.1)') 
     .      'Normal matrix ill-conditioned, rcond0 = ',rcond0
        call report_stat('STATUS','ORBFIT','orbits/norm_solve',' '
     .                 ,message,0)
      endif
 
        
      call report_stat('STATUS','ORBFIT','orbits/norm_solve',' '
     .                 ,'Normal equations solve',0)
      return
      end

