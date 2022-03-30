      subroutine svant ( satcrd,iprn,xhat_t,yhat_t,zhat_t,dx
     .                 , svantpart )
C
C     Correct the satellite vector wrt center-of-earth for the
c     offset between the phase center and center-of-mass

C     R. W. King   27 February 1991 - last modified rwk 050209
C
c         satcrd(6,2)    : SV wrt earth (inertial system) L1, L2  (km)
c         isat           : index of SV in ITSAT array
c         itsat(maxsat)  : PRN #'s of SVs on T-file
c         xhat_t,yhat_t,zhat_t : unit vectors of satellite  
c         dx(3,2)        : SV antenna wrt SV center of mass L1 L2 (m)



      implicit none

      include '../includes/dimpar.h'
       
      real*8 satcrd(6,2),xhat_t(3),yhat_t(3),zhat_t(3),dx(3,2)
     .     , svantpart(3,3)

      integer iprn,i,j
                    
      logical debug/.false./

      if( debug ) then 
        print *,'SVANT satcrd ',satcrd
        print *,'SVANT dx xhat yhat zhat '
     .   ,dx,xhat_t,yhat_t,zhat_t
        print *,'iprn dx  ',iprn,dx
      endif 
      do j=1,2
        satcrd(1,j)= satcrd(1,j) +
     .       dx(1,j)*1.d-3*xhat_t(1) + dx(2,j)*1.d-3*yhat_t(1) 
     .     + dx(3,j)*1.d-3*zhat_t(1)
        satcrd(2,j)= satcrd(2,j) +
     .       dx(1,j)*1.d-3*xhat_t(2) + dx(2,j)*1.d-3*yhat_t(2) 
     .     + dx(3,j)*1.d-3*zhat_t(2)
        satcrd(3,j)= satcrd(3,j) +
     .       dx(1,j)*1.d-3*xhat_t(3) + dx(2,j)*1.d-3*yhat_t(3) 
     .     + dx(3,j)*1.d-3*zhat_t(3)
      enddo 
   

c       Partial derivatives of the inertial vector wrt x, y, z offsets  
c       (units of partials are meters)

      do i=1,3
       svantpart(i,1) = xhat_t(i)
       svantpart(i,2) = yhat_t(i)
       svantpart(i,3) = zhat_t(i)
      enddo
                    
      return
      end
