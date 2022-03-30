Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
      subroutine sbout
c
c     output tabular ephemeris
c     Rick Abbot - November 1984
c
      implicit none

      include '../includes/dimpar.h'
      include '../includes/arc.h'

      integer*4 i,j,k,l,iunit
     
      character*1 upperc

      real*8 outvec(maxyt2)
                   
      logical debug
      data debug/.false./

      npstep=0
      nrec=nrec+1

      iunit=isunit 
      if (nsign.lt.0) iunit=3        

cd      print *,'SBOUT nsign iunit ',nsign,iunit     

      if (kount.eq.0) go to 150

c     fill output vector with positions and partials wrt position,
c     skip over velocities and partials wrt velocity......

      k = 0          
      do i = 1,neq,6 
        if ( upperc(apar).ne."V") then  
          do j = 0,2
            k = k+1     
            outvec(k) = y(i+j,4) 
          enddo
        else 
          do j = 0,5
            k = k+1     
            outvec(k) = y(i+j,4)
          enddo 
        endif
      enddo
c** t not now defined
c      if( debug ) write(*,*) t,(outvec(j),j=1,k)  
      write (iunit) (outvec(j),j=1,k) 
cd      print *,'  wrote on iunit outvec(1-k) ',(outvec(j),j=1,k)
      goto 200
c old way!!!!!!
c      write (iunit)
c     $            (y(l,4),l=1,3),  (y(l,4),l=7,9),  (y(l,4),l=13,15),
c     $            (y(l,4),l=19,21),(y(l,4),l=25,27),(y(l,4),l=31,33),
c     $            (y(l,4),l=37,39),(y(l,4),l=43,45),(y(l,4),l=49,51),
c     $            (y(l,4),l=55,57)
c      go to 200
  150 write (iunit) (y(l,4),l=1,3) 
cd      print *,'SBOUT wrote on iunit y(1-3,4) ',(y(l,4),l=1,3)
  200 continue   
      return
      end
