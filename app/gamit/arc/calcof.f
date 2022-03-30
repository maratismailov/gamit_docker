Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
      subroutine calcof
c
c     subroutine to compute the coefficients for the adams-moulton integration
c     Rick Abbot - October 1984
c     reference: numerical intergration of differential equations by three
c     methods - equations for, description of computer routines, discussion of
c     results, W.B.Smith, Lincoln Lab Technical Note 1968-31, August 21,1968

      implicit none

      include '../includes/dimpar.h'   
      include '../includes/arc.h'

      integer*4 np,i,k,l

      real*8 a,b
      dimension a(13),b(13)  

c adam,calcof
c      common/coeff/cadams(13),dadams(13)     
c      real*8 cadams,dadams

       logical debug/.false./

c.....in the technical note the subscript 0 corresponds to array index 1 in
c.....this subroutine (so that the array indices are offset by 1 from the
c.....subscripts in the note)
c
      np=npredt+1
c
      a(1) = 1.d0
      a(2) = 1.d0 - 0.5d0
      a(3) = 1.d0 - 1.d0/3.d0 - 0.5d0*a(2)
      if (np.le.3) go to 500
      do 200 i=3,npredt
      a(i+1) = 1.d0
      do 100 l=2,i+1
  100 a(i+1) = -1.d0/(i+1-l+2)*a(l-1) + a(i+1)
  200 continue
c
  500 b(1) = 1.d0
      b(2) = -0.5d0
      b(3) = -1.d0/3.d0 - 0.5d0*b(2)
      if (np.le.3) go to 1000
      do 700 i=3,ncoret
      b(i+1) = 0.d0
      do 600 l=2,i+1
  600 b(i+1) = -1.d0/(i+1-l+2)*b(l-1) + b(i+1)
  700 continue
c
 1000 do 1200 k=1,npredt
      cadams(k) = 0.d0
      do 1100 i=k,npredt
 1100 cadams(k) = cadams(k) + (-1)**(k-1)*binom(i,k)*a(i)
 1200 continue
c
      do 1700 k=1,ncoret
      dadams(k) = 0.d0
      do 1600 i=k,ncoret
 1600 dadams(k) = dadams(k) + (-1)**(k-1)*binom(i,k)*b(i)
 1700 continue
c            
      if(debug) print *,'CALCOF a b cadams dadmas 1 & npredt '
     .             ,           a(1),b(1),cadams(1),dadams(1)  
     .             ,  a(npredt),b(npredt),cadams(npredt),dadams(npredt)

      return
      end
