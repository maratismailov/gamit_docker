Copyright 1993 Massachusetts Institute of Technology and The Regents of the
c              University of California, San Diego.  All Rights Reserved.
c
      subroutine weop
c
c     Written by: Simon McClusky 5/4/1994
c
c     Add earth orientation weights to normal matrix
c
      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'
      include 'parameters.h'

      integer maxc,indx,ind,jj,i,j,k   

      parameter(maxc=6*(6+1)/2)

      real*8 coveop(maxc) 

c      write(6,10)
c10    format('Weighting earth orientation parameters')
c
c  assumption is that earth orientation parameters are between islot1=8001
c  and islot1=8006 (directly after orbits).
      indx = 0
      do 20 i=1,ntpart
        if(islot1(i).ge.80001.and.islot1(i).le.80006) goto 25
        indx=indx+1
20    continue
25    continue
c
c  created lower triangle of eop weight matrix.
       do 200 j=1,6
         JJ=(J*J-J)/2
         DO 200 K=1,J
           IND=K+JJ
           coveop(IND)=0.D0
           IF(J.NE.K) GO TO 200
           coveop(IND)=eop_apr(j)**2
c           print*, 'j ind coveop(ind)',j,ind,coveop(ind)
c  invert to get weight
           if(coveop(ind).gt.0.d0) then
             coveop(IND)=1.D0/coveop(IND)
           else
             coveop(ind)=0.d0
           endif
200    continue
c
c  Add weight matrix to lower triangle normal matrix

      call addwgt( indx,maxc,6,coveop )

c      do 300 j=1,6
c        JJ1=(J*J-J)/2
c        JJ2=J+INDX
c        JJ3=(JJ2*JJ2-JJ2)/2
c        DO 300 K=1,J
c          IND1=K+JJ1
c          IND2=K+INDX+JJ3
c          a(IND2)=a(IND2)+coveop(IND1)
c          ALC(IND2)=ALC(IND2)+coveop(IND1)
c          WRITE(6,*) IND1,IND2,coveop(IND1)
c300   continue
c
      return
      end

