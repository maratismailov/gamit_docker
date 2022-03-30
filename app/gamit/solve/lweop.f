Copyright 1993 Massachusetts Institute of Technology and The Regents of the
c              University of California, San Diego.  All Rights Reserved.
c
      subroutine lweop

c     Written by: S. McClusky 5/4/1994
c
c     Add station weight constraints to normal matrix for loose solution
c
      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'
      include 'parameters.h'

      integer maxc,indx,ind,jj,i,j,k

      parameter(maxc=6*(6+1)/2)

      real*8 coveop(maxc),temp

c     debug
c      integer*4 jj1,jj2,jj3,ind1,ind2


c      if( logprt ) write(6,10)
c10    format(' loosely weighting earth orientation parameters')
c
c   print and write to Q-file the loose-solution constraints   

c**     If a priori's all zero, assume that not adjusted.  This is a bit of a 
c       kluge.  It would be better to leave the values out of the normal equations
c       altogether and set ieopw(1) = 0 
      if ( eop_apr2(1).eq.0.d0 .and. eop_apr2(2).eq.0.d0 .and.
     .      eop_apr2(3).eq.0.d0 .and. eop_apr2(4).eq.0.d0 ) then        
            if( logprt ) write(6,159) 
            write(10,159)
159         format(/,1x,'Pole position not adjusted')
      else
        if( logprt ) write(6,160)
        write(10,160)
160     format(/,
     1    2x,'A priori pole position errors in arcs and arcs/day',/,
     2    '   Xp   Xp_rate   Yp    Yp_rate ',/)
          if( logprt ) write( 6,161) (eop_apr2(j),j=1,4)
          write(10,161) (eop_apr2(j),j=1,4)
161       format(1x,4(f5.3,3x))   
       endif
      if ( eop_apr2(5).eq.0.d0 .and. eop_apr2(6).eq.0.d0 ) then
           if( logprt ) write(6,179) 
           write(10,179)
179        format(/,1x,'UT1 and rate not adjusted')
      else
        if( logprt ) write(6,180)
        write(10,180)
180     format(/,
     1    2x,'A priori earth rotation errors in sec and sec/day',/,
     2    '  ut1   ut1_rate',/)
          if( logprt ) write( 6,182) (eop_apr2(j),j=5,6)
          write(10,182) (eop_apr2(j),j=5,6)
182       format(1x,2(f5.3,3x))   
      endif
c
c  assumption is that earth orientation parameters are between islot1=80001
c  and islot1=80006 (directly after orbits).
      indx = 0
      do 20 i=1,ntpart
        if(islot1(i).ge.80001.and.islot1(i).le.80006) goto 25
        indx=indx+1
20    continue
25    continue
c
c
c   created lower triangle of eop weight matrix.
       do 200 j=1,6
         JJ=(J*J-J)/2
         DO 200 K=1,J
           IND=K+JJ
           coveop(IND)=0.D0
           IF(J.NE.K) GO TO 200
           coveop(IND)=eop_apr2(j)**2
c           print*, 'j ind coveop(ind)',j,ind,coveop(ind)
c   temp = (new constraint)/(old constraint)
c   (remember that weight = 1.0/(constraint**2) )
c   weight increment = (new weight) - (old weight)
c                    = 1.0/(new variance) - 1.0/(old variance)
c                    = (1.0 - temp*temp)/(new variance)
c   new weight = old weight + weight increment
c
           if(eop_apr(j).gt.0.d0) then
             temp = eop_apr2(j)/eop_apr(j)
           else
C   no eop constraints in constrained solutions
             temp=1.d0
           endif
C   invert to get weight
           if(coveop(ind).gt.0.d0) then
             coveop(IND)=(1.D0-TEMP*TEMP)/coveop(IND)
           else
             coveop(ind)=0.0d0
           endif
200    continue

c   Add weight matrix to lower triangle normal matrix

      call addwgt( indx,maxc,6,coveop )
 
c      do 300 j=1,6
c        JJ1=(J*J-J)/2
c        JJ2=J+INDX
c        JJ3=(JJ2*JJ2-JJ2)/2
c        DO 300 K=1,J
c          IND1=K+JJ1
c          IND2=K+INDX+JJ3  
c          print *,'ADDWGT: data  ind2 alc ',ind2,alc(ind2)
c          a(IND2)=a(IND2)+coveop(IND1)
c          ALC(IND2)=ALC(IND2)+coveop(IND1) 
c         print*,'ADDWGT: indx ind1 ind2 alc ',indx,ind1,ind2,alc(ind2)
c        print *,'LWEOP ind1 ind2 coveop ',IND1,IND2,coveop(IND1)
c300   continue
c
      RETURN
      END
