Copyright (c) Massachusetts Institute of Technology and the University of
California at San Diego, 1995.  All rights reserved.
c
      Subroutine COVST2( covkep )
c
c     convert state-vectors to keplerian elements including covariances
c
      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'   
      include 'parameters.h'

      integer indx1,indxs1,igood,isat,ielem,ind,ij,i,j,k

      real*8 cove(21),temp1(6),orbsig(6),covkep(maxcov),sigwm
     .     , factor,temp

c     find first orbital element in entire parameter list (indx1)
c
      indx1=0
      do 10 i=1,ntpart
         if(islot1(i).gt.500.and.islot1(i).le.1100) go to 11
         indx1=indx1+1
   10 continue
   11 continue

c  loop over all satellites
      igood=0

      if( logprt ) write(     6,220)
      if (iqflag.eq.1) write(10,220)
  220 format(/,1x,
     .' Post-fit Keplerian orbital errors in parts in 10**7',//,
     .'          a         e         I       Node   Perigee',
     .'    M Anom.     w+M')

      do 100 isat=1,nsat
         indxs1=indx1+(norb)*(isat-1)
c        skip unadjusted satellite (assuming all elements are fixed or free)
         if(free(indxs1+1).eq.0) go to 100
         igood=igood+1
         ielem=21*(igood-1)
         call copy1d(ielem+1,ielem+21,-ielem,covkep,cove)
c
c      compute error of sum of perigee and mean anomaly
       sigwm=dsqrt(cove(15)+cove(21)+2.d0*cove(20))
c
      factor=1.d7
      do 211 j=1,6
c        convert semi-major axis error units to parts in 10**7
         ind=(j*j-j)/2+j
         orbsig(j)=dsqrt(cove(ind))
c
c     no factor 3, no rms scaling
c
c        temp=3.d0*sclerr*orbsig(j)
         temp = orbsig(j)
         if(j.eq.1) temp1(j)=temp*factor/26000.d0
c        convert other element error units to parts in 10**7
         if(j.ge.2) temp1(j)=temp*factor
  211 continue
c
      sigwm=sigwm*factor
c
      if( logprt) write( 6,221) isprn(isat),(temp1(k),k=1,6),sigwm
      if (iqflag.eq.1)
     .   write(10,221) isprn(isat),(temp1(k),k=1,6),sigwm
  221 format(1x,'PRN',1x,i2,1x,7(f8.3,2x))

c     compute error correlation matrix, and write to o-file only
        
      if( ioflag.eq.1 ) then         
        write(15,6000)
 6000   FORMAT (/,1X,'ERROR CORRELATION MATRIX:',/)
        ij=0
        do 600 i=1,6
          do 500 j=1,i
            ij=ij+1
            cove(ij)=cove(ij)/(orbsig(i)*orbsig(j))
  500     continue
          write(15,7000) i,(cove(j),j=ij-i+1,ij)
  600   continue
 7000   format (1x,i3,'. ',(6x,6f8.4))
      endif

  100 continue
c
      return
      end
