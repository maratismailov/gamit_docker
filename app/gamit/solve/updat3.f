Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
C
      Subroutine UPDAT3( covkep )
c
c     inquiring if updating m-, s-, g- files

      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'  
      include 'parameters.h'

      integer ind

      real*8 covkep(maxcov)

c     update merge file  
      if (ipfil(1).eq.1) call update(mfiln)
c     update l-file with station coordinates
      if(ipfil(2).eq.1) call upl
c     update i-file with clock parameters
C This is incomplete for now
      if(ipfil(3).eq.1) call upi
c     check if satellite partials are present
c     if not skip question
      do 60 ind=1,ntpart
 60   if((islot1(ind).gt.500).and.(islot1(ind).le.2400)) go to 70
      go to 100

c     propagate state-vector to keplerian covariances
 70   call covst2(covkep)
      if( ipfil(4).eq.1 )
     .     call uporb 
c
 100  continue
      return
      end
