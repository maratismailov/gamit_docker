      Subroutine LWZEN    

c     Add station weight constraints to normal matrix for loose solution
c     See wzen.f for information about on the gauss markov weighting system.

c     Written by R. King from Y. Bock routines WZEN and LWSTAT  24 Sept 1993 
c     and modified from S. McClusky additions to WZEN for Markov weights Jun 1994  
c     Modified for average zenith delays and new Markov process by R King May 2004.


      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'
      include 'parameters.h'

      integer*4 maxc
      parameter(maxc=maxatm*(maxatm+1)/2)   
     
               
      character*256 message

      integer*4 indx,i

      real*8  covatm(maxc)



c     WRITE(6,999)
c 999 FORMAT(' WEIGHTING ZENITH DELAY PARAMETERS')
             

c  Print and write to Q-file the loose-solution constraints

        if( logprt ) write(6,130)
     .    (i,rlabel((i-1)*3+1)(1:4),sitnam(i),zen_apr2(i)
     .    ,zen_mar2(i),zen_tau2(i),i=1,nsite)
        write(10,130)
     .    (i,rlabel((i-1)*3+1)(1:4),sitnam(i),zen_apr2(i)
     .    ,zen_mar2(i),zen_tau2(i),i=1,nsite)
  130   format(//
     .        ,' A priori zenith-delay errors in meters',/
     .        , 100(i3,2x,a4,1x,a12,2x,3f10.3/))    

c  Test to make sure the Markov constraints are the same as the
c  tight solution

      do i=1,nsite
        if(zen_mar2(i).ne.zen_mar(i) .or. zen_tau2(i).ne.zen_tau(i) )
     .    then 
            write(message,'(a,i3,a)') 
     .        'Loose zenith delay Markov for site ',i
     .        ,' not equal tight solution value'
            call report_stat('FATAL','SOLVE','lwzen',' ',message,0)   
        endif
      enddo
                

c  First weight the average values  
 
c     new weight = old weight + weight increment
c     weight increment = (new weight) - (old weight)
c                      = 1.0/(new variance) - 1.0/(old variance)

c    atmospheric parameters follow the station coordinates
      indx = 3*nsite
      do i=1,nsite  
        if( zen_apr(i).gt.0.d0) then  
            covatm(1) = 1.d0/zen_apr2(i)**2 - 1.d0/zen_apr(i)**2 
        else
          covatm(1) = 1.d0/zen_apr2(i)**2
        endif
        call addwgt(indx,1,1,covatm)
        indx = indx + 1
      enddo      
         
c     With new average zenith delays now passed to the h-file, we always keep
c     the same weighting of the piecewise linear parameters in the tight and
c     loose solutions, so no need to add or adjust weights in the normal matrix
  
      return
      end

