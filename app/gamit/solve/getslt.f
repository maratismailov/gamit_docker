c
      Subroutine GETSLT( ierinv )
c
c     solve normal equation to get least square solution
c     and postfit goodness
c
      implicit none

      include '../includes/dimpar.h'
      include 'solve.h' 
      include 'parameters.h' 
c
      integer ierinv,nded,ifast  
c*** For debug
c      integer*4 i,j,jel

      real*8 goodft,rnew

      logical zeroad

c
      zeroad = .true.
      nded = ntpart-nlive
c
c     to form subrix of live parameters and to invert
c     that matrix
c                 
c      print *,'GETSLT before SOLVE1 '
c      call printa( 1440 )    


      call solve1( nded,ierinv)

c      print *,'ntpart nlive nded ',ntpart,nlive,nded     

c      print *,'GETSLT after SOLVE1 '
c      call printa( 1440 )
      if(ierinv.ne.0) return
c
c     derive solutions, l2flag = 1 : LC mode
c
      if (l2flag.ne.1) then

         call solve2( adjust,nded,rnew,zeroad ) 
c         print *,'GETSLT after SOLVE2 '
c         call printa( 1440 )
c         call printan22( 490 )
      else

         call solvlc( adjust,nded,chi2,rnew,zeroad) 
c         print *,'GETSLT after SOLVLC '
c         call printa( 1440)  
      endif
c
c     in the case of x2 = 0, all results have been obtained
c
       if ( zeroad ) goto 20
c
c     ifast = 0 means do both chi2 and params
c
      ifast = 0
      call solve3(adjust,nded,rnew,ifast)      
c      print *,'GETSLT after SOLVE3 '
c      call printa( 1440 ) 

c
c     compute postfit goodness of fit
c
 20   goodft = chi2/dble(nobs-nlive)
      goodft = dabs(goodft)
      sclerr = dsqrt(goodft)
c      print *,'GETSLT nobs nlive chi2 goodft sclerr '
c    .       ,        nobs,nlive,chi2,goodft,sclerr
c
      return
      end
