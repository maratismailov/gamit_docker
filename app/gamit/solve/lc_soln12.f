      Subroutine LC_SOLN12
                 
c     Do the LC solution after separate L1 and L2
c     This is a renamed version of mode=2 for MODELC, separated
c     from mode=2 (WL bias fixing) for clarity.  Since I'm not
c     at all sure that the present scheme for L1+L2 then LC buys
c     us anything, this whole sequence might go away. (What we need,
c     if anything, is the original L1+L2 with constraint, without LC 
c     following.

      implicit none
     
c     global common blocks

      include '../includes/dimpar.h'
      include 'solve.h'
      include 'parameters.h'
            
   
      integer ii,nded,ierinv,ifast

      real*8 goodft,rnew
              
      logical zeroad


      write (6,190)
      write(10,190)
 190  format(//5x,'  ionospheric corrected observations (lc)',/,
     *    5x,'  after l1,l2-l1 separate mode solution')
      call wlmode(free,1)
      nobs = nobs/2

      do ii = 1,nsat
         iseen(ii) = iseen(ii)/2
      enddo

      nded = ntpart-nlive
      call lcnorm(ierinv,1)
      call qhead1(2)
c     print *,'LC_SOLUTIONC  ntpart,nded,nlive: ',ntpart,nded,nlive
c     write(6,'(5i20,/)') (free(i),i=1,ntpart)
c     write(6,'(5d20.10,/)') (a(i),i=1,ntpart)
      call report_stat('STATUS','SOLVE','lc_solution',' '
     .      , 'Solving LC normal equations after L1/L2 separate',0)
      call lcnorm(ierinv,2)
      if( ierinv.ne.0 )  call report_stat('FATAL','SOLVE','ls_solution'
     .  ,' ','Bad inversion of LC normal equations after L1/L2 solution'
     .      ,0)
      call solvlc(adjust,nded,rnew,zeroad)      
c
c     in the case of (x2 .ne. 0), we have to calculate more
c
      if ( .not. zeroad ) then
         ifast = 0
         call solve3(adjust,nded,rnew,ifast) 
      endif


c     compute postfit goodness of fit

      goodft = chi2/dble(nobs-nlive)
      goodft = dabs(goodft)
      sclerr = dsqrt(goodft)

      return
      end


