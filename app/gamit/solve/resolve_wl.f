      Subroutine RESOLVE_WL( mfix ) 
                                                                
c     Resolve the wide-lane ambiguities (only option currently is decision-function)
c     Called from GET_WIDELANE; calls ROUND_BIAS.  R. King 9 January 2004

c     Code from D. Dong's routine GETWL code for IJOB=2 (statement 100ff).   


      implicit none     

c-- Argument

c       return number of biases fixed
      integer*4 mfix
         
c-- Global common blocks

      include '../includes/dimpar.h'
      include 'solve.h'
      include 'parameters.h'
       
c-- Local variables

      integer lane,mopt,idiag,ifix,nlive0,mlive0,mlive,iter
     .      , jelf,isig,i1,i,j,j2,lw,lwl

      real*8 bka,usig,ssig,dummy 
c**   temporary to test TAH code
      character*1 symbl    
      character*256 message

                 
c-- Functions

   
c** Flag for debug printout
      integer*4 idebug
c      idebug = 0  : none 
c               1  : print to q- and log file decision values for fixed and not-fixed biases
c               2  : full print of logic 
      idebug = 1    
                      

c-- Print the criteria for wide-lane ambiguity resolution

      if( logprt ) write (6,10)
      write (10,10)
 10   format (
     . /,7x,'==== estimate and fix wide-lane (L2-L1) ambiguities =+==',
     . //,2x,'Algorithm: LC_HELP',/,
     .    2x,'L2-L1 biases estimated from phases and ionospheric ',
     .   'constraint using the decision',/,2x,'function',
     .   ' and chi-square search.  If a bias is not fixed using ',
     .   'the decision function',/,
     .   2x,'and P2 pseudo-range (PR) is available, the PR wide-lane ',
     .   'estimate is used with the ',/
     .   ,2x,'PR decision function criteria if ',/,
     .   7x,'|PR estimate - phase estimate| < 0.4 cycle',/,
     .   2x,'Finally, for biases fixed from phases,',
     .   ' the PR estimate overrides if',/,
     .   7x,'PR estimate is a factor of two closer to an integer than',
     .   ' the phase estimate.',/
     .   2x,'For codeless receivers, half-integer values are allowed.')   
      write (10,20)  
      if( logprt ) write (6,20)
 20   format(/,'---For resolution (only), WL and NL bias sigmas scaled '
     .  ,'by nrms of deviations from integer for each baseline',/)    

  
c-- Set the options for round_bias 
          
c     this routine called now only for WL with LC_HELP or LC_AUTCLN
       mopt = 2
c     this must be WL
       lane = 2
c     only decision function supported now
       noptin = 6

c-- Fix coordinate, atmospheric, and orbital parameters
  
      ifix = 1
      do i= 1,lpart
        if( free(i).ne.0 ) then
c         currently hard-wired to fix only 1 bias at a time
          nfix(ifix) = 1
          bdev(ifix) = adorg(i)-adjust(i)
          adjust(i) = adorg(i)
          free(i) = 0
          nlive0 = nlive
          nlive = nlive-1
          call bnew( nlive0,ifix )
        endif
      enddo     
          

c-- Compute chi2 from the normalized rms

      chi2=sclerr**2*dble(nobs-nlive)  


c-- Print the WL criteria and nrms   
                 
      bka = bkappa*1.d6   
      write(10,30) bka,sclerr,wldev,wlsig 
      if( logprt) write(6,30) bka,sclerr,wldev,wlsig
 30   format(/,'  LC_HELP  ion constraint ',f6.2,' ppm    nrms ',f6.2
     .    ,'    wldev ',f4.2,'    wlsig ',f4.2)
      if( logprt ) write(6,20) bka,sclerr,wldev,wlsig
      write (message,'(a,d9.2)')
     .     'Phase widelane nrms = ',sclerr
      call report_stat('STATUS','SOLVE','get_widelane',' ',message,0)   
      if( sclerr.gt.10.d0 ) 
     .   call report_stat('WARNING','SOLVE','get_widelane',' '
     .           ,'Phase widelane nrms > 10, relax ion constraint ',0)      

c-- Get the bias sigmas for computing scale factors and printout
         
      i = 0   
      j= lpart
      if(idebug.gt.1) print *,'RESOLVE_WL lpart l1bias ',lpart,l1bias
      lwl = l1bias
      do j2 = 1,lwl*iband    
        j = j + 1
        if (free(j).gt.0) then
          i = i+1
          idiag = i*(i+1)/2
          sigma(i) = dsign(dsqrt(dabs(a(idiag))),a(idiag))
        endif   
      enddo

      if( idebug.gt.1 ) print *,'RESOLVE_WL bias  sigma(1-50) '
     .    ,(sigma(i),i=1,50) 

c     Get and print the scale factors for sigmas, computed as the normalized rms of 
c     the deviation from an integer compared to the sigma for each baseline      

      call get_bias_scale( mopt )
                

c-- Print the bias estimates (WL mainly of interest but include L1 for completeness)

      if( logprt ) write(6,40)
      write(10,40)
 40   format(
     .  /,8x,'Label',18x,'Estimate',5x,'Uncertainty',13x,'Pseudorange',/
     .  ,31x,'(cycles)   Unscaled   Scaled',6x,'Estimate  Uncertainty')
      i = 0
      j = lpart
      lw = 0         
      if(idebug.ge.2) print *,'RESOLVE_WL lpart l1bias ',lpart,l1bias  
       do j2 = 1,lwl*iband
         j = j + 1
         if (j2.gt.lwl) lw = lw + 1
         if (free(j).eq.0) then
           symbl = ' '
           call qhead4(j2,j,dummy,dummy,symbl,lwl,lw,1)
         else  
           symbl = '*'
           i = i+1
           usig = sigma(i)   
           ssig = usig * bscale(i)                  
           if(idebug.gt.1) print *,'Bef qhead4 i bscale usig ssig '
     .          ,i,bscale(i),usig,ssig
           call qhead4(j2,j,usig,ssig,symbl,lwl,lw,2)
         endif   
      enddo


c-- Copy the fixed parameter adjustments into a temporary array

      i1 = lpart+1
      call copy1d(i1,ntpart,0,adjust,adorg)  

c---  Set the counter for live biases (msig= # live non-bias parameters + 1 )
c                                     ( see also BFWORK< GLOHED)
      msig = 1
                                                     
   
c-- Resolve the wide-lane biases     

      iter = 0    
c **  Temporary for using/testing TAH shortcut for fixing WL biases.
c **  Call NBIAS which calls TAH-modified version of NRBIAS.
c     --- loop over all biases, returning to statement 100 after each --
               
 100  mlive0 = nlive
      if(idebug.gt.1) then 
         print *,'RESOLVE_WL call ROUND_BIAS mopt lane nlive,nfix '
     .   ,mopt,lane,nlive
        write(*,'(10i5)') (nfix(i),i=1,5)   
      endif
      call round_bias( lane,mopt )  
      if( idebug.gt.1) print *,'  after ROUND_BIAS nlive ',nlive
*     MOD TAH 031202: iteration counter to report number of biases fixed
      iter = iter + 1
      mlive = nlive  

c     update solution with fixed integer biases
      ifix = mlive0 - mlive
* MOD TAH 031202: Only execute code below if ifix is > 0
*     (Added because goto 284 commented out).      
      if( idebug.gt.1 ) print *,'mlive0 mlive ifix ',mlive0,mlive,ifix
      if( ifix.gt.0 ) then 
         call bnew( mlive0,ifix,free,bdev,nfix )   
         do i = 1,mlive
*           MOD TAH 980310: Replaced calc by jelf call for HP compiler bug
            idiag = jelf(i+1)            !  = I*(I+1)/2
            sigma(i) = dsign(dsqrt(dabs(a(idiag))),a(idiag)) 
         enddo

c     set isigma array with slot of live parameters
c      print *,'RESOLVE_WL 3 ntpart isigma ',ntpart
c      write(*,'(10i7)') (isigma(i),i=1,ntpart)
      isig = 0
      do i=1,ntpart
        if( free(i).ne.0 ) then
           isig = isig+1
           isigma(isig) = i 
        endif
      enddo
c      print *,'RESOLVE_WL 4 ntpart isigma ',ntpart
c      write(*,'(10i7)') (isigma(i),i=1,ntpart)
      endif
*      MOD TAH 980310: Removed the ipuse test in goto so that
*     algorithm will be iterated whenever new biases are being fixed.
c     IF (MLIVE.LT.MLIVE0.and.ipuse.le.1) GOTO 123
      if( mlive.lt.mlive0 ) goto 100

c    --- end of loop on biases -----     
    
c-- Report the number of biases fixed and number of live/dead parameters
                                        

      if( logprt ) write(6,115) l1bias
      write(10,115) l1bias
  115 format(i5,' Phase ambiguities in solution')
                                           
c     Note: 'sclerr' is nrms but no longer used to scale the uncertainties'
      if( logprt ) write(6,120) iter-1,sclerr
      write(10,120) iter-1,sclerr
  120 format(i5,' WL ambiguities resolved (ion nrms = ',f6.2,')')
      if( logprt ) write(6,130)
      write(10,130)
  130 format(/,' Wide-lane bias-fixing complete')  
  
c     set number of biases fixed to 'iter' from TAH logic
      mfix = iter     

      return
      end


