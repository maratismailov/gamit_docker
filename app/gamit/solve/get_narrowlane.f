
      Subroutine GET_NARROWLANE ( mopt )

c     R. King from D. Dong's original BISOPT and NBIAS.  Eliminates all bias-fixing
c     options except rounding with decision function.  Called from SOLVE.
c     

c     Narrow_lane ambiguity resolution for L1_ONLY, L2_ONLY, LC_HELP, or LC_AUTCLN  mopt=1
c     Narrow_lane ambiguity resolution for L1L2_INDEPENDENT - mopt=3
          
c     Control for searching, set according to method:
c      mopt = 1  search L1 (or L2) biases only
c           = 3  search both L2-L1 and L1 biases (not used??)


      implicit none
                     
c-- Global parameters in common

      include '../includes/dimpar.h'
      include 'solve.h'  
      include 'parameters.h'
                
c-- Local variables

      integer*4 mopt,lane,idiag,mlive0,isig,ifix,mlive,iter,i
            

c-- External function    

      integer jelf
 

      
c-- Report what we're doing and the criteria used

      call report_stat('STATUS','SOLVE','get_narrowlane',' '
     .          , 'Resolving narrow-lane ambiguities',0)
      if( logprt ) write(6,10) nldev,nlsig,nlcut,nlrat,nldmax
      write(10,10) nldev,nlsig,nlcut,nlrat,nldmax
   10 format(/,' Narrow-lane bias-fixing criteria: ',
     .       ' deviation : ',f5.2,'  sigma : ',f5.2,
     .       '  decision func. : ',f8.1,' ratio : ',f6.1,/,
     .        36x,'maximum distance :',f8.1)    
c     narrowlane
      lane = 1
                 
c       print *,'GET_NARROWLANE idxb '
c       write(*,'(10i7)') (idxb(i),i=1,100)


c-- Get the scale factors for sigmas, computed as the normalized rms of 
c   the deviation from an integer compared to the sigma for each baseline

      call get_bias_scale( mopt ) 


c     --- loop over all biases, stopping when the next cannot be fixed -----

      iter = 0                      

c     initialize the array indicating to LCLOOS and LOOS12 whether
c     a bias was resolved or was fixed because it was dependent or no obs
      do i=1,maxdd
        nlres(i) = 0
      enddo
       
 100  mlive0 = nlive
c      print *,'GET_NARROWLANE call ROUND_BIAS lane nlive,nfix '
c     .    ,lane,nlive
c      write(*,'(10i5)') (nfix(i),i=1,5)   
c**RWK: this superseded by ROUND_BIAS:      CALL NBIASR(lane,SIGMA,FREE,MOPT,NOPTIN,NFIX,BDEV)
      call round_bias( lane,mopt )  
c      print *,'  after ROUND_BIAS nlive ',nlive
      iter = iter + 1
      mlive = nlive  

c     update the solution with the fixed integer bias
      ifix = mlive0 - mlive
*     MOD TAH 031202: Only execute code below if ifix is > 0
*     (Added because goto 284 commented out).
      if( ifix.gt.0 ) then 
        call bnew( mlive0,ifix )
        do i = 1,mlive  
*         MOD TAH 980310: Replaced calc by jelf call for HP compiler bug
          idiag = jelf(i+1)            !  = I*(I+1)/2
          sigma(i) = dsign(dsqrt(dabs(a(idiag))),a(idiag))
        enddo

c     set the isigma array with indices of the live parameters
c      print *,'GET_NARROWLANE 1 ntpart isigma ',ntpart
c      write(*,'(10i7)') (isigma(i),i=1,ntpart)
      isig = 0
      do i = 1,ntpart
        if( free(i).ne.0 ) then 
           isig = isig + 1
           isigma(isig) = i
        endif
      enddo
c      print *,'GET_NARROWLANE 2 ntpart isigma ',ntpart
c      write(*,'(10i7)') (isigma(i),i=1,ntpart)
      endif
*     MOD TAH 980310: Removed the ipuse test in goto so that algorithm will 
*     be iterated while ever new biases are being fixed.      
c      ipuse = 0
C     IF (MLIVE.LT.MLIVE0.and.ipuse.le.1) GOTO 123
      if( mlive.lt.mlive0 ) goto 100

c     --- end of loop over biases -------------------------


c-- Report the number of biases fixed and number of live/dead parameters

      if( logprt ) write(6,120) iter-1
      write(10,120) iter-1
  120 format(i5,' NL ambiguities resolved')
      nlive = 0
      do i = 1,ntpart
        if(free(i).eq.1) nlive = nlive + 1
      enddo
      if( logprt ) write(6,130)
      write(10,130)
  130 format(/,' Narrow-lane bias-fixing complete')  
      call qhead1 (2)

      return
      end
