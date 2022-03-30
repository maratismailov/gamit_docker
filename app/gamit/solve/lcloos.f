      Subroutine LCLOOS( iop )
C
C     LC solution with loose sit and sat constraints
C     iop = 3: LC mode, bias free, loose constraint solution
C     iop = 4: LC mode, bias fixed, loose constraint solution
c
      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'
      include 'parameters.h'
c
      integer fre2(maxprm),iop
      logical zeroad


      equivalence (isigma,fre2)
              
      integer*4 nded,nlive0,ifast,ierinv,k0,k1,k2,k3
     .        , ia1,ia2,ifix,imb,imb1,iarray,i

      real*8 rnew,goodft
         
      logical debug/.false./

c       Biases-free solution

      if( iop.eq.3 ) then
 
c---- write a header to screen  and  Q-file

      call report_stat('STATUS','SOLVE','lcloos',' '
     .          , 'Performing LC biases-free loose solution',0)  
      write(10,'(/,a,/)') 'Performing LC biases-free loose solution'
                        
c---- new station constraints
      if ( sitest ) call lwstat

c---- new satellite constraints
      if ( satest ) call lwsat

c----- new zenith-delay contraints
      if( zenest ) call lwzen 

c----- atmospheric gradient constraints. iflag = 1 N/S gradient, = 2 E/W gradient
      if( gradest) then
          call lwgrad 
      endif
                      
c----- new earth orientation constraints
      if( eopest ) call lweop

c---- copy LC mode normal matrix
      call lcnorm(ierinv,1)

c---- replace free array (bias part)        
c      print *,'LCLOOS 3 idxb(1) ',idxb(1)
      k0 = idxb(1)
      if (k0.lt.0) k0 = -k0
      k1 = k0
      k3 = 0
c---- check L1 biases                  
         k2 = l1bias        
         ia1 = k1
         ia2 = k1+k2-1       
c         print *,'k2 ia1 ia2 free(211-270) ',k2,ia1,ia2
c         write(*,'(10i7)') (free(i),i=211,270)
         do i = ia1,ia2
            fre2(i) = free(i)
            adorg(i) = adjust(i)
            adjust(i) = 0.0d0
            k3 = k3+1
            if (free(i).eq.0) then   
c             if the parameter number is in the nlres array, the bias
c             was resolved (rather than fixed as dependent)
              if( iarray(i,nlres,k2).ne.0 ) then
c                print *,'found in iarray i k2 ',i,k2
                nlive = nlive + 1
                idxb(k3) = iabs(idxb(k3))  
                free(i) = 1 
              endif
            endif
         enddo
         k3 = k3+k2
c      print *,'k3 idxb free(211-270) ',k3
c      write(*,'(10i7)') (idxb(i),i=1,80) 
c      write(*,'(10i7)') (free(i),i=211,270)    
c      print *,' fre2 '
c      write(*,'(10i7)') (fre2(i),i=211,270)    
      nded = ntpart-nlive

      call lcnorm(ierinv,2)
      if( ierinv.ne.0 ) call report_stat('FATAL','SOLVE','lcloos',' '
     .                           ,'Inversion error in LCNORM(2)',0)

c---- get LC loose constraint solution
      call solvlc(adjust,nded,rnew,zeroad)
c---- in the case of x2=0, all results have been obtained in solvlc.
      ifast=0
      if ( .not.zeroad ) call solve3(adjust,nded,rnew,ifast)
      goodft = chi2/dble(nobs-nlive)
      goodft = dabs(goodft)
      sclerr = dsqrt(goodft)
      if(debug) print *,'LCLOOS nobs nlive chi2 goodft sclerr '
     .                        , nobs,nlive,chi2,goodft,sclerr 

c      Biases-fixed solution

      elseif ( iop.eq. 4 ) then
      
c --- write a header to the screen and q-file
      call report_stat('STATUS','SOLVE','lcloos',' '
     .          , 'Performing biases-fixed loose solution',0)   
      write(10,'(/,a)') 'Performing LC biases-fixed loose solution'

c---- replace biases by original fixed value   
c      print *,'LCLOOS 4 idxb(1) ',idxb(1)
      k0 = idxb(1)
      if (k0.lt.0) k0 = -k0
      k1 = k0
      ifix = 1
      k3 = 0
         k2 = l1bias        
         ia1 = k1
         ia2 = k1+k2-1
c         print *,'k2 ia1 ia2 free(211-270) ',k2,ia1,ia2
c         write(*,'(10i7)') (free(i),i=211,270)
         imb1 = 0
         do i = ia1,ia2
            k3 = k3+1
            if (free(i).eq.1.and.fre2(i).eq.0) then
               imb = msig+imb1
               nfix(ifix) = imb
               bdev(ifix) = adorg(i)-adjust(i)
               adjust(i) = adorg(i)
               free(i) = 0
               idxb(k3) = -iabs(idxb(k3))
               nlive0 = nlive
               nlive = nlive-1
               call bnew( nlive0,ifix )
            endif
            if (free(i).eq.1) imb1 = imb1 + 1
         enddo
         k3 = k3+k2
c      print *,'k3 idxb free(211-270) ',k3
c      write(*,'(10i7)') (idxb(i),i=1,80) 
c      write(*,'(10i7)') (free(i),i=211,270)


      chi2 = sclerr**2*dble(nobs-nlive)

      endif
c-------------------------------------------------------------
c
      return
      end


