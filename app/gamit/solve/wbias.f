      Subroutine wbias( ia1,last_bias )

c       ia1 :  Index in full normal equation matrix (A) of last non-bias parameter
c       ndimb:  Dimension of bias matrix (set = maxwm1 (dimpar.h) in calling routine)
          
      include '../includes/dimpar.h'                                       
      include 'solve.h'

c        INTEGER MAXWM1
c        PARAMETER (MAXWM1=(MAXOBS*(MAXOBS+1)/2))
                  
      integer*4 ia1,last_bias,j,jj,ind,numbias,k

      real*8 wgtbias(maxwm1)

              
c     get the size of the bias matrix
      numbias=(last_bias-ia1)/iband   
     
c      print *,'WBIAS maxwm1 ia1 last_bias iband numbias '
c     .        , maxwm1,ia1,last_bias,iband,numbias

c The constraint now set in batch file and stored in common

c  Fill the weight matrix
      
      do j=1,numbias
        jj=(j*j-j)/2
        do k=1,j
          ind=k+jj    
          if( ind.gt.maxwm1 ) then
           print *,'j jj k ind ',j,jj,k,ind    
c            stop 4
          endif
          wgtbias(ind)=0.d0
          if( j.eq.k) then
            wgtbias(ind) = 1.d0/bias_apr**2
          endif
        enddo
      enddo   
c      print *,'Bias weights 1-6 ',(wgtbias(j),j=1,6)

c  Add weight matrix to normal matrix 
           
c    A matrix (ia1+1  to  numbias ) and
c    AN22 matrix (1 to numbias)
c      call addwgt( ia1,maxwm1,numbias,wgtbias )  
c      print *,'ia1 maxwm1 numbias ',ia1,maxwm1,numbias 
                              
c      AN22 matrix ( 1 to numbias ) 
      call addwgtb( ia1,maxwm1,numbias,wgtbias )  
     
      return
      end
        

 
