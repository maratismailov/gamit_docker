Copyright 1993 Massachusetts Institute of Technology and The Regents of the
c              University of California, San Diego.  All Rights Reserved.

      SUBROUTINE WZEN 

c     Add zenith delay weights to normal matrix.

c     Written R King 1993
c     Modified by Simon McClusky for Gauss Markov weighting June 1994 
c     Modified by R King and T Herring for average zenith delay and a new 
c       formulation of the Gauss Markov process  May 2004
c
      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'
      include 'parameters.h'
           
      character*3 upperc     

      integer maxc,indx,ind    

      integer*4 i,k,n,icnt,ier  ! Removed j, iadd

 
      parameter(maxc=maxatm*(maxatm+1)/2)

      real*8 covatm(maxc),be(maxc),u(maxsit),bterm(maxsit)
     .     , delta(maxsit),dtime(maxsit),lterm(maxsit),sum,rcond
             
c     debug
c     integer*4 ict
c     real*8 temp1(maxc),temp(maxc),


* MOD TAH 040406: Added quanities for Kalman filter to set mean
*     of fogm process to zero.
      real*8 kgain(maxatm)    ! Kalman Gain vector
     .,      atva             ! Variance term ATVA+bterm
     .,      covzmean(maxc)   ! Covariance matrix with near zero variance of mean
      integer*4 jel           ! Function to return index in lower triangular matrix 


c     debug
c     integer*4 ict
c     real*8 temp1(maxc),temp(maxc),
c     WRITE(6,999)
c 999 FORMAT(' WEIGHTING ZENITH DELAY PARAMETERS')

     
c First weight the average values
                                                   
c     atmospheric parameters follow the station coordinates
      indx =  3*nsite  
      do i=1,nsite  
         covatm(1) = 1.d0/zen_apr(i)**2    
         call addwgt(indx,1,1,covatm)   
         indx = indx + 1
      enddo
      
c     Now weight the piecewise-linear function if used (# assumed the 
c     same for every station)

      if( nzen.gt.1 ) then

         icnt=0
         do i=1,nsite
            nzen=nzen
            indx= 4*nsite+icnt

c           Directly generate the covariance matrix and then apply the 
c           condition that the mean of the process should be zero with
c           a small uncertainity (1 mm sigma for the mean is used to
c           avoid a singular covariance matrix when it is inverted).

c           The covariance function for a first-order Gauss Markov process
c           is given by:
c           R(t) = LT_var*exp(-(t/tau))
c           where LT_var is long term variance (lterm(i) in code)
c                 t is absolute time difference
c                 tau is correlation time
c           Since the time increments are constant, we can replace t with
c           (n-k)*delta_t and since mulitplication in the exp term is the
c           same has taking exp(delta_t/tau) to the (n-k) power, the final
c           direct evaluation becomes:
c           R(t) = R((n-k)*delta-t) = lterm(i)*u(i)**(n-k)
c           where u(i) = exp(delat_t/tau(i))
c
            ind = 0
            sum = 0.d0
c
c           compute time span of each atmospheric parameter
            if (upperc(zenmod).eq.'PWL') then
               dtime(i)=(((iend+nzen-2)-(istart-1))/
     .                        (nzen-1))*inter/3600.d0
            else
               dtime(i)=(((iend+nzen-1)-(istart-1))/nzen)*i
     .                                             nter/3600.d0
            endif
c
c           compute point to point variances (delta) from m/sqrt(hr) rms.
c           move overall rms to variance
            delta(i) = (zen_mar(i)**2)*dtime(i)  

c           compute u(i) for particular station
            if ( zen_tau(i) .lt. 0.001 ) then
                u(i) = 0.d0
            else
                u(i) = exp(-dtime(i)/zen_tau(i))
            endif
c
c           compute long term variance lterm(i)
            lterm(i) = delta(i)/(1 - u(i)**2)

c           print *, 'zen_apr,zen_mar,zen_tau'
c     .          ,zen_apr(i),zen_mar(i),zen_tau(i)
c           print *, 'delta,dtime,u',delta(i), dtime(i),u(i)
c           print *, 'lterm',lterm(i)

* MOD TAH 040406: Directly generate the covariance matrix.  Also compute
*           the sum of all elements.  This value is needed to apply the 
*           condition that the mean of the process be near-zero
            ind = 0
            sum = 0
            do n = 1, nzen
              do k = 1,n
                  ind = ind + 1
*                 Covariance matrrix is lterm(i)*exp(-dt*(n-k)/tau)
                  covatm(ind) = lterm(i)*u(i)**(n-k)
                  if( k.eq.n ) then
                      sum = sum + covatm(ind)
                  else     ! For off-diagonal elements need to add twice
                      sum = sum + 2*covatm(ind)
                 endif  
              end do
            end do

****        Now apply the constraint that the mean is nearly zero 
*           variance (+- 1mm)
*           To apply zero mean constraint we use a Kalman filter observation
*           equation that says sum of estimates should be zero +- sqrt(bterm(i)).
*           In Matlab: the code is. K is Kalman gain vector; 
*           bterm = 1e-6;
*           K = sum(covatm)'/(sum(sum(covatm))+bterm);
*           covzmean is the covariance matrix with near zero variance of its
*           mean.  (The bterm set how small the variance is.
*           covzmean = covatm - K*sum(covatm)
*           Fortran code below implements the Matlab code above.
            bterm(i) = 1.d-6    ! Mean +- 1 mm
            atva = sum + bterm(i)
*           Form the Kalman gain vector (vector in this case because only one "obs")
            do n = 1,nzen
               sum = 0
               do k = 1,nzen
                  ind = jel(n,k)
                  sum = sum + covatm(ind)
              end do
              kgain(n) = sum/atva
            end do
*****       Now generate the near-zero variance mean covariance matrix
            do n = 1,nzen
*               form the sum of the columns of the covariance matrix
               sum = 0
               do k = 1, nzen
                 ind = jel(k,n)
                 sum = sum + covatm(ind)
               end do
*              Now compute the new covariance matrix
               do k = 1,n
                  ind = jel(k,n)
                  covzmean(ind) = covatm(ind) - kgain(k)*sum
               end do
            end do
****        Finally, copy the zero-mean covariance back to covatam
            ind = 0
            do n = 1,nzen
               do k = 1, n
                  ind = ind + 1
                  covatm(ind) = covzmean(ind)
               end do
            end do
c
c            WRITE(6,210)
c  210       FORMAT(/,'Covariance matrix for zenith delays')
c            ICT=0
c            DO 240 II=1,nzen
c              DO 220 JJ=1,II
c                 ICT=ICT+1
c                 **ucommenet declaraion if this used
c                 TEMP(JJ)=COVATM(ICT)
c  220         CONTINUE
c              WRITE(6,230) (TEMP(IJ),IJ=1,II)
c  230         FORMAT(1X,9D12.5)
c  240       CONTINUE
C
C            invert to get weight
C
            CALL INVER2(COVATM,BE,1,nzen,rcond,IER)
            if( ier.ne.0 )call report_stat('WARNING','SOLVE','wzen',' '
     .,      'Bad inverse of a priori zenith-delay covariance matrix',0)
c
c            WRITE(6,250)
c  250       FORMAT(/,'Weight matrix for zenith delays')
c            ICT=0
c            DO 280 II=1,nzen
c               DO 260 JJ=1,II
c                  ICT=ICT+1
c                  **uncomment declaration if this used
c                  TEMP1(JJ)=COVATM(ICT)
c  260         CONTINUE
c              WRITE(6,270) (TEMP1(IJ),IJ=1,II)
c  270         FORMAT(1X,9D12.5)
c  280       CONTINUE
C
C           Add weight matrix to normal matrix                 

            call addwgt( indx,maxc,nzen,covatm )

c            do j=1,nzen
c               JJ1=(J*J-J)/2
c               JJ2=J+INDX
c               JJ3=(JJ2*JJ2-JJ2)/2
c               DO K=1,J
c                  IND1=K+JJ1
c                  IND2=K+INDX+JJ3
c                  a(IND2)=a(IND2)+COVATM(IND1)
c                 ALC(IND2)=ALC(IND2)+COVATM(IND1) 
c                 print *,'indx ind1 ind2 alc ',indx,ind1,ind2,alc(ind2)
c               end do
c            end do         
C
             icnt=icnt+nzen

         end do
c     end if on nzen > 1
      endif

      RETURN
      END
