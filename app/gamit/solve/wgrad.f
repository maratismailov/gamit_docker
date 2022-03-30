Copyright 1998 Massachusetts Institute of Technology and The Regents of the
c              University of California, San Diego.  All Rights Reserved.

      Subroutine WGRAD
      
c     Add a priori atmospheric gradient weights to normal matrix 
c     Use full VCV of weights to represent a first order Gauss-Markov 
c     process to model between parameter variances.

c     Written by : R. King 1993, S. McClusky 940611; R. King 980917 

      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'
      include 'parameters.h'
             
      integer*4 maxc
      parameter(maxc=(maxgrad/2)*((maxgrad/2)+1)/2)
       
      integer*4 indx,ind,ind2,iadd,itype
     .        , ncoord,navgzen,jgrad,jzen,n,ier,icnt,i,j,k
c     debug                 
c     .        , ict

      real*8  covgrad(maxc),be(maxc),u(maxsit),bterm(maxsit)
     .     ,  delta(maxsit),dtime(maxsit),lterm(maxsit),sum,rcond
c     debug
c    .     , temp1(maxc),temp(maxc)

      character*3 upperc

c      print *,'WGRAD: Weighting gradient delay parameters'
        
c  Skip over station coordinates and zenith-delay parameters
      ncoord = 3*nsite  
      navgzen = nsite    
      jzen = nsite*nzen
c      print *,' ncoord navgzen jzen idtgrad '
c      ,ncoord,navgzen,jzen,idtgrad
                                 
c  Loop over N/S and E/W gradients

      do 1000 itype = 1, 2

c  If this is E/W gradients, skip over N/S gradients
      jgrad = 0              
      if( itype.eq.2 ) jgrad = ngrad * nsite

c  Loop over all gradient parameters

      icnt=0
      do 500 i=1,nsite
        indx=ncoord+navgzen+jzen+jgrad+icnt
c        print *
c        ,'WGRAD: ngrad,ncoord,navagzen,jgrad,icnt'
c     .             ,ngrad,ncoord,navgzen,jgrad,icnt

c       Constant (single) gradient in each direction

        if( ngrad.eq.1 ) then
c         compute the diagonal covariance elements 
          covgrad(1) = grad_apr(i,itype)**2
c         invert to get weight
          if(covgrad(1).gt.0.d0) then
            covgrad(1) = 1.d0/covgrad(1)
          else
            covgrad(1) = 1.d0
          endif
c         print*,'covgrad: ',covgrad(1)
c         add gradient weight matrix to lower triangle normal matrix
c         indx is last element (used for multiple gradients); need current one here
          indx = indx + 1    
          ind2 = indx+(indx*indx-indx)/2
          a(ind2)=a(ind2)+covgrad(1)
          alc(ind2)=alc(ind2)+covgrad(1)
c          print *,'indx,ind2,covgrad',indx,ind2,covgrad(1)

c       Piecewise linear model with stochastic constraints

        else

c       A first-order Gauss Markov process is used to calculate VCV for gradient
c       parameters.  The process is defined as Ai+1 = Ai*e**(-dtime/tau)+Wi
c       where: A     = first order Gauss Markov process
c             dtime  = time between gradient estimates
c             tau    = correlation time of the process
c             Wi     = white noise
c
c      Letting: U = e**(-dtime/tau)+Wi
c
c      Therefore: A =  --                          --
c                      | 1     0     0     0    ... |
c                      | U     1     0     0    ... |
c                      | U**2  U     1     0    ... |
c                      | U**3  U**2  U     1    ... |
c                      |  .     .    .     .    ... |
c                      |  .     .    .     .    ... |
c                      |  .     .    .     .    ... |
c                      --                          --
c      Now: Apriori gradient variance matrix Q is given by
c                 Q =  --                          --
c                      | L     0     0     0    ... |
c                      | 0     D     0     0    ... |
c                      | 0     0     D     0    ... |
c                      | 0     0     0     D    ... |
c                      |  .     .    .     .    ... |
c                      |  .     .    .     .    ... |
c                      |  .     .    .     .    ... |
c                      --                          --
c      Where: L     = the long term variance of the gradients
c             D     = the point to point variance of the gradients
c
c                 B =  --                          --
c                      | B     B     B     B    ... |
c                      | B     B     B     B    ... |
c                      | B     B     B     B    ... |
c                      | B     B     B     B    ... |
c                      |  .     .    .     .    ... |
c                      |  .     .    .     .    ... |
c                      |  .     .    .     .    ... |
c                      --                          --
c
c     Where: B = the overall level of the gradients
c
c    The apriori gradient VCV matrix is given by: B+A*Q*A_T ( A_T is the transpose of A )
c      Where: A*Q*A_T = --                                                            --
c                           | L          LU          LU**2           LU**3             ... |
c                           | LU      LU**2+D      LU**3+DU        LU**4+DU**2         ... |
c                           | LU**2   LU**3+DU    LU**4+DU**2+D   LU**5+DU**3+DU       ... |
c                           | LU**3  LU**4+DU**2  LU**5+DU**3+DU  LU**6+DU**4+DU**2+D  ... |
c                           | .          .             .               .               ... |
c                           | .          .             .               .               ... |
c                           | .          .             .               .               ... |
c                           --                                                            --
c
c    This matrix can be expressed as:
c                                                              k-1
c                                                              ---
c                         A*Q*A_T(n,k) = LU**((n-1)+(k-1)) + D \  U**((n-k)+(0,2,4,6,....) + B
c                                                              /
c                                                              ---
c                                                              j=1
c   Where n = row number
c         k = column number
c
c
        ind = 0
        sum = 0.d0

c       Compute time span of each gradient parameter
        if (upperc(gradmod).eq.'PWL') then
        dtime(i)=(((iend+ngrad-2)-(istart-1))/(ngrad-1))
     .             *inter/3600.d0
        else
          dtime(i)=(((iend+ngrad-1)-(istart-1))/ngrad)
     .             *inter/3600.d0
        endif

c       Compute point to point variances (delta) from m/sqrt(hr) rms.
c       Move overall rms to variance
        delta(i) = (grad_mar(i,itype)**2)*dtime(i)
        bterm(i) = (grad_apr(i,itype)**2)

c       Compute u(i) for particular station
        if ( grad_tau(i,itype) .lt. 0.001 ) then
          u(i) = 0.d0
        else
          u(i) = exp(-dtime(i)/grad_tau(i,itype))
        endif

c       Compute long term variance lterm(i)
        lterm(i) = delta(i)/(1 - u(i)**2)

c        print *,'grad_apr,grad_mar,tau'
c     .        ,grad_apr(i,itype),grad_mar(i,itype),grad_tau(i,itype)
c        print *,'bterm,delta,dtime,u',bterm(i),delta(i),dtime(i),u(i)
c        print *, 'lterm',lterm(i)

        do  n = 1,ngrad
          do  k = 1,n
            iadd = 0
            sum = 0.d0
            do  j = 1,(k-1)
               sum = sum + delta(i)*u(i)**((n-k)+iadd)
               iadd = iadd + 2  
            enddo
               ind = ind+1
               covgrad(ind) = (lterm(i)*(u(i)**((n-1)+(k-1)))) + sum
     .                        + bterm(i)
c            print *,'n,k,ind,covgrad(ind)',n,k,ind,covgrad(ind)
          enddo
        enddo
c
c       WRITE(6,210)
c  210  FORMAT(/,'Covariance matrix for gradith delays')
c       ICT=0
c       DO 240 II=1,ngrad
c         DO 220 JJ=1,II
c           ICT=ICT+1
c           **ucommenet declaration if this used
c           TEMP(JJ)=COVGRAD(ICT)
c  220    CONTINUE
c         WRITE(6,230) (TEMP(IJ),IJ=1,II)
c  230    FORMAT(1X,9D12.5)
c  240  CONTINUE


C      Invert to get weight
       call inver2(covgrad,be,1,ngrad,rcond,ier)
       if( ier.ne.0 )call report_stat('WARNING','SOLVE','wgrad',' '
     .    , 'Bad inverse of a priori gradith-delay covariance matrix',0)
c       print *,'weight ',(covgrad(k),k=1,6)
c
c       WRITE(6,250)
c  250  FORMAT(/,'Weight matrix for gradient delays')
c       ICT=0
c       DO 280 II=1,ngrad
c         DO 260 JJ=1,II
c           ICT=ICT+1
c           **uncomment declaration if this used
c           TEMP1(JJ)=COVGRAD(ICT)
c  260    CONTINUE
c         WRITE(6,270) (TEMP1(IJ),IJ=1,II)
c  270    FORMAT(1X,9D12.5)
c  280  CONTINUE
C
          
c      Add weight to normal matrix

      call addwgt( indx,maxc,ngrad,covgrad )

c      do 300 j=1,ngrad
c      JJ1=(J*J-J)/2
c      JJ2=J+INDX
c      JJ3=(JJ2*JJ2-JJ2)/2
c      DO 300 K=1,J
c      IND1=K+JJ1
c      IND2=K+INDX+JJ3
c      a(IND2)=a(IND2)+COVGRAD(IND1)
c      ALC(IND2)=ALC(IND2)+COVGRAD(IND1)
c      print *,'indx ind1 ind2 alc ',indx,ind1,ind2,alc(ind2)
C     WRITE(6,*) IND1,IND2,COVGRAD(IND1)
c  300 CONTINUE
        
c-----end if on constant or piecewise linear
      endif    

      icnt=icnt+ngrad

c----- end loop on stations
  500 continue

c-----end loop on N/S and E/W gradients
 1000 continue

      return
      end

