Copyright 1993 Massachusetts Institute of Technology and The Regents of the
c              University of California, San Diego.  All Rights Reserved.
c
      Subroutine LWGRAD 
c
c     Written by: Simon McClusky 9/30/1996; modified for multiple gradients R. King 980917
c
c     Add loose gradient weight increment to normal matrix (see WGRAD for model)

      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'
      include 'parameters.h'
c
      integer*4 itype,indx,ind2,i,j,ncoord,navgzen,icnt,jzen,jgrad
     .         , iadd,ind,ict,ii,jj,k,n,ier
      character*3 upperc
                    
      integer*4 maxc
      parameter(maxc=(maxgrad/2)*((maxgrad/2)+1)/2)   

      real*8 u(maxsit),bterm(maxsit),delta(maxsit),dtime(maxsit),sum
     .     , u2(maxsit),bterm2(maxsit),delta2(maxsit),sum2,lterm(maxsit)
     .     , lterm2(maxsit)
     .     , covgrad(maxc),covgradt(maxc),covgrad2(maxc),be(maxc)
     .     , temp(maxc)
c    debug:
c     .     , temp1(maxc),temp2(maxc),temp3(maxc)

      real*8 m2cyc,at10deg,rcond

c     initialization to avoid compiler warning
      indx = 0
          
c     compute the scaling constant
      m2cyc = 1/(299792458.d0/1.57542d9)
c Changed as the G77 compiler does not support intrinsic trig functions working with degrees 
c      at10deg = 1/(sind(10.d0)*tand(10.d0)+0.003)
      at10deg = 1/(sin(10.d0/180.d0*pi)*tan(10.d0/180.d0*pi)+0.003)
 

c     print and write to the Q-file the loose-solution constraints
c     -- use different formats for single and multiple gradients
         
      if( ngrad.gt.1 ) then
        if( logprt ) write(6,160)
        write(10,160)
      else
        if( logprt ) write(6,162)
        write(10,162)  
      endif
      do i=1,nsite   
        if( ngrad.gt.1 ) then
          if( logprt )write( 6,161) i,rlabel((i-1)*3+1)(1:4),sitnam(i)
     .      ,ngrad
     .      ,( grad_apr2(i,j)*at10deg/m2cyc,grad_mar2(i,j)*at10deg/m2cyc
     .        ,grad_tau2(i,j),j=1,2)
          write(10,161) i,rlabel((i-1)*3+1)(1:4),sitnam(i),ngrad
     .      ,( grad_apr2(i,j)*at10deg/m2cyc,grad_mar2(i,j)*at10deg/m2cyc
     .         ,grad_tau2(i,j),j=1,2)  
        else
          if( logprt ) write( 6,163) i,rlabel((i-1)*3+1)(1:4),sitnam(i)
     .                  ,(grad_apr2(i,j)*at10deg/m2cyc,j=1,2)
          write( 10,163) i,rlabel((i-1)*3+1)(1:4),sitnam(i)
     .                  ,(grad_apr2(i,j)*at10deg/m2cyc,j=1,2)
         endif
      enddo
  160 format(/,
     .    3x,'A priori atmospheric gradient error at 10 degrees'
     .   ,' elevation angle',/,
     .   'Station                 # N/S: A priori(m) '
     .   ,' Mar(m/sqrt(hr)) Correl(hr)   E/W: A priori(m) '
     .,  ' Mar(m/sqrt(hr)) Correl(hr)',/)         
  161     format(i3,1x,a4,2x,a12,1x,i3,2(6x,f8.5,f16.5,f12.1))
  162 format(/,
     .    3x,'A priori atmospheric gradient error in meters at ',
     .       '10 degrees elevation angle',/,
     .       'Station                  North-South  East-West',/)
  163     format(i3,1x,a4,2x,a12,3x,2(3x,f8.5)) 
             
 
c  assumption is that atmospheric gradient parameters are directly after 
c  multiple zenith delays in the parameter list

c  Skip over station coordinates and zenith-delay parameters
      ncoord = 3*nsite 
      navgzen = nsite         
      jzen = nsite*nzen
                            
c  Loop over N/S and E/W gradients

      do 1000 itype=1,2

c  If this is E/W gradients, skip over N/S gradients
      jgrad = 0
      if( itype.eq.2 ) then 
         jgrad = nsite*ngrad
      endif

c  Loop over stations  
      icnt=0
      do 500 i = 1,nsite 
        indx=ncoord+navgzen+ngrad+icnt
c        print *
c     . ,'LWGRAD: ngrad,ncoord,navgzen,mzen,ngrad,icnt'
c     .             ,ngrad,ncoord,navgzen,nzen,mrad,icnt
                
c       Constant gradient in each direction

        if( ngrad.eq.1 ) then

c         compute the diagonal covariance elements 
          covgrad(1) = grad_apr2(i,itype)**2   
c          print *, 'LWGRAD grad_apr2 ',covgrad(1)
c         invert to get weight
          if(covgrad(1).gt.0.d0) then
            temp(1) = grad_apr2(i,itype)/grad_apr(i,itype)
          else
            temp(1) = 1.d0
          endif    
c          print *, 'LWGRAD grad_apr grad_apr2 ',covgrad(1)
          

c          weight increment = (new weight) - (old weight)
c                           = 1.0/(new variance) - 1.0/(old variance)
c                           = (1.0 - temp*temp)/(new variance)
c          new weight = old weight + weight increment 

          covgrad(1) = (1.d0 - temp(1)*temp(1))/covgrad(1)

c         add gradient weight matrix to lower triangle normal matrix 
c         indx is last element (used for multiple gradients); need current one here
          indx = indx + 1   
          ind2 = indx+(indx*indx-indx)/2
          a(ind2)=a(ind2)+covgrad(1)
          alc(ind2)=alc(ind2)+covgrad(1)
c          print *,'indx,ind2,covgrad',indx,ind2,covgrad(1)

c       Piecewise linear model

        else 

c         Compute time span of each gradient parameter 
c          print*,'gradmod iend ngrad istart inter '
c     .          , gradmod,iend,ngrad,istart,inter
          if (upperc(gradmod).eq.'PWL') then
            dtime(i)=(((iend+ngrad-2)-(istart-1))/(ngrad-1))
     .             *inter/3600.d0
          else
            dtime(i)=(((iend+ngrad-1)-(istart-1))/ngrad)
     .             *inter/3600.d0
          endif

c         Compute point to point variances (delta) from long period sigma.
c         Move overall sigma to variance
          delta(i) = (grad_mar(i,itype)**2)*dtime(i)
          delta2(i) = (grad_mar2(i,itype)**2)*dtime(i)
          bterm(i) = grad_apr(i,itype)**2
          bterm2(i) = grad_apr2(i,itype)**2

c         Compute u(i) for particular station
          if ( grad_tau(i,itype) .lt. 0.001) then
            u(i) = 0.d0
          else
            u(i) = exp(-dtime(i)/grad_tau(i,itype))
          endif
          if ( grad_tau2(i,itype) .lt. 0.001) then
             u2(i)=0.d0
          else
             u2(i) = exp(-dtime(i)/grad_tau2(i,itype))
          endif

c         Compute long term variance lterm(i)
          lterm(i) = delta(i)/(1 - u(i)**2)
          lterm2(i) = delta2(i)/(1 - u2(i)**2)
             

c          print *,'itype grad_apr grad_mar tau '
c     .     ,itype,grad_apr(i,itype),grad_mar(i,itype),grad_tau(i,itype)
c          print *,'delta,bterm,dtime,u ',delta(i),bterm2(i),u(i)
c          print *,'lterm ',lterm(i)
c          print *,'itype grad_apr2,grad_mar2,tau2'
c     .  ,itype,grad_apr2(i,itype),grad_mar2(i,itype),grad_tau2(i,itype)
c          print *,'delta2,bterm2,dtime,u2 '
c     .        ,delta2(i),bterm2(i),dtime(i),u2(i)
c          print *,'lterm2',lterm2(i)
          ind = 0
          do  n = 1,ngrad
            do  k = 1,n
              iadd = 0
              sum = 0.d0
              sum2 = 0.d0
              do  j = 1,(k-1)
                 sum = sum + delta(i)*u(i)**((n-k)+iadd)
                 sum2 = sum2 + delta2(i)*u2(i)**((n-k)+iadd)
                 iadd = iadd + 2
              enddo
               ind = ind+1
               covgradt(ind) = (lterm(i)*(u(i)**((n-1)+(k-1)))) + sum
     .                        + bterm(i)
               covgrad2(ind) = (lterm2(i)*(u2(i)**((n-1)+(k-1)))) + sum2
     .                        + bterm2(i)
c               print *,'n,k,ind,covgradt(ind)',n,k,ind,covgradt(ind)
c               print *,'n,k,ind,covgrad2(ind)',n,k,ind,covgrad2(ind)
            enddo
          enddo

c        WRITE(6,210)
c  210   FORMAT(/,'Covariance matrix for loose gradients')
        ICT=0
        DO 240 II=1,ngrad
         DO 220 JJ=1,II
           ICT=ICT+1
           TEMP(JJ)=COVgrad2(ICT)
  220    CONTINUE
c         WRITE(6,230) (TEMP(IJ),IJ=1,II)
c  230    FORMAT(1X,9D12.5)
  240  CONTINUE

C       Invert to get weight
         CALL INVER2(COVgradT,BE,1,ngrad,rcond,IER)
         if( ier.ne.0 )call report_stat('WARNING','SOLVE','lwgrad',' '
     .     , 'Bad inverse of tight-soln gradient covariance matrix',0)
         CALL INVER2(COVgrad2,BE,1,ngrad,rcond,IER)
         if( ier.ne.0 )call report_stat('WARNING','SOLVE','lwgrad',' '
     .     , 'Bad inverse of loose-soln gradient covariance matrix',0)
c         print *,'LWGRAD ind ',ind      
c         print *,'1 covgradt ',(covgradt(ii),ii=1,6)   
c         print *,'1 covgrad2 ',(covgrad2(ii),ii=1,6)
c         WRITE(6,250)
c  250    FORMAT(/,'Weight matrix for loose gradient delays')
c         ICT=0
c         DO 280 II=1,ngrad
c           DO 260 JJ=1,II
c             ICT=ICT+1
c           **uncomment declaration if used
c            TEMP1(JJ)=COVgrad2(ICT)
c  260      CONTINUE
c           WRITE(6,270) (TEMP1(IJ),IJ=1,II)
c  270      FORMAT(1X,9D12.5)
c  280    CONTINUE

c       WRITE(6,282)
c  282  FORMAT(/,'Weight matrix for tight gradients')
c       ICT=0
c       DO 288 II=1,ngrad
c         DO 284 JJ=1,II
c           ICT=ICT+1
c           **uncomment declaration if used
c          TEMP2(JJ)=COVgradT(ICT)
c 284    CONTINUE
c         WRITE(6,286) (TEMP2(IJ),IJ=1,II)
c  286    FORMAT(1X,9D12.5)
c  288  CONTINUE
C
c      (remember that weight = 1.0/(constraint**2) )
c      weight increment = (new weight) - (old weight)
c      new weight = old weight + weight increment
c
       do  k=1,ind
         covgrad(k)=covgrad2(k)-covgradt(k)
       enddo
              
c       print *,' ind  ',ind
c       print *,'2 covgradt',(covgradt(ii),ii=1,6)   
c       print *,'2 covgrad2 ',(covgrad2(ii),ii=1,6)
c       print *,'2 covgrad ',(covgrad(ii),ii=1,6)   


c       WRITE(6,292)
c  292  FORMAT(/,'Weight matrix increment for loose gradient delays')
c       ICT=0
c       DO 298 II=1,ngrad
c         DO 294 JJ=1,II
c           ICT=ICT+1
c           TEMP3(JJ)=COVgrad(ICT)
c  294    CONTINUE
c         WRITE(6,296) (TEMP3(IJ),IJ=1,II)
c  296    FORMAT(1X,9D12.5)
c  298  CONTINUE

c  Add weight increment to normal matrix
                   
c      print *,'LWGRAD 1 indx ngrad maxc ',indx,ngrad,maxc 
c      print *,'  ALC ',(alc(ii),ii=1,5)                       

      call addwgt( indx,maxc,ngrad,covgrad )

c      print *,'LWGRAD 2 indx ngrad maxc ',indx,ngrad,maxc
c      print *,'  ALC ',(alc(ii),ii=1,5)
 
c      do 300 j=1,ngrad
c      JJ1=(J*J-J)/2
c      JJ2=J+INDX
c      JJ3=(JJ2*JJ2-JJ2)/2
c      DO 300 K=1,J
c      IND1=K+JJ1
c      IND2=K+INDX+JJ3
c      if (l2flag.eq.1.or.l2flag.eq.3.or.l2flag.eq.4 ) then 
c         print *,'data  ind2 alc ',ind2,alc(ind2)
c         ALC(IND2)=ALC(IND2)+COVgrad(IND1)  
c      else
c         A(IND2)=A(IND2)+COVgrad(IND1)
c      endif
c      print*,'LWGRAD: indx ind1 ind2 alc ',indx,ind1,ind2,alc(ind2)
c  300 CONTINUE
            
         
c-----end if on constant of piecewise linear model
      endif

      icnt=icnt+ngrad
                                       
c-----end loop on stations
500   continue

c-----end loop on N/S and E/W gradients
1000  continue

      return
      end

