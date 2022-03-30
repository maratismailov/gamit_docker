      Subroutine AVGLOAD( icall, istat,jsat,elev,iphi )  

c     Form the normal equations to estimate the effective atmospheric and
c     hydrological loading corrections to be passed to GLOBK on the h-file
c     The inversion is done in glohed.  R. King 050224

c     Called only if atmlmod or hydrolmod non-blank    

c       icall = 1 : increment the normal equations for this epoch/site/sat
c       icall = 2 : invert the normal equations for all sites

      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'    
      include 'models.h'
      include 'parameters.h'
                       
c  Calling arguments

      integer*4 istat,jsat,iphi(maxsit,maxsat)
      real*8 elev(maxsit,maxsat)       

c  Local
      logical first,debug
      integer*4 icall,i,j,k
      real*8 omc,eradkm,dN,dE,dU,dzen,ainc(5,5),binc(5)
     .     , adjload(5,maxsit),w,w1(3),sumw1(3,maxsit),latd
     .     , testload 
      data eradkm/6378.137d0/,first/.true./,debug/.false./
      save first,w,w1,sumw1,adjload
       
      if( icall.eq. 2 ) go to 100

c  Initialize all the arrays on first entry
        
      if( first ) then 
        do k=1,nsite  
          do i=1,3
            sumatmload(i,k) = 0.d0
            sum1atmload(i,k) = 0.d0        
            sumhydload(i,k) = 0.d0
            sum1hydload(i,k) = 0.d0    
            sumw1(i,k) = 0.d0
          enddo
          count(k) = 0.d0  
          do i=1,5   
            batm(i,k) = 0.d0 
            bhyd(i,k) = 0.d0    
            adjload(i,k) = 0.d0
            do j=1,5
              amatl(i,j,k) = 0.d0
            enddo
          enddo     
          arithavg(k) = .false.
c         If the site is within 5 degree of the pole, use a weighted arithmetic mean
c         instead of a least squares solution to avoid the problem if an ill-define East
          latd = preval( 3*(k-1)+1) /convd           
          if(dabs(latd).gt.85.d0 ) arithavg(k) = .true.
          if( debug ) 
     .       print *,'k,latd arithavg ',istat,latd,arithavg(k)
          enddo 
        first = .false.     
      endif
           

c  See if this observations was included in the solution
                                   
      if( iphi(istat,jsat).eq.0 ) return

                                 
c  Get the partial derivatives wrt N, E, U, and zenith delay

c     C-file partials are cycles per radian or km: convert to cyc/mm
      dN = tpart(1,istat,jsat) / eradkm / 1.d6
      dE = tpart(2,istat,jsat) / eradkm / 1.d6
      dU = tpart(3,istat,jsat) / 1.d6
c     Zen delay partial is cycles per m: convert to cyc/mm
      dzen = tpart(4,istat,jsat) / 1.d3    
             
c  Get the weights

c     weights (inverse variance in mm**2; units cancel in normal equations)
      w = sit_err(istat)**2 + sat_err(jsat)**2
      if( err_mod(istat).eq.'elevation') 
     .     w = w + sit_elv(istat)**2/(dsin(elev(istat,jsat)))**2
      w = 1.d0/w*1.d-6                      


c  Form the increments to the coefficient matrix, same for atml and hydro 
          
      if( .not.arithavg(istat).or.debug ) then 
        ainc(1,1) = dN**2
        ainc(1,2) = dN*dE
        ainc(1,3) = dN*dU
        ainc(1,4) = dN*dzen  
        ainc(1,5) = dN
        ainc(2,1) = dE*dN
        ainc(2,2) = dE**2
        ainc(2,3) = dE*dU
        ainc(2,4) = dE*dzen  
        ainc(2,5) = dE
        ainc(3,1) = dU*dN
        ainc(3,2) = dU*dE         
        ainc(3,3) = dU**2
        ainc(3,4) = dU*dzen   
        ainc(3,5) = dU
        ainc(4,1) = dzen*dN
        ainc(4,2) = dzen*dE
        ainc(4,3) = dzen*dU
        ainc(4,4) = dzen**2  
        ainc(4,5) = dzen
        ainc(5,1) = dN
        ainc(5,2) = dE
        ainc(5,3) = dU
        ainc(5,4) = dzen
        ainc(5,5) = 1.d0   
c       increment the left- and right-hand sides of the normal equations
        do i=1,5
          do j=1,5
            amatl(i,j,istat) = amatl(i,j,istat) + w * ainc(i,j) 
          enddo       
        enddo   
      endif              
              
c     get the weights for the arithmetic average method
      if( arithavg(istat).or.debug ) then
c       weight by the observation error only
        sumweight(istat) = sumweight(istat) + w        
c       weight by the partial  
        w1(1) = w*dN**2
        w1(2) = w*dE**2
        w1(3) = w*dU**2                           
        do i=1,3
          sumw1(i,istat) = sumw1(i,istat) + w1(i)
        enddo
        count(istat) = count(istat) + 1.d0
      endif

c  See if atmospheric loading was applied

      if( atmlmod.ne.'        ' ) then                                 
        if( .not.arithavg(istat).or.debug)  then
c         compute the change in phase from the loading applied
          omc =  dN*atmload(1,istat) + dE*atmload(2,istat) 
     .         + dU*atmload(3,istat)
cd        print *,'omc ',omc
          binc(1) = dN * omc
          binc(2) = dE * omc
          binc(3) = dU * omc
          binc(4) = dzen * omc
          binc(5) = omc 
cd        print *,'binc ',binc
c         increment the right-hand side 
          do i=1,5
            batm(i,istat) = batm(i,istat) + w * binc(i)
          enddo               
        endif
        if( arithavg(istat).or.debug ) then
          do i=1,3
            sumatmload(i,istat) = sumatmload(i,istat)+atmload(i,istat)*w
            sum1atmload(i,istat) = 
     .         sum1atmload(i,istat) + atmload(i,istat)*w1(i)
          enddo 
        endif
      endif
                               

c  See if hydrological loading was applied

      if( hydrolmod.ne.'        ' ) then  
        if( .not.arithavg(istat).or.debug)  then
c         compute the change in phase from the loading applied
          omc =  dN*hydload(1,istat) + dE*hydload(2,istat) 
     .         + dU*hydload(3,istat)
          binc(1) = dN * omc
          binc(2) = dE * omc
          binc(3) = dU * omc
          binc(4) = dzen * omc
          binc(5) = omc 
c         increment the right-hand side 
          do i=1,5
            bhyd(i,istat) = bhyd(i,istat) + w * binc(i)
          enddo  
        endif
        if( arithavg(istat).or.debug ) then
          do i=1,3
            sumhydload(i,istat) = sumhydload(i,istat)+hydload(i,istat)*w
            sum1hydload(i,istat) = 
     .         sum1hydload(i,istat) + hydload(i,istat)*w1(i)
          enddo 
        endif
      endif
                        
      return
                                        
c -----------------end of code for incrementing (icall = 1)
     
100   continue
c     Calculate the averages--by inverting or taking the weighted mean
c     (called once by NORMD so loop over all stations)
          
      do k=1,nsite    
        if( atmlmod.ne.'        ' ) then                                 
          if( .not.arithavg(k).or.debug ) then
            call solveatm(5,amatl(1,1,k),batm(1,k),adjload(1,k))
            if( debug ) then 
              print *,'Atm load estimate for site ',k
              print *,(adjload(i,k),i=1,5) 
            endif
            if( .not.arithavg(k) ) then
              do i=1,3
                atmlavg(k,i) = adjload(i,k)
              enddo   
            endif
            if( debug ) then
              write(*,'(a,5f12.5)') 'Sigma ',(dsqrt(amatl(i,i,k)),i=1,5)
              print *,'Covariance matrix '
              do i=1,5
                write(*,'(5d14.5)') (amatl(i,j,k),j=1,5)
              enddo
              print *,'Correlation matrix'
              do i=1,5  
                write(*,'(5d14.5)') 
     .           ( amatl(i,j,k)/dsqrt(amatl(i,i,k)*amatl(j,j,k)),j=1,5)
              enddo    
            endif
          endif               
          if( arithavg(k).or.debug ) then
            do i=1,3
              sumatmload(i,k) = sumatmload(i,k)/sumweight(k) 
              sum1atmload(i,k) = sum1atmload(i,k)/sumw1(i,k)
c             use the obs-weighted value since partials may be bad close to the pole
              if( arithavg(k) ) atmlavg(k,i) = sumatmload(i,k)
            enddo 
            if( debug ) then 
              write(*,'(a,f10.0,3f8.2)') 
     .           'Atmload Obs-count Obs-wghted mean  '
     .           , count(k),(sumatmload(i,k),i=1,3)     
              write(*,'(a,10x,3f8.2)')   
     .          '               Obs-part-wgtd mean  '
     .         , (sum1atmload(i,k),i=1,3)
            endif
          endif
        else        
          do i=1,3
            atmlavg(k,i) = 0.d0
          enddo
        endif   

        if( hydrolmod.ne.'        ' ) then  
          if( .not.arithavg(k).or.debug ) then
            call solveatm(5,amatl(1,1,k),bhyd(1,k),adjload(1,k))
            if( debug ) then 
              print *,'Hydro load estimate for site ',k
              print *,(adjload(i,k),i=1,5) 
            endif
            do i=1,3
              hydrolavg(k,i) = adjload(i,k)
            enddo   
            if( debug ) then
              write(*,'(a,5f12.5)') 'Sigma ',(dsqrt(amatl(i,i,k)),i=1,5)  
              print *,'Covariance matrix '
              do i=1,5
                write(*,'(5d14.5)') (amatl(i,j,k),j=1,5)
              enddo
              print *,'Correlation matrix'
              do i=1,5  
                write(*,'(5d14.5)') 
     .           ( amatl(i,j,k)/dsqrt(amatl(i,i,k)*amatl(j,j,k)),j=1,5)
              enddo    
            endif
          endif
          if( arithavg(k).or.debug ) then
            do i=1,3
              sumhydload(i,k) = sumhydload(i,k)/sumweight(k) 
              sum1hydload(i,k) = sum1hydload(i,k)/sumw1(i,k)
c             use the obs-weighted average since partials may be bad close to pole
              hydrolavg(k,i) = sumhydload(i,k)
            enddo 
            if( debug ) then 
              write(*,'(a,f10.0,3f8.2)') 
     .           'Hydro load Obs-count Obs-wghted mean  '
     .           , count(k),(sumhydload(i,k),i=1,3)     
              write(*,'(a,10x,3f8.2)')   
     .          '               Obs-part-wgtd mean  '
     .            , (sum1hydload(i,k),i=1,3)
            endif
          endif
        else        
          do i=1,3
            hydrolavg(k,i) = 0.d0
          enddo
        endif 
                                      
c       check for reasonableness 
        testload=dsqrt(atmlavg(k,1)**2+atmlavg(k,2)**2+atmlavg(k,3)**2)
        if( testload.gt.100.d0 ) call report_stat('FATAL','SOLVE'
     .    ,'avgload',' ','Atmospheric load > 10 cm, something wrong',0)
        testload = 
     .      dsqrt(hydrolavg(k,1)**2+hydrolavg(k,2)**2+hydrolavg(k,3)**2)
        if( testload.gt.100.d0 ) call report_stat('FATAL','SOLVE'
     .    ,'avgload',' ','Hydrologic lLoad > 10 cm, something wrong',0) 

      enddo
 
      return
      end


c--------------------------------------------------------------

      Subroutine SOLVEATM(nparam,amat,bvec,adjust )

c     Solve the normal equations   

      implicit none
     
      integer*4  mxoprm,job,info,nparam,i,j 
                              
      real*8 amat(5,5),bvec(5),adjust(5)

      real*8 scale(5),rcond0,z(5),dterm2(2)

      character*256 message      

      parameter(mxoprm=5)
                

c     Invert the normal equations using linpack (blas 1) routines
c     assumng the that the matrix is positive definite (and thus
c     symmetric).  
           
c*    call report_stat('STATUS','SOLVE','avgload',' '
c*   .                 ,'Solving the normal equations',0)
 
c     first scale the matrix
      do i=1,nparam
        scale(i) = 1.d0/dsqrt(amat(i,i))
      enddo
      do i=1,nparam
        do j=1,nparam
          amat(i,j) = scale(i)*scale(j)*amat(i,j)
        enddo 
      enddo
 
 
c     factor the matrix
      call dpoco( amat,mxoprm,nparam,rcond0,z,info )
      
      if ( info.eq.0 .and.rcond0.ge.1.d-16 ) then      

c       do the inversion (job = 1 means no determinant)
        job = 1                                   
c       dpodi takes a full matrix but returns the inverse in upper triangle
        call dpodi( amat,mxoprm,nparam,dterm2,job)
c       fill in the lower half and rescale 
        do i=1,nparam
          do j = 1, i
            amat(j,i) = scale(i)*scale(j)*amat(j,i)
c*              amat(j,i) = amat(j,i) / (scale(i)*scale(j))
            amat(i,j) = amat(j,i)
          enddo
        enddo

c       multiply inv(amat) * bvec to get the adjustment vector 
        call dgemv ('N', nparam, nparam, 1.d0, amat, mxoprm, bvec
     .             , 1, 0.d0, adjust, 1)

                                     
      else
                         
        write(message,'(a,d7.1)') 
     .      'Normal matrix ill-conditioned, rcond0 = ',rcond0
        call report_stat('STATUS','SOLVE','avgload',' '
     .                 ,message,0)
      endif
 
        
c      call report_stat('STATUS','SOLVE','avgload',' '
c     .                 ,'Normal equations solved for average loading',0)
      return
      end


      


         


              


