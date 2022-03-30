      Subroutine EPHTRP( pjd,ibody, ivel, rvec )

c     Perform an Everett interpolation to get coordinates at the requested epoch

c     R. King from older GAMIT (and originally PEP) routines solred and 
c     lunred.  See Ash(1972) for the formulas.   10 January 2018 

c     pjd     r*8  PEP JD in coordinate time

c     ibody   i*4  Body to interpolate
c                  3 Earth
c                  2 Venus
c                  5 Jupiter
c                 10 Moon

c     ivel    i*4  0 position only
c                  1 position and velocity

c     rvec(6  r*8 position (i=1,3) and velocity i=4,6) for requested body
  

      implicit none
                            
      real*8 pjd,mpjd,rvec(6)
      integer*4 ibody,ivel                                                 

c  Coordinate arrays and tabular intervals 
      common/ephtab/earth(6,15),moon(6,120),venus(6,15),jupiter(6,15)
     .  ,tsun(15),tmoon(120),jdstart,intsun,intmoon,intvenus,intjupiter
      integer*4 jdstart,intsun,intmoon,intvenus,intjupiter
      real*8 earth,moon,venus,jupiter,tsun,tmoon

c  Coordinate array for interpolation (10 points centered around the current time
      real*8 yy(10,3)

c  Local time-spacing and relative times
      real*8 dint,tymid,p
                   
c  Function and variables for report_stat call
      integer*4 rcpar,len
      character*80 prog_name,message

c  Local indices
      integer*4 kmid,i,j,k

      logical debug/.false./

c  Get the calling module name for report_stat
      len = rcpar(0,prog_name)

      
c  Fill the interpolation arrays and compute the indices for the requested body
                      
      if( ibody.eq.3 ) then  
        dint = dfloat(intsun)  
        mpjd = pjd - 2400000.d0 
c       find and index the middle record
        do i=1,15             
          if(debug)
     .   print *,'dint mpjd i tsun tsun1 ',dint,mpjd,i,tsun(i),tsun(i+1)
          if( tsun(i).le.mpjd.and.(tsun(i)-mpjd).le.dint) then
            kmid = i 
            tymid = tsun(i)   
          endif
        enddo             
c       set the first record and populate the y-vectors    
        k = kmid - 5
        do i=1,10
          do j=1,3
            yy(i,j) = earth(j,k+i)
          enddo     
        enddo    
        if(debug) then
          print *,'EPHTRP Earth kmid k  ',kmid,k 
          write(*,'(a)') 'yy-vectors 1-3 :'
          do i=1,10 
            write(*,'(3f16.2)') (yy(i,j),j=1,3)
          enddo 
        endif
      elseif( ibody.eq.2 ) then
        dint = dfloat(intvenus)  
        mpjd = pjd - 2400000.d0 
c       find and index the middle record
        do i=1,15             
         if(debug) 
     .   print *,'dint mpjd i tsun tsun1 ',dint,mpjd,i,tsun(i),tsun(i+1)
         if( tsun(i).le.mpjd.and.(tsun(i)-mpjd).le.dint) then
          kmid = i 
          tymid = tsun(i)   
          endif
        enddo  
        if(debug) print *,'Sun kmid tymid ',kmid,tymid 
                    
c       set the first record and populate the y-vectors    
        k = kmid - 5
        do i=1,10
          do j=1,3
            yy(i,j) = venus(j,k+i)
          enddo     
        enddo    
      elseif( ibody.eq.5 ) then
        dint = dfloat(intjupiter)  
        mpjd = pjd - 2400000.d0 
c       find and index the middle record
        do i=1,15             
          if(debug)
     .   print *,'dint mpjd i tsun tsun1 ',dint,mpjd,i,tsun(i),tsun(i+1)
          if( tsun(i).le.mpjd.and.(tsun(i+1)-mpjd).le.dint) then
            kmid = i 
            tymid = tsun(i)   
          endif
        enddo             
c       set the first record and populate the y-vectors    
        k = kmid - 5
        do i=1,10
          do j=1,3
            yy(i,j) = jupiter(j,k+i)
          enddo     
        enddo                   
      elseif( ibody.eq.10 ) then
        dint = 2.d0**intmoon
        mpjd = pjd - 2400000.d0 
c       find and index the middle record
        do i=1,120             
          if(debug) print *,'dint mpjd i tmoon tmoon1 '
     .          ,dint,mpjd,i,tmoon(i),tmoon(i+1)
          if( tmoon(i).le.mpjd.and.(tmoon(i)-mpjd).le.dint) then
            kmid = i 
            tymid = tmoon(i)   
          endif
        enddo    
       if(debug) print *,'Moon kmid tymid ',kmid,tymid 
c       set the first record and populate the y-vectors    
        k = kmid - 5
        do i=1,10
          do j=1,3
            yy(i,j) = moon(j,k+i)
          enddo     
        enddo             
      else   
        write(message,'(a,i3,a)') 'ibody=',ibody,' not coded'
        call report_stat('FATAL',prog_name,'lib/ephtrp',' ', message,0)
      endif

c  Reset kmid for the y-vectors 
      kmid = 5 
      if( debug ) then
        do i=1,10
          write(*,'(a,i3,3f16.2)') 'yy ',i,(yy(i,j),j=1,3)
        enddo     
      endif 

c  Fractional times within the tabular interval 
      p = (mpjd-tymid)/dint                                                                                      
      if(debug) write(*,'(a,2f12.4,f8.6)') 'mpjd tymid p',mpjd,tymid,p
      if( debug ) then
        print *,'EPHTRP pjd kmid dint tymid p ',pjd,kmid,dint,tymid,p
        print *,'y-vectors for ibody ',ibody
        do i=1,10
          write(*,'(3f16.1)') (yy(i,j),j=1,3)
        enddo
      endif  

c  Interpolate the requested value 
      call evtrp(yy,kmid,p,dint,ivel,rvec)

      if( debug ) write(*,'(a,i2,f12.3,6f16.1)')   
     .     'ivel pjd rvec ',ivel,pjd,(rvec(i),i=1,6)
      return
      end

c---------------------------------------------------------------------------
                   
      Subroutine EVTRP( yy,kmid,p,dint,ivel,rvec )

c     Interpolate coordinates given the tabular values, y-vectors
c     R, King from the PEP coding of Ash [1972]

      implicit none
                     
c     Input                
c       yy(10,3)  position values to be interpolated  
c       kmid      index of the yy vector nearest below the interpolation epoch (usually 5)  
c       p         fractional time of the requested time within the interval
c       dint      tabular interval (4. for Earth, Venus, Jupiter), 0.4 for Moon)
c       ivel      0 if position only; 1 if position and velocity
      integer*4 ivel,kmid
      real*8 yy(10,3),p,dint,rvec(6) 

c      Output
c        rvec(6)  position and velocity 
            
c  Everett interpolation coefficients (constant) 
      real*8 evcf(5,5),fact(9),dfact(4)
      DATA evcf/
     .    1.7873015873015873D0,
     .   -0.4960317460317460D0, 0.1206349206349206D0,
     .   -0.1984126984126984D-01, 0.1587301587301587D-02,
     .   -0.9359567901234568D0,
     .    0.6057098765432098D0,-0.1632716049382716D0,
     .    0.2779982363315696D-01,-0.2259700176366843D-02,
     .    0.1582175925925926D0,
     .   -0.1171296296296296D0, 0.4606481481481481D-01,
     .   -0.8796296296296296D-02, 0.7523148148148148D-03,
     .   -0.9755291005291005D-02,
     .    0.7605820105820106D-02,-0.3505291005291005D-02,
     .    0.8597883597883598D-03,-0.8267195767195767D-04,
     .    0.1929012345679012D-03,
     .   -0.1543209876543210D-03, 0.7716049382716048D-04,
     .   -0.2204585537918871D-04, 0.2755731922398589D-05/
      DATA fact
     .   /1.0D0,3.0D0,5.0D0,7.0D0,9.0D0,2.0D0,4.0D0,6.0D0,18.0D0/
     . ,   dfact/6.d0,20.d0,42.d0,72.d0/


       integer*4 j10,j1l,ijk0,nr,n,i,j,k,kk
       real*8 ytrp(5,2,3),p2,q,q2,f(4),s1,s2
       logical debug/.false./
                                        
      if(debug) then 
        write(*,'(a,i3,2f10.5,i3)') 'EVTRP kmid p dint  '
     .                         ,  kmid,p,dint
        write(*,'(a)') 'yy-vectors 1-3 :'
        do i=1,10 
          write(*,'(3f16.2)') (yy(i,j),j=1,3)
        enddo 
      endif 
                  
c  Compute the interpolation y-vectors 
      do k=1,3  
        do j=1,2   
          if(j.eq.1) nr = kmid 
          if(j.eq.2) nr = kmid + 1     
          if(debug.and.k.eq.1) write(*,'(a,3i3)') 'kmid j nr ',kmid,j,nr
          do n=1,4
            f(n) = yy(nr+n,k) + yy(nr-n,k)      
            if(debug.and.k.eq.1) 
     .        write(*,'(a,2i3,3f16.2)') 
     .            'For k=1 nr n yy-dwn yy-up fn '
     .              ,nr,n,yy(nr-n,k),yy(nr+n,k),f(n)
          enddo                  
          do i = 1,5
            ytrp(i,j,k) = evcf(1,i)*yy(nr,k) +
     .                 (evcf(2,i)*f(1) +
     .                 (evcf(3,i)*f(2) +
     .                 (evcf(4,i)*f(3) +
     .                 (evcf(5,i)*f(4)))))
          if(debug.and.k.eq.1 ) 
     .      write(*,'(a,i3,2f16.2)') 
     .       ' i evcf yy-0 ',i,evcf(1,i),yy(nr,k)
          enddo
        enddo
      enddo
      if(debug) then
        do i=1,5
         write(*,'(a,i2,3f16.2)') 'j=1 i ytrp ',i,(ytrp(i,1,k),k=1,3)
        enddo    
        do i=1,5
         write(*,'(a,i2,3f16.2)') 'j=2 i ytrp ',i,(ytrp(i,2,k),k=1,3)
        enddo
      endif
         
c  Get the auxiliary time quantities

      q = 1.d0 -p
      p2 = p**2
      q2 = q**2
      if(debug) write(*,'(a,4f12.4)') 'p q p2 q2 ',p,q,p2,q2

c  Do the interpolation
                           
      do i=1,3
        rvec(i) = q*(ytrp(1,1,i) 
     .      + q2*(ytrp(2,1,i)+q2*(ytrp(3,1,i)+q2*(ytrp(4,1,i)
     .                                               +q2*ytrp(5,1,i)))))
     .      +     p*(ytrp(1,2,i) 
     .      + p2*(ytrp(2,2,i)+p2*(ytrp(3,2,i)+p2*(ytrp(4,2,i)
     .                                               +p2*ytrp(5,2,i)))))
      enddo   
      if(debug) write(*,'(a,3f16.6)') 'EVTRP pos  ',(rvec(i),i=1,3) 

                                                            
c  Interpolate positions to get velocities if required

      if(ivel.eq.1 ) then 
        do i=1,3 
          s1 = 0.d0
          s2 = 0.d0
          do k=5,2,-1
            s1 = s1*p2 + ytrp(k,2,i)*fact(k)    
            s2 = s2*q2 + ytrp(k,1,i)*fact(k)
          enddo                                   
          rvec(i+3) = (s1*p2+ytrp(1,2,i)-s2*q2-ytrp(1,1,i))/dint
        enddo  
        if(debug) write(*,'(a,3f16.10)') 'EVTRP rvec vel '
     .      ,(rvec(i),i=4,6)   
      endif
            
      if(debug) then
        print *,'STOP in EVTRP' 
        stop 1
      endif

      return
      end 



