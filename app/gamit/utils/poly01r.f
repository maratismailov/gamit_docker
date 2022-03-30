***** poly01 <weigh=yes/no>  <order=0/1>   <fake=yes/no>  <res=yes/no> <n_sigma>
* MOD TAH 990303: Added additional argument of sigma editing.
*                 Output also changed to add sigma and an edit flag (E) at
*                 the ends of the edited data.
* MOD TAH 010705: Changes routine to use rcpar rather than getarg because
*                 of different calling numbers in getarg on HP. 
* MOD TAH 020304: Added removing a refernce x value so that year plots over
*                 short time scales will correctly yield RMS.  Also cleaned up
*                 bounds for fake_x values.
* MOD TAH 020418: Re'Mod'd the boundaries for the lines.
* MOD TAH 040714: Increased significant digits in header
*****
*****

****** comparison with svdlin. estimates are ok.
***** Uncertainties are not. Murray estimates uncertainties as follows:
***** The uncertainties are obtained from the squareroot of (multiplying the covariance-relevant element) with the chi**2)
***** dof a : N_obs - 1
***** dof b : N_obs - 2


* calculates (weighted)mean and fits  a 1st order polynomial(y=a+bx)
** input x(km) y(mm) s(mm)


* n_max : maximum data points
      integer*4 n_max
      parameter (n_max = 9999)

* x data, y data, s sigmas, covariance-weights
      real*8 x(n_max),y(n_max),s(n_max),w(n_max)

* MOD TAH 990303: New variables
* n_sigma - Sigma limit to edit data
* sres    -- Scaled residual
* max_res -- Maximum scaled residual
* avg     -- Avergae with the slope set to zero
      real*8 n_sigma, sres, max_res, avg
      
* run_ns  - Character string with sigma limit
      character*16 run_ns
      
* edflag  - Flag set to E for edited points, blank for good
      character*1  edflag
      
* ierr    - IOSTAT error reading n_sigma from runstring 
* i_with_max -- Largest outlier point
* n_act      -- Actual number of remaining observations
      integer*4 ierr, i_with_max, n_act
            
* still_edits -- Logical to indicate that we are still editing
      logical still_edits      

* loops
      integer*4 i
* number of observations
      integer*4 n_obs
*     account for weigths
      character*10 weigh
* intermediate parameters  sumyw,sumw,sums,sumxw,sumxyw,sumxxw,det
      real*8 sumyw,sumw,sums,sumxw,sumxyw,sumxxw,det
* arbitrary x-values for nicer plots
      real*8 fake_x
* reference value of x to avoid rounding error with year scales
*     for small amounts of data,  Added min_x and max_x values
      real*8 ref_x, min_x, max_x

* (weighted) mean
      real*8 wmean,unc_mean
* a  + b X     constant :a     slope : b
      real*8 constant,slope,unc_constant,unc_slope
* order  0: 0th order polynomial
* order  1: 1th order polynomial
      character*1 order
*     chisquare
      real*8 chisqr
* normalized rms, weighted rms
      real*8 nrms, wrms
* write fake start and stop x
      character*10 fake, res
* number of command line arguments
      integer*4 n_arg, iargc 
* zero 
      real*8 zero
***** poly01 <weigh=yes/no>  <order=0/1>   <fake=yes/no>  <res=yes/no>
*****
****


      call report_stat('CLEAR','POLY01R',' ',' ', ' ',0)
      call report_stat('CLEAR','LIB',' ',' ', ' ',0)
 

      n_arg = iargc()
      weigh = "no"
      order = "0"    
      fake  = "no"
      res   = "no"
      zero  = 0.d0
      if ( n_arg .ge. 1 )  call  rcpar(1,weigh)
      if ( n_arg .ge. 2 )  call  rcpar(2,order)
      if ( n_arg .ge. 3 )  call  rcpar(3,fake)
      if ( n_arg .ge. 4 )  call  rcpar(4,res)
*     MOD TAH 990303: see if new  argument
      if ( n_arg .ge. 5 )  then
          call  rcpar(5,run_ns)
          read(run_ns,*,iostat=ierr) n_sigma
          if( ierr.ne.0 ) then 
             write(*,80) ierr, run_ns
 80          format('IOSTAT Error ',i4,' occurred decoding ',
     .              a)
             n_sigma = 0.d0
          endif
      else
          n_sigma = 0.d0
      endif        
c      write(*,*)  n_arg,weigh,order,fake,res,n_sigma
c     initialize these to avoid compiler warning of undefined variable at 600
      wrms = 0.d0
      nrms = 0.d0
      unc_slope = 0.d0
      unc_constant = 0.d0
      slope = 0.d0
      unc_mean = 0.d0
      wmean = 0.d0
      min_x = 1.d6 
      max_x = -1.d6
      
* MOD RWK 010103: Put n_obs inside loop to avoid use of loop variable i
      n_obs = 0
* MOD TAH 020418: Set ref_x to 0 value, will be updated below when the data
*     is read.  (Ensure no strange values if data not read correctly).
      ref_x = 0
      do i=1,n_max
        w(i) = 1.d0
        s(i) = 1.d0
        if(weigh(1:3).eq."yes") then
          read(*,*, end = 900 ) x(i),y(i),s(i)

          if(s(i).ne.0.d0) w(i)=1.d0/(s(i)*s(i))
          n_obs = n_obs + 1
        else
          read(*,*, end = 900 ) x(i),y(i)
          n_obs = n_obs + 1
        endif
* MOD TAH 020304: Remove reference value from x.  Also keep track of
*       min and max x values
        if( i.eq.1 ) ref_x = x(1)
        x(i) = x(i) - ref_x
        if( x(i).lt.min_x ) min_x = x(i)
        if( x(i).gt.max_x ) max_x = x(i)

      enddo
c*  900  n_obs = i -1
  900 continue

c      write(*,'(i5,3f10.3)') (i, x(i),y(i),w(i),i=1,n_obs)

* MOD TAH 990303: Loop added here to edit data:

      still_edits = .true.
      n_act = n_obs
      do while ( still_edits )
        
******   WMean **********
*        SUM y*w
*        SUM   w
         sumyw = 0.d0
         sumw  = 0.d0
         sums  = 0.d0
         n_act = 0
         do i=1,n_obs

*          Check to see if sigma is negative indicating an edited
*          data point         
           if( s(i).gt.0 ) then
              n_act = n_act + 1
              sumyw  = sumyw + y(i)*w(i)
              sumw   = sumw  +      w(i)
              sums   = sums  + s(i)*s(i)
           end if
         enddo
         wmean = sumyw / sumw
c        where did the original theory for unc mean come from?
c        current unc_mean from King/McClusky 950904.
c        unc_mean     = sqrt ( sums )
         unc_mean     = sqrt ( 1/sumw )

         if ( order .eq. "1" ) then
*******     y = A + B*x
            sumxw  = 0.d0
            sumxyw = 0.d0
            sumxxw = 0.d0
            n_act = 0
            do i=1,n_obs

*              Check for edited data
               if( s(i).gt.0 ) then  
                  n_act = n_act + 1          
                  sumxxw  = sumxxw + x(i)*x(i)*w(i)
                  sumxyw  = sumxyw + x(i)*y(i)*w(i)
                  sumxw   = sumxw  + x(i)     *w(i)
               endif
            enddo
            det = sumw* sumxxw - sumxw * sumxw
            constant = ( sumxxw * sumyw - sumxw * sumxyw ) / det
            slope    = (-sumxw  * sumyw + sumw  * sumxyw ) / det
            avg      = sumyw/sumw
            unc_constant = sqrt ( sumxxw / det )
            unc_slope    = sqrt ( sumw / det )
         endif


*        order 0 :
         if ( order .eq. "0" ) then
            n_act = 0 
            chisqr  = 0.0d0
            do i=1,n_obs
               if( s(i).gt.0 ) then 
                   n_act = n_act + 1
                  chisqr = chisqr + (( y(i) - wmean ) / s(i))**2
               endif
            enddo

C           if ( n_obs .gt. 1 ) then
C              nrms = sqrt ( chisqr / (n_obs - 1. ) )
C              wrms = sqrt ( ( n_obs /  (n_obs - 1. ) ) *  chisqr / sumw )
            if ( n_act .gt. 1 ) then
               nrms = sqrt (chisqr/(n_act - 1))
               wrms = sqrt(float(n_act)/(float(n_act)-1.)*chisqr/sumw)   
            else
               nrms = 999.99d0
               wrms = 999.99d0
            endif

         elseif ( order .eq. "1" ) then
*           order 1 :
           n_act = 0 
           chisqr = 0.d0
           do i=1,n_obs
              if( s(i).gt.0 ) then 
                 n_act = n_act + 1
                 chisqr = chisqr + ((y(i)-constant-slope*x(i))/s(i))**2
              endif
            enddo
            if ( n_act .gt. 2 ) then
               nrms = sqrt (chisqr/(n_act-2))     
               wrms = sqrt(float(n_act)/(float(n_act)-2.)*chisqr/sumw)   
            else
               nrms = 999.99d0
               wrms = 999.99d0
            endif             

         else    
            call report_stat('FATAL','POLY01R','poly01r',' '
     .         ,'Invalid order',0) 
            stop
         endif
         
* MOD TAH 990303: loop over data and edit the largest outlier if needed
         still_edits = .false.
         if( n_sigma.gt.0 ) then 
            max_res = 0.d0
            i_with_max = 0
            do i = 1, n_obs
               if( order .eq. "0" ) then 
                  sres = (y(i)-wmean)/(s(i)*nrms)
               else
                  sres = (y(i)-constant-slope*x(i))/(s(i)*nrms)
               endif
               if( (abs(sres).gt.n_sigma .and. 
     .             abs(sres).gt.max_res) .and.
     .             s(i).gt.0                   ) then 
                   still_edits = .true.
                   max_res = abs(sres)
                   i_with_max = i
               end if
            end do

*           See if we edited any data
            if( i_with_max.gt.0 ) then
                s(i_with_max) = -s(i_with_max)
            end if
         end if   
                  
*
* MOD TAH 990303: End loop over editing data         
      enddo 

*     Now output the results
      if (weigh(1:3).eq."yes") then
          write(*,600) wmean,unc_mean,
     >        constant,unc_constant,slope,unc_slope, n_act, nrms,wrms,
     >        n_obs
 600      format("  wMEAN:",2x,2f15.4,
     >           " constant:",2x,2f15.4,
     >           "    slope:",2x,2f20.10,
     >           "      obs:",2x,i5,
     >           "     nrms:",2x,f10.4,
     >           "     wrms:",2x,f10.4,
     >           "  all obs:",2x,i5)
      else
          write(*,601) wmean,unc_mean,
     >        constant,unc_constant,slope,unc_slope, n_act, nrms,wrms,
     >        n_obs
 601      format("  wMEAN:",2x,2f15.4,
     >           " constant:",2x,2f15.4,
     >           "    slope:",2x,2f20.10,
     >           "      obs:",2x,i5,
     >           "     nrms:",2x,f10.4,
     >           "     wrms:",2x,f10.4,
     >           "  all obs:",2x,i5)
      endif

      if ( order .eq. "1" ) then
***      this is for repeatability plots
C        if ( x(1) .le.  x(2) ) then
C            fake_x = x(1) / 2. 
C        else
C            fake_x = x(1) * 2. 
C        endif
         fake_x = min_x - (max_x-min_x)/5.d0
* MOD TAH 020418: Changed fake_x to -2000
         fake_x = -2000
         if ( res(1:3) .eq. "yes" ) then
             if(fake(1:3) .eq."yes") write(*,610)fake_x+ref_x,zero
         else
             if(fake(1:3).eq."yes") 
     .              write(*,610)fake_x+ref_x,constant+slope*fake_x
         endif
*        write increasing order
         do i=1,n_obs
            if( s(i).gt.0 ) then
               edflag = " "
            else
               edflag = "E"
            endif 
            if ( res(1:3) .eq. "yes" ) then
                write(*,620) x(i)+ref_x,zero,y(i)-(constant+slope*x(i)),
     .                       abs(s(i)), edflag                
            else
                write(*,620) x(i)+ref_x,constant+slope*x(i), 
     .                       y(i), abs(s(i)), edflag
            endif
         enddo
C        if ( x(n_obs) .gt.  x(n_obs-1) ) then
C           fake_x = x(1) * 2. 
C        else
C           fake_x = x(1) / 2.
C        endif
         fake_x = max_x + (max_x-min_x)/5.d0
* MID TAH 020418: Changes fake_x to + 3000
         fake_x = 3000
         if ( res(1:3) .eq. "yes" ) then
            if(fake(1:3) .eq."yes") write(*,610)fake_x+ref_x,zero
         else
            if(fake(1:3).eq."yes") 
     .             write(*,610)fake_x+ref_x,constant+slope*fake_x
         endif

      else  if ( order .eq. "0" ) then
***      this is for repeatability plots
         if( fake(1:3) .eq."yes")  then
             fake_x = min_x - (max_x-min_x)/5.d0
             fake_x = -2000
             if ( res(1:3).eq."yes") then
                write(*,610)fake_x+ref_x,zero
             else
                write(*,610)fake_x+ref_x,wmean
             endif
*            write increasing order
             do i=1,n_obs
                if( s(i).gt.0 ) then
                   edflag = " "
                else
                   edflag = "E"
                endif 
                if ( res(1:3) .eq. "yes" ) then
                   write(*,620) x(i)+ref_x,wmean,y(i)-wmean, 
     .                                      abs(s(i)),edflag
                else
                   write(*,620) x(i)+ref_x,wmean,y(i)-wmean, 
     .                                      abs(s(i)),edflag
                endif   
             enddo
             fake_x = max_x + (max_x-min_x)/5.d0
             fake_x = 3000
             if ( res(1:3) .eq. "yes" ) then
                write(*,610)fake_x+ref_x,zero 
             else
                write(*,610)fake_x+ref_x,wmean 
             endif
         endif
      endif
 610  format(f20.4,5x,f20.4,1x,a1)
 620  format(f20.4,5x,f20.4,5x,f20.4,1x,f20.4,1x,a1)

      stop
      end
