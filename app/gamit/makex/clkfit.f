      Subroutine clkfit( nepc,jd,t,clock,use
     .                 , nseg,npseg,goodseg,jd0,t0,clkoff,rate
     .                 , rms,iprn,debug,iprndb )

c     Fit a 2nd order polynomial to one or more segments of SV clock values,
c     taken either from an SV3 file or a navigation file

c     R. King 28 July 2015; based on K. Feigl routine gamit/lib/clkera.f.
c     Unlike clkera, designed for erratic ground receivers using crystal
c     clocks. For SV clocks, I have assumed that if the first point (or
c     any point) is bogus, it will be marked as such in the SP3 file, 
c     and that any two points that generate a frequency greater than 
c     1.d-9 is a jump.
                                                   
      implicit none
             
      include '../includes/makex.h'
c      Dimensions in makex.h:
c        maxepc - input clock array values (dimension shared with observaton epochs)
c        maxseg - segments in clock fit 

c     Input:
c       nepc          - number of epochs in the input array
c       jd(maxepc)    - PEP Julian day of the input epochs
c       t(maxepc)     - seconds-of-day of the input epochs
c       clock(maxepc) - clock values (sec) 
c        use(maxepc)  - T if useable, F if not usualbe for fit

      integer*4 nepc,jd(maxepc)
      real*8 t(maxepc),clock(maxepc)
      logical use(maxepc)

c     Output:
c       nseg           - number of segments (between jumps) 
c       npseg(maxseg)  - number of points in each segment 
c       goodseg(maxseg)- T or F for whether a fit was obtained
c       jd0(maxseg)    - JD of first epoch for each segment
c       t0(maxseg)     - seconds-of-day of first epoch for each segment
c       clkoff(maxseg) - clock offset at Jd0/t0 for each segment
c       rate(maxseg)   - clock rate for each segment
           
      integer*4 nseg,jd0(maxseg)
      real*8 t0(maxseg),clkoff(maxseg),rate(maxseg)
      logical goodseg(maxseg)
                                                               
c     partial derivatives
      real*8 pp(2),pp1(2)
                                
c     number of points and starting and stopping points for each segment
      integer*4 npseg(maxseg),startseg(maxseg),stopseg(maxseg)

c     needed for poly01c
      integer*4 n_obs
      real*8 x(maxepc),y(maxepc),s(maxepc)
c     needed for residuals
      real*8 resid,sum2,rms(maxseg)
                                                    
c     other local  
      integer*4 i,j,k
      real*8 df,dt,oldtmt,newtmt,oldclk,newclk,clockref
      logical new_segment,end_segment
             
c     satellite to debug
      integer*4 iprn,iprndb
      logical  debug

c     Outer loop over segments (usually only 1)
                           
      new_segment = .true. 
      end_segment = .false.
      nseg = 1                                 

c     loop over all the data
      do i=1,nepc 
            
        if(new_segment.and.use(i)) then
          npseg(nseg) = 1      
          startseg(nseg) = i
          jd0(nseg) = jd(i)
          t0(nseg) = t(i)     
          if(debug.and.iprn.eq.iprndb) print *
     .     ,'new_segment nepc i nseg startseg npseg jd0 t0 '
     .     , nepc,i,nseg,startseg(nseg),npseg(nseg),jd0(nseg),t0(nseg)
c         initialize quantities for allen variance
          oldtmt = 0.d0
          oldclk = clock(i)
          new_segment = .false.

        else 
c         old segment: test for jump or end of list
          if( use(i) ) then             
            npseg(nseg) = npseg(nseg) + 1 
            if( debug.and.iprn.eq.iprndb ) print *
     .       ,'old segment nepc i nseg startseg npseg jd0 t0 t'
     .       , nepc,i,nseg,startseg(nseg),npseg(nseg),jd0(nseg),t0(nseg)
     .       , t(i)
            newtmt = (jd(i)-jd0(nseg))*86400.d0 + (t(i) - t0(nseg))
            newclk = clock(i)  
            df = (newclk-oldclk)/(newtmt-oldtmt) 
            if( dabs(df).gt.1.d-9 ) then
c             declare a jump detected if the frequency implied  is greater than 1 part in 10**9 
              new_segment = .true. 
              end_segment = .true.
              if(debug.and.iprn.eq.iprndb) 
     .          print *,'jump newclk oldclk newtmt oldtmt df'
     .                 , newclk,oldclk,newtmt,oldtmt,df
                npseg(nseg) = npseg(nseg) - 1
                stopseg(nseg) = i - 1  
              endif
            endif  
          endif 
          if( i.eq.nepc ) then
            stopseg(nseg) = i
            end_segment = .true.
          endif
c         compute the offset and rate for this segment 
          if(end_segment ) then
            n_obs = 0
            do j=startseg(nseg),stopseg(nseg)
              if( use(j) ) then
                n_obs = n_obs + 1 
                x(n_obs)=jd(jd(j)-jd0(nseg))*86400.d0 + (t(j)-t0(nseg))
                y(n_obs) = clock(j)
                s(n_obs) = 1.d-6    
               endif
            enddo   
            if( debug.and.iprn.eq.iprndb ) 
     .         print *,'Filled the arrays, startseg stopseg n_obs'
     .           , startseg(nseg),stopseg(nseg),n_obs
            if( n_obs.ge.3 ) then
              if(debug.and.iprn.eq.iprndb ) then
                do k=1,n_obs
                   print *,'k t clk s ',k,x(k),y(k),s(k)
                enddo
              endif
              call poly01c( maxepc,npseg(nseg),x,y,s
     .                    , clkoff(nseg),rate(nseg)
     .                    , iprn,debug,iprndb )
              if(debug.and.iprn.eq.iprndb) 
     .             print *,'nseg npseg clkoff rate '
     .                     ,nseg,npseg(nseg),clkoff(nseg),rate(nseg)
              goodseg(nseg) = .true.
            else
              goodseg(nseg) = .false.
            endif
c         endif on new_segment
          endif

c     end of loop over all points  
      enddo
               
c     Now compute the residuals 
      do k=1,nseg               
       if(debug.and.iprn.eq.iprndb) 
     .    print *,'residuals nseg npseg startseg stopseg '
     .           ,k,npseg(k),startseg(k),stopseg(k)
        sum2 = 0.d0
        do i=startseg(k),stopseg(k)
          if(use(i) ) then
            dt = (jd(i)-jd0(nseg))*86400.d0 + (t(i)-t0(nseg))
            resid = clock(i) - (clkoff(nseg)+rate(nseg)*dt)
            sum2 = sum2 + resid**2
            if(debug.and.iprn.eq.iprndb ) 
     .         print *,'i dt clock clkoff rate resid '
     .                ,i,dt,clock(i),clkoff(k),rate(k),resid
          endif
        enddo
        rms(k) = dsqrt(sum2/npseg(k))
      enddo                                

      return
      end

c----------------------------------------------------------------------------
 
      Subroutine poly01c( n_max, n_obs ,x,y,s, constant,slope
     .                  , iprn,debug,iprndb ) 

      implicit none
* calculates (weighted)mean and fits  a 1st order polynomial(y=a+bx)

* n_max dimension of arrays
* n_obs number of observations
      integer*4 n_max,n_obs
* x data, y data, s sigmas, covariance-weights
      real*8 x(n_max),y(n_max),s(n_max),w(n_max)

* sres    -- Scaled residual
* max_res -- Maximum scaled residual
* avg     -- Avergae with the slope set to zero
      real*8  sres, max_res, avg
       
* loops
      integer*4 i
*     account for weigths
      character*10 weigh
* intermediate parameters  sumyw,sumw,sums,sumxw,sumxyw,sumxxw,det
      real*8 sumyw,sumw,sums,sumxw,sumxyw,sumxxw,det
* reference value of x to avoid rounding error with year scales
*     for small amounts of data,  Added min_x and max_x values
      real*8 ref_x, min_x, max_x

* (weighted) mean
      real*8 wmean,unc_mean
* a  + b X     constant :a     slope : b
      real*8 constant,slope,unc_constant,unc_slope
*     chisquare
      real*8 chisqr
* normalized rms, weighted rms
      real*8 nrms, wrms
* zero 
      real*8 zero              

* for debug                
      logical debug 
      integer*4 iprn,iprndb  

      weigh = "no"  
      zero  = 0.d0
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
      do i=1,n_obs
        w(i) = 1.d0
        s(i) = 1.d0
        if(weigh(1:3).eq."yes") then
          if(s(i).ne.0.d0) w(i)=1.d0/(s(i)*s(i))
        endif
* MOD TAH 020304: Remove reference value from x.  Also keep track of
*       min and max x values
        if( i.eq.1 ) ref_x = x(1)
        x(i) = x(i) - ref_x
        if( x(i).lt.min_x ) min_x = x(i)
        if( x(i).gt.max_x ) max_x = x(i)

      enddo

c      write(*,'(i5,3f10.3)') (i, x(i),y(i),w(i),i=1,n_obs)

cd      print *,'POLY01C n_obs ',n_obs 
       
******   WMean **********
*        SUM y*w
*        SUM   w
         sumyw = 0.d0
         sumw  = 0.d0
         sums  = 0.d0
         do i=1,n_obs
              sumyw  = sumyw + y(i)*w(i)
              sumw   = sumw  +      w(i)
              sums   = sums  + s(i)*s(i)
         enddo
         wmean = sumyw / sumw
c        where did the original theory for unc mean come from?
c        current unc_mean from King/McClusky 950904.
c        unc_mean     = sqrt ( sums )
         unc_mean     = sqrt ( 1/sumw )

*******  Offset and slope
*******  y = A + B*x
            sumxw  = 0.d0
            sumxyw = 0.d0
            sumxxw = 0.d0
            do i=1,n_obs        
                  sumxxw  = sumxxw + x(i)*x(i)*w(i)
                  sumxyw  = sumxyw + x(i)*y(i)*w(i)
                  sumxw   = sumxw  + x(i)     *w(i)
            enddo
            det = sumw* sumxxw - sumxw * sumxw
            constant = ( sumxxw * sumyw - sumxw * sumxyw ) / det
            slope    = (-sumxw  * sumyw + sumw  * sumxyw ) / det
            avg      = sumyw/sumw
            unc_constant = sqrt ( sumxxw / det )
            unc_slope    = sqrt ( sumw / det )
	    
********  Statistics
           chisqr = 0.d0
           do i=1,n_obs
                 chisqr = chisqr + ((y(i)-constant-slope*x(i))/s(i))**2
            enddo
            if ( n_obs .gt. 2 ) then
               nrms = sqrt (chisqr/(n_obs-2))     
               wrms = sqrt(float(n_obs)/(float(n_obs)-2.)*chisqr/sumw)   
            else
               nrms = 999.99d0
               wrms = 999.99d0
            endif             

*     Now output the results
      if( debug.and.iprn.eq.iprndb ) then        
         print *,'CLKFIT for PRN ',iprn,' n_obs ',n_obs
        if (weigh(1:3).eq."yes") then
          write(*,600) wmean,unc_mean,
     >        constant,unc_constant,slope,unc_slope, n_obs, nrms,wrms,
     >        n_obs
 600      format("  wMEAN:",2x,2f15.9,
     >           " constant:",2x,2f15.9,
     >           "    slope:",2x,2d16.8,
     >           "      obs:",2x,i5,
     >           "     nrms:",2x,f10.4,
     >           "     wrms:",2x,f10.4,
     >           "  all obs:",2x,i5)       
        else
          write(*,601) wmean,unc_mean,
     >        constant,unc_constant,slope,unc_slope, n_obs, nrms,wrms,
     >        n_obs
 601      format("  wMEAN:",2x,2f15.9,
     >           " constant:",2x,2f15.9,
     >           "    slope:",2x,2d16.8,
     >           "      obs:",2x,i5,
     >           "     nrms:",2x,f10.4,
     >           "     wrms:",2x,f10.4,
     >           "  all obs:",2x,i5)       
        endif
      endif

      end

