      Subroutine check_met(ndim,nmet,code,met_val,ngood)

c      Routine to check a single or array of met values for reasonableness.
                                                                            
      integer*4 ndim,nmet,ngood
      real*4 met_val(ndim)  
      character*1 code
            
c  Input
c     ntm  et  :  number of (time) values
c     code     : single-character code for values
c                 'Z'  zenith hydrostatic delay (mm)
c                 'W'  zenith set delay  (m)
c                 'P'  pressure (hPa)
c                 'T'  temperature (C)
c                 'H'  relative humidity (%)
c     met_val  :  values
c     
c   Output 
c       ngood =  0   no time with all good values
c       ngood =  1   first (or only) value ok
c       ngood >  1   first value not ok, use the 'ngood'th values
                          
      logical good
                   

      ngood = 1
      good = .false.
      do while (.not.good)
 
        if( code.eq.'Z' ) then     
          if( ngood.gt.nmet ) then
            ngood = 0
            return  
          endif
          if( met_val(ngood).gt.1000.d0 .and. 
     .        met_val(ngood).lt.2500.d0 ) then
            good = .true.   
            return
          endif

        elseif( code.eq.'W' ) then  
          if( ngood.gt.nmet ) then
            ngood = 0
            return
          elseif( met_val(ngood).gt.-1.d0 .and. 
     .        met_val(ngood).lt.0.5d0 ) then
            good = .true.  
            return  
          endif

        elseif( code.eq.'P' ) then   
          if( ngood.gt.nmet ) then
            ngood = 0
            return
          elseif( met_val(ngood).gt.600.d0 .and. 
     .        met_val(ngood).lt.1100.d0 ) then
            good = .true.    
            return
          endif

        elseif( code.eq.'T' ) then   
          if( ngood.gt.nmet ) then
            ngood = 0
            return
          elseif( met_val(ngood).gt.-100.d0 .and. 
     .       met_val(ngood).lt.60.d0 ) then
            good = .true.
            return
          endif

        elseif( code.eq.'H' ) then  
          if( ngood.gt.nmet ) then
            ngood = 0
            return
          elseif( met_val(ngood).gt.-1.d0 .and. 
     .        met_val(ngood).lt.100.d0 ) then
            good = .true.
            return
          endif                                         
           
        endif  

        ngood = ngood + 1 

      enddo

      return
      end

          

