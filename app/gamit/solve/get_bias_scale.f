      Subroutine GET_BIAS_SCALE( mopt ) 

c     Calculate a scale factor for each L1 or L2-L1 bias sigma, common 
c     for each baseline.   R. King 070221
  
      implicit none 

      include '../includes/dimpar.h'
      include 'solve.h'
      include 'parameters.h'


c   Input                      
 
      integer*4  mopt

c      mopt = 1  NL (L1 or L2) biases
c      mopt = 2  WL L2-L1 biases
c      mopt = 3  L1 and L2 biases (independent)

c   Output      

c      real*8 bscale(maxbis) in /bbii/ of solve.h, dimension is maxbis, # of L1 +  L2-L1 (or L2) biases
c        value is scale factor for each bias, same for all biases on a given baseline

                   
c   Common for wavelength factors for biases 
c      integer half(maxbis),lwave
c      common/lhalf/lwave(maxsit,maxsat,2),half
 
c  Adjustments (from solve.h:
c        adjust(maxprm) : adjustments    

c  Function to determine is WL for corresponding NL has been resolved
c        ---appended to bfwork.f
      logical wl_fixed
                    
c  Local
      integer*4 maxbas
      parameter (maxbas=maxsit*(maxsit-1)/2)
      integer*4 numbias(maxbas),bindx(maxbis),blhalf(maxbas)
     .      ,nbl,live,ia1,ia2,la1,la2,ih,i
      real*8 sum2dev(maxbas),sum2rat(maxbas),blscale(maxbas)
     .     , bldev(maxbas),dadj,dev,adjwl
      character*2 alane 
      character*9 bllabel(maxbas),lastlabel  
      character*256 message
      logical use_dev

c** Flag for debug printout
      integer*4 idebug
c      idebug = 0  : none 
c               1  : print to q- and log file decision values for fixed and not-fixed biases
c               2  : print to log the fixing criteria
c               3  : full print of logic 
       idebug = 1


c Initialize the arrays

      do i=1,maxbas  
        numbias(i) = 0
        sum2dev(i) = 0.d0 
        sum2rat(i) = 0.d0 
        blscale(i) = 1.d0
      enddo            
      do i=1,maxbis
        bscale(i) = 1.d0 
      enddo

c Loop through the parameter array to accumalate the normaized 
c root-sum-squares for each baseline
               
      if( idebug.ge.3) print *,'GET_BIAS_SCALE mopt lpart l1bias nlive '
     .             ,mopt,lpart,l1bias,nlive     
      lastlabel = ' '  
c     set indices for bias loop, both within the absolute parameter array
c     (for adjust values) and array of free parameters (for sigmas)
c     lpart is # of non-bias parameters; l1bias is the # of L1 (or L2-L1) biases 
      if( mopt.eq.1 .or. mopt.eq.3 ) then
c       NL:  non-bias prameters may be fixed or free 
c       L1 L2 independent: use the L1 scale factors for L2  
        ia1 = lpart + 1
        ia2 = ia1 -1 + l1bias
        la1 = 1
        la2 = lpart
        live = 0 
        do i=la1,la2
          if(free(i).gt.0) live = live + 1 
        enddo 
        if( idebug.ge.3) print *,'NL ia1 ia2 la1 la2 live '
     .     ,ia1,ia2,la1,la2,live
      elseif( mopt.eq.2 ) then   
c       WL:  all parameters fixed except biases   
        ia1 = lpart + l1bias + 1               
        ia2 = ia1 -1 + l1bias         
        la1 = lpart+1
        la2 = lpart + l1bias           
        live = 0 
        do i=la1,la2
         if(free(i).gt.0 ) live = live + 1
        enddo
        if(idebug.ge.3) print *,'WL ia1 ia2 la1 la2 live '
     .        ,ia1,ia2,la1,la2,live     
      else    
c       without this branch, live,ia1,ia2 may be logically uninitialized 
        print *,'SOLVE/get_bias_scale: mopt must be 1,2, or 3'
        stop
      endif
      ibias = 0   
      nbl = 0 
      do i=ia1,ia2 
        if( free(i).gt.0 ) then
          live = live + 1  
          ibias = ibias + 1 
          if( mopt.eq.2 .or.
     .        (mopt.eq.3.and.iband.eq.2) ) then
            ih = half(ibias + l1bias)
          else
            ih = half(ibias)
          endif
          if(idebug.ge.3) print *,'i rlabel ',i,rlabel(i)
          if( rlabel(i)(6:14).ne.lastlabel ) then
            nbl = nbl + 1           
            if( nbl.gt.maxbas ) then
              write(message,'(a,i4,a)') '# baselines > maxbas ('
     .           ,maxbas,')'    
              call report_stat('STATUS','SOLVE','get_bias_scale',' '
     .                   ,message, 0)
            endif
            bllabel(nbl) = rlabel(i)(6:14)   
            if(idebug.ge.3) print *,'i nbl bllabel lastlabel '
     .                    ,i,nbl,lastlabel
            lastlabel = bllabel(nbl)
            blhalf(nbl) = ih   
          endif            
c         set the baseline pointer for each bias
          bindx(ibias) = nbl         

c         see if bias to be included in calculation:
c            exclude if sigma > 2   (catch weakly determined biases; 
c                    the limit >1 is ok because the a priori scaling might be wrong
c            exclude if NL with WL not fixed
          use_dev = .true.
          if( sigma(live).ge.2.d0 ) use_dev = .false.
          if( mopt.eq.1.and.iband.eq.2 ) then
             adjwl = adjust(i+l1bias)
             if( .not.wl_fixed(ih,adjwl) ) use_dev = .false. 
          endif 
          if(idebug.ge.3) then
            print *,' NLi WLi live sigmaNL sigmaWL freeNL freeWL '
     .             ,   i,i+l1bias,live,sigma(i),sigma(i+l1bias)
     .             ,   free(i),free(i+l1bias)
          endif  

          if( use_dev ) then  
            numbias(nbl) = numbias(nbl) + 1  
            if(idebug.ge.3) print *,'i nbl numbias adjust live sigma '
     .         , i,nbl,numbias(nbl),adjust(i),live,sigma(live) 
            dadj = dabs( dnint( adjust(i)*dble(ih) ) )
            dev = ( dabs(adjust(i)*dble(ih)) - dadj )
     .         / dble(ih)
            sum2dev(nbl) = sum2dev(nbl) + dev**2
            sum2rat(nbl) = sum2rat(nbl) + (dev/sigma(live))**2
            if(idebug.ge.3) 
     .         print *,'nbl ibias,ih,dadj dev sum2dev sum2rat '
     .        ,nbl,ibias,ih,dadj,dev,sum2dev(nbl),sum2rat(nbl)  
          endif
        endif
      enddo
                          
c Calculate the nrms of deviations for each baseline 
c (if number of biases for a baseline less than 5, skip
c and use scale for previous baseline)
       
      do i=1,nbl
        if( numbias(i).ge.5 ) then  
          bldev(i) = dsqrt(sum2dev(i)/dfloat(numbias(i)))
          blscale(i) = dsqrt(sum2rat(i)/dfloat(numbias(i))) 
          if(idebug.ge.3) 
     .       print *,'i numbias bldev blscale '
     .          ,i,numbias(i),bldev(i),blscale(i)
         else 
           if( i.ge.1 ) then   
             if( numbias(i).gt.0 ) then
               bldev(i) = dsqrt(sum2dev(i)/dfloat(numbias(i))) 
             else
               bldev(i) = 0.d0
             endif
             blscale(i) = blscale(i-1)     
             if(idebug.ge.3) 
     .         print *,'i numbias bldev blscale '
     .             ,i,numbias(i),bldev(i),blscale(i)
           endif
         endif
      enddo

c Print the scale factors into the q-file 
                     
      alane = 'NL'
      if( mopt.eq.2 ) alane = 'WL' 
      write(10,1) alane 
   1  format(/,2x,a2
     .    ,' Baseline   Wavelength  # Biases  RMS Deviation '
     .    ,' Sigma Scale Factor')
      if( logprt ) write(6,1)
      do i=1,nbl       
        write(10,'(3x,a9,10x,i1,6x,i3,4x,f10.3,10x,f6.3)') 
     .      bllabel(i),blhalf(i),numbias(i),bldev(i),blscale(i)
        if( logprt ) write(6,'(3x,a9,6x,f6.3,12x,f6.3)') 
     .      bllabel(i),blhalf(i),numbias(i),bldev(i),blscale(i)
      enddo  
      write(10,'(1x)')
          
c Assign the scale factors to each bias parameter 

c    When resolving WLs, use the WL scale factor also for
c    the NL biases for consistency in the printout

c    When resolving L1 + L2 independent (together), use the
c    L1 scale factors also for L2.

      if(idebug.ge.2 )print *,'End of GET_BIAS_SCALE ibias ',ibias
      do i=1,ibias
        bscale(i) = blscale(bindx(i))     
        if(idebug.ge.2) print *,'i bscale ',i,bscale(i)
        if( mopt.eq.2 .or. mopt.eq.3 ) then
          bscale(i+ibias) = bscale(i)  
          if(idebug.ge.2)  print *,'i+ibias bscale '
     .       ,i+ibias,bscale(i+ibias)
        endif
      enddo     

c Hidden feature to turn off NL rescaling:
      if( mopt.eq.1 .and. .not.nlscale ) then
        write(10,'(/,a,/)') '** NL rescaling not used **'  
        if( logprt ) write(6,'(a)') '** NL rescaling not used **'
        do i=1,ibias
          bscale(i) = 1.d0
        enddo
      endif

      return
      end
      
