      Subroutine CHECK_ADJUST

c     Check the magnitudes of the adjustments to coordinates, and
c     write into the q-file any that are greater than than the
c     update tolerance or greater than twice the a priori constraints.

c     R. King 17 January 2007

c     The adjustments larger than the tolerance (normally 10-30 cm) 
c     indicate possible non-linearlities and are used updat3/upl
c     and sh_gamit to update the l-file for the next iteration.

c     The adjustments larger than twice the a priori constraints
c     cause a bad solution.  These may be used by sh_sigelv to 
c     loosen the constraints for the next iteration.

                             

      implicit none     

      include '../includes/dimpar.h'
      include 'solve.h'               
      include 'parameters.h'
                    
      integer*4 istat,iparm,ilive,nlarge 
            
      real*8  erad,adjm,sigm  

      character*120 message

      data erad/6378137.d0/


c  Check for adjustments larger than the sestbl tolerance for updating the l-file

      write(message,'(a,f5.2)') 
     .  'Adjustments larger than l-file tolerance:',coord_upd_tol
      write(10,'(/,a)') message
      if( logprt ) write(6,'(/,a)') message
      ilive = 0                
      nlarge = 0
      do istat = 1, nsite
c         latitude (radians) 
        iparm = 3*(istat-1) + 1     
        if( free(iparm).gt.0 ) then
          ilive = ilive + 1  
          adjm = adjust(iparm)*erad 
          sigm = sigma(ilive)*erad
c          print *,'istat iparm ilive adjm sigm ',istat,iparm,ilive
c     .           , adjm,sigm  
          if( dabs(adjm).gt.coord_upd_tol .and. sigm.lt.10.d0 ) then  
            nlarge = nlarge + 1
            write(message,'(a5,a,1x,a14,2x,f8.3)') 
     .           keyword(12),'LFTOL', rlabel(iparm)(1:14),adjm
            write(10,'(a)') message
            if( logprt ) write(10,'(a)') message
          endif
        endif
c         longitude (radians)
        iparm = 3*(istat-1) + 2
        if( free(iparm).gt.0 ) then
          ilive = ilive + 1            
          adjm = adjust(iparm)*dcos(postvl(iparm-1))*erad 
          sigm = sigma(ilive)*dcos(postvl(iparm-1))*erad  
c          print *,'istat iparm ilive adjm sigm ',istat,iparm,ilive
c     .           , adjm,sigm  
          if( dabs(adjm).gt.coord_upd_tol .and. sigm.lt.10.d0 ) then  
            nlarge = nlarge + 1  
            write(message,'(a5,a,1x,a14,2x,f8.3)')  
     .           keyword(12),'LFTOL', rlabel(iparm)(1:14),adjm 
            write(10,'(a)') message
            if( logprt ) write(6,'(a)') message
          endif
        endif
c         radius  
        iparm = 3*(istat-1) + 3
        if( free(iparm).gt.0 ) then  
         ilive = ilive + 1    
         adjm = adjust(iparm)*1.d3
         sigm = sigma(ilive)*1.d3      
c          print *,'istat iparm ilive adjm sigm ',istat,iparm,ilive
c     .           , adjm,sigm  
          if( dabs(adjm).gt.coord_upd_tol .and.sigm.lt.10.d0 ) then       
            nlarge = nlarge + 1
            write(message,'(a5,a,1x,a14,2x,f8.3)')  
     .            keyword(12),'LFTOL', rlabel(iparm)(1:14),adjm    
            write(10,'(a)') message
            if( logprt ) write(6,'(a)') message 
          endif
        endif
      enddo 
      if( nlarge.eq.0 ) write(10,'(a)') '  NONE'


c  Check for adjustments larger than twice the a priori constraint
                

      write(message,'(a)') 
     .   'Adjustments larger than twice the a priori constraint:'
      write(10,'(/,a)') message
      if( logprt ) write(6,'(/,a)') message
      ilive = 0   
      nlarge = 0
      do istat = 1, nsite
c         latitude (radians) 
        iparm = 3*(istat-1) + 1     
        if( free(iparm).gt.0 ) then
          ilive = ilive + 1  
          adjm = adjust(iparm)*erad 
          sigm = sigma(ilive)*erad  
c          print *,'istat iparm ilive adjm sigm ',istat,iparm,ilive
c     .           , adjm,sigm  
          if( dabs(adjm).gt.2.d0*1.d3*stat_apr(istat,1) ) then  
            nlarge = nlarge + 1 
            write(message,'(a5,a,1x,a14,2x,f8.3)')
     .         keyword(12),'APTOL', rlabel(iparm)(1:14),adjm
            write(10,'(a)') message
            if( logprt ) write(6,'(a)') message
          endif
        endif
c         longitude (radians)
        iparm = 3*(istat-1) + 2
        if( free(iparm).gt.0 ) then
          ilive = ilive + 1   
          adjm = adjust(iparm)*dcos(postvl(iparm-1))*erad 
          sigm = sigma(ilive)*dcos(postvl(iparm-1))*erad  
c          print *,'istat iparm ilive adjm sigm ',istat,iparm,ilive
c     .           , adjm,sigm  
         if( dabs(adjm).gt.2.d0*1.d3*stat_apr(istat,2) ) then    
            nlarge = nlarge + 1
            write(message,'(a5,a,1x,a14,2x,f8.3)')  
     .           keyword(12),'APTOL', rlabel(iparm)(1:14),adjm 
            write(10,'(a)') message
            if( logprt ) write(6,'(a)') message
          endif
        endif
c         radius  
        iparm = 3*(istat-1) + 3
        if( free(iparm).gt.0 ) then  
         ilive = ilive + 1    
         adjm = adjust(iparm)*1.d3
         sigm = sigma(ilive)*1.d3  
c          print *,'istat iparm ilive adjm sigm ',istat,iparm,ilive
c     .           , adjm,sigm  
         if( dabs(adjm).gt.2.d0*1.d3*stat_apr(istat,3) ) then  
            nlarge = nlarge + 1
            write(message,'(a5,a,1x,a14,2x,f8.3)')  
     .           keyword(12),'APTOL', rlabel(iparm)(1:14),adjm 
            write(10,'(a)') message
            if( logprt ) write(6,'(a)') message
          endif
        endif
      enddo 
      if( nlarge.eq.0 ) write(10,'(a)') '  NONE'

      
      return
      end
