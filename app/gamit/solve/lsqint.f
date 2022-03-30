C
      Subroutine LSQINT

c    This routine is now almost a shell.  It call LSQIO to convert
c    coordinates to cartesian, and prints out the observable options,
c    sites, and satellites

      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'                     
      include 'parameters.h'


C Convert site coordinates to Cartesian

      call lsqio

c Print the observable 
                              
      if( l2flag.eq.0 ) 
     .     write(10,'(a)') ' L1 receiver - only L1 observations'
      if( l2flag.eq.-1 ) 
     .     write(10,'(a)') ' L1 observations only' 
      if( l2flag.eq.-2 ) write(10,'(a)') ' L2 observations only'
      if( l2flag.eq.1 ) write(10,'(a)') 
     .    ' Ionospheric corrected observations (LC)'
      if( l2flag.eq.2 ) write(10,'(a)') 
     .    ' Separate L1 and L2 observations'    
      if( l2flag.eq.3 ) write(10,'(a)') 
     .    ' LC solution with ionosphere constraint bias-fixing'
      if(l2flag.eq.4) then
        write(10,'(a)') '  LC solution with AUTCLN bias-fixing'
        if( bias_apr.gt.0.d0 ) write(10,'(a,f6.0,a)')  
     .    ' --Bias constraints = ',bias_apr,' cycles'    
      endif
      if( l2flag.eq.5 ) write(10,'(a)') 
     .    ' L1 and L2 independent observations'
                        
c Print the station information

c** rwk 100210: This output is superfluous, so remove it as we add rcvr/ant info further down 
c               (rcvr/ant info not available at this point in the run)
c      if( logprt ) write(6,'(/,a,7x,a)') ' Tracking stations'
c     .   ,'Receiver            SwVer   Antenna'
c      write(10,'(/,a,7x,a)') ' Tracking stations'
c     .  ,'Receiver            SwVer   Antenna'
c      i1 = 1
c      do is=1,nsite
c        if( logprt ) 
c     .       write(6,'(i3,2x,a4,2x,a12,2x,a20,1x,f5.2,2x,a20)') 
c     ,       is,rlabel(i1)(1:4),sitnam(is)
c     .      ,rcvr_type(is),rcvr_swver(is),ant_type(is)
c        write(10,'(i3,2x,a4,2x,a12,2x,a20,1x,f5.2,2x,a20)') 
c     .       is,rlabel(i1)(1:4),sitnam(is) 
c     .      ,rcvr_type(is),rcvr_swver(is),ant_type(is)
c        i1 = i1 + 3
c      enddo

c Print the satellite     
c** rwk 100210: This summary also superfluous, remove in favor of a horizontal listing of
c               PRN (as well as channel)  in the 'Satellites used' section
c      if( logprt ) write(6,'(/,a,/,(i3,a,i2))') ' Satellites observed'
c     .    , (is,'.. PRN',isat(is),is=1,nsat)
c      write(10,'(/,a,/,(i3,a,i2))') ' Satellites observed'
c     .    , (is,'.. PRN',isat(is),is=1,nsat)

      RETURN
      END
