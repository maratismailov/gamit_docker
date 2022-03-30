      Subroutine otlcmc( jd, t, otidemod, notl, icall, dx )

c     Calculate corrections for ocean tidal loading to transform between a terrestrial, 
c     solid Earth (CE) frame and a joint (E + oceans) center-of-mass (CM) frame.

      implicit none              

      integer*4 jd,notl,icall
      real*8 t,dx(3)
      character*8 otidemod
                                
c  Input                  
c      jd       : PEP Julian day
c      t        : UT1 seconds-of-day  (UTC ok as approximation)    
c      otidemod : Ocean tide model 
c      notl     : number of tidal components (54 for NAO, 11 for others so far)
c      icall    : =1 read in the correction coefficients; =2 compute the corrections

c  Output
c     dx  : Cartesian offset of the CE frame with respect to the CM frame
     
      logical unit_ok,found
      integer*4 len,rcpar,ioerr,i,j
      character*80 prog_name     
      character*120 line                        
      character*256 message
      real*8 xjd,xmjd,oangle(54),cmccof(54,6)   
     

      save cmccof 
        
c  Get the program name for report_stat calls

c     get calling program name and X-file name for report_stat
      len = rcpar(0,prog_name)
  
      if( icall.eq.1 ) then        

c  Read the CMC correction coefficients      
               

c     make sure the unit number is not taken
      inquire(unit=46,exist=unit_ok,iostat=ioerr) 
      if( unit_ok ) then
        open(file='otlcmc.dat',unit=46,iostat=ioerr,status='old')
        if( ioerr.ne.0 )   call report_stat('FATAL',prog_name
     .            ,'lib/otlcmc','otlcmc.dat'
     .        ,'Error opening file for COM ocean loading correction',0)
      else
        call report_stat('FATAL',prog_name,'lib/otlcmc',' '
     .    ,' Unit 46 not available to open otlcmc.dat',ioerr)
      endif 
      found = .false.
      do while( .not.found ) 
        read(46,'(a)',iostat=ioerr) line 
        if( ioerr.eq.-1 ) then    
          write(message,'(a,a7,a)') 
     .      'EOF on otlcmc.dat, model ',otidemod(1:7),' not found'
          call report_stat('FATAL',prog_name,'lib/otlcmc'
     .      ,' ',message,ioerr)
        elseif( ioerr.ne.0 ) then
           call report_stat('FATAL',prog_name,'lib/otlcmc','otlcmc.dat'
     .        , 'Error reading OTL model name',ioerr)
        else   
          if( line(1:5).eq.'MODEL' ) then     
             if(line(7:13).eq.otidemod(1:7)) then      
               found = .true.
               do i=1,notl
                 read(46,'(41x,3(2x,2e12.4))',iostat=ioerr) 
     .              (cmccof(i,j),j=1,6)   
                 if( ioerr.ne.0)  then
                   write(message,'(a,a8,a,i2,a)') 
     .               'Error reading OTL CMC values for model ',otidemod
     .               ,' (notl=',notl,')'
                   call report_stat('FATAL',prog_name,'lib/otlcmc'
     .                ,'otlcmc.dat',message,ioerr)
                 endif
               enddo 
             else
              continue
             endif  
          else
            continue
          endif
        endif
      enddo      
      close( unit=46 )         
      return
      
c-------------end of initial call to read and save coefficients

      elseif( icall.eq.2 ) then

c  Get the angular arguments

      if( notl.eq.11 ) then 

c       Scherneck routine expects true Julian date (= PEP_JD -0.5)
        xjd = jd - 0.5d0
        call ocearg ( t,xjd,oangle ) 
                             
c       if# tidal components = 54, assume NAO
      elseif (notl.eq.54 ) then

c        Matsumoto routine expects Modified Julian date 
c           (= JD - 2400000 = PEP_JD - 2400001 )  
        xmjd = jd - 2400001  + t/86400.d0
        call ocearg2 ( xmjd, oangle )  
  
      endif
          

c  Compute the corrections
      
      dx(1) = 0.d0
      dx(2) = 0.d0
      dx(3) = 0.d0
      do i=1,notl       
c       Note: Scherneck tabulates coefficients in order Z X Y 
         dx(1) = dx(1) + cmccof(i,3)*dcos(oangle(i)) 
     .                 + cmccof(i,4)*dsin(oangle(i))
         dx(2) = dx(2) + cmccof(i,5)*dcos(oangle(i))
     .                 + cmccof(i,6)*dsin(oangle(i))
         dx(3) = dx(3) + cmccof(i,1)*dcos(oangle(i)) 
     .                 + cmccof(i,2)*dsin(oangle(i))
      enddo
               
c      print *,'OTLCMC xjd, oangle ',xjd,(oangle(i),i=1,notl) 
c     endif on icall
      endif
                 
      return
      end

                

                                 

      
