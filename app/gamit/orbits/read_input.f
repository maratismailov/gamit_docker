c
      Subroutine read_input

c     R King 2 December 2000, based on S. McClusky routine of same name

c-----Read batch file with key word style               


      implicit none

      include '../includes/dimpar.h'
      include '../includes/orbits.h'
      include 'orbfit.h'

      character*16 code   
      character*256 wcmd            
      character*256 message

      integer nsvs_excl,excl_prn(maxsat)
      integer ib,lift_arg,count_arg
      integer*4 iarray
          
      integer ip(3),iv(3),ir(maxorb),i,j,jj
      integer lcmd 
      integer ierr

      logical debug/.false./

c----- Read all input controls

c------------ Check if any satellites to be excluded and set list-------------
c
      lcmd = 0
      call getcmd(lucmd,'exclude',wcmd,lcmd,1)   
      nsvs_excl = 0
      if (lcmd .ne. 0 .or. nbad_sat .ne. 0 ) then
        nsvs_excl = count_arg(wcmd)
        if (nsvs_excl.ne.0 .or. nbad_sat .ne. 0 ) then
          do  i = 1,nsvs_excl
            ib = lift_arg(wcmd,code,i)
            read(code(:ib),*,iostat=ierr) excl_prn(i)
            if( ierr.ne.0 ) then
              print *,'IOSTST ierr ',ierr, i, ' ib ',ib, 
     .                ' code ',code(3:ib),'|'
              excl_prn(i) = 9
            endif
          enddo	  
	  if ( nbad_sat .gt. 0 ) then
	    if ( ibad_sat(nbad_sat) .gt. 0 ) then
	      do i = 1,nbad_sat
	        excl_prn(nsvs_excl + i) = ibad_sat(i)
	      enddo
	      nsvs_excl = nsvs_excl + nbad_sat
	    endif        
       if(debug) print *,'nsvs_excl,excl_prn ',nsvs_excl,excl_prn 
	  endif
          write(message,'(a,45i3)') 'Excluding PNs '
     .                       ,(excl_prn(i),i=1,nsvs_excl)
          call report_stat('STATUS','ORBFIT','orbits/orbfit',' '
     .         ,message,0)
        endif
      endif        
      nsat = 0
      do i=1,nprn(1)
        if( iarray(itsat(i,1),excl_prn(1),nsvs_excl).eq.0 ) then
          nsat = nsat + 1
          isat(nsat) = itsat(i,1)
          satnam(nsat) = satnamt(i,1)
        endif
      enddo
cd     if(debug) print *,'iter nbad_sat ',iter,nbad_sat(iter) 
      
     
c------------ Get the max misfit tolerance value (if given) -------------
c
      lcmd = 0
      max_fit_tol = 0.0
      call getcmd(lucmd,'max_fit_tol',wcmd,lcmd,1)   
      if( lcmd.gt.0 )  then
        read(wcmd,*) max_fit_tol
        write(message,'(a,f6.4)') 'Max orbit misfit tolerance is: '
     .                       ,max_fit_tol
        call report_stat('STATUS','ORBFIT','orbits/orbfit',' '
     .         ,message,0)
      endif
             
c------------------- Set parameters to be estimated -----------------------

      do i=1,mxoprm
       islot(i) = 0
       prmnam(i) = ' ' 
      enddo         
      nparam = 0
                 
c     Pointers to parameters (islot values)
c       1-3      translation
c       4-6      inertial rotation
c       7        scale
c       8-10     terrestrial rotation
c       11-13    terrestrial rotation rate   
c       n01-n19  satellite parameters (partials interpolated from ephemeris)
c                  where n is the satellite index on the reference T-file
c                  (n=1-nsat)

c-----Translation (default is to estimate)
                         
      call zero1i( 1,3,ip )
      call getcmd(lucmd,'trans',wcmd,lcmd,1)    
      if( lcmd.gt.0 )  read(wcmd,*) ip(1),ip(2),ip(3)
      do i = 1,3  
        if (lcmd.le.0 .or. ip(i).gt.0 ) then 
          nparam = nparam + 1
          islot(nparam) = i  
          apr_prm(nparam) = 0.d0
          if( i.eq.1 ) write(prmnam(nparam),'(i3,a)') nparam
     .       ,'. X-TRANS (m)'   
          if( i.eq.2 ) write(prmnam(nparam),'(i3,a)') nparam
     .       ,'. Y-TRANS (m)'  
          if( i.eq.3 ) write(prmnam(nparam),'(i3,a)') nparam
     .       ,'. Z-TRANS (m)'  
        endif
      enddo
            
c-----Inertial rotation (default is not to estimate)
          
      call zero1i( 1,3,ip )
      call getcmd(lucmd,'i_rot',wcmd,lcmd,1)    
      if( lcmd.gt.0 ) read(wcmd,*) ip(1),ip(2),ip(3)
      do i = 1,3  
        if (ip(i).gt.0 ) then 
          nparam = nparam + 1
          islot(nparam) = 3 + i   
          apr_prm(nparam) = 0.d0
          if(i.eq.1) write(prmnam(nparam),'(i3,a)') nparam
     .       ,'. INERT X-ROT (mas)'
          if(i.eq.2) write(prmnam(nparam),'(i3,a)') nparam
     .       ,'. INERT Y-ROT (mas)'
          if(i.eq.3) write(prmnam(nparam),'(i3,a)') nparam
     .       ,'. INERT Z-ROT (mas)'
        endif
      enddo  
 
             
c-----Scale (default is to estimate)
          
      call zero1i( 1,3,ip )
      call getcmd(lucmd,'scale',wcmd,lcmd,1)    
      if( lcmd.gt.0 )  read(wcmd,*) ip(1)
      if (lcmd.le.0 .or. ip(1).gt.0 ) then 
         nparam = nparam + 1
        islot(nparam) = 7 
        apr_prm(nparam) = 0.d0
        write(prmnam(nparam),'(i3,a)') nparam,'. SCALE (ppb)'
      endif

c-----Terrestrial rotations (default is to estimate)
          
      call zero1i( 1,3,ip )
      call getcmd(lucmd,'t_rot',wcmd,lcmd,1)    
      if( lcmd.gt.0 ) read(wcmd,*) ip(1),ip(2),ip(3) 
      do i = 1,3  
        if (ip(i).gt.0 ) then 
          nparam = nparam + 1
          islot(nparam) = 7 + i  
          apr_prm(nparam) = 0.d0
          if(i.eq.1) write(prmnam(nparam),'(i3,a)') nparam
     .        ,'. TERR X-ROT (mas)'
          if(i.eq.2) write(prmnam(nparam),'(i3,a)') nparam
     .        ,'. TERR Y-ROT (mas)'
          if(i.eq.3) write(prmnam(nparam),'(i3,a)') nparam
     .        ,'. TERR Z-ROT (mas)'
        endif
      enddo
    
c-----Terrestrial rotation rates (default is not to estimate)
               
      call zero1i( 1,3,ip )
      call getcmd(lucmd,'t_rat',wcmd,lcmd,1)
      if( lcmd.gt.0 ) read(wcmd,*) ip(1),ip(2),ip(3)
      do i = 1,3  
       if (ip(i).gt.0 ) then 
        nparam = nparam + 1
        islot(nparam) = 10 + i    
        apr_prm(nparam) = 0.d0
* MOD TAH 200608: Fixed Y and Z labels (all said X)
        if(i.eq.1) write(prmnam(nparam),'(i3,a)') nparam
     .     ,'. TERR X-ROTD (mas/d)'
        if(i.eq.2) write(prmnam(nparam),'(i3,a)') nparam
     .     ,'. TERR Y-ROTD (mas/d)'
        if(i.eq.3) write(prmnam(nparam),'(i3,a)') nparam
     .     ,'. TERR Z-ROTD (mas/d)'
       endif
      enddo

c------------------------ Set Satellite Parameters ---------------------
       
c                       (defaults are not to estimate)
c     position
      call getcmd(lucmd,'pos',wcmd,lcmd,1)
      if (lcmd.gt.0) read(wcmd,*) ip(1),ip(2),ip(3) 
c     velocity
      call getcmd(lucmd,'vel',wcmd,lcmd,1)  
      if (lcmd.gt.0) read(wcmd,*) iv(1),iv(2),iv(3)  
c     radiation parameters
      call getcmd(lucmd,'srad',wcmd,lcmd,1)     
      if (lcmd.gt.0) read(wcmd,*) (ir(i),i=1,13)
      do j=1,nsat    
        jj = iarray(isat(j),itsat(1,1),nprn(1) )
        do i=1,3
          if (ip(i).gt.0 ) then
             nparam = nparam + 1
             islot(nparam) = j*100 + i   
             apr_prm(nparam) = satics(i,jj)
             write(prmnam(nparam),'(i3,a,i2,a,a4,a)') nparam
     .          ,'. PN',isat(j),' ',icsnam(i),' (m)'
          endif
        enddo
        do i=1,3
          if (iv(i).gt.0 ) then
             nparam = nparam + 1
             islot(nparam) = j*100 + 3 + i  
             apr_prm(nparam) = satics(i+3,jj) 
             write(prmnam(nparam),'(i3,a,i2,a,a4,a)') nparam
     .          ,'. PN',isat(j),' ',icsnam(i+3),' (m/s)'
          endif
        enddo
       do i=1,13
          if (ir(i).gt.0 ) then
            nparam = nparam + 1
            islot(nparam) = j*100 + 6 + i      
            apr_prm(nparam) = satics(i+6,jj)
             write(prmnam(nparam),'(i3,a,i2,a,a4)') nparam
     .          ,'. PN',isat(j),' ',icsnam(i+6)
          endif
        enddo
      enddo
               
      if(debug) then                
        write(*,'(a,i4,a)') 'NPARAM=',nparam,' ISLOT:'
        write(*,*) (islot(i),i=1,nparam)
      endif 

C DEBUG               
c      write(iscrn,'(a,35i4)') 'nprn1,itsat',nprn(1)
c     .      ,(itsat(i,1),i=1,nprn(1))
c      write(iscrn,'(a,35i4)') 'nsat,isat',nsat,(isat(i),i=1,nsat)
c      write(iscrn,'(a,i4,a)') 'NPARAM=',nparam,' ISLOT:'
c      write(iscrn,*) (islot(i),i=1,nparam)
      return
      end

      subroutine zero1i(istart,iend,irray)
c     initialize 1-d integer array

      implicit none
      integer iend,irray,istart,i
      dimension irray(iend)

      do i = istart,iend
         irray(i) = 0
      enddo

      return
      end

