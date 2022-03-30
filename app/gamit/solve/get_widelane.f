      Subroutine GET_WIDELANE  
                          
c**   MOVE last_nonbias to common eventually

c     Get the widelane ambiguities from the autcln file or try to resolve them
c     Replaces old routines modelc (mode 1) and getwl; called by SOLVE.   
c     R. King 8 Dec 2003
   

      implicit none
         
c  Nearly global common blocks

      include '../includes/dimpar.h'
      include 'solve.h'
      include 'parameters.h'
          

c  Local variables
                                   
      integer*4 last_nonbias
      integer*4 lbias1dd,lwl,lw,ifirst,isn,ioerr
     .        , ifix,nlive0,mfix,i1,i2,j2,jb,jp
     .        , idum,i,j,k

      logical dep_bias ! Set true when bias is dependent and already
                       ! eliminated. Widelane is set to zero in this case.

      real*8 fac,dummy,usig,ssig     

      character*1 symbl 
          
c* This now set in gethed.f from fL1, fL2 on c-file
c      data gear/0.779220779220779d+00/
      data dummy/0.0d0/
 
c*** test an assumption /block8/ lpart = last_nonbias 
c      print *,'GET_WIDELANE lpart ',lpart
      last_nonbias = lpart          
      
c  Set iband for bias resolution
      iband = 2 

c  Set the local variable for DD biases from that saved from DOPT    
       lbias1dd = l1bias

                                             
c-- Move the non-bias adjusts into a temporary vector

      call copy1d(1,lpart,0,adjust,adorg) 
c      print *,'GET_WIDELANE adj->adorg lpart ',lpart
c      print *,'GET_WL 1 nlive lpart adjust,adorg 1,210,244,245 '
c     .     ,nlive,lpart,adjust(1),adorg(1),adjust(210),adorg(210)
c     .     ,adjust(244),adorg(244),adjust(245),adorg(245)


c-- Read back the L1/L2 solution from unit 27
c**  skip this for LC_AUTCLN?
      
c      call flush(27)
* MOD TAH 030628: See comments above about closing and
*     opening files.
      close(27,iostat=ioerr)
c     see if a local disk rather than the processing directory to be used
      open(unit=27,file=ftmp27, form='unformatted',status='unknown'
     .    , iostat=ioerr)     
      if( ioerr.ne.0 ) call report_stat('FATAL','SOLVE','get_widelane'
     .  ,ftmp27,'Error re-opening scratch unit',ioerr)

C      REWIND 27
      READ (27,iostat=ioerr) R2SUM,SCLERR,NLIVE,CHI2 
      if(ioerr.ne.0) call report_stat('FATAL','SOLVE','get_widelane',' '
     .  ,'Error reading first record of temp file 27',ioerr)
      i1=nlive*(nlive+1)/2
      call tapdrv(27,1,nlive,i1,a,5)
      call tapdrv(27,1,nlive,nlive,sigma,4)
      read (27,iostat=ioerr) (FREE(J),J=1,NTPART)
      if(ioerr.ne.0) call report_stat('FATAL','SOLVE','get_widelane',' '
     .  , 'Error reading second record of temp file 27',ioerr)
      call tapdrv(27,1,ntpart,ntpart,adjust,4)
      READ (27,iostat=ioerr) (IDXB(J),J = 1,MBIAS) 
      if(ioerr.ne.0) call report_stat('FATAL','SOLVE','get_widelane',' '
     . ,'Error reading third record of temp file 27',ioerr)
      NOBS = NOBS*2     
c      print *,'GET_WIDELANE red 27 ntpart nlive adjust 1 210 244 245 '
c     .     ,ntpart,nlive,adjust(1),adjust(210),adjust(244),adjust(245)


c  Use the WL biases read from the AUTCLN file or resolve them here
         
      if( l2flag.eq.4 ) then 

c       LC_AUTCLN  
                   
        call report_stat('STATUS','SOLVE','get_widelane',' '
     .          , 'Fixing wide-lane ambiguities from AUTCLN N-file',0)  
        if( logprt ) write(6,10)
        write(10,10)
 10    format (/,7x,'====Wide-lane ambiguities from the '
     .        , 'N-file (LC_AUTCLN) =====',/)    
        jp = last_nonbias + l1bias
        do jb = 1,l1bias
          jp = jp + 1
          symbl = rorx(jb)
          call qhead4(jb,jp,dummy,dummy,symbl,idum,idum,3)
        enddo       

c       print *,'GET_WL 2 nlive adjust adorg (1,210,244,245) '
c     .     ,nlive,adjust(1),adorg(1),adjust(210),adorg(210)
c     .     ,adjust(244),adorg(244),adjust(245),adorg(245)
                

c       values read in by read_biases, called by dopt, and stored in common /acbiases/
        mfix = 0   
        k = 0                               
c        print *,'GET_WIDELANE wlval ',(wlval(i),i=1,35)
c        print *,' free(211-270) '
c        write(*,'(10i7)') (free(i),i=211,270)
        do i=1, l1bias
            i1 = last_nonbias+lbias1dd+i 

* MOD TAH 080609: Check current status of parameter.  If free(i1) is zero
*           this is a dependent bias and we will force back to zero.  If
*           this is not done then narrow lane is non-integer.
            dep_bias = .false.
            if( free(i1).eq.0 ) then
                dep_bias = .true.
            endif
  
            if( rorx(i).eq.'X' ) then 
              if( .not.dep_bias ) then 
                 adjust(i1) = wlval(i)
              else
                 adjust(i1) = 0.d0
                 write(10,220) i1, wlval(i)
 220             format('Reset Fixed Dependent bias ',i5,' from ',
     .                   F8.1,' to  0.0')
              end if

              mfix = mfix + 1
              ifix = 1
              free(i1) = 0  
              nfix(1) = 1 
              bdev(1) = 0.d0
              nlive0 = nlive
              nlive = nlive-1    
c             print *,'GET+WIDELANE call BNEW i i1 adjust(1 210 244 245)'
c     .          ,i,i1,adjust(1),adjust(210),adjust(244),adjust(245)
              call bnew( nlive0,ifix )                   
c             print *,'  after BNEW adjust(1 210 244 245'
c     .            ,adjust(1),adjust(210),adjust(244),adjust(245)
            else
* MOD TAH 080527: Save the LC-Autcln value into adjust.  Need to make
*           a little different from an integer or solve thinks that it
*           have been fixed
* MOD TAH 080609: Reset to zero if dependent value
              if( .not.dep_bias ) then 
                 adjust(i1) = wlval(i) + sign(1.d-4,wlval(i))
              else
                 adjust(i1) = 0.d0
                 write(10,240) i1, wlval(i)
 240             format('Reset Free  Dependent bias ',i5,' from ',
     .                   F8.1,' to  0.0')
              end if
            endif   
c***          endif
        enddo      
c        print *,'GET_WIDELANE numwl k lbias1dd mfix nlive '
c     .      ,numwl,k,lbias1dd,mfix,nlive    
c        print *,' free(211-270) '
c        write(*,'(10i7)') (free(i),i=211,270)  
c        print *,' adjust(211-270) '
c        write(*,'(10f7.2)') (adjust(i),i=211,270)    
        if( logprt ) write(6,15) l1bias
         write(10,15) l1bias
   15   format(i5,' Phase ambiguities in solution')
        if( logprt ) write(6,20) mfix
        write(10,20) mfix
   20   format(i5,' WL ambiguities resolved by AUTCLN')


c  Autcln file not available, resolve the WL ambiguities the old way
     
      elseif (l2flag.eq.3 ) then   

c       LC_HELP

        call report_stat('STATUS','SOLVE','get_widelane',' '
     .          , 'Resolving wide-lane ambiguities ',0)  
        call resolve_wl(mfix)  

c       write out the summary of biases after fixing the WLs
        if( logprt ) write(6,30) 
        write(10,30)   
  30  format(/,' Summary of biases after fixing WLs '
     .       ,' (Non-bias parameters fixed, uncertainties still '
     .       ,'scaled by WL baseline deviation nrms)')  
      if( logprt ) write(6,40)
      write(10,40)
  40   format(/
     .   ,8x,'Label',18x,'Estimate',5x,'Uncertainty',13x,'Pseudorange',/
     .  ,31x,'(cycles)   Unscaled   Scaled',6x,'Estimate  Uncertainty')
        i1=0
        j = lpart 
        lw = 0
          symbl='*'
          lwl = l1bias
          do j2 = 1,lwl*iband
            j = j + 1
            if (j2.gt.Lwl) lw = lw + 1
            if (free(j).eq.0) then
              symbl = ' '
              call qhead4(j2,j,dummy,dummy,symbl,lwl,lw,1)
            else
              i1=i1+1
              symbl='*'
              usig = sigma(i1)
              ssig = bscale(i1)*usig 
              call qhead4(j2,j,usig,ssig,symbl,lwl,lw,2)
            endif   
          enddo   
         
      else 
        call report_stat('FATAL','SOLVE','get_widelane',' '
     .          , 'Invalid l2flag in WL resolution',0) 
      endif
       
c      print *,'GET_WIDELANE after resolution lbias1dd nlive free'
c     .     ,lbias1dd,nlive
c       do i=1,lbias1dd
c         j= last_nonbias+lbias1dd+i 
c         write(*,'(i4,f7.1,i4)') j,adjust(j),free(j)
c       enddo 

          
c-- Count the remaining live parameters

      nlive = 0    
      do i= 1,ntpart
        if(free(i).eq.1) nlive = nlive+1 
      enddo    
c      print *,'GET_WIDELANE after recount nlive ',nlive



c-- Copy the WL-fixed  adjustments into the temporary array

      i1 = lpart+1
      call copy1d(i1,ntpart,0,adjust,adorg)  
c      print *,'GET_WIDELANE adj->adorg ntpart ',ntpart
     

c-- Read the LC bias free solution back into storage from unit 29

*     MOD TAH 030627: added flush call back
c      print *,'Closing unit 29'
      close(29,iostat=ioerr)
c      print *,' Error closing 29 ',ioerr
      open(unit=29,file=ftmp29, form='unformatted',status='unknown'
     .     ,iostat=ioerr)  
      if( ioerr.ne.0 ) call report_stat('FATAL','SOLVE','get_widelane'
     .   ,ftmp29,'Error re-opening scratch unit',ioerr)   
C      call flush(29)
C     REWIND 29
      read (29,iostat=ioerr) r2sum,sclerr,nlive,chi2,msig  
      if(ioerr.ne.0) call report_stat('FATAL','SOLVE','get_widelane',' '
     . ,'Error reading first record of temp file 29',ioerr)
      i1=nlive*(nlive+1)/2
      call tapdrv(29,1,nlive,i1,a,5)
      call tapdrv(29,1,nlive,nlive,sigma,4)
      read (29,iostat=ioerr) (FREE(J),J = 1,NTPART)   
      if(ioerr.ne.0) call report_stat('FATAL','SOLVE','get_widelane',' '
     .  ,'Error reading second record of temp file 29',ioerr)
      call tapdrv(29,1,ntpart,ntpart,b,4)
      call tapdrv(29,1,ntpart,ntpart,borg,4)
      call tapdrv(29,1,ntpart,ntpart,adjust,4)
      read (29,iostat=ioerr) (idxb(j),j = 1,mbias)
      if(ioerr.ne.0) call report_stat('FATAL','SOLVE','get_widelane',' '
     .  ,'Error reading third record of temp file 29',ioerr)
      nobs = nobs/2     

c      print *,'GET_WIDELANE read 29 r2sum ntpart nlive '
c     .       ,' adjust 1,210,244,245 '
c     . ,r2sum,ntpart,nlive,adjust(1),adjust(210),adjust(244),adjust(245)
c      print *,'a 1 22155 ',a(1),a(22155)


c-- Update the WL biases

      fac = gear/(1.0d0-gear)
      ifirst = idxb(1)
      if (ifirst.lt.0) ifirst = -ifirst
      i1 = ifirst
      i2 = i1-1
        isn = l1bias
        i1 = i1+isn
        i2 = i2+isn*2 
c        print *,'GET_WIDELANE updating biases i1 i2 isn ',i1,i2,isn
        do i = i1,i2
          adjust(i-isn) = adjust(i-isn)+fac*adorg(i)
          adjust(i) = adorg(i)
        enddo
        i1 = i1+isn      
                   
c      print *,'GET_WIDELANE new adjust(1 210 244 245) '
c     .       , nlive,adjust(1),adjust(210),adjust(244),adjust(245)

         
c-- Print the number of live and dead parameters

      call qhead1 (2)

      return
      end



