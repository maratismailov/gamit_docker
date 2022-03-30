      Subroutine LC_SOLUTION  
                          
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
      integer*4 idiag,nded,ioerr
     .        , ierinv,ifast,i,i1,j 
      real*8 goodft,rnew  

     
      logical zeroad

c  Function
   
      integer*4 jelf

c*** test an assumption /block8/ lpart = last_nonbias 
c      print *,'GET_WIDELANE lpart ',lpart
      last_nonbias = lpart          


c  Get the sigmas from the solution   (do this earlier??)
             
      do i=1,nlive
        idiag=i*(i+1)/2
        idiag = jelf(i+1)
        sigma(i) = dsign(dsqrt(dabs(a(idiag))),a(idiag))
      enddo
      

c  Store the full L1/L2 solution  (do this earlier?)

c*      flush(27)
      rewind 27
      write (27) r2sum,sclerr,nlive,chi2
      i1=nlive*(nlive+1)/2
      call tapdrv(27,1,nlive,i1,a,2)
      call tapdrv(27,1,nlive,nlive,sigma,1)
      write (27) (free(j),j=1,ntpart)
      call tapdrv(27,1,ntpart,ntpart,adjust,1)
      write (27) (idxb(J),j=1,mbias)    
c      print *,'LC_SOLUTION write 27 L1/L2 soln nlive ntpart '
c     .   ,' adjust 1 210 244 245 '
c     .     ,nlive,adjust(1),adjust(210),adjust(244),adjust(245)


c  Remove the bias adjustments and set the pointers to indicate WL biases are fixed                                

      
      i1=lpart+1
      call zero1d(i1,ntpart,adjust)  
      call wlmode(free,1)


c  Halve the observation counters--total and by satellite
      nobs = nobs/2
      do i = 1,nsat
         iseen(i) = iseen(i)/2
      enddo

c  Read back the WL free normal equations

      nded = ntpart-nlive
      call lcnorm(ierinv,1)
      call qhead1(2)
c      print *,'LC_SOLUTION read 28 LC NEs r2sum  ntpart,nded,nlive: '
c     .    ,r2sum,ntpart,nded,nlive
c     write(6,'(5i20,/)') (free(i),i=1,ntpart)
c     write(6,'(5d20.10,/)') (a(i),i=1,ntpart)  

c  Solve the LC normal equations

      call report_stat('STATUS','SOLVE','lc_solution',' '
     .      , 'Solving LC normal equations after L1/L2 separate',0)
      call lcnorm(ierinv,2)
c      if( ierinv.ne.0 )  call report_stat('FATAL','SOLVE','lc_solution'
c     .  ,' ','Bad inversion of LC normal equations after L1/L2 solution'
c     .  ,0)
      call solvlc(adjust,nded,rnew,zeroad)      
c     in the case of (x2 .ne. 0), we have to calculate more
      if ( .not. zeroad ) then
         ifast = 0
         call solve3(adjust,nded,rnew,ifast) 
      endif

c  Compute the postfit goodness of fit

      goodft = chi2/dble(nobs-nlive)
      goodft = dabs(goodft)
      sclerr = dsqrt(goodft)    

c  Recompute the sigmas for saving on unit 29 (they aren't actually
c  needed before they are recalculated from the diagonal elements of 'a'
c  in LSQERR, but consistency helps the debugging
             
      do i=1,nlive
        idiag=i*(i+1)/2
        idiag = jelf(i+1)
        sigma(i) = dsign(dsqrt(dabs(a(idiag))),a(idiag))
      enddo
      
c-- Save the LC-mode bias-free normal equations on unit 29

* MOD TAH 030627: Uncommented flush call to see if this fixes
*     endfile: truncation failed in endfile
*     apparent state: unit 29 named tmp_net2.29
*     error 
c      call flush(29)
* MOD TAH 030628: Adding flush did not work, so closed and
*     reopened both units 27 and 29.  Note: Need to regenerate
*     file names here.
      close(29,iostat=ioerr) 
      open(unit=29,file=ftmp29, form='unformatted',status='unknown'
     .     ,iostat=ioerr)  
      if( ioerr.ne.0 ) call report_stat('FATAL','SOLVE','lc_solution'
     .  ,ftmp29,'Error re-opening scratch unit',ioerr)   
c      REWIND 29
      WRITE (29) R2SUM,SCLERR,NLIVE,CHI2,MSIG
      i1=nlive*(nlive+1)/2
      call tapdrv(29,1,nlive,i1,a,2) 
c** RWK 040203: The sigmas written to unit 29 here are the sigmas from the L1/L2
c        solution, not LC; i.e. they are not consistent with the diagonal
c        elements of 'a'.  However, they are not used again before
c        they are recalculated from 'a' in LSQERR for printout. 
      call tapdrv(29,1,nlive,nlive,sigma,1)
      WRITE (29) (FREE(J),J = 1,NTPART)
      call tapdrv(29,1,ntpart,ntpart,b,1)
      call tapdrv(29,1,ntpart,ntpart,borg,1)
      call tapdrv(29,1,ntpart,ntpart,adjust,1)
      WRITE (29) (IDXB(J),J = 1,MBIAS)   
c      print *,'LC_SOLUTION wrte 29 ntpart nlive adjust 1 210 244 245'
c     .     ,ntpart,nlive,adjust(1),adjust(210),adjust(244),adjust(245)  
c      print *,'a 1 22155 ',a(1),a(22155)

      call report_stat('STATUS','SOLVE','lc_solution',' '
     .          , 'LC solution complete',0)

      return
      end



