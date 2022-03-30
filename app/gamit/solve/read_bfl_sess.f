      Subroutine READ_BFL_SESS
c
c     read session part of batch file with key word style
c            Dong 930716

      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'
      include 'parameters.h'

      character*256 message
      character*120 wcmd
      integer match_name,lift_arg,count_arg
      character*4 snam,upperc,sites(maxsit)
      character*16 lowerc,code
      character*3 alabel
      character*9 scmd
      character*21 form1
c      character*14 form2 
      character*16 form2
      character*1  suse(maxsit)
c
      integer*4 issn,ioerr
      integer lcmd,i,ic,ib1,j,ib,k,iii

      real*8 error1,error2

c     find session commands
                              
c     this now dummy since multi-session not supported
      issn = 1
      write(scmd,'(a8,i1)') 'session_',issn
      call getcmd(5,scmd,wcmd,lcmd,2)       

c     site and satellite selections

c     default: use all sites and satellites
      do i = 1,nsite
         sites(i) = cfiln(i)(2:5)
         call uppers(sites(i))
         iusest(i) = 1
      enddo
      do i = 1,nsat
         iusesa(i) = 1
      enddo
c     screen output
      if( logprt ) write (6,'(/,'' Site and satellite selections:'')')
c     select sites and satellites
      call getcmd(5,'includ',wcmd,lcmd,3)    
      if (lcmd.le.0) goto 30
      ic = count_arg(wcmd)
      if (ic.le.0) goto 30
c     read site or satellite name and match their index
      do 10 k = 1,ic
         ib = lift_arg(wcmd,code,k)
         if (ib.le.0) goto 10
c        all sites
         if (lowerc(code(1:7)).eq.'all_sit') then
            do i = 1,nsite
               iusest(i) = 1
            enddo
            goto 10
         else
            snam = upperc(code(1:4))
            ib1 = match_name(nsite,4,sites,snam)
            if (ib1.gt.0.and.ib1.le.nsite) iusest(ib1) = 1
         endif
c        all satellites
         if (lowerc(code(1:7)).eq.'all_sat') then
            do i = 1,nsat
               iusesa(i) = 1
            enddo
            goto 10
         else
            snam = upperc(code(1:4))
            if (snam(1:2).eq.'PN') then
               read (snam,'(2x,i2)') ib1
               do j = 1,nsat
                  if (ib1.eq.isprn(j)) iusesa(j) = 1
               enddo
            endif
         endif
 10   continue

c     remove sites and satellites
 30   call getcmd(5,scmd,wcmd,lcmd,2)
      call getcmd(5,'exclude',wcmd,lcmd,3)   
      if (lcmd.le.0) goto 50
      ic = count_arg(wcmd)
      if (ic.le.0) goto 50
c     read site or satellite name and match their index
      do 40 k = 1,ic
         ib = lift_arg(wcmd,code,k)
         if (ib.le.0) goto 40
c        all sites
         if (lowerc(code(1:7)).eq.'all_sit') then
            do i = 1,nsite
               iusest(i) = 0
            enddo
            goto 40
         else
            snam = upperc(code(1:4))
            ib1 = match_name(nsite,4,sites,snam)
            if (ib1.gt.0.and.ib1.le.nsite) iusest(ib1) = 0
         endif
c        all satellites
         if (lowerc(code(1:7)).eq.'all_sat') then
            do i = 1,nsat
               iusesa(i) = 0
            enddo
            goto 40
         else
            snam = upperc(code(1:4))
            if (snam(1:2).eq.'PN') then
               read (snam,'(2x,i2)') ib1
               do j = 1,nsat
                  if (ib1.eq.isprn(j)) iusesa(j) = 0
               enddo
            endif
         endif
 40   continue

 50   continue   
      write (alabel,'(i3)') maxsit
      form1 = '(/,a,/,1X,' //alabel// '(2x,i3))'
      form2 = '(1X,' //alabel// '(3x,a1))  '
      do i = 1,nsite
         suse(i) = 'Y'
         if (iusest(i).le.0) suse(i) = 'N'
      enddo
      write(10,form1) ' Stations used',(iii,iii=1,nsite)
      write(10,form2) (suse(iii),iii=1,nsite)
      do i = 1,nsat
         suse(i) = 'Y'
         if (iusesa(i).le.0) suse(i) = 'N'
      enddo     
      
      write(10,'(/,a)') ' Satellites used ( Channel / PRN )'
      write(10,'(1x,32(2x,i2))') (i,i=1,nsat)
      write(10,'(1x,32(2x,i2))') (isprn(i),i=1,nsat)
      write(10,'(1x,32(3x,a1))') (suse(i),i=1,nsat)
             
c     data weighting
      error1 = 0.d0
      error2 = 0.d0   
      call getcmd(5,scmd,wcmd,lcmd,2)
      call getcmd(5,'error',wcmd,lcmd,3)  
      if (lcmd.le.0) then 
c        no argument for 'error model' means new-style input: call get_err_apr
         call get_err_apr( scmd,error1,error2,1 )
      else                            
c        non-blank argument for 'error model' means old-style input; no longer supported
         call report_stat('FATAL','SOLVE','read_bfl_sess',' ',
     .    'None or old-style data wgts input: rerun FIXDRV',0)
      endif   
c     if a noise-file has been input, read it and override the solve batch-file values
      if( nfiln(1:1).ne.' ') then  
         call get_err_apr( scmd,error1,error2,2 )
      endif
              
c------ Set weighting factors

      aphi=dble(error1)
      bphi=dble(error2) 
c      print *,'aphi bphi ',aphi,bphi
c convert mm to L1 cycles
      aphi=aphi/190.29367d0 
c     leave distance-proportional term dimensionaless (convert from ppm)
      bphi = bphi*1.d-6
      do i = 1, nsite
            sit_err(i) = sit_err(i)/190.29367d0
            sit_elv(i) = sit_elv(i)/190.29367d0
         enddo
         do i = 1,nsat
            sat_err(i) = sat_err(i)/190.29367d0
         enddo
c distance proportional weighting
      iwght=1
c uniform weighting
      if(error2.eq.0.) iwght=0
c      print *,'READ_BFL_SESS sit_err sit_elv ',sit_err(1),sit_elv(1)

c     ionospheric and atmospheric constraints
           
      call getcmd(5,scmd,wcmd,lcmd,2)  
      if( l2flag.ne.4 ) then               
         akappa = 0.
         bkappa = 8.0
         call getcmd(5,'ionos',wcmd,lcmd,3)
         if (lcmd.gt.0) then
           read(wcmd(1:lcmd),*,iostat=ioerr) akappa,bkappa
           if(ioerr.ne.0) call report_stat('FATAL','SOLVE'
     .       ,'read_bfl_sess'
     .       ,' ','Error reading ionsopheric constraints',ioerr)
         endif  
         if( logprt ) write(6,75) akappa,bkappa
         write(10,75) akappa,bkappa
 75      format(/,' Assumed ionosphere error'/,
     1           ' constant : ',f7.1,'  ppm : ',f7.2)
c        convert mm to L1 cycles
         akappa = akappa/190.29367d0     
c        leave proportional term unitless (convert from ppm)
         bkappa = bkappa*1.d-6     
c         print *,'  akappa bkappa ',akappa,bkappa
      endif
      iatcon = 2
      call getcmd(5,scmd,wcmd,lcmd,2)
      call getcmd(5,'atmos',wcmd,lcmd,3)
      if (lcmd.gt.0) then
        if (wcmd(1:1).eq.'Y'.or.wcmd(1:1).eq.'y') iatcon = 1
      endif
      if(iatcon.eq.1) then
        if(logprt) write(6,'(/,a,/)') ' Tropospheric constraint applied'
        write(10,'(/,a,/)') ' Tropospheric constraint applied'
       endif
c      Note:  current code does not allow the tropospheric spatial
c             constraint to be used with multiple zenith delays.
c             If numzen>1, iatcon is reset in LSQUAR.

c     wide-lane bias-fixing criteria

      call getcmd(5,scmd,wcmd,lcmd,2)
      call getcmd(5,'wide',wcmd,lcmd,3)
      if (lcmd.le.0) then
         wldev = 0.15
         wlsig = 0.15
         wlcut = 1000.0
         wlrat = 10.0
         wldmax = 500.
      else
         read (wcmd,*) wldev,wlsig,wlcut,wlrat,wldmax
         if( wldmax.le.0.d0 ) then
            wldmax = 500.d0
            write(message,'(2a,1x,f4.0,a)') 'Batch file missing value'
     .          ,' for maximum BL length for WL, using ',wldmax,' km'
     .    ,' maximum BL length for WL, using ',wldmax,' km'
            call report_stat('WARNING','SOLVE','read_bfl_sess',' '
     .                      , message,0)
         endif
      endif
      call getcmd(5,scmd,wcmd,lcmd,2)
      call getcmd(5,'pseudorange',wcmd,lcmd,3)
      if (lcmd.le.0) then
         prdev = 0.15
         prsig = 0.15
         prcut = 1000.0
      else
         read (wcmd,*) prdev,prsig,prcut
      endif

      return
      end
