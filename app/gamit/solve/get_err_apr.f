c Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
c
      Subroutine GET_ERR_APR( scmd,error1,error2,job )
         
c     Read a solve batch file to get the a priori data weights; called by read_bfl.  S. McClusky  August 1995

c     Now called by read_bfl_sess rather than read_bfl to allow session-by-session weights,
c     and called twice--first time to read main solve batch file (job=1), then a second time
c     to read station- and satellite-dependent values written into the N-file (same format as 
c     batch file) by autcln (or manually) (job=2). -- R. King July 1998 / January 2007

      implicit none 

      include '../includes/dimpar.h'
      include 'solve.h'
      include 'parameters.h'

      logical bflag,print_flag

      character*4 snam,upperc,buf4
      character*9 scmd
      character*10 tmp_model
      character*16 code
      character*120 wcmd
      character*256 message 

      integer i,j,i1,i0,ib,ic,type,ioerr,in,job
      integer lcmd,count_arg,lift_arg,mchkey
      integer indorb,ival,indx

      real*8 tmp_error,tmp_elv,error1,error2
                    
c     if this is the second call, the N-file name will be non-blank; in this case
c     open the N-file and don't reinitialize the data arrays
c
c       print_flag :  true if first call and no second call to be made (no N-file)
c                          or if second call with N-file
c                     false if first call with second call to be made
  
                                          
      if( job.eq.1 ) then 
c       set input to main batch file and initialise data arrays
        in = 5
        call zero1d(1,maxsit,sit_err) 
        do i=1,maxsit
          sit_err(i) = 10.d0
          err_mod(i) = 'uniform'
        enddo
        call zero1d(1,maxsit,sit_elv)
        call zero1d(1,maxsat,sat_err)   
        print_flag = .true.
c       see if N-file to be read in second call
        call getcmd(in,scmd,wcmd,lcmd,2)       
        call getcmd(in,'noise',wcmd,lcmd,3)
        if( lcmd.gt.0 ) then
           nfiln = wcmd(1:lcmd)
           print_flag = .false.
        endif
      elseif (job.eq.2 ) then 
        in = 13
        print_flag = .true.   
      else                
        print_flag = .true.
        call report_stat('FATAL','SOLVE','get_err_apr',' ',
     .   'Illegal value of ijob',0)
      endif   
    
      bflag = .false.
    
c     Get a priori data error for sites

      call getcmd(in,scmd,wcmd,lcmd,2)  
      call getcmd(in,'error',wcmd,lcmd,3)   
      if( lcmd.lt.0 ) then
        call report_stat('WARNING','SOLVE','get_err_apr',' ',
     .'A priori measurement error not found, default used, rerun FIXDRV'
     .     ,0) 
      endif
      do 150 i0 = 1,1000
         call getcmd(in,'stn_error',wcmd,lcmd,3)
         if (lcmd.le.0) goto 155
c        decompose the command line
         ic = count_arg(wcmd)
c        not enough arguments
         if (ic.le.1) goto 150
c        pointer to
         ib = lift_arg(wcmd,code,1)
         if (ib.le.0) goto 150
         type = 1
         if (upperc(code(1:4)).eq.'ALL_') type = 2
         if (type.eq.1) snam = upperc(code(1:4))
c PT960102: remove free format character read
         indx = ib+1
         call read_line(wcmd,indx,'CH',ioerr,ival,tmp_model)

c PT960129: If there is no elevation specified then don't read it but
c           set it to zero (to avoid DEC problems!)
         if(ic.eq.3)then
           tmp_elv = 0.d0
           read (wcmd(indx:),*,iostat=ioerr)tmp_error
         else
           read (wcmd(indx:),*,iostat=ioerr)tmp_error,tmp_elv
         endif

         if( ioerr .gt. 0 ) then  
           call report_stat('FATAL','SOLVE','get_err_apr',' '
     .       ,'Error reading error model values',0)
         endif
         do 140 i=1,nsite
            i1=3*(i-1)+1
            if (type.eq.1) then
               buf4=upperc(rlabel(i1)(1:4))
               if (buf4.ne.snam) goto 140
            endif
            err_mod(i)=tmp_model
            sit_err(i)=tmp_error
            sit_elv(i)=tmp_elv
            if ( err_mod(i) .eq. 'baseline' ) then
              error1 = sit_err(i)
              error2 = sit_elv(i)
              bflag = .true.
            endif
c  Cannot mix with uniform and elevation error models with baseline model,
c  or specifiy individual station weights with the baseline model.
            if (type .eq. 1 .and. (err_mod(i) .eq. 'baseline'
     .          .or. bflag )) then    
                write(message,'(2a)') 
     .             'Can only use baseline error model in all_sites mode'
     .          ,', not with site-dependent elevation or uniform models'
                call report_stat('FATAL','SOLVE','get_err_apr',' '
     .                          ,message,0)
            endif
  140    continue
  150 continue
  155 if ( print_flag ) then
        if ( err_mod(1) .ne. 'baseline') then
          if( logprt ) write(6,160)
          write(10,160)
  160     format(/,4x,'A priori receiver measurement error models and',
     .' std devs in mm',/
     .  ,'Station                     Model       Std dev    Elev   ',/)
        else
          if( logprt ) write(6,161)
          write(10,161)
  161     format(/,4x,'A priori receiver measurement error models and',
     .' std devs in mm',/
     .  ,'Station                     Model       Std dev    ppm    ',/)
        endif
        do 165 i=1,nsite    
          if( logprt ) write( 6,162) i,rlabel((i-1)*3+1)(1:4),sitnam(i)
     .              ,err_mod(i),sit_err(i), sit_elv(i)
          write(10,162) i,rlabel((i-1)*3+1)(1:4),sitnam(i),err_mod(i)
     .                 ,sit_err(i), sit_elv(i)
  162     format(i3,2x,a4,1x,a12,5x,a10,1x,2(f7.2,3x))
c         SOLVE will fail unpredictably if data noise is zero, so check
          if(sit_err(i).eq.0.d0 .and. sit_elv(i).eq.0.d0 ) then
            write(message,163) i
  163       format('a priori error for site ',i2
     .              ,' is zero')
            write(10,'(a)') message
           call report_stat('FATAL','SOLVE','get_err_apr',' ',message,0)
          endif
  165   continue
      endif

c     Get a priori for satellite data error 

c       find the index of the first orbit parameter
      indorb = 0
      do  i=1,ntpart
c         this should work just testing on islot1 = 501, I think --rwk
          if((islot1(i).gt.500).and.(islot1(i).le.2400)) then
            indorb = i
            go to 170   
          endif
      enddo     
  170 call getcmd(in,scmd,wcmd,lcmd,2)
      call getcmd(in,'error',wcmd,lcmd,3)
      do 250 i0 = 1,1000
         call getcmd(in,'sat_error',wcmd,lcmd,3)
         if (lcmd.le.0) goto 255
c        decompose command line
         ic = count_arg(wcmd)
c        not enough arguments
         if (ic.le.1) goto 250
c        pointer to
         ib = lift_arg(wcmd,code,1)
         if (ib.le.0) goto 250
         type = 1
         if (upperc(code(1:4)).eq.'ALL_') type = 2
         if (type.eq.1) snam = upperc(code(1:4))  
         read (wcmd(ib+1:lcmd),*) tmp_error
         if( tmp_error.ne.0.d0 .and. indorb.eq.0 ) then
           call report_stat('WARNING','SOLVE','get_err_apr'
     .   ,' ','Code bug: cannot get sat meas error in BASELINE mode',0)
           do i=1,nsat
             sat_err(i) = 0.d0
           enddo
         else
           do 240 i=1,nsat
              i1=norb*(i-1)+indorb
              if (type.eq.1) then
                 j = mchkey(rlabel(i1),snam,20,4)
                 if (j.le.0) goto 240
              endif
              sat_err(i) = tmp_error
 240       continue  
         endif
 250  continue
 255  if( print_flag ) then
       if( logprt ) write(6,260)
       write(10,260)
 260   format(/,
     1    4x,'A priori satellite measurement error std devs in mm',/,
     2    'Satellite           Std dev   ',/)
           do 265 i=1,nsat
           if( logprt ) write( 6,262) i,isprn(i),sat_err(i)
           write(10,262) i,isprn(i),sat_err(i)
  262      format(i3,2x,'PRN ',i3,5x,f7.2)
  265      continue
      endif
        
      return
      end
