Copyright 1996 Massachusetts Institute of Technology and The Regents of the
c              University of California, San Diego.  All Rights Reserved.

       Subroutine READ_BFL( icall,constraints,free_fix )

c     Read batch file with key word style
c     Written by Dong Jul 93; modified by King and McClusky Aug 93-Jun 94

      implicit none
      
      include '../includes/dimpar.h'
      include 'solve.h'
      include 'parameters.h'

      integer*4 lift_arg,count_arg,icall,ilen,ioerr
     .        , nblen,indx,norbm,indorb,ntzen,ntgrad,indzen,indgrad,iyr2 
     .        , nskip,julday,jdstart,jdref,iyrref,idoyref,monref,idayref
     .        , ihrref,minref,ntfile,nrfile,idble,isngle

      integer*4 ii  ! ! Loop counter for DEBUG

c     integer*4 mchkey

      real*8 shft(3),shftdot(3),grot(3),grotdot(3),scale,scaledot,gepoch
     .     , temp,fixdrv_vers
     .     , stawgts(maxsit),satwgts(maxsat)
     .     , sodstart,sodref,delt2,refyr,decyrs,secref

      character*1 upperc
      character*2 ayr
      character*4 free_fix
      character*5 constraints,datum
      character*6 keyw
      character*9 stone,code
      character*16 lowerc,uname,tfilnm(maxtfl)
      character*19 dattim 
      character*57 soln_epoch
      character*120 wcmd
      character*256 message
                        
      logical fcheck
      dimension temp(maxprm)

      integer lcmd,i,ic,ib1

      logical debug/.false./

c      icall= 1: read controls for initial solution
c      icall= 2: read file-update options after solution

      goto (10,300) icall

c------ Initial solution - read all input controls   (icall = 1)

c------ Read the FIXDRV version number  
       
 10   call getcmd(5,'FIXDRV version',wcmd,lcmd,1)
      if (lcmd.le.0 ) then
        fixdrv_vers = 0.0
      else
        read(wcmd,*) fixdrv_vers
      endif


c------ Read the owner

      call getcmd(5,'owner',wcmd,lcmd,1)
      if (lcmd.le.0) then
         call report_stat('FATAL','SOLVE','read_bfl',' '
     .    ,'Missing Owner name in batch file',0)
      else
         owner = wcmd(1:lcmd)
      endif  

              
c------ Read the scratch file directory

       call getcmd(5,'scratch',wcmd,lcmd,1)
       if( lcmd.gt.0 ) then
          scratch_dir = wcmd(1:lcmd) 
        else
          scratch_dir = ' '
        endif

c-----  Read the log-file option

      call getcmd(5,'log',wcmd,lcmd,1)
      if( lcmd.le.0 ) then   
c       default is print on for backward compatibility
        logprt = .true.
      else   
        if( wcmd(1:1).eq.'N' ) then
          logprt = .false.
        else
          logprt = .true.
        endif
      endif   

c------ Read the control for skipping the loose solutions

      call getcmd(5,'skip',wcmd,lcmd,1)     
      if( lcmd.le.0 ) then
c       default is to do the loose solutions
        do_loose = .true.
      else   
        if( wcmd(1:1).eq.'Y' ) then
          do_loose = .false.
        else
          do_loose = .true.
        endif
      endif   

c ----- Read the bias-weighting option (experimental)

      call getcmd(5,'bias_apr',wcmd,lcmd,1)
      if( lcmd.le.0 ) then
        bias_apr = 1000.d0
      else
        read(wcmd,*) bias_apr
      endif          

c ----- Read the r-condition number for removing a dependent bias (experimental)

      call getcmd(5,'bias_rcond',wcmd,lcmd,1)
      if( lcmd.le.0 ) then 
c*** set the default rcond (ier 104 in inversion) high so that too many dependent
c    biases are removed (i.e., an independent bias not estimated)
        bias_rcond = 1.d4
      else
        read(wcmd,*) bias_rcond
      endif          
                     
c------ Read the debug flag for bias-fixing  (change later to add levels?)
                     
      call getcmd(5,'bias_debug',wcmd,lcmd,1)
      if( lcmd.le.0 ) then   
        bias_debug = .false.
      else   
        if( wcmd(1:1).eq.'Y' ) then
          bias_debug = .true.
        else
          bias_debug = .false.
        endif
      endif     

c------ Read the Q-file name and open the q-file and o-file

 50   call getcmd(5,'Q-file',wcmd,lcmd,1)
      if (lcmd.le.0) then  
           call report_stat('FATAL','SOLVE','read_bfl',' '
     .    ,'Missing Q-file name in batch file',0)
      else
         qfiln = wcmd(1:lcmd)
      endif              
      if( logprt ) write (6,60) qfiln
 60   format (/,' Q-file name: ',2x,a)
      open (unit=10,file=lowerc(qfiln),form='formatted',err=50,
     .      status='unknown')     
      call newfil(qfiln,ofiln,'O')
      call lowers (ofiln)
      open(unit=15,file=ofiln,status='unknown')


c------ M-file name

      call getcmd(5,'M-file',wcmd,lcmd,1)
      if (lcmd.le.0) then
         call report_stat('FATAL','SOLVE','read_bfl',' '
     .    ,'Missing M-file name in batch file',0)
      else
         mfiln = wcmd(1:lcmd)
      endif
      if( logprt ) write (6,'(/,a,2x,a)') ' M-file name: ',mfiln
      if (fcheck(mfiln)) then
         call mopens (mfiln,'old',11,ioerr)
      else
         call report_stat('FATAL','SOLVE','read_bfl',mfiln
     .                      ,'Could not find M-file:',0)
      endif
      if (ioerr.ne.0) then
         call report_stat('FATAL','SOLVE','read_bfl',mfiln
     .                   ,'Could not find M-file:',ioerr)
      endif     


c------ Read the m-file headers to get times and parameter flags    
c       --subsequent records for each c-file read in gethed

c       ntpart (stored in common/block2/ is the total number of parameters
c         (all sites, clocks, zenith delays, gradients, orbital parameters,
c         and biases). It is set initially here from the value on record 1 
c         of the m-file and corresponds to the non-zero elements of islot1 
c         (see cfmrg/fills1.f), including one-way L1 and L2-L1 biases.  It 
c         is redefined in lsquar after independent double difference biases 
c         are formed. 
c       npartm(nsite) (stored in common/parts/) is the number of parameter
c         partials plus one-way biases on each c-file and corresponds to the 
c         non-zero  elements of islot2 (see cfmrg/fills2.f).  It is read from
c         record 3 of the m-file. In the current code, it is assumed to be the 
c         same for every c-file.  
c       npartc (stored in common/parts/) is the number of partial types (no
c         biases) on the c-file.  Its value is the same for all c-file record 5s, 
c         but the partial in the tpart array may be zero if the satellite is not 
c         present at a particular epoch.  
             
      call readm1 (11,
     .             nepoch,ntpart,
     .             rlabel,idms,islot1,preval,temp,
     .             nsat,isprn,
     .             nsite,sitnam,
     .             nrfile,obfiln,
     .             ntfile,tfilnm,norbm)  
c     only one t-file now that multi-session is not supported
      tfiln = tfilnm(1)
      if(debug) then   
* MOD TAH 190609: Only ouput first tfile names used with only one sessiom.  Long
*     term mod is to set (MAXTFL=1 in dimpar.h simce only one sesssion.
         print *,'READ_BFL ntfile tfilnm ',ntfile,tfilnm(1:ntfile)
         print *,'READ_BFL tfiln ',tfiln
         print *,'READ_BFL ntpart nrfile rlabel ',ntpart,nrfile
         print *,'READ_BFL obfiln(1) ',obfiln(1)
         do i=1,ntpart,100
            print *,'READ_BFL rlabel ',i,rlabel(i) 
         enddo
      endif      
      call readm2 (11,
     .             it0,t00,
     .             nepoch,inter,nskip,
     .             nsite,cfiln,stawgts,
     .             nsat,satwgts )     
         if( nsite.gt.maxsit .or. nsat.gt.maxsat ) then
             write(message,'(2a,4i4)') 'Dimension mismatch on M-file'
     .               ,'--nsite, maxsit,nsat,maxsat:'
     .               ,   nsite,maxsit,nsat,maxsat
             call report_stat('FATAL','SOLVE','read_bfl',' ',message,0)
          else
             call report_stat('STATUS','SOLVE','read_bfl',' '
     .                       , 'Read records 1 and 2 of the M-file',0)
         endif    


c----- Read N-file name and open the file
                     
      nfiln = ' '
      call getcmd(5,'session',wcmd,lcmd,2)       
      call getcmd(5,'noise',wcmd,lcmd,3)
      if( lcmd.gt.0 ) then
         nfiln = wcmd(1:lcmd)
      endif
      if( nfiln(1:1).ne.' ') then
        if( fcheck(nfiln) ) then
          open(unit=13,file=nfiln,iostat=ioerr,form='formatted'
     .               ,status='old') 
          if( ioerr.ne.0 )  call report_stat('FATAL','SOLVE'
     .            ,'read_bfl',nfiln,'Error opening N-file',ioerr)
        else
          call report_stat('FATAL','SOLVE','read_bfl',nfiln
     .            ,'Cannot find specified N-file',ioerr)
        endif  
      endif


c------ H-file mode

      call getcmd(5,'H-file',wcmd,lcmd,1)
      if (lcmd.le.0) then
         ihmode = 0
      else
         read(wcmd,*) ihmode
      endif
      gloprt = .true.
c     construct M-file and H-file names
      ilen = nblen(qfiln)
      mfiln= '                '
      hfiln= '                '
      mfiln(1:ilen)=qfiln(1:ilen)
      mfiln(1:1)='m'        
c*    Lengthen name of h-file (only) to include 2-digit year prior to day: hxxxxa.yyddd
c*    Interim implementation is to read the m-file here (above),
c*    but longer-term we will change the q- and m-file names too,
c*    so can get the full h-file name from these.
      hfiln(1:ilen)=qfiln(1:ilen)
      hfiln(1:1) = 'h'
      iyr2 = mod(it0(3),100)
      write(ayr,'(i2.2)') iyr2
      hfiln(8:9) = ayr
      hfiln(10:12) = qfiln(8:10) 
      if (gloprt .and.  
c        standard mode: write h-file only for final solutions
     .  ( (ihmode.eq.0.and.hfiln(6:6).eq.'a')
c        second mode: write h-file for all  solutions
     .   .or. ihmode.ge.1 ) )  then
        open(unit=21,file=hfiln,form='formatted',status='unknown')
c**  The following needed for DEC but prohibited on Sun:
c**     .       , recl=250)   
      endif


c------ Write the program header to the Q-file and H-file
                                   

      call lsqprt(vers,owner)


c------ Select the datum and get its parameters for the baseline computations
c       (coordinates themselves are now always spherical geocentric)

      call getcmd(5,'datum',wcmd,lcmd,1)
      if (lcmd.le.0) then
         idatum = 0   
         datum = '     ' 
      else
         read(wcmd,*) idatum 
         if( idatum.eq.0 ) then
           datum = '     '  
         else
           call report_stat('FATAL','SOLVE','read_bfl',' '
     .        ,'Geodetic coordinates no longer supported',0)
         endif
      endif
      call read_gdatum( 20,datum,semi,finv,shft,shftdot
     .                , grot,grotdot,scale,scaledot,gepoch)  
c     convert the semi-axis to kilometers
      semi = semi*1.d-3

c------ phase difference options

      call getcmd(5,'differen',wcmd,lcmd,1)
      if (lcmd.le.0) then
         call report_stat('FATAL','SOLVE','read_bfl',' '
     .                   ,'Missing differencing mode',0)
      else
         if (upperc(wcmd(1:1)).eq.'S') then
            isngle = 1
            idble = 0
            indx = 1
            call report_stat('FATAL','SOLVE','read_bfl',' '
     .                   ,'Single-differences not supported',0)
         else
            isngle = 0
            idble = 1
            indx = 2
         endif
      endif
      if( logprt )
     .    write(6,'(/,'' Phase differencing mode: '',a)') wcmd(1:lcmd)
 

c------ Epochs interval

      call getcmd(5,'epoch',wcmd,lcmd,1)
      if (lcmd.le.0) then
         istart = 1 
c        set default to 24hrs at 30s interval, but maybe ok to make this very large?
c        (was 225 !!) until 980217 --rwk
         iend = 2880
         idecim = 1
      else
         ic = count_arg(wcmd)
         if( ic.eq.2 ) then
           read(wcmd,*) istart,iend
c          set decimation to 1 to maintain backwards compatibility
           idecim = 1
         elseif ( ic.eq.3 ) then
           read(wcmd,*) istart,iend,idecim
c          set decimation to 1 to maintain backwards compatibility
           if ( idecim .eq. 0 ) idecim = 1
         else
            call report_stat('FATAL','SOLVE','read_bfl',' '
     .                      ,'Bad entry for epochs line',0)
         endif
      endif      

c------- Write session info to q-file and o-file
              
c     get the midpoint epoch and write it out in exactly the same form as
c     in GLOBK print files, to facilitate extraction of information for plotting  
      jdstart =julday( it0(1),it0(2),it0(3) ) 
      sodstart=  3600.d0*t00(1) + 60.d0*t00(2) + t00(3)
      jdref = jdstart
      sodref = sodstart       
      delt2 = ( dfloat(nepoch)/2.d0 ) * inter
      call timinc(jdref,sodref, delt2) 
      call dayjul(jdref,iyrref,idoyref)        
      call monday(idoyref,monref,idayref,iyrref)  
      refyr = decyrs( iyrref,idoyref,sodref)
      call ds2hms(iyrref,idoyref,sodref,ihrref,minref,secref)      
c     create an epoch line matching exactly the GLOBK lines for grep'ing
      write(soln_epoch,'(a,i4,a,i2,a,i2,i3,a,i2,a,f9.4,a)') 
     .   " Solution refers to     : ",iyrref,"/",monref,"/",idayref
     .     ,ihrref,":",minref,"    (",refyr,")"
      write (10,95)  soln_epoch,istart,iend,idecim  
      write (15,95)  soln_epoch,istart,iend,idecim
      if( logprt ) write (6,95) soln_epoch,istart,iend,idecim
  95  format (/,a,//,' Epoch interval: ',i5,1x,' - ',i5,//,
     .          ' Decimation interval: ',i3)
                                   
c------ Combination mode

      call getcmd(5,'combin',wcmd,lcmd,1)
      if (lcmd.le.0) then
         l2flag = 1
      else
         do i = 1,6
            keyw(i:i) = upperc(wcmd(i:i))
         enddo
         if (keyw.eq.'L2_ONL') then
           l2flag =-2
         elseif (keyw.eq.'L1_ONL') then
           l2flag =-1
         elseif (keyw.eq.'L1_REC') then
           call report_stat('FATAL','SOLVE','read_bfl',' '
     .         ,'L1_RECEIVER option no longer supported: use L1_ONLY',0)
c           l2flag = 0
         elseif (keyw.eq.'LC_ONL' .or. keyw.eq.'LC    ') then
           l2flag = 1
         elseif (keyw.eq.'L1&L2 ') then
           l2flag = 2
         elseif (keyw.eq.'LC_HEL') then
           l2flag = 3
         elseif (keyw.eq.'LC_AUT') then
           l2flag = 4
         elseif (keyw.eq.'L1,L2_') then
           l2flag = 5
         else
           write(message,'(a,a6,a)') 'Combination keyword = ',keyw
     .                          , ' not recognized'
           call report_stat('FATAL','SOLVE','read_bfl',' ',message,0)
         endif
      endif
      if(logprt) write (6,'(/,'' Combination mode: '',a)') wcmd(1:lcmd)


c------ Bias-search approach

            call getcmd(5,'approa',wcmd,lcmd,2)
            if (lcmd.le.0) then
               noptin = 6
            else   
               call report_stat('FATAL','SOLVE','read_bfl',' ',
     .'Only decision-function approach to bias searching now allowed',0)
c              ib1 = mchkey(wcmd,'neares',lcmd,6)
c              if (ib1.gt.0) noptin = 1
c              ib1 = mchkey(wcmd,'aid',lcmd,3)
c              if (ib1.gt.0) noptin = 2
c              ib1 = mchkey(wcmd,'cumu',lcmd,4)
c              if (ib1.gt.0) noptin = 3
c              ib1 = mchkey(wcmd,'expec',lcmd,5)
c              if (ib1.gt.0) noptin = 4
c              ib1 = mchkey(wcmd,'searc',lcmd,5)
c              if (ib1.gt.0) noptin = 5
c              ib1 = mchkey(wcmd,'decis',lcmd,5)
c              if (ib1.gt.0) noptin = 6
c              ib1 = mchkey(wcmd,'prior',lcmd,5)
c              if (ib1.gt.0) noptin = 7
            endif
            if (noptin.ge.4) then
               call getcmd(5,'approa',wcmd,lcmd,2)
               call getcmd(5,'criter',wcmd,lcmd,3)
               if (lcmd.le.0) then
                  nldev = 0.15
                  nlsig = 0.15
                  nlcut = 1000.0
                  nlrat = 10.0
                  nldmax = 500.  
               else
                  read (wcmd,*) nldev,nlsig,nlcut,nlrat ,nldmax
                 if( nldmax.le.0.d0 ) then
                   nldmax = 500.d0 
                   write(message,'(2a,f4.0,a)') 'Batch file missing'
     .              ,' value for maximum BL length for NL, using '
     .              , nldmax,' km'
                   call report_stat('WARNING','SOLVE','read_bfl',' '
     .                             ,message,0)

                  endif
               endif   
               call getcmd(5,'nlscale',wcmd,lcmd,3)    
               if( wcmd(1:1).eq.'N' ) then
                 nlscale = .false.
               else
                 nlscale = .true.
               endif
            endif


c------ Solution mode 

c**RWK 061007: 'quick' solution and implicit biases no longer supported, so
c              suppress echo and trap old batch files

      call getcmd(5,'quick',wcmd,lcmd,1)
      if (lcmd.le.0) then
         lquick = 0
      else
         if (upperc(wcmd(1:1)).eq.'F') lquick=0
         if (upperc(wcmd(1:1)).eq.'Q') lquick=1
      endif
c**      write (6,'(/,'' Quick solution choice:  '',a)') wcmd(1:lcmd)
      call getcmd(5,'biases',wcmd,lcmd,1)
      if (upperc(wcmd(1:1)).eq.'I'.and.lquick.eq.0) then
        lquick=2
c**        write (6,'(/,'' Implicit biases'')')
      endif
c**061007
      if( lquick.ne.0 ) call report_stat('FATAL','SOLVE','read_bfl',' '
     .   ,'Quick solution and implicit biases no longer supported',0)


               
c------ See what partials are available on the C-files      

c      Note that ntpart, the total number of parameters, always
c      includes both L1 and L2-L1 biases, even if the receivers
c      are single frequency.  This makes the SOLVE logic easier.
c      See also cfmrg/fills1.f and solve/bcheck.f. 
      sitpar =  .false.
      satpar =  .false.
      zenpar =  .false.
      gradpar = .false.
      eoppar =  .false.
      do i = 1,ntpart
        if(  islot1(i).ge.1.and.islot1(i).le.300 ) sitpar = .true.
c RWK 190429: Upper limit changed to include more SRPs and also SVANTs. 
        if(  islot1(i).ge.501.and.islot1(i).le.2700 ) satpar = .true.
        if( (islot1(i).ge.301.and.islot1(i).le.400) ) zenpar = .true.
        if(  islot1(i).ge.24001.and.islot1(i).le.29000 ) gradpar= .true.
        if(  islot1(i).ge.80001.and.islot1(i).le.80006 ) eoppar = .true.
      enddo   

c------- Count the zenith-delay parameters according to the M-file values
c         ntzen (local) is the total number of zenith delays (all sites)
c         nzen (solve.h) is the number of zenith delays per site
        
      if( zenpar ) then
        ntzen = 0
        do i=1,ntpart  
          if( (islot1(i).gt.300.and.islot1(i).le.400) .or.
     .       (islot1(i).gt.21500).and.(islot1(i).le.24000) ) then
c               save index for first atm parameter (avg if #/site=1, multiple if #>1)
                if( ntzen.eq.0 ) indzen = i
                ntzen = ntzen + 1
           endif
        enddo  
        nzen = ntzen/nsite 
c       if nzen > 1, we need to subtract off the number of avg ZD parameters
        if( nzen.gt.1 ) nzen = nzen - 1 
c       set nzen in common /atmprms/ for rmjunk.f 
c       (checked after call to get_zen_apr if zenith delay actually estimated)
      endif   

c------- Count the gradient parameters (NS or EW) according to the M-file values

      if( gradpar ) then
        ntgrad = 0
        do i=1,ntpart
          if(islot1(i).gt.24000.and.islot1(i).le.29000) ntgrad=ntgrad+1
        enddo
        if( ntgrad.gt.0 ) gradpar = .true.
        ngrad = ntgrad/2/nsite   
c       set ngrad in common /atmprms/  rmjunk.f 
c       (checked after call to get_grad_apr if gradients actually estimated)
                     


c------- Find the index of the first orbit parameter

        do  i=1,ntpart
c           this should work just testing on islot1 = 501, I think --rwk
            if((islot1(i).gt.500).and.(islot1(i).le.2400)) then
              indorb = i
              go to 100   
            endif
         enddo   
  100   continue 
        if(debug) print *,'READ_BFL indorb ',indorb 
   
                       
c       Number of orbital parameters to be estimated  - default is M-file number

c** temporary--don't allow a mismatch between m-file and SOLVE
        norb = norbm
c       call getcmd(5,'orbit_param',wcmd,lcmd,1)
c       if (lcmd.gt.0) then
c         read(wcmd,*) norb
c       endif
c       if( norb.gt.norbm ) then
c         write(*,'(a,i2,a,i2)') 'READ_BFL: norb=',norb,' > norbm=',normb
c         stop
c       endif


c------ Convert coordinates and printout station, satellites, and analysis options
           
      call lsqint

c------ Set the 'free' array for parameter estimation

c     parameters commanded to be estimated
      call getcmd(5,'paramete',wcmd,lcmd,2)
      do  i = 1,1000
         call getcmd(5,'estima',wcmd,lcmd,3)
         if (lcmd.le.0) goto 120   
         if(debug) print *,'READ_BFL wcmd call set_par 1'
     .     ,wcmd(1:lcmd)
         call set_para(wcmd(1:lcmd),1)
      enddo     

c     parameters commanded to be fixed
 120  call getcmd(5,'paramete',wcmd,lcmd,2)
      do  i = 1,1000
         call getcmd(5,'fix',wcmd,lcmd,3)
         if (lcmd.le.0) goto 130
         if(debug) print *,'READ_BFL wcmd call set_par -1'
     .     ,wcmd(1:lcmd)
         call set_para(wcmd(1:lcmd),-1)
      enddo                                                 

c     convert overwrite protected free matrix back to 0's and 1's.
 130  do i = 1,ntpart
         if(free(i).le.0) free(i)=0
         if(free(i).gt.0) free(i)=1 
      enddo
   
c     set the group flags for estimatation  
      sitest =  .false.
      zenest =  .false.
      gradest = .false. 
      satest =  .false.
      eopest =  .false. 
      svantest = .false.
      do i = 1,ntpart
        if( islot1(i).ge.1 .and. islot1(i).le.300 .and.
     .      free(i).eq.1 ) sitest = .true. 
        if(  ( (islot1(i).ge.301 .and. islot1(i).le.400) .or.
     .       (islot1(i).ge.21501.and.islot1(i).le.24000) ) .and.
     .        free(i).eq.1 ) zenest = .true.   
        if( (islot1(i).ge.24001.and.islot1(i).le.29000) .and.
     .        free(i).eq.1 ) gradest = .true.
        if( (islot1(i).ge.501.and.islot1(i).le.2400) .and.
     .        free(i).eq.1 ) satest = .true.              
        if( (islot1(i).ge.2401.and.islot1(i).le.2700) .and.
     .        free(i).eq.1 ) svantest = .true.              
        if( (islot1(i).ge.80001.and.islot1(i).le.80006) .and.
     .        free(i).eq.1 ) eopest = .true.
      enddo
      if(debug) print *
     .  ,'READ_BFL sitest zenest gradest satest svantest eopest '
     .  ,          sitest,zenest,gradest,satest,svantest,eopest 

c------ For quick solution, fix biases always (ignore analyst mistake in batch file)
                          
c** quick no longer supported
      if( lquick.ge.1 ) then  
        call report_stat('FATAL','SOLVE','read_bfl',' ' 
     .       , ' Quick algorithm no longer suppored',0)
      endif
      nlive=0
      do  i=1,ntpart
        if(free(i).eq.1) nlive=nlive+1
      enddo      
      if(debug) then 
        print *,'READ_BFL ntpart nlive free ',ntpart,nlive 
        write(6,'("READ_BFL FREE ",100(i5,1x,i1))')  
     .                           (ii,free(ii),ii=1,ntpart)
      endif 
              

c------ Solve cutoff angle

         call get_solve_cut
              

c-----  Minimum signal-to-noise (in common /cutoff/ but not used

         minsnr = 0 

c------ A priori constraints on stations, satellites, and zenith-delays

c     Change to force a priori constraints on zenith delays, atmospheric 
c     gradients, Earth-orientation parameters, and SV antenna offsets.
c     Stations and satellites can still be left completely free in the 'tight' 
c     solution.  There are no longer separate flags for tight and loose
c     since in applying the loose constraints we only need to know if
c     tight ones need to be first removed. --rwk 961230
c     Add clock aprioris to prevent singularity w/ missing data.---rwk 061005
      sitwgt =  .false.
      satwgt =  .false.
      zenwgt =  .true.
      gradwgt = .true.
      eopwgt =  .true.  
      svantwgt = .true.      
      clkwgt = .true.
c     station coordinates      
      if (sitest ) call get_sit_apr     
c     clocks --- not coded
c     zenith delays
      if ( zenest ) then   
         call get_zen_apr( indzen )
      endif
c     atmospheric gradients   
      if( nzen.eq.1 ) then
        indgrad = indzen + nzen*nsite 
      else
        indgrad = indzen + (nzen+1)*nsite 
      endif
      if ( gradest )  then
        call get_grad_apr(indgrad) 
      endif 

c     satellite parameters              
      if( satest )   
     .   call get_sat_apr(indorb )   
c     earth orientation parameters
      if( eopest ) call get_eop_apr
      endif  
      goto 800  

c ------------------------------------------------------------------


c---File update options after solutions (repeated for tight and loose)

c       Q-file output option

 300  if ( constraints.eq.'tight' .and. free_fix.eq.'free' )
     .    stone = 'tight_fre'
      if ( constraints.eq.'tight' .and. free_fix.eq.'fixd' )
     .    stone = 'tight_fix'
      if ( constraints.eq.'loose' .and. free_fix.eq.'free' )
     .    stone = 'loose_fre'
      if ( constraints.eq.'loose' .and. free_fix.eq.'fixd' )
     .    stone = 'loose_fix'
      call getcmd(5,stone,wcmd,lcmd,2)
c     default: not print, not update
      iqflag = 2
      ioflag = 2
      do i = 1,4
         ipfil(i) = 2
      enddo

c       Print- out options

      call getcmd(5,'print',wcmd,lcmd,3)
      if (lcmd.le.0) goto 320
      ic = count_arg(wcmd)
      if (ic.le.0) goto 320
      do 310 i = 1,ic
         ib1 = lift_arg(wcmd,code,i)
         if (ib1.le.0) goto 310
         if (upperc(code(1:1)).eq.'Q') iqflag = 1
         if (upperc(code(1:1)).eq.'O') ioflag = 1
 310  continue          
c     set default threshold correlation for printing to old value for consistency
c     (though FIXDRV will use the new default of 0.99999 [usually none])
 320  call getcmd(5,stone,wcmd,lcmd,2)
      correl_prt = 0.98
      call getcmd(5,'correl_prt',wcmd,lcmd,3)    
c     get the  argument and count the number of characters 
      if (lcmd.gt.0 ) then     
        ib1 = lift_arg(wcmd,code,1) 
        if( ib1.le.0 ) call report_stat('FATAL','SOLVE','read_bfl',' '
     .    ,'Missing argument for correlation print threshold',0)
        read ( code(1:lcmd),*) correl_prt
      endif 

c        File-update options

      call getcmd(5,stone,wcmd,lcmd,2)
      call getcmd(5,'update',wcmd,lcmd,3)
      if (lcmd.le.0) goto 400
      ic = count_arg(wcmd)
      if (ic.le.0) goto 400
      do 340 i = 1,ic
         ib1 = lift_arg(wcmd,code,i)
         if (ib1.le.0) goto 340
         if (upperc(code(1:1)).eq.'M') ipfil(1) = 1
         if (upperc(code(1:1)).eq.'L') ipfil(2) = 1
         if (upperc(code(1:1)).eq.'I') ipfil(3) = 1
         if (upperc(code(1:1)).eq.'G') ipfil(4) = 1
 340  continue 
c     Adjustment tolerance (meters) for updating L-file coordinates (default 30 cm)
      call getcmd(5,stone,wcmd,lcmd,2)
      coord_upd_tol = .30
      call getcmd(5,'coord_upd_tol',wcmd,lcmd,3)
c     get the  argument and count the number of characters 
      if (lcmd.gt.0 ) then     
        ib1 = lift_arg(wcmd,code,1) 
        if( ib1.le.0 ) call report_stat('FATAL','SOLVE','read_bfl',' '
     .    ,'Missing argument for coordinate-update threshold',0)
        read ( code(1:lcmd),*) coord_upd_tol 
      endif 

      if (ipfil(1).eq.1) then
c        read the input m-file name
         call getcmd(5,stone,wcmd,lcmd,2)
         call getcmd(5,'input_m',wcmd,lcmd,3)
         if (lcmd.le.0 ) then
           call report_stat('FATAL','SOLVE','read_bfl',' '
     .                     ,'Input M-file name not found',0)
         else
           minf = wcmd(1:lcmd)
         endif
c        read the output m-file name
         call getcmd(5,stone,wcmd,lcmd,2)
         call getcmd(5,'output_m',wcmd,lcmd,3)
         if (lcmd.le.0 ) then
           call report_stat('WARNING','SOLVE','read_bfl',' '
     .          , 'Output M-file name not found; using mfile.out',0)
           moutf = 'mfile.out       '
         else
           moutf = wcmd(1:lcmd)
         endif
      endif
      if (ipfil(2).eq.1) then
c        read the input l-file name
         call getcmd(5,stone,wcmd,lcmd,2)
         call getcmd(5,'input_l',wcmd,lcmd,3)
         if (lcmd.le.0 ) then
           write(*,*) '** Input L-file name not found'
         else
           linf = wcmd(1:lcmd)
         endif
         linf = wcmd(1:lcmd)
c        read the output l-file name
         call getcmd(5,stone,wcmd,lcmd,2)
         call getcmd(5,'output_l',wcmd,lcmd,3)
         if (lcmd.le.0 ) then
           write(*,*) '** Output L-file name not found; using lfile.out'
           loutf = 'lfile.out       '
         else
           loutf = wcmd(1:lcmd)
         endif
      endif
      if (ipfil(3).eq.1) then
c        read the input i-file name
         call getcmd(5,stone,wcmd,lcmd,2)
         call getcmd(5,'input_i',wcmd,lcmd,3)
         if (lcmd.le.0 ) then
           write(*,*) '** Input I-file name not found'
         else
           iinf = wcmd(1:lcmd)
         endif
c        read the output i-file name
         call getcmd(5,stone,wcmd,lcmd,2)
         call getcmd(5,'output_i',wcmd,lcmd,3)
         if (lcmd.le.0 ) then
           write(*,*) '** Output I-file name not found; using ifile.out'
           ioutf = 'ifile.out       '
         else
           ioutf = wcmd(1:lcmd)
         endif
      endif
      if (ipfil(4).eq.1) then
c        read the input g-file name
         call getcmd(5,stone,wcmd,lcmd,2)
         call getcmd(5,'input_g',wcmd,lcmd,3)
         if (lcmd.le.0 ) then
           write(*,*) '** Input G-file name not found'
         else
           ginf = wcmd(1:lcmd)
         endif
c        read the output g-file name
         call getcmd(5,stone,wcmd,lcmd,2)
         call getcmd(5,'output_g',wcmd,lcmd,3)
         if (lcmd.le.0 ) then
           write(*,*) '** Output G-file name not found; using gfile.out'
           goutf = 'gfile.out       '
         else
           goutf = wcmd(1:lcmd)
         endif
c**        read one-line G-file description
c**        no longer read from input file - create here
c        call getcmd(5,stone,wcmd,lcmd,2)
c        call getcmd(5,'g_line',wcmd,lcmd,3)
c        gline = wcmd(1:lcmd)
         call getusr(uname)
         call runtim(dattim)
         write(gline,'(a,a16,1x,a16,1x,a19)')
     .         'By SOLVE ',qfiln,uname,dattim
      endif

 400  continue
c.....
 800  return
      end
