      subroutine somake (
     .       mfile,qfile,lfile,gfile,nfile,xfile,bfil2
     .    ,  lsess,lsitt,lsesfo,lxfil,lcfil
     .    ,  xorc,kstns,snames,nsat,norb,observ,idec
     .    ,  zenest,nzen,zenmod,gradest,ngrad,gradmod
     .    ,  iyr,idoy,idatum,proc,ierr_sestbl,ierr_sittbl
     .   ,   fixdrv_vers,scratch_dir )
C
C     Write a new- (FONDA-) style batch file for SOLVE
C    `   
c     ---Input---
c            * Files names *
C            MFILE   M-file
C            LFILE   L-file
C            GFILE   Current G-file name. Updated in this routine
c            NFILE   Optional auxillary solve batch file with stn & sat dependent noise
c            XFILE   X-files in solution
c            BFIL2   secondary B-file
c            * Unit numbers
c            LSESS   sestbl.
c            LSITT   sittbl.
c            LSESFO  session.info
c            LXFIL   X-file
c            LCFIL   C-file
c            * Solution controls *
c            XORC    input is 'X' or 'C' file
c            SNAMES  4-character site codes 
c            NSAT            number of satellites on T-file
c            IYR     year of session
c            IDOY    day-of-year of session
c            IDATUM  datum for coordinate converstions
c
c    ---Output---
c            IERR_SESTB, IERR_SITTB   return code for sestbl. and sittbl. reads
c                    =0 if ok; = -1 if entry not found, unless REQD = .false.
c                    = -2 if format error; =-3 if entry improper

      implicit none

      include '../../gamit/includes/dimpar.h'

      integer*4 lsess,lsitt,lxfil,lcfil,lsesfo
     .        , nsat,isn
     .        , kstns,ksats,nxepoc,ischan1(maxsat),isessn
     .        , istart,iend,check_flg
     .        , nblen,lchar,ichar1,ichar2,ichar3,nout
     .        , iyr,idoy,iwkn,int1,jsat,norb,numarg,count_arg
     .        , ihmode,idatum,ieop,narg,idec,mssn
     .        , ill,ierr_sestbl,ierr_sittbl,i0,i,j
     .        , ival,err,indx,nzen,ngrad,ioerr,nskip,nsamp,mchkey

      real*4 error_a,error_b,akappa,bkappa
     .     , satcv(maxsat,maxorb),statcv(maxsit,3)
     .     , czen_all,zenvar_all,zentau_all
     .     , czen(maxsit),vzen(maxsit),tzen(maxsit)
     .     , cgrad_all,gradvar_all,gradtau_all
     .     , cgrad(maxsit),vgrad(maxsit),tgrad(maxsit)
     .     , wldev,wlsig,wlcut,wlrat,wldmax
     .     , nldev,nlsig,nlcut,nlrat,nldmax
     .     , cwob,cwobr,cut1,cut1r,prdev,prsig,prcut
     .     , allcut,cutoff(maxsit)
     .     , correl_prt,coord_upd_tol,p1relweight,biasapr,brcond

      real*8 sow

      character*1  yes,trocon,sitcon,orbcon,xorc
     .          ,  svantcon, printhfile, zenithdelay,ans 
     .          ,  nlscale
      character*2  buf2,pn,aeop
      character*3  fixcrd(maxsit),lowerc,zenmod,gradmod,proc,type_cut
     .          ,  afix
      character*4  snames(maxsit),satnam(maxsat)
     .          ,  bconstrain
      character*5  quick,exprmt,fixdrv_vers,position
      character*6  aczen,acgrad,buf6,bias,biassearch
      character*8  biases
      character*9  error_model
      character*10 hmode,avzen,avgrad
      character*11 biaspath
      character*15 observ,observp
      character*16 lfile,mfile,qfile,gfile,gfile1
     .           , uname,lfile1,xfile(maxsit),bfil2,nfile
      character*17 decision
      character*20 acwob,acut1,pweight
      character*21  string
      character*80 scratch_dir
      character*256 message,args,line,line2,blnk256

      logical reqd,fcheck,debug,found,xsats,xstns,kbit,in_sestbl
     .      , eof,zenest,gradest


      blnk256 = ' '
      yes = 'y'
                         
c  Add the SOLVE run to the experiment bat file and open the secondary file

      write(17,10) bfil2
   10 format('solve <  ',A16)
      open( 20, file=bfil2, form='formatted',status='unknown' )

C Get user name for Q-, H-file headers

      call getusr(uname)

c  Write the header comments, version, mode, and owner

      write(20,101)
  101 format( 
     . '*',
     . '-------------------------------------------------------------',/
     .,'*   << key-word-controlled format >>                        *',/
     .,'* symbol ":" must exist in command lines as separator       *',/
     .,'* any non-blank character at first column means comment line ',/
     .,'* empty after ":" means comment line too                    *',/
     .,'-------------------------------------------------------------',/
     .,'*',/
     .,'-------------- Part 1 -- Files and Global Controls')

      write(20,'(a,a5)') ' FIXDRV version:     ',fixdrv_vers  
      write(20,'(a)')    ' operation mode:    batch'
      write(20,'(a,a4)') ' owner:            ',proc  
      if( scratch_dir(1:1).ne.' ' ) then
        write(20,'(2a)') ' scratch directory: ',scratch_dir
      endif

                                    
c  Set the quick/full control: quick no longer supported, so hardwire
 
c     if ( irun.eq.1 ) quick = 'quick'
      quick = 'full '

c  Write the Q-, H-, and M-file names  (**reorder and make datum, h-mode explicit)

c     Control for writing constrained-solution into H-file
c         ihmode = 0 : default, write loose solutions only
c         ihmode = 1 : write tight solutions only (for quick solution)
c                = 2 : write loose and constrained solutions     
c         Note:  In autcln post-fit mode, no h-file will be written for the
c                prefit solution if imode = 0 
       ihmode = 0
       reqd = .false.
       call rdsest( 16, 'H-file solutions',10, hmode, lsess, reqd, ill )
       if( ill.ne.0 ) ierr_sestbl = ill
       if (quick.eq.'quick') ihmode=1
       if( hmode(1:3).eq.'ALL' .and. quick.eq.'full') ihmode=2
      write(20,110) qfile,ihmode,idatum,mfile
  110 format(
     .  ' Q-file name:       ',a16,/
     ., ' H-file mode:       ',i1,/
     ., ' datum code:        ',i1,/
     ., ' M-file name:       ',a16 )

c  Choice of experiment : BASEL RELAX ORBIT 

      reqd = .false.
      exprmt = 'BASEL'
      call rdsest( 20,'Choice of Experiment',5,exprmt,lsess,reqd,ill )
      if( ill.ne.0 ) ierr_sestbl = ill

c  Bias search options --- only decision function now allowed
      reqd = .false.
      call rdsest( 11,'Bias Search',6,bias,lsess,reqd,ill )
      if( bias.eq.'      ') bias = 'DEC_FN'
      if( ill .ne. 0 ) ierr_sestbl = ill

C  Bias constraint
      reqd = .false. 
      call rdsest( 14,'Bias Constrain',4,bconstrain,lsess,reqd,ill )
      if( bconstrain.eq.'    ' ) bconstrain = 'NONE'
      if( ill .ne. 0 ) ierr_sestbl = ill    
                
c  NL sigma rescale (default yes, hidden feature to turn off)
      reqd = .false. 
      call rdsest( 8,'NL scale',1,nlscale,lsess,reqd,ill )  
      if( nlscale.eq.' ' ) nlscale = 'Y'
      if( ill .ne. 0 ) ierr_sestbl = ill    

C Check for validity of previous two options
      if(biaspath(1:3).eq.'ALL' .and. biassearch(1:3).eq.'SUB') then  
        write(message,'(a,a11,a6)') 
     .        'Invalid biaspath/biassearch option pair in sestbl.:'
     .       , biaspath, biassearch 
        call report_stat('WARNING','FIXDRV','somake',' ',message,0)
        ierr_sestbl = -3
      endif

c  Write observable options

      write(20,120) observ
  120 format(
     . ' phase difference options: double difference',/
     .,' combination mode:         ',a15)
c       ambiguity resolution 
        if( observ(1:7).ne.'LC_ONLY') then 
c         defaults:
          nldev = 0.15
          nlsig = 0.15
          nlcut = 1000.
          nlrat = 10.
          nldmax = 500.
          reqd = .false.
          call rdsest( 23, 'Ambiguity resolution NL', 30, line
     .               , lsess, reqd, ill )
          if( ill.ne.0 ) ierr_sestbl = ill
          i0 = index( line, '=' ) + 1
c         LINE will be blank if nothing found
          if( line.ne.blnk256 ) then
c           need to count arguments to avoid ugly abort for out-of-date sestbl.
            args = ' '
            args = line(i0:nblen(line))
            numarg= count_arg(line)
            if( numarg.eq.5 ) then
               read( line(i0:nblen(line)), *)
     .          nldev,nlsig,nlcut,nlrat,nldmax
            else
              read( line(i0:nblen(line)), *) nldev,nlsig,nlcut,nlrat 
              write(message,'(2a)') 
     .             'Missing max BL in sestbl. Ambiguity resolution NL,'
     .           , ' set to default 500 km'
             call report_stat('WARNING','FIXDRV','somake',' ',message,0)
            endif
          endif

c         Write the ambiguity resolution approach

          if(bias.eq.'DEC_FN') then   
c           Dong and Bock decision function approach
            decision = 'decision function'
            write(20,125) decision,nldev,nlsig,nlcut,nlrat,nldmax
  125      format('    bias search approach:  ',a17,/
     .,'    search path:           narrow lane',/
     .,'    search criteria:   ', 4f8.2,f9.1 ) 
c           write this only if you want to turn off rescaling 
            if( nlscale.eq.'N') then
              write(20,'(a)') '    nlscale: N'
            endif  
          else
            call report_stat('FATAL','FIXDRV','somake',' '
     .   ,'Only decision function supported for ambiguity resolution',0)
          endif
      endif


c  Write start/stop epochs and decimation factor for SOLVE

      reqd = .false.
      istart = 1
c     get epochs from X-file or session.info
      if( fcheck(xfile(1)) ) then
          iend = nxepoc( xfile(1), xorc, lxfil, lcfil )
      else if( fcheck('session.info') ) then
         debug = .false.
         check_flg = 2.
c        this a kluge; isessn not really known here
         isessn=1    
c*         print *,'SOMAKE iyr doy ',iyr,idoy
         call rsesfo(lsesfo,debug,check_flg,iyr,idoy,isessn,ksats
     .              ,ischan1,iwkn,sow,int1,iend,found )    
      else
          write(message,'(a,a16,a)') ' Neither X-file (',xfile(1)
     .         , ') nor session.info available--cannot get # epochs' 
          call report_stat('FATAL','FIXDRV','somake',' ',message,0)
      endif
      reqd = .false.
      call rdsest(13, 'Select Epochs', 30, line, lsess, reqd, ill )
      if( ill.ne.0 ) ierr_sestbl = ill
c     line will be blank if nothing found
      i0 = index( line, '=' ) + 1
      if( line.ne.blnk256 ) read(line(i0:nblen(line)),*) istart,iend
      write(20,130) istart,iend,idec
  130 format(' start and end epochs:   ',3(1x,i4 ))

c  Write cutoff angle

c      default now 0 to use autcln cutoffs
       allcut = 0.
       type_cut = 'all'
       in_sestbl = .false.
c      read cutoff for all sites from sestbl.
       reqd = .false.        
       call rdsest(15,'Elevation Cutoff',6,buf6,lsess,reqd,ill )
       if( ill.ne.0 ) ierr_sestbl = ill 
       if( buf6(1:4).ne.'    ' ) then
          read(buf6,'(f6.0)',iostat=ioerr) allcut    
          if( ioerr.ne.0 )  call report_stat('FATAL','FIXDRV','somake'
     .      ,' ','Error reading cutoff from sestbl.',ioerr)
          in_sestbl = .true.  
       endif           
c      override with sittbl. values if present
       reqd = .false.                  
       do i = 1,kstns
         cutoff(i) = allcut
         call rdsitt(snames(i),6,'CUTOFF',nout,line,lsitt,reqd,ill)
         if( ill.ne.0 ) ierr_sittbl = ill 
         if( line(1:6).ne.'      ') then
            read( line(1:6),'(a6)',iostat=ioerr) buf6 
            if( ioerr.ne.0 ) then
                 write(message,'(a,a4)') 
     .              'Error reading cutoff from sittbl.; site=',snames(i)
                     call report_stat('FATAL','FIXDRV','somake',' '
     .                              , message,ioerr)
            endif   
            if( buf6(1:6).ne.'      ' ) then 
c              this leads to a misread if no decimal point - use f6.0 or *
c               read(buf6,'(f6.2)',iostat=ioerr) cutoff(i) 
                read(buf6,'(f6.0)',iostat=ioerr) cutoff(i)
               if( ioerr.ne.0 ) then
                 write(message,'(a,a4)') 
     .              'Error reading cutoff from sittbl.; site=',snames(i)
                     call report_stat('FATAL','FIXDRV','somake',' '
     .                              , message,ioerr)
               endif
               if( cutoff(i).ne.allcut ) then
                  type_cut = 'stn'
                  if( in_sestbl ) then
                      write(message,'(a,f6.2,a,f6.2,a,a4)')
     .                   '  Warning:  sestbl elevation cutoff (',allcut
     .               ,') overriden by sittbl entry (',cutoff(i),') for '
     .               , snames(i)
                      call report_stat('WARNING','FIXDRV','somake',' '
     .                               , message,0)
                  endif
               endif
            endif
         endif
       enddo 
      write(20,'(a)') ' set cutoff_elevation:'
      if( type_cut.eq.'all' ) then
         write(20,'(a,f6.2)') '   cutoff: all_sites',allcut
      elseif( type_cut.eq.'stn' ) then
         do i = 1,kstns
           write(20,'(a,1x,a4,1x,f6.2)')
     .     '   cutoff:',snames(i),cutoff(i)
         enddo
      else
         call report_stat('FATAL','FIXDRV','somake',' '
     .                   ,'Bad logic for elevation cutoff',0)
      endif
                        
c Write bias-parameter apr (global, possibly for debug only)
         
      call rdsest( 12,'Bias apriori',12,line,lsess,reqd,ill)
      if( ill.ne.0 ) ierr_sestbl = ill
      i0 = index( line, '=' ) + 1
      if( line.ne.blnk256 ) then
        read(line(i0:nblen(line)),*) biasapr 
      else
        biasapr = 1000.  
      endif 
      write(20,'(a,f7.2)') ' bias_apr: ',biasapr   

c  Condition-ratio for inversion to detect dependent biases    

      call rdsest( 10,'Bias rcond',12,line,lsess,reqd,ill)
      if( ill.ne.0 ) ierr_sestbl = ill
      i0 = index( line, '=' ) + 1
      if( line.ne.blnk256 ) then
        read(line(i0:nblen(line)),*) brcond
      else  
c** default set high to avoid missing an independent bias
c   tah/scm/rwk 070501
        brcond = 10000.
      endif 
      write(20,'(a,g12.3)') ' bias_rcond: ',brcond            
      call rdsest(10,'Bias_debug',1,ans,lsess,reqd,ill )
      if( ans.eq.'Y'.or. ans.eq.'y') then
        write(20,'(a)') ' bias_debug: Y'      
      else
         write(20,'(a)') ' bias_debug: N' 
      endif


c Write control for screen print
                       
      call rdsest( 11,'SOLVE print',1,ans,lsess,reqd,ill ) 
      if( ill.ne.0 ) ierr_sestbl = ill
      if( ans.eq.' '.or.ans.eq.'N' ) then    
        write(20,'(a)') ' log print: N'      
      else
        write(20,'(a)') ' log print: Y' 
      endif
                          
c Allow skipping the loose solution

      call rdsest( 10,'Skip loose',1,ans,lsess,reqd,ill )
      if( ans.eq.' '.or.ans.eq.'N' ) then    
        write(20,'(a)') ' skip loose: N'      
      else
        write(20,'(a)') ' skip loose: Y' 
      endif

c------------------------------------------------------------------------


c  Write parameter options

      write(20,230)
  230 format('*',/,'-------------- Part 2 -- Parameters',/
     .          ,' set parameters:')

C  Station-dependent parameters

       write(20,231)
  231  format('    estimate:   all_sites all_parameters')
       write(20,232)
  232  format('    fix:        all_sites clock')
       if( observ(1:3).eq.'L1_' ) write(20,233)
  233  format('    fix:        all_sites l2-bias')
       if( observ(1:3).eq.'L2_' ) write(20,234)
  234  format('    fix:        all_sites l2-bias',/
     .      ,'#               l1-biases used for l2-only observations')
       if( quick.eq.'quick' ) write(20,235)
  235  format('    fix:        all_sites l1-bias l2-bias',/
     .      ,'#               no biases for quick solution')

c      Zenith-delay estimation
       if( .not.zenest ) write(20,236)
  236  format('    fix:        all_sites zenith') 
       if( .not.gradest ) write(20,237)
  237  format('    fix:        all_sites grad')
       if( exprmt.eq.'ORBIT' ) write(20,238)
  238  format('    fix:        all_sites coord')


c      loop over stations
       do 260 i = 1, kstns

c        Check if any coordinates to be fixed :  FIXCRD = DDD / YYY / NNN
         reqd = .false.
         fixcrd(i) = 'NNN'   
         call rdsitt( snames(i),3,'FIX',nout,afix,lsitt,reqd,ill )   
         if( ill.ne.0 ) ierr_sittbl = ill  
         if( afix(1:1).ne.' ' ) then 
            read(afix,'(a3)') fixcrd(i)
         endif  
         if( fixcrd(i).eq.'YYY' ) then
            write(20,240) snames(i)
  240       format('   fix:  ',a4,' coordinate')
         endif
               
C        Set clock parameters
         reqd = .false.
         call rdsitt( snames(i),5,'CLKFT',nout,line,lsitt,reqd,ill ) 
         if( ill.ne.0 ) ierr_sittbl = ill
         if( line(1:1).ne.' ' ) then
           if( line(1:1).eq.'Y')
     .        write(30,251) snames(i)
  251         format(' estimate: ',a4,' clock_offset')
           if ( line(2:2).eq.'Y' .or. line(3:3).eq.'Y' ) 
     .         call report_stat('FATAL','FIXDRV','somake',' '
     .   ,'Can no longer estimate station clock rate or acceleration',0)
         endif
  260 continue


c  Satellite parameters

c     Set satellite parameters
                                       
      if( exprmt.eq.'RELAX' .or. exprmt.eq.'ORBIT' ) then
         write(20,'(a)') '    estimate:   all_sats all_parameters'
         if( ill.ne.0 ) ierr_sestbl = ill
         reqd = .false.
         call rdsest(14,'SV antenna off',1,svantcon,lsess,reqd,ill) 
         if( svantcon.ne.'Y') then
             svantcon = 'N'
             write(20,'(a)') '    fix:        all_sats svant'
         else
             if( norb+3.gt.maxorb ) then
               write(message,'(a,i2,a,i2)') 
     .            '# orbit + svant parameter constraints =',norb+3
     .           ,' exceeds MAXORB= ',maxorb
               call report_stat('FATAL','FIXDRV','somake',' ',message,0)
              endif 
 
         endif
       endif

c  **No changes in satellite parameters estimated currently allowed by sestbl.**
c    If added, the code at satellite constraints to get PN numbers must be
c    moved up.


c  Global (earth-orientation) parameters
          
c     set default:  All if estimating orbital parameters, none if coordinates only  
      if( exprmt.eq.'BASELINE' ) then
        ieop = 0
      else
        ieop = 15   
      endif
      reqd = .false.
      call rdsest( 12,'Estimate EOP',2,aeop,lsess,reqd,ill )
      if( ill.ne.0 ) ierr_sestbl = ill
      if( aeop(1:1).ne.' ' ) read(aeop,'(i2)') ieop
      if( ieop.eq.15 ) then
c**       write(20,'(a)') '    estimate:   global all_parameters'
          write(20,'(a)')
     .       '    estimate:   global wob ut1 wob_rate ut1_rate'
      else
        write(20,'(a)') '    fix:        global all_parameters'
        if(kbit(ieop,1)) write(20,'(a)') '    estimate:   global wob'
        if(kbit(ieop,2)) write(20,'(a)') '    estimate:   global ut1'
        if(kbit(ieop,3)) write(20,'(a)')
     .                              '    estimate:   global wob_rate'
        if(kbit(ieop,4)) write(20,'(a)')
     .                              '    estimate:   global ut1_rate'
      endif

      write(20,290)
  290 format(' exit set:' )

c----------------------------------------------------------------------------


c  Write a priori constraints

      write(20,300)
  300 format('*',/
     . '-------------- Part 3 -- A priori Constraints',/
     .,' set apriori constraints:')


c  Station-coordinate constraints

      reqd = .false.
      call rdsest( 18,'Station Constraint',1,ans,lsess,reqd,ill ) 
      if( ans.eq.' ') then
         sitcon = 'Y'
      else
         sitcon = ans
      endif
      if( ill.ne.0 ) ierr_sestbl = ill
      if( sitcon.eq.'Y') then
c       loop over stations
        do  i = 1, kstns
          reqd = .true.
          call rdsitt( snames(i),13,'COORD.CONSTR.'
     .               , nout,line,lsitt,reqd,ill ) 
c          print *,'MNET i sname nout line ill '
c     .           , kstns,i,sname,nout,line,ill
          read( line(1:nout),*,iostat=ioerr ) (statcv(i,j),j=1,3)
          if( ioerr.ne.0 ) call report_stat('FATAL','FIXDRV','somake'
     .      , ' ','Error reading station a priori constraints',ioerr)
c rwk 990324  The following changed for simplicity; maybe not necessary or robust
c          if( ill.ne.0 ) ierr_sestbl = ill
c          i1 = index( line( 1:nout), '.' )
c          i2 = index( line(i1:nout), ' ' ) + i1 - 1
c          i3 = index( line(i2:nout), '.' ) + i2 - 1
c          i4 = index( line(I3:nout), ' ' ) + I3 - 1
c          read( line(   1:i2  ), '(f5.3)',iostat=ioerr,err=305 )  
c     .         statcv(i,1)
c          read( line(i2+1:i4-1), '(f5.3)',iostat=ioerr,err=305 ) 
c     .         statcv(i,2)
c          read( line(i4+1:nout), '(f5.3)',iostat=ioerr,err=305 ) 
c     .          statcv(i,3)    
c  305     if( ioerr.ne.0 ) then
           if( ioerr.ne.0 ) then
            write(message,'(a,a4)') 
     .       'Error reading coord constraints from sittbl.; site='
     .        ,snames(i)
            call report_stat('FATAL','FIXDRV','somake',' ', message
     .                      ,ioerr)
          endif
c         check for zero values -- not allowed
          do j=1,3
            if( statcv(i,j).eq.0.d0 ) then
               write(message,'(a,i2,a)')  'A priori sigma for site '
     .                                     , i,' is zero in sittbl.'
               call report_stat('WARNING','FIXDRV','somake',' '
     .                         , message,0)
               ierr_sittbl = -1
            endif
          enddo
          write(20,310) snames(i),(statcv(i,j),j=1,3)
  310     format('    tight_apr_coord:  ',a4,3(1x,f8.3) )   
c          print *,'wrote tight_apr_coord '
        enddo
      endif       
c rwk 070724: Tighten the loose constraints for short baselines to avoid
c             numerical problems; assume baselines are short for single frequency
      if( observ(1:2).eq.'L1'.or.observ(1:2).eq.'L2') then
        write(20,'(a)') '    loose_apr_coord:  all_  1.  1.  1.'
        call report_stat('WARNING','FIXDRV','somake',' '
     .     ,'Loose-solution constraints set to 1 m for L1 and/or L2',0)
      else
        write(20,'(a)') '    loose_apr_coord:  all_  10.  10.  10.'
      endif

c  Zenith-delay numbers and constraints
                  
      if( zenest ) then
c       model and number of zenith delays read in bmake for m-file; cannot change here.
        czen_all = 0.50
        zenvar_all = .02
        zentau_all = 100.
        reqd = .false.
        call rdsest( 10,'Zenith Con',6,aczen,lsess,reqd,ill )
        if( ill.ne.0 ) ierr_sestbl = ill
        call rdsest( 10,'Zenith Var',10,avzen,lsess,reqd,ill )
        if( ill.ne.0 ) ierr_sestbl = ill
        if( aczen(1:1).ne.' ' ) read(aczen,'(f6.0)') czen_all
        if( avzen(1:1).ne.' ' ) read(avzen,*) zenvar_all,zentau_all
        write(20,320) nzen,zenmod
  320   format('    zenith delays: all_sites ',i3,1x,a3)
C       override constraints with sittbl. values 
        do  i = 1, kstns
          czen(i) = czen_all
          vzen(i) = zenvar_all
          tzen(i) = zentau_all
          call rdsest( 10,'Zenith Con',6,aczen,lsess,reqd,ill )
          if( ill.eq.0 ) then
            call rdsitt( snames(i),6,'ZCNSTR',nout,line,lsitt
     .               , reqd,ill )
            if( ill.ne.0 ) ierr_sittbl = ill
            if( line(1:6).ne.'      ') 
     .       read( line(1:nout),'(f6.0)',iostat=ioerr,err=325) czen(i)
            call rdsitt( snames(i),6,'ZENVAR',nout,line,lsitt
     .               , reqd,ill )
            if( ill.ne.0 ) ierr_sittbl = ill
            if( line(1:6).ne.'      ') 
     .        read( line(1:nout),'(f6.0)',iostat=ioerr,err=325) vzen(i)
            call rdsitt( snames(i),6,'ZENTAU',nout,line,lsitt
     .                 , reqd,ill )
            if( ill.ne.0 ) ierr_sittbl = ill
            if( line(1:6).ne.'      ') 
     .         read( line(1:nout),'(f6.0)',iostat=ioerr,err=325) tzen(i) 
  325       if( ioerr.ne.0 ) then
              write(message,'(a,a4)') 
     .       'Error reading atmospheric constraints from sittbl.; site='
     .            ,snames(i)
              call report_stat('FATAL','FIXDRV','somake',' ', message
     .                        ,ioerr)
            endif
          else
             ierr_sittbl = ill
          endif
          write(20,330) snames(i),czen(i),vzen(i),tzen(i)
  330     format('    tight_apr_zenith:  ',a4,f6.3,2x,f6.3,2x,f5.1 )
c         loose should be the same as tight since GLOBK does not estimate them
          write(20,331) snames(i),czen(i),vzen(i),tzen(i)
  331     format('    loose_apr_zenith:  ',a4,f6.3,2x,f6.3,2x,f5.1 )
        enddo                                                       
      endif

c  Atmospheric gradient constraints     

      if( gradest ) then 
c       (For now assign equal gradients to NS and EW directions) 
c       (Number of gradients     set in cmmake for m-file.  Cannot change here.)
        cgrad_all = 0.03
        gradvar_all = .01
        gradtau_all = 100.
        reqd = .false.
        call rdsest( 10,'Gradient Con',6,acgrad,lsess,reqd,ill )
        if( ill.ne.0 ) ierr_sestbl = ill
        call rdsest( 10,'Gradient Var',10,avgrad,lsess,reqd,ill )
        if( ill.ne.0 ) ierr_sestbl = ill
        if( acgrad(1:1).ne.' ' ) read(acgrad,'(f6.0)') cgrad_all
        if( avgrad(1:1).ne.' ' ) read(avgrad,*) gradvar_all,gradtau_all
        write(20,340) ngrad,gradmod
  340   format('    gradients    : all_sites ',i2,1x,a3)     
C       override with sittbl. values - ngrad must be the same for now
        do  i = 1, kstns
          cgrad(i) = cgrad_all
          vgrad(i) = gradvar_all
          tgrad(i) = gradtau_all
          call rdsest( 12,'Gradient Con',6,acgrad,lsess,reqd,ill )
          if( ill.eq.0 ) then
            call rdsitt( snames(i),6,'GCNSTR',nout,line,lsitt
     .               , reqd,ill )
            if( ill.ne.0 ) ierr_sittbl = ill
            if( line(1:6).ne.'      ') 
     .       read( line(1:nout),'(f6.0)',iostat=ioerr,err=345) cgrad(i)
            call rdsitt( snames(i),6,'GRDVAR',nout,line,lsitt
     .               , reqd,ill )
            if( ill.ne.0 ) ierr_sittbl = ill
            if( line(1:6).ne.'      ') 
     .       read( line(1:nout),'(f6.0)',iostat=ioerr,err=345) vgrad(i)
            call rdsitt( snames(i),6,'GRDTAU',nout,line,lsitt
     .               , reqd,ill )
            if( ill.ne.0 ) ierr_sittbl = ill
            if( line(1:6).ne.'      ') 
     .        read( line(1:nout),'(f6.0)',iostat=ioerr,err=345) tgrad(i) 
  345       if( ioerr.ne.0 ) then
              write(message,'(a,a4)') 
     .       'Error reading gradient constraints from sittbl.; site='
     .            ,snames(i)
             call report_stat('FATAL','FIXDRV','somake',' ', message
     .                      ,ioerr)
            endif
          else
           ierr_sittbl = ill
          endif
          write(20,346) snames(i),cgrad(i),vgrad(i),tgrad(i)
  346     format('    tight_apr_gradient:  ',a4,f6.3,2x,f6.3,2x,f5.1 )
c         loose should be the same as tight since GLOBK does not estimate them
          write(20,347) snames(i),cgrad(i),vgrad(i),tgrad(i)
  347      format('    loose_apr_gradient:  ',a4,f6.3,2x,f6.3,2x,f5.1 )
        enddo 
      endif 


c  Satellite constraints (integrated parameters + SV antenna offsets)

      reqd = .false.
      call rdsest( 20,'Satellite Constraint',1, orbcon,lsess,reqd, ill )
      if( ill.ne.0 ) ierr_sestbl = ill  
      if( orbcon.eq.' ')  orbcon = 'N'
      if( orbcon.eq.'Y' ) then
c        header line must follow - read it    
         read( lsess,'(a256)') line 
c        Note:  norb is the number of (integrated) orbital parameters on the T-file
c               narg is the number of orbit parameters in the sestbl. entry, which
c                  can be less (if an oversight--rest will be filled in) or more
c                  if SV antenna offsets are estimated
c               maxorb is the dimension in dimpar.h, and must be large enough to 
c                  include both the integrated parameters and 3 SV ant offsets

c        new-style: one line for all satellites
         if( index(line,'all').ne.0 ) then
              read( lsess, '(a256)',iostat=ioerr) line
              if( ioerr.ne.0 ) call report_stat('FATAL','FIXDRV'
     .        ,'somake',' ','Error reading satellite constraints',ioerr)
c             check number of orbital parameters
              narg = count_arg(line)   
              if( narg.gt.maxorb ) then       
               write(message,'(i2,a,i2)') narg
     .           , 'orbit constraints is sestbl.but dimension MAXORB= '
     .           ,maxorb
               call report_stat('FATAL','FIXDRV','somake',' ',message,0)
              endif 
              if( narg.lt.norb ) then  
                 write(message,'(i2,a,i2,a)') norb
     .                ,' orbit parameters but only ',narg
     .                ,' constraints in the sestbl.'
                 call report_stat('WARNING','FIXDRV','somake',' '
     .                           ,message,0)   
                 call report_stat('WARNING','FIXDRV','somake',' '
     .                    ,'  Missing values set to last one',0)
                 if( narg.lt.9 ) ierr_sestbl = -3
              endif
              do i=1,nsat 
                read( line,*,iostat=ioerr)  ( satcv(i,j), j=1,narg )  
               if( ioerr.ne.0 ) call report_stat('FATAL','FIXDRV'
     .        ,'somake',' ','Error reading satellite constraints',ioerr)
                if( narg.lt.norb ) then
c                 fill in the rest of the radiation-pressure parameters from those there
                  do j=narg+1,norb
                    satcv(i,j) = satcv(i,narg)
                  enddo  
                endif  
                if( svantcon.eq.'Y' .and. narg.lt.(norb+3) ) then
c                 set SV antenna constraints to default 10 cm
                  do j=norb+1,norb+3
                    satcv(i,j) = 0.1d0
                  enddo
                endif
              enddo 
c             check for zero values -- not allowed
              do j=1,norb
                if( satcv(1,j).eq.0.d0 )
     .            call report_stat('WARNING','FIXDRV','somake',' '
     .           ,'A priori sigma for satellites from sestbl is zero',0)
              enddo
              write(20,'(2a)')
     .        '*     units are ppm for elements, percent for rad parms'
     .              ,', m for SV antenna offsets'   
              write(20,'(a,6(1pe8.1),13(1pe8.1))') 
     .            '    tight_apr_orbit:  all_',(satcv(1,j),j=1,norb)
              if( svantcon.eq.'Y' ) write(20,'(a,3f8.3)') 
     .            '    tight_apr_svant:  all_',(satcv(1,j+norb),j=1,3)

c        old-style: one line for every satellite
         elseif( index(line,'ch#').ne.0 ) then
c           get PN #'s from session.info if available
            if( fcheck('session.info')) then
c             this a kluge; isessn not really known here
              isessn=1
              debug = .false.
              check_flg = 2 
              call rsesfo( lsesfo,debug,check_flg,iyr,idoy,isessn
     .                   , ksats,ischan1,iwkn,sow,int1,iend,found )    
            endif
c           nsat from X-, C-, T-, or G-file in BMAKE; current SOLVE cannot handle
c           mismatched numbers of satellites among sessions
c           now read constraints for each satellite
            do 350 i = 1,ksats
              read( lsess, '(a256)',iostat=ioerr ) line
              if( ioerr.ne.0 ) call report_stat('FATAL','FIXDRV'
     .        ,'somake',' ','Error reading satellite constraints',ioerr)
              i0 = index( line, 'NO.' ) + 3
              if( i0 .eq. 3 )  then
                 call report_stat('WARNING','FIXDRV','somake',' '
     .    , '# satellites in sestbl. constraint list <  # satellites',0)
                 ierr_sestbl = -3
              endif
              if( line(i0+1:i0+1) .EQ. ' ' )  then
                 read( line(i0:i0), '(i1)' )  jsat
              else
                 read( line(i0:i0+1), '(i2)' )  jsat
              endif
c             much more user friendly if this is free format
              i0 = index( line, '=' ) + 1
              i0 = max(i0,i0) 
   
c             check number of orbital parameters 
              line2 = blnk256
              line2(1:256-i0) = line(i0+1:256)
              narg = count_arg(line2)     
              if( narg.gt.maxorb ) then
                write(message,'(i2,a,i2)') narg
     .            , 'orbit constraints is sestbl.but dimension MAXORB= '
     .           ,maxorb
               call report_stat('FATAL','FIXDRV','somake',' ',message,0)
              endif 
              if( narg.lt.norb ) then 
                 write(message,'(i2,a,i2,a)') norb
     .                ,' orbit parameters but only ',narg
     .                ,' constraints in the sestbl.'
                 call report_stat('WARNING','FIXDRV','somake',' '
     .                          , message,0) 
                 call report_stat('WARNING','FIXDRV','somake',' '
     .                    ,'  Missing values set to last one',0)
                 if( narg.lt.9 ) ierr_sestbl = -3
              endif
              read( line(i0:nblen(line)), *)  ( satcv(i,j), j=1,narg )  
              if( narg.lt.norb ) then
c                fill in the rest of the radiation-pressure parameters from those there
                 do j=narg+1,norb
                    satcv(i,j) = satcv(i,narg)
                 enddo  
              endif   
c             check for zero values -- not allowed
              do j=1,norb
                if( satcv(i,j).eq.0.d0 ) then
                  write(message,'(a,i2,a)') 
     .               'A priori sigma for satellite ', jsat,' is zero'
                  call report_stat('WARNING','FIXDRV','somake',' '
     .               ,message,0)
                endif
              enddo    
              if( svantcon.eq.'Y' .and. narg.lt.(norb+3) ) then
c               set SV antenna constraints to default 10 cm
                do j=norb+1,norb+3
                  satcv(i,j) = 0.1d0
                enddo
              endif
c             allow different constraints only if PN #s available from
c             session.info.  Can avoid this by putting PN #s in the sestbl.
c             or reading them from the X-file, if available.
              if( fcheck('session.info') ) then
c               get PN #s from session.info
                write(buf2,'(i2)') ischan1(i)
                read(buf2,'(a2)') pn
                if( pn(1:1).eq.' ' ) pn(1:1)='0'
                satnam(i) = 'PN'//pn
                if( i.eq.1 ) write(20,'(a)')
     .         '*     units are ppm for elements, percent for rad parms' 

                write(20,'(a,a4,6(1pe8.1),9(1pe8.1))') 
     .          '    tight_apr_orbit:  ',satnam(i),(satcv(i,j),j=1,norb)
                if( svantcon.eq.'Y' ) write(20,'(a,a4,3f8.3)') 
     .        '    tight_apr_svant:  ',satnam(i),(satcv(1,j+norb),j=1,3)
              endif
c          end of loop on satellites
  350      continue
           if( .not.fcheck('session.info')) then
              call report_stat('WARNING','FIXDRV','somake',' '
     .  , 'No session.info so making all satellite constraints equal',0)
              write(20,'(a,6(1pe8.1),9(1pe8.1))') 
     .           '    tight_apr_orbit:  all_',(satcv(1,j),j=1,norb)   
           endif

         else
            call report_stat('WARNING','FIXDRV','somake',' '
     .          ,'Satellite constraint header line missing',0 )
         endif
c     end of condition on orbital constraints
      endif
          
      if( exprmt.eq.'RELAX' .or. exprmt.eq.'ORBIT' )  then
          if( norb.eq.9 ) then
            write(20,360)
  360     format(
     .'    loose_apr_orbit:  all_ 1. 1. 1. 1. 1. 1.',
     .' 1000. 1000. 1000.' )
          elseif( norb.eq.12 ) then
            write(20,361)
  361     format(
     .'    loose_apr_orbit:  all_ 1. 1. 1. 1. 1. 1.',
     .' 100. 100. 100. 100. 100. 100.' )
          elseif( norb.eq.15 ) then
            write(20,362)
  362     format(
     .'    loose_apr_orbit:  all_ 1. 1. 1. 1. 1. 1.',
     .' 100. 100. 100. 100. 100. 100. 100. 100. 100.' )
          elseif( norb.eq.19 )  then
            write(20,363)
  363     format(
     .'    loose_apr_orbit:  all_ 1. 1. 1. 1. 1. 1.',
     .' 100. 100. 100. 100. 100. 100. 100. 100. 100.',
     .' 100. 100. 100. 100.' )
          endif       
          if( svantcon.eq.'Y' ) then
c           set loose defaults = 5 m   ok?
            write(20,'(a)') '    loose_apr_svant:  all_ 5. 5. 5.'
          endif
      endif

c  Global (earth-orientation) constraints

      if( ieop.gt.0 ) then
c        set defaults = 10 m and 1 m/day (in arcsec and sec) or read from sestbl.
         cwob = 0.3
         cwobr = 0.03
c        cut1 = 0.2
c        set UT1 tight (.02 ms=1 cm) for tight solution
         cut1 = 2.e-5
         cut1r = 0.02
         reqd = .false.
c        read them all if any are turned on
         call rdsest( 9,'Wobble Con',20,acwob,lsess,reqd,ill )
         if( ill.lt.0 ) return
         if( acwob(1:1).ne.' ' ) read (acwob,*) cwob,cwobr
         call rdsest( 7,'UT1 Con',20,acut1,lsess,reqd,ill )
         if( ill.lt.0 ) return
         if( acwob(1:1).ne.' ' ) read (acut1,*) cut1,cut1r
          write(20,'(a)')
     .         '*     units are s, s/d for UT1, arcs arcs/d for wobble'
         if( kbit(ieop,1).or.kbit(ieop,3) ) then
            write(20,'(a,2f8.5)') '    tight_apr_wob:  ',cwob,cwobr
            cwob = 3.
            cwobr = 0.3
            write(20,'(a,2f8.5)') '    loose_apr_wob:  ',cwob,cwobr
         endif
         if( kbit(ieop,2).or.kbit(ieop,4) ) then
            write(20,'(a,2f8.5)') '    tight_apr_ut1:  ',cut1,cut1r
            cut1 = 0.02
            cut1r = 0.02
            write(20,'(a,2f8.5)') '    loose_apr_ut1:  ',cut1,cut1r
         endif
      endif

      write(20,390)
  390 format(' exit set:')

c-----------------------------------------------------------------------


c  Session-dependent Options

      write(20,400)
  400 format('*',/
     . '---------------Part 4 -- Session Options' )



c   Site selections  
c** rwk 060818: hardwire to 1 session
         mssn = 1

         write(20,410) 
  410    format(' set session_1 options:',/
     .          '    include:       all_sites all_sats' )
         reqd = .false.
         call rdsest( 14, 'Station Number', 5, line, lsess
     .              , reqd, ill)   
         if( line(1:1).eq.' ' ) line(1:1) = '*'
         if (ill .ne. 0 )  return
         if( index(line(1:2),'*') .ne. 0 )  then
c           if station number is '*' read no more and include all stations
            do i=1,kstns
              line(i:i) = 'y'
            enddo
         else
            string = 'Stations Session     '
            if( mssn.lt.10 ) then
              write( string(18:18), '(i1)' ) isn
            else
               write( string(17:18), '(i2)' ) isn
            endif
            call rdsest( 18, string, kstns, line, lsess, reqd, ill )
            if( ill.ne.0 ) ierr_sestbl = ill
         endif
c        loop for excluding sites, but requires a session-dependent
c        array of site names--easier to do this manually for now, and later
c        to change sestbl. to read the names directly
         xstns= .false.
         do i=1,kstns
           if( line(i:i).ne.'y' .and. line(i:i).ne.'Y' )  xstns=.true.
         enddo
         if( xstns) then
              call report_stat('WARNING','FIXDRV','somake',' '
     .                 , 'Excluding sites not allowed by FIXDRV:',0)
              call report_stat('WARNING','FIXDRV','somake',' '
     .                 , '  edit SOLVE batch file manually',0)
              call report_stat('WARNING','FIXDRV','somake',' '
     .                 , ' Format: ''exclude:  MOJA'' ',0)
              ierr_sestbl = -3
         endif

c   Satellite selections

      ksats = nsat       
      do 420 isn = 1,mssn
         call rdsest( 16,'Satellite Number',5,line,lsess,reqd,ill )
         if( line(1:1).eq.' ' ) line(1:1) = '*'
         if( ill .ne. 0 )  ierr_sestbl = ill
         if( index(line(1:2),'*') .ne. 0 )  then
c           if satellite number is '*' read no more and include all sats
            do i=1,ksats
              line(i:i) = 'y'
            enddo
         else
            string = 'Satellites Session   '
            if( mssn.lt.10 ) then
              write( string(20:20), '(i1)' ) isn
            else
               write( string(19:20), '(i2)' ) isn
            endif
            call rdsest( 18,string(1:18),ksats,line,lsess,reqd,ill )
            if( ill.ne.0 ) ierr_sestbl = ill
         endif
         xsats =.false.
         do i=1,ksats
           if( line(i:i).ne.'y' .and. line(i:i).ne.'Y' ) xsats=.true.
         enddo
         if( xsats ) then
             call report_stat('WARNING','FIXDRV','somake',' '
     .                , 'Excluding satellites not allowed by FIXDRV:',0)
              call report_stat('WARNING','FIXDRV','somake',' '
     .                 , '  edit SOLVE batch file manually',0)
              call report_stat('WARNING','FIXDRV','somake',' '
     .                 , ' Format: ''exclude:  PN03'' ',0)
              ierr_sestbl = -3
         endif
  420  continue
              
   
c  Error model (read also in bmake to write out the sh_sigelv command)

c        sestbl. input allows baseline- or elevation-dependent weights 
c        N-file allows station- or satellite-dependent weights
       
c     set defaults
      error_model = 'uniform  '
      error_a = 10.
      error_b = 0.
c     --- old entry  : commented out with trap by rwk 050503
      reqd = .false.
      call rdsest( 17, 'Measurement Error', 30, line, lsess, reqd, ill )
      if( ill.ne.0 ) ierr_sestbl = ill
      if( line(1:5).ne.'     ') then      
        call report_stat('FATAL','FIXDRV','somake',' '
     .         ,'Measurement Error no longer supported in sestbl. ',0)
c        ichar1 = index( line, 'MM' ) - 1
c        ichar2 = index( line, '+' ) + 1
c        ichar3 = index( line, 'PPM' ) - 1
c        if ( ichar1 .GE. 10 )  then
c           call report_stat('WARNING','FIXDRV','somake',' '
c    .          , 'Session table "Measurement Error" wrong format.',0)
c           ierr_sestbl = -1
c        endif
c        write( format(3:3), '(I1)' )  ichar1
c        read( line(1:ichar1), format )  error_a
c        lchar = ichar3 - ichar2 + 1
c        if( lchar .ge. 10 )  then
c           call report_stat('WARNING','FIXDRV','somake',' '
c    .          , 'Session table "Measurement Error" wrong format.',0)
c           ierr_sestbl = -1
c        endif
c        write( format(3:3), '(I1)' )  lchar
c        read( line(ichar2:ichar3),format )  error_b
      endif
c     ---new entries via sestbl.
      call rdsest( 13, 'Station Error', 20, line, lsess, reqd, ill )
      if( line(1:5).ne.'    ' ) then
         indx = 1
         call read_line(line,indx,'CH',err,ival,error_model)
         call lowers(error_model)
         if( error_model.eq.'uniform  ' ) then
            read(line(indx:),*) error_a
         else if( error_model.eq.'baseline ' .or.
     .            error_model.eq.'elevation' ) then
            read(line(indx:),*) error_a,error_b
         endif
       endif
      if( error_model.eq.'baseline '.and.error_b .eq. 0. ) then
         error_model='uniform'
        write(20,430) error_model,error_a
  430   format('    error model:',/,'      stn_error: all_sites ',a9
     .        ,f5.1,/,'       sat_error: all_sats  0.0')
      else
        write(20,431) error_model,error_a,error_b
  431   format('    error model:',/,'      stn_error: all_sites ',a9
     .         ,2f6.1)
c       satellite weighting not yet allowed for baseline mode
        if( error_model.ne.'baseline' ) then
           write(20,432)
  432      format('      sat_error: all_sats  0.0')
        endif
      endif
c     N-file with station- or satellite-dependent weights from autcln (blank if none written)   
      write(20,'(a,a10)') '      noise file name:        ',nfile  


c  Tropospheric constraints

      reqd = .false.
      call rdsest( 24, 'Tropospheric Constraints', 1, trocon
     .           , lsess, reqd, ill )
      if( ill.ne.0 )  ierr_sestbl = ill
      if( trocon.eq.' ') trocon = 'N'
      write(20,440) trocon
  440 format('    atmosphere constraint:     ',a)

                  
c  Wide-lane ambiguity options 
             
      if( observ(1:7).eq.'L1&L2  ' .or.
     .    observ(1:7).eq.'L1,L2_I' .or.
     .    observ(1:7).eq.'LC_HELP' .or.
     .    observ(1:7).eq.'LC_AUTC' ) then
c          ionospheric constraints:  AKAPPA, BKAPPA
c          AKAPPA = constant term (mm),  BKAPPA = linear term(ppm)
        reqd = .false. 
        akappa = 0.d0 
        bkappa = 8.0d0
        call rdsest(23,'Ionospheric Constraints',30,line,lsess,reqd,ill)
        if( ill.ne.0 ) ierr_sestbl = ill
        ichar1 = index( line, 'MM' ) - 1  
        ichar2 = index( line, '+' ) + 1
        ichar3 = index( line, 'PPM' ) - 1
        if( ichar1 .ge. 10 )  then
          call report_stat('WARNING','FIXDRV','somake',' '
     .          , 'Ionsopheric Constraints in sestbl. wrong format',0)
         ierr_sestbl = -3
        endif 
        read( line(1:ichar1),'(f4.0)' )  akappa
        lchar = ichar3 - ichar2 + 1
        if( lchar .ge. 10 )  then
           call report_stat('WARNING','FIXDRV','somake',' '
     .          , 'Ionsopheric Constraints in sestbl. wrong format',0)
           ierr_sestbl = -3
        endif
        read( line(ichar2:ichar3), '(f4.0)' )  bkappa    
C       Check if consistent with type of solution
        if(     observ .eq. 'L1,L2_INDEPEND.'
     .   .or. observ .eq. 'L1_ONLY        '
     .   .or. observ .eq. 'L2_ONLY        '
     .   .or. observ .eq. 'L1_RECEIVER    ' ) then
          bkappa = 0.
          write(message,'(a15,a)') observ
     .           , ' Warning: Ionospheric constraint set to zero'
          call report_stat('WARNING','FIXDRV','somake',' ',message,0)
        endif
        write(20,445) akappa,bkappa   
c       Note: this format does not allow for constraints < 0.1 ppm (usually ok)
  445   format('    ionosphere constraint:   ',2f6.1)
C       wide-lane ambiguity resolution parameters
c       new defaults:
        wldev = 0.15
        wlsig = 0.15
        wlcut = 1000.
        wlrat = 10.
        wldmax = 500.
        reqd = .false.
        call rdsest( 23, 'Ambiguity resolution WL', 30, line
     .             , lsess, reqd, ill )
        if( ill.ne.0 ) ierr_sestbl = ill
        i0 = index( line, '=' ) + 1
        if( line.ne.blnk256 ) then
c         need to count arguments to avoid ugly abort for out-of-date sestbl.
          args = ' '
          args = line(i0:nblen(line))
          numarg= count_arg(line)
          if( numarg.eq.5 ) then
            read( line(i0:nblen(line)), *)
     .          wldev,wlsig,wlcut,wlrat,wldmax
          else
            read( line(i0:nblen(line)), *) wldev,wlsig,wlcut,wlrat
            write(message,'(2a)') 
     .           'Missing max BL in sestbl. Ambiguity resolution WL,'
     .         , ' set to default 500 km'
            call report_stat('WARNING','FIXDRV','somake',' ',message,0)
            wldmax = 500.
          endif
        endif
        write(20,450) wldev,wlsig,wlcut,wlrat,wldmax
  450   format
     .   ( '    wide lane ambiguity criteria: ',2f6.2,f8.1,f6.2,f9.1 )
C       wide-lane pseudorange ambiguity resolution parameters
c       conservative defaults:
        prdev = 0.05
        prsig = 0.05
        prcut = 1000.
        reqd = .false.
        call rdsest( 23, 'Ambiguity resolution PR', 30, line
     .             , lsess, reqd, ill )
        if( ill.ne.0 ) ierr_sestbl = ill
        i0 = index( line, '=' ) + 1
c       LINE will be blank if nothing found
        if( line.ne.blnk256 )
     .    read( line(i0:nblen(line)), *)  prdev,prsig,prcut
        write(20,460) prdev,prsig,prcut
  460   format( '    pseudorange ambiguity criteria: ',2f6.2,f8.1,f6.2 )  
      endif

      write(20,495)
  495 format(' exit set:')


c----------------------------------------------------------------------------

c  Solution Options

      write(20,500)
  500 format('*',/
     . '-------------- Part 5 -- Solution Options' )
   
c  Get threshhold for printing correlations
               
      correl_prt = 0.9999
      call rdsest( 17, 'Correlation print', 30, line
     .           , lsess, reqd, ill )
      if( ill.ne.0 ) ierr_sestbl = ill
      i0 = index( line, '=' ) + 1
c     LINE will be blank if nothing found
      if( line.ne.blnk256 )
     .  read( line(i0:nblen(line)), *)  correl_prt

c  Get threshold for updating L-file coordinates   
 
      coord_upd_tol = .3
      call rdsest( 16, 'Update tolerance', 30, line
     .           , lsess, reqd, ill )
      if( ill.ne.0 ) ierr_sestbl = ill
      i0 = index( line, '=' ) + 1
c     LINE will be blank if nothing found
      if( line.ne.blnk256 )
     .  read( line(i0:nblen(line)), *)  coord_upd_tol  

c  Write options for the tight, bias-free solution

      write(20,510) 
  510 format(' set tight_free solution option:')
      write(20,511) mfile,mfile,lfile 
  511 format(
     .  '    print out solution:    q-file ofile',/
     ., '    update file option:    m-file l-file g-file',/
     ., '    input_m file name:     ',a16,/
     ., '    output_m file name:    ',a16,/
     ., '    input_l file name:     ',a16)
c     create an updated L-file unless experiment choice is ORBIT-only
c       (the input lfile name is updated optionally in BMAKE)
      if ( exprmt.ne.'ORBIT' ) then
        call upnam1(lfile,lfile1)
        call lowers(lfile1)
        write(20,512) lfile1,coord_upd_tol
  512   format('    output_l file name:    ',a16,/
     .        ,'    coord_upd_tol:      ',f8.3 )
      endif                                
c     adjustment tolerance for updating L-file coordinates
c     --default is 30 cm; to get all, set = 0.
      
c     create an updated-G-file unless experiment choice is BASELINE
      if ( exprmt.ne.'BASELINE' ) then
        call upnam1 (gfile,gfile1)
        call lowers(gfile1)
        write(20,513) gfile,gfile1
  513   format('    input_g file name:     ',a16,/
     .         '    output_g file name:    ',a16 )
      endif
      write(20,514) correl_prt
  514   format('    correl_prt:            ',f8.6,/
     ., ' exit set:')


c  Write options for the tight, bias-fixed solution
      
      if( observ(1:7).ne.'LC_ONLY') then 
        write(20,520) 
  520   format(' set tight_fix solution option:') 
        write(20,511) mfile,mfile,lfile
c       Create an updated L-file unless experiment choice is ORBIT-only
c       (the input lfile name is updated optionally in BMAKE)
        if ( exprmt.ne.'ORBIT' ) then
          call upnam1(lfile,lfile1)
          call lowers(lfile1)
          write(20,512) lfile1,coord_upd_tol
        endif   
c       create an updated-G-file unless experiment choice is BASELINE
        if ( exprmt.ne.'BASELINE' ) then
         call upnam1 (gfile,gfile1)
         call lowers(gfile1)
         write(20,513) gfile,gfile1
        endif
        write(20,514) correl_prt
      endif

c  Write options for the loose, bias-free solution

      write(20,530)
  530 format(' set loose_free solution option:',/
     .,'    update file option: ',/
     .,' exit set:' )


c  Write options for the loose, bias-fixed solution
        
      if( observ(1:7).ne.'LC_ONLY') then 
        write(20,540)
  540   format(' set loose_fix solution option:',/
     .,'    print out solution: ',/
     .,' exit set:' )
      endif

      return
      end
