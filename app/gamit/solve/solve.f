Copyright 1995 Massachusetts Institute of Technology and The Regents of the
c              University of California, San Diego.  All Rights Reserved.

      PROGRAM SOLVE

c GAMIT least square estimation package.

c Principal authors:  Y. Bock, D. Dong, R. King

c     catalog of files:
c        C-file : data file               binary   unit: 31-130   input
c        Q-file : solution (verbose)      ascii    unit: 10       output
c        M-file : initial merging file    binary   unit: 11       input
c        M-file : updated merging file    binary   unit: 12       output
c        N-file : autcln noise/bias file  ascii    unit: 13       input
c        O-file : solution (terse)        ascii    unit: 15       output
c        I-file : initial site file       ascii    unit: 16       input
c        I-file : updated site file       ascii    unit: 18       output
c        G-file : initial orbit parameter ascii    unit: 16       input
c        G-file : updated orbit parameter ascii    unit: 17       output
c        T-file : tabulated orbit file    binary   unit:          input
c        H-file : covariance matrices     ascii    unit: 21       output
c
c     auxillary units:
c        unit  5 :   terminal                                     input
c        unit  6 :   screen                                       output
c        unit 26 :   for debugging                                output
c        unit 27 :   total normal matrix storage                  scratch
c        unit 28 :   LC mode normal matrix storage (wl free)      scratch
c        unit 29 :   LC mode normal matrix storage (wl fixed)     scratch
c
c     options after obtaining one solution:
c        ambiguity searching
c        after L1, L2 separate mode, run LC mode
c        after LC mode tight constraints, run LC mode loose
c            constraints
c        exit program
c
c     old options no longer supported:
c        old-style batch files
c        interactive input
c        single differences
c        repeat solution with different parameters
c        repeat solution with new data
c
c     Program flow is controlled by the input choice of observable, which is
c     used to set the internal flag 'l2flag' in common /flags/:
c
c         L2_ONLY          l2flag =-2
c         L1_ONLY          l2flag =-1
c         L1_RECEIVER      l2flag = 0
c         LC_ONLY          l2flag = 1
c         L1&L2            l2flag = 2  (with ion constraint)
c         LC_HELP          l2flag = 3
c         LC_AUTCLN        l2flag = 4
c         L1L2_INDEPENDENT l2flag = 5
c
c     Flag 'lquick' in /flags/ controls whether regular solution with explicit biases (=0),
c     quick solution (=1) or regular solution with implicit biases (=2).  The last (newest)
c     option is just quick but with an extra bias parameter only if the data point is
c     explicit flagged
c
c     Solution characteristics and printout for each of these options are set using

c          char*5  constraints :  'tight'  or  'loose'
c          char*4  free_fix    :  'free'   or  'fixd'
c          char*4  phase_obs   :  'LC'  'L1'  'L2'   'L1L2'  'L12I'

      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'
      include 'models.h'
      include 'parameters.h'
c      
      logical fcheck 
      integer ioerr,trimlen
      character*1 lowerc,nbatch
      character*4 free_fix,phase_obs
      character*5 constraints      
                         
c     define some constants
      pi    = 4.d0*datan(1.d0)
      convd = pi/180.d0  
c     These are WGS84 values from gdetic.dat
      semi = 6378.137d0 
      finv = 298.257223563d0 

c
c==========================================================================
                                                       
c---Stop if a previous step has failed
  
      if( fcheck('GAMIT.fatal') )
     .  call report_stat('FATAL','SOLVE','solve',' '
     .                  ,'GAMIT.fatal exists: SOLVE not executed',0)
  
 

c---Write to screen and status file the program name and version

      call sversn(vers)

c---Check for an incorrect batch file

      read(5,'(a1)') nbatch
      if(lowerc(nbatch).eq.'b') call report_stat('FATAL','SOLVE','solve'
     .                                      ,' ','Oldstyle batchfile',0)
      rewind 5
c     open unit 26 for possible debug
c     open(26,file=lowerc('lsqx.prt'))

c---Get input options and scratch directory, open the Q-file, read the M-file,
c   set parameter controls, and print out the analysis options

      constraints = 'tight'
      free_fix = 'free'    
      call read_bfl ( 1,constraints,free_fix )    
                        
c --Open normal matrix scratch files
                    
c     file names generated from q-file to make the scratch files unique 
c     in the tmp directory when running multiple jobs 
c     change by Robert Pickle to make the file names unique across years:
cc      ftmp27 = 'tmp_' // qfiln(1:10) // '.27'
cc      ftmp28 = 'tmp_' // qfiln(1:10) // '.28'
cc      ftmp29 = 'tmp_' // qfiln(1:10) // '.29'      
      ftmp27 = 'tmp_' // hfiln(1:12) // '.27'
      ftmp28 = 'tmp_' // hfiln(1:12) // '.28'
      ftmp29 = 'tmp_' // hfiln(1:12) // '.29'

c     see if a local disk rather than the processing directory to be used
      if( scratch_dir(1:1).ne.' ' ) then   
        ftmp27 = scratch_dir(1:trimlen(scratch_dir)) // '/' // ftmp27
        ftmp28 = scratch_dir(1:trimlen(scratch_dir)) // '/' //ftmp28
        ftmp29 = scratch_dir(1:trimlen(scratch_dir)) // '/' //ftmp29   
      endif 
      open(unit=27,file=ftmp27, form='unformatted',status='unknown'
     .    , iostat=ioerr)     
      if( ioerr.ne.0 ) call report_stat('FATAL','SOLVE','solve',ftmp27
     .     ,'Error opening scratch unit',ioerr)
      open(unit=28,file=ftmp28, form='unformatted',status='unknown'
     .     ,iostat=ioerr)       
      if( ioerr.ne.0 ) call report_stat('FATAL','SOLVE','solve',ftmp28
     .     ,'Error opening scratch unit',ioerr)
      open(unit=29,file=ftmp29, form='unformatted',status='unknown'
     .     ,iostat=ioerr)  
      if( ioerr.ne.0 ) call report_stat('FATAL','SOLVE','solve',ftmp29
     .     ,'Error opening scratch unit',ioerr)   

c---Read in data and perform initial biases-free solution
                              
      call lsquar     

c---Now branch according to choice of observable

c  ** Single-frequency with ambiguity resolution
      if( l2flag.le.0) then

c       read output/update controls and output biases-free solution
        if( l2flag.eq.0 .or. l2flag.eq.-1) phase_obs = 'L1  '
        if( l2flag.eq.-2) phase_obs = 'L2  '
        call read_bfl ( 2,constraints,free_fix )
        call write_soln (constraints,free_fix,phase_obs)
        if(l2flag.le.0) call bisopt (1)
        free_fix = 'fixd'
c       output the biases-fixed solution
        call read_bfl ( 2,constraints,free_fix )
        call write_soln (constraints,free_fix,phase_obs)

c       do the loose solutions : first biases free, then fixed
c       this not yet tested                                   
        if( do_loose ) then
          constraints = 'loose'
          free_fix = 'free'
          call loos12( 3 )
          call read_bfl ( 2,constraints,free_fix )
          call write_soln (constraints,free_fix,phase_obs)
          free_fix = 'fixd'
          call loos12( 4 )
          call read_bfl ( 2,constraints,free_fix )
          call write_soln (constraints,free_fix,phase_obs)
        endif

c  ** LC with no ambiguity resolution

      else if( l2flag.eq.1 ) then

c       read output/update controls and output constrained the solution
        call read_bfl ( 2,constraints,free_fix ) 
        phase_obs = 'LC  '
cd        print *,'SOLVE calling WRITE_SOLN '
cd         stop 1 
        call write_soln (constraints,free_fix,phase_obs)
cd        print *,'SOLVE after WRITE_SOLN do_loose ',do_loose 
c       do the loose solution
        if( do_loose ) then
          constraints = 'loose'
          free_fix = 'free'
cd          print *,'SOLVE calling LCLOOS '
          call lcloos( 3 )               
cd          print *,'SOLVE aft LCLOOS sclerr ',sclerr 
          call read_bfl ( 2,constraints,free_fix )  
cd          print *,'SOLVE aft READ_BFL 2 sclerr ',sclerr
          call write_soln (constraints,free_fix,phase_obs)
        endif


c **  L1 & L2 with ionospheric constraint

      else if( l2flag.eq.2 ) then

c       read output/update controls and output biases-free L1/L2 solution
        phase_obs = 'L1L2'
        call read_bfl ( 2,constraints,free_fix )
        call write_soln (constraints,free_fix,phase_obs)

c       Here we resolve widelanes
        call bisopt( 2)
        free_fix = 'fixd'
c       output the biases-fixed L1/L2 solution
        call read_bfl ( 2,constraints,free_fix )
        call write_soln (constraints,free_fix,phase_obs)

c       do LC after L1/L2 separate 
        call lc_solution
c       output the biases-free LC solution
        free_fix = 'free'
        phase_obs = 'LC  '
        call read_bfl ( 2,constraints,free_fix )
        call write_soln (constraints,free_fix,phase_obs)
c       resolve the narrowlane biases using LC
        call bisopt(1)
c       output the biases-fixed LC solution
        free_fix = 'fixd'
        call read_bfl ( 2,constraints,free_fix )
        call write_soln (constraints,free_fix,phase_obs)

c       do the L1/L2 loose solutions : first biases free, then fixed
        constraints = 'loose'
        if( do_loose ) then
          free_fix = 'free'
          phase_obs = 'L1L2'
          call loos12( 3 )
          call read_bfl ( 2,constraints,free_fix )
          call write_soln (constraints,free_fix,phase_obs)
          free_fix = 'fixd'
          call loos12( 4 )
          call read_bfl ( 2,constraints,free_fix )
          call write_soln (constraints,free_fix,phase_obs)
        endif


c **  LC_HELP or LC_AUTCLN

      else if ( l2flag.eq.3 .or. l2flag.eq.4 ) then
                 
        phase_obs = 'LC  '
c       get the LC, biases-free solution
        call lc_solution       
c       get the wide-lane ambiguities from AUTCLN or try to resolve them
        call get_widelane                        
c       read output/update controls and output the LC-free solution
        call read_bfl ( 2,constraints,free_fix ) 
        call write_soln (constraints,free_fix,phase_obs)
c       now try to resolve the narrow-lane biases   
c**      rwk sequence
        call get_narrowlane( 1 ) 
c       output the biases-fixed solution
        free_fix = 'fixd'
        call read_bfl ( 2,constraints,free_fix )
        call write_soln (constraints,free_fix,phase_obs)

c     do the loose solutions : first biases free, then fixed
             
        if( do_loose ) then
          constraints = 'loose'
          free_fix = 'free'
          call lcloos( 3 )
          call read_bfl ( 2,constraints,free_fix )
          call write_soln (constraints,free_fix,phase_obs)
          free_fix = 'fixd'
          call lcloos( 4 )
          call read_bfl ( 2,constraints,free_fix )
          call write_soln (constraints,free_fix,phase_obs)
        endif

c  ** L1,L2 independent solutions

      elseif(l2flag.eq.5) then

c       read output/update controls and output biases-free solution
        phase_obs = 'L12I'
        call read_bfl ( 2,constraints,free_fix )
        call write_soln (constraints,free_fix,phase_obs)
        call bisopt (3)
        free_fix = 'fixd'
c       output the biases-fixed solution
        call read_bfl ( 2,constraints,free_fix )
        call write_soln (constraints,free_fix,phase_obs)

c       do the loose solutions : first biases free, then fixed
c       this not yet tested
        if( do_loose ) then
          constraints = 'loose'
          free_fix = 'free'
          call loos12( 3 )
          call read_bfl ( 2,constraints,free_fix )
          call write_soln (constraints,free_fix,phase_obs)
          free_fix = 'fixd'
          call loos12( 4 )
          call read_bfl ( 2,constraints,free_fix )
          call write_soln (constraints,free_fix,phase_obs)
        endif

      else
         call report_stat('FATAL','SOLVE','solve',' '
     .     ,'l2flag ambiguous',0)
      endif

c---- end of run - close output files

c     close debug file
c     close(unit = 26)
      if( logprt) write(6, '(a)') 'Normal stop in SOLVE'
      write(10,'(a)') 'Normal stop in SOLVE'
      close(unit = 10)

* MOD TAH 971111: Close the temp files with delete option.
      close(27,status='delete')
      close(28,status='delete')
      close(29,status='delete')
      if(ioflag.eq.1) close(unit = 15)
c**   clear IEEE underflow and imprecise exceptions
c**   i = ieee_flags('clear','all','',out)

      call report_stat('STATUS','SOLVE','solve',' ', 'Normal stop',0)
      stop
      end
