 
      program htoglb

      implicit none
 
*     This program will convert an ascii Q file output from the
*     GPS propgram SOLVE to a global file capable of being read by
*     GLOBK.
*
*     In this first version, only station coordinates and satellite
*     orbit information are copied to the global files.  All other
*     parameters are ignored.
*
*     Runstring:
*     % htoglb <dir> [inputs]
*     Where <dir> is the name of the directory for output, and
*           [inputs ..] are the names of the input files.
*
*     The output files are created as
*     dir/hyymmddhhX.glb
*     where X is the last letter of the hfile name before the doy
*             extent.
*     If there are more than one covariance matrix in the q file the
*     the sequence goes .glb .glc .gld . .. etc.
*
* MOD TAH 920812: Added features which allow htoglb to read h-files
*     generated from solvem in fonda.  (These are terrestrial data
*     types typcially).
* MOD TAH 930208: Added feature to read the new h-file type entry from
*     the keys: line of the h-file.  Only implememted for GPS h-files.
* MOD TAH 940517: Added polar motion and UT1 parameter estimates to the
*     allowed parameters.
* MOD TAH 941005: Version 3.21 introduced that treats bug in 9.28 of
*     solve; and allows JPL satcov files to be read.
* MOD TAH 950206: Version 3.24 Added GSFC VLBI covariance files, UT SLR
*     files (no covarinance) and upgraded to read 9.3 Solve hfiles.
* MOD TAH 950729: Updated to handle additional radiation parameters in
*     Gamit H-files.
* MOD TAH 980122: Added dynamic memory allocation
* MOD TAH 110202: Moved the opening of the ephemeris file so that it does
*     not keep getting re-opened.  (Problem with some compilers)
* MOD TAH 200723: Added rm_param option for normal equations.

      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
 
*   ierr    - IOSTAT error on file open
*   rcpar   - Reads runstring
*   nq      - NUmber of q file being processed
*   nt      - Actual number of file (excludes count on options)
*   len_run - Length of runstring
*   num_in_file - NUmber of covariance matrix in file
*   trimlen - Length of string
*   ephem_unit  - Unit number for ephemeris file
*   ephem_err   - IOSTAT error opening ephemeris file.
 
      integer*4 ierr, rcpar, nq, len_run, num_in_file, trimlen
      integer*4 ephem_unit, ephem_err, indx, nt

 
*
*   vma_data(1) - Array to hold the covariance
*           - matrix and solution vector. (Matrix is contained
*           - as a full square matrix)
*           - Declare as only 1, set with malloc later.
 
      integer*4 vma_data(1)

      character*12 hsver, cdum
 
 
      common / vma_area / vma_data
 
****  Start, decode the directory name for the output

      write(*,50) hsver(htoglb_version)
 50   format(/' ++++++++++++++++++++++++++++++++++++++',
     .       /' + HTOGLB                 Version ',a,
     .       /' ++++++++++++++++++++++++++++++++++++++',/)
 
      len_run = rcpar(1, glb_dir )
      if( len_run.le.0 ) then
          call proper_runstring('htoglb.hlp','htoglb', 1)
      end if
      if( glb_dir(1:len_run).eq.'-' ) glb_dir = ' '

      len_run = rcpar(2, ephem_file )

*     Now open empheris file (this is appended
*     to any existing ones
      if( trimlen(ephem_file).gt.0 ) then
          ephem_unit = 200
          call open_lu( ephem_unit, ephem_file, ephem_err,
     .                  'append')
          call report_error('IOSTAT',ephem_err,'open', 
     .                      ephem_file, 0,'htoglb')
      else
          ephem_err = -1
      endif
      len_run = 1   ! Set to 1 to make sure loop activates below
                    ! even when no ephem_file is given (TAH 110208)

*     Now loop over the inputs in the runstring and make the
*     global files


      call report_stat('status','htoglb','main',' ','Started',0)
 
      nq = 0
      nt = 0
      usr_vma_space = 2*max_vma_space
      comp_eigen = .false.
      igs_ptname = .false.
      expt_code = ' '
      auto = .true.
      decon_str = ' '
      gamit_str = 'RT'
      sngday = .true.
      rm_param = .false.
      rm_site = .true.    ! When sites are used, set to false.
      
      do while ( len_run.gt.0 )
          nq = nq + 1
          len_run = rcpar(nq+2, hfile )
          
*         Check to see if constraint level passed
          if( hfile(1:3).eq.'-C=' .or. hfile(1:3).eq.'-c=' ) then
              indx = 4
              call read_line(hfile,indx,'R8',ierr,snx_constraint,
     .                       cdum)
              call report_error('IOSTAT',ierr,'decod',hfile,0,
     .                          'htoglb')
              write(*,80) snx_constraint
 80           format('SINEX Constraint set to ',f12.4,' m')

          else if( hfile(1:2).eq.'-e' .or. hfile(1:2).eq.'-E' ) then
              comp_eigen = .false.

          else if( hfile(1:2).eq.'-v' .or. hfile(1:2).eq.'-V' ) then
              comp_eigen = .true.

          else if( hfile(1:2).eq.'-s' .or. hfile(1:2).eq.'-S' ) then
              igs_ptname = .true.

          else if( hfile(1:3).eq.'-N=' .or. hfile(1:3).eq.'-n=' ) then
*             Decode the style of name desired
              name_format = hfile(4:)
              call casefold(name_format)

* MOD TAH 200504: Changes string test from 1:3 to 1.2 (former only -a
*         option worked.
          else if( hfile(1:2).eq.'-a' .or. hfile(1:2).eq.'-A' ) then
*             Turn off automatic update of BASELINE mode GAMIT hiles
              auto = .false.
              gamit_str = ''   ! Set the decon_str (could also be
*                                turned off with -d=X)

          else if( hfile(1:3).eq.'-F=' .or. hfile(1:3).eq.'-f=' ) then
*             Decode the style of name desired
              expt_code = hfile(4:)

*         See if the rotation/translation deconstrained is passed.
*         OPtions are R - rotation (10mas added), T - translation (1m added)
          else if( hfile(1:3).eq.'-D=' .or. hfile(1:3).eq.'-d=' ) then
              decon_str = hfile(4:)
              call casefold(decon_str)
              gamit_str = decon_str     ! Allows control TR options
 

          else if( hfile(1:3).eq.'-M=' .or. hfile(1:3).eq.'-m=' ) then

*             Memory being set.  See if not already done.
              if( nt.gt.0 ) then
                  write(*,85)
 85               format('**WARNING** Memory can only be set before'
     .                  ,' first file read')
              else
                  indx = 4
                  call read_line(hfile,indx,'I4',ierr,usr_vma_space,
     .                       cdum)
                  call report_error('IOSTAT',ierr,'Mem decod',hfile,0,
     .                          'htoglb')
*                 Convert mbytes to I*4 words     
                  usr_vma_space = usr_vma_space*(1024.)**2/4
              end if
* MOD TAH 200727: See if -MD/-md option passed.
          else if( hfile(1:3).eq.'-MD' .or. hfile(1:3).eq.'-md' ) then
               sngday = .false.   ! Keep all days in CODE SINEX files.
               write(*,'(a)' ) 'Keeping CODE multiday positions'
              
          else if( len_run.gt.0 ) then 
          
*             Assign memory if this is the first file
              if( nt.eq.0 ) then
                  write(*,90) usr_vma_space*4/(1024.)**2
 90               Format(' Allocating ',F8.3,' Mbytes for run')
                  call htog_mem( vma_data, usr_vma_space, istart_vma )
              end if        
 
*             If we have an entry procress
              call report_stat('status','htoglb','main',hfile,
     .                         'Processing file ',0)
              nt = nt + 1
              write(*,100) nt, hfile(1:len_run)
 100          format(/,64('-'),/,' Processing file ',i3,' h-file ',a)
 
*             See if we can open
              open(100, file=hfile, status='old', iostat=ierr)
              call report_error('IOSTAT',ierr,'open',hfile,0,
     .                        'htoglb')
 
*             If open OK, continue.  Otherwize get next file
              if( ierr.eq.0 ) then

                  call init_hfr

                  call make_glb(100, num_in_file, vma_data,
     .                     usr_vma_space )
                  write(*,150) num_in_file, hfile(1:len_run)
 150              format(2x,i3,' Solutions extracted from ',a)

                  if( ephem_err.eq.0 ) then
                      call write_ephem( ephem_unit )
                  end if

              end if
*                     ! Runstring has length
          end if
*                     ! Looping over the runstring
      end do
 
****  Thats all
      end

CTITLE htog_mem
 
      subroutine htog_mem( vma_data, usr_vma_space, istart_vma )
 
      implicit none

*     Routine to assign the memory for the htoglb  run.
  
* PASSED variables
 
*   vma_data(1)     - Location in memory where we
*                   - start the addressing from
*   usr_vma_space   - Number of I*4 words to allocate
 
 
      integer*4 vma_data(1), usr_vma_space
 
* Local Variables
 
*   vma_i4wrd       - Numner of bytes needed for this
*                   - run.
*   vma_start       - Start of the vma_data array.
*   mallocg          - Assigns memory
*   memloc          - Address that starts a big enough area
*                   - of memory
*   memassign       - Integer*8 function to return the location in
*                     vma_data array where there is space available
 
      integer*4 vma_i4wrd

      integer*8 istart_vma, memassign


*     Allocate the memory for the run.  Number of I*4 words is 
*     passed 
      vma_i4wrd = usr_vma_space
c      memloc = mallocg(vma_i4wrd)
* MOD 040801: Replaced mallocg call with memassign that is 
*     64-bit compatable.
      istart_vma = memassign( vma_i4wrd, 1, loc(vma_data)) 

      if( istart_vma.eq.0 ) then
          write(*,120) vma_i4wrd/(1024.0**2)*4
 120      format('*** DISASTER *** Not enough memory.  Try a larger',
     .                ' computer ',F8.2,' Mbytes needed')
          stop 'HTOGLB: Not enough memory'
      end if
 
****  Thats all
      return
      end
 
 
 

