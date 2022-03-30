
      program glbtosnx

      implicit none
     
*     This program will read Globk ver 1.0 and greater binary
*     hfiles and write out sinex files.
*
*     Runstring is:
*     glbtosnx <dir> <command file> <input binary hfile> <output file name>

      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
      include '../includes/sinf_def.h'
*
*   vma_data(max_vma_space) - Array to hold the covariance
*           - matrix and solution vector. (Matrix is contained
*           - as a full square matrix)
 
      integer*4 vma_data(1)
      real*8 vma_r8(1)
      
      integer*4 kerr
      integer*8 icov_parm, isol_parm  
      
      equivalence ( vma_r8, vma_data )  
 
****  Decode the runstring   
 
      call decode_gts_runs
      
***** Now read in the binary hfile

      call fmpopen(cglb_dcb,kerr,glb_file,'ro',0)
      call report_error('FMPopen',kerr,'open',glb_file,1,
     .                  'glbtosnx')
     
*     Read the header of the file     
      call rw_glb_header('R', kerr)

* MOD TAH: 970610: If kerr is 2001, then file needs to be swapped
      if( kerr.ne.0 ) then
          if( kerr.eq.2001 ) then
              write(*,120) 
 120          format('**WRONG BYTE ORDER** Run swaph') 
              call report_error('RW_GLB_HEADR',kerr,'read',glb_file,
     .                          1,'glbtosnx') 
          else
              call report_error('RW_GLB_HEADR',kerr,'read',glb_file,
     .                          1,'glbtosnx') 
          end if
      end if 
      
*     Assign the memory for reading the covariance matrix
      call gts_mem(vma_data(1), icov_parm, isol_parm)
      
*     Now read all the binary file contents
      call read_glb(vma_data(icov_parm), vma_data(isol_parm))

*     Now create and write the sinex file
      call create_snx(cnum_parn, vma_data(icov_parm),
     .                vma_data(isol_parm))
      
*     Thats all
      end
      
CTITLE DECODE_GTS_RUNS

      subroutine decode_gts_runs
      
      implicit none

*     Routine to read the runstring for glbtosnx.
      
      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'

* LOCAL VARIABLES
*   rcpar   - Reads runstring
*   len_run - Length of runstring

      integer*4 rcpar, len_run, nr
      integer*4 jerr   ! IOSTAT error
      logical finished
      
*   runstring - Runstring read from file

      character*256 runstring


***** Get the directory for the output files
      len_run = rcpar(1, glb_dir )
      if( len_run.le.0 ) then 
          call proper_runstring('glbtosnx.hlp','glbtosnx/glb_dir',1)
      end if
      
*     Now loop over the command_file/input commands
      finished = .false.
      Var_Limit = 20.d0  ! Initialize variance limit Same value 
                         ! 20 m^2 in <=2.08 version)
      nr = 1
      glbtosnx_opt = 0
      snx_comfile = ' '
      igs_ptname = .false.
      do while ( .not.finished)
         nr = nr + 1
         len_run = rcpar(nr, runstring)
         if( len_run.gt.0 ) then

*            For the moment ignore.  Later will be options
*            command file.  
* MOD TAH 210112: Allowed upper and lower case versions..
             if( runstring(1:2).eq.'-h' .or. 
     .           runstring(1:2).eq.'-H' ) then
                 call sbit(glbtosnx_opt,1,1)
             elseif (runstring(1:2).eq.'-s' .or.
     .               runstring(1:2).eq.'-S' ) then
                 igs_ptname = .true.
             elseif (runstring(1:3).eq.'-V=' .or.
     .               runstring(1:3).eq.'-v=') then
                 read(runstring(4:),*,iostat=jerr) Var_Limit
                 call report_error('IOSTAT',jerr,'decoding',
     .               runstring,0,'GLBTOSNX/DECODE_GTS_RUNS')
                 if( jerr.ne.0 ) Var_Limit = 20.d0 
             else
                 snx_comfile = runstring
                 finished = .true.
                 nr = nr + 1
             end if
         else
             finished = .true.
             nr = nr + 1
         end if
      end do
      
*     Now get the input file
      len_run = rcpar(nr, glb_file)
      if( len_run.eq.0 ) 
     .    call proper_runstring('glbtosnx.hlp','glbtosnx/glb_file',1)
     
*     Now see if name of snxfile or code passed/
      nr = nr + 1
      hfile = ' '
      len_run = rcpar(nr, hfile)
      
***** Thats all
      return
      end
      
CTITLE GTS_MEM
        
      subroutine gts_mem(vma_data, icov_parm, isol_parm)     

      implicit none

*     Routine to assign the memory for the htosnx run

      include '../includes/kalman_param.h'
      include '../includes/glb_hdr_def.h'
      include 'htoglb_comm.h'
      
* PASSED VARIABLES    

* vma_data -- Array to hold data
* icov_parm, isol_parm -- Starting positions for arrays

      integer*4 vma_data(*)
      integer*8 icov_parm, isol_parm

* LOCAL variables
  
*   vma_i4wrd       - Number of bytes needed (I*8)
*   memassign       - Integer*8 function to return the location in
*                     vma_data array where there is space available
 
 
      integer*8 vma_i8wrd
      integer*8 memassign8
    
* MOD TAH 190511: Introduced when number of parameters > 32767.
      integer*8 i8  !  i8 value of 1 to force I8 calculations


****  Start by computing the number of I*4 words needed
      i8 = 1
     
      vma_i8wrd = (cnum_parn*(cnum_parn+i8))*2 + 128

* MOD TAH 980209: Set minimum number of bytes to allow all
*     the apriori orbit parameters to be read.
      if( vma_i8wrd.lt.4096 ) vma_i8wrd = 4096

* MOD TAH 040801: Replaced with memassign which is 64-bit compatable.            
      istart_vma = memassign8(vma_i8wrd,1, loc(vma_data)) 
      
      if( istart_vma.eq.0 ) then
          write(*,120) vma_i8wrd/(1024.0**2)*4
 120      format('*** DISASTER *** Not enough memory.  Try a larger',
     .                ' computer ',F8.2,' Mbytes needed')
          stop 'GLBTOSNX: Not enough memory'
      end if
      
      icov_parm = istart_vma
      isol_parm = icov_parm + (i8*cnum_parn)*cnum_parn*2
      
****  Thats all
      return
      end
