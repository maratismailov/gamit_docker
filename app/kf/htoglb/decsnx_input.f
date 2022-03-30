CTITLE DECSNX_INPUT
 
      subroutine decsnx_input (unit, line )
 
      implicit none

*     Routine to decode the input  blocks from SINEX:
 
      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
 
* PASSED VARIABLES
 
*   unit        - Unit number
 
      integer*4 unit
 
*   line        - Line read from input file
 
      character*(*) line
 
* LOCAL VARIABLES
 

*   indx        - Pointer in string
*   trimlen     - Length of string
 
      integer*4 indx, trimlen
 
*   block_found - True is block found
 
      logical block_found

****  Start decoding the types of site blocks
      block_found = .false.
 
      indx = index(line,'INPUT/ACKNOW')
      if( indx.gt.0 ) then
          block_found = .true.
          call snx_finbl(unit)
      end if

      indx = index(line,'INPUT/HIST')
      if( indx.gt.0 ) then
          block_found = .true.
          call dec_input_hist( unit, line  )
      end if
 
      indx = index(line,'INPUT/FILES')
      if( indx.gt.0 ) then
          block_found = .true.
          call dec_input_files( unit, line  )
      end if
 
 
****  Thats all the blocks so far
      if( .not.block_found ) then
          write(*,500) line(1:trimlen(line))
500       format('DECSNX_INPUT: Unknown block type',/a)
          call snx_finbl(unit)
      end if
 
****  Thats all
      return
      end
 
CTITLE DEC_INPUT_FILES  

      subroutine dec_input_files(unit, line )
 
      implicit none

*     Routine to read the input files block
 
      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
 
* PASSED VARIABLES
 
*   unit        - Unit number of sinex file
 
      integer*4 unit
 
*   line        - Line read from file
 
      character*(*) line
 
* LOCAL VARIABLES
 
*   ierr        - IOSTAT error
*   ns          - Solution number to which the entry applies
*   cr_yr, cr_doy, cr_sec - Creation yr, doy and seconds of day
 
      integer*4 ierr, ns, cr_yr, cr_doy, cr_sec, i, trimlen,
     .          jerr

*   cr_epoch    - Creation JD

      real*8 cr_epoch

*   filename   - Original data file name
*   descript   - Description of solution
*   creator    - Creator of file

      character*4  creator
      character*32 filename, descript
 
*   end_block   - Set true at end of block
 
      logical end_block
 
****  Start decoding
 
      end_block = .false.
      ierr = 0
      do i = 1, cnum_soln_recs
         qskalobs_file(i) = 'No_Entry_available'
         qsexpt_title(i)  = 'No_Entry_available'
      end do

      do while ( ierr.eq.0 .and. .not.end_block )
 
          read(unit,'(a)', iostat=ierr) line
          if( ierr.eq.0 .and. line(1:1).eq.' ' ) then

*****         Read the files line
              read(line,120, iostat=jerr) creator, cr_yr, cr_doy, 
     .                       cr_sec, filename,
     .                       descript
 120          format(1x, a3, 1x,i2,1x,i3,1x,i5,1x, a29,1x,a32)
              call yds_to_jd( cr_yr,cr_doy,cr_sec, cr_epoch)

*             Now find the corresponding history record
              ns = 0
              do i = 1, cnum_soln_recs
                 if( qscreator(i)(1:3).eq.creator(1:3) .and.
     .               qsrun_time(i).eq.cr_epoch ) then
                     if( ns.eq.0 ) then
                         ns = i
                         qskalobs_file(i) = filename
                         qsexpt_title(i) = descript
                     else
cd                        write(*,210) line(1:trimlen(line)), ns, i
cd 210                     format('**WARNING** Ambiguous INPUT/FILE',
cd     .                          ' entry',/,a,/,
cd     .                          ' Matches INPUT/HISTORY entries ',i3,
cd     .                          ' and ',i3)
                     end if
                 end if
              end do
              qexpt_title = qsexpt_title(cnum_soln_recs)
              
              if( ns.eq.0 ) then
                  write(*,220) line(1:trimlen(line))
 220              format('**WARNING**  No history entry found for ',
     .                   'INPUT/FILES entry:',a)
              end if

          end if
          if( line(1:1).eq.'-' .or. ierr.ne.0 ) then
              end_block = .true.
          end if
      end do
 
****  Thats all
      return
      end
 
CTITLE DEC_INPUT_HIST

      subroutine dec_input_hist(unit, line )
 
      implicit none

*     Routine to read the input/history block
 
      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
 
* PASSED VARIABLES
 
*   unit        - Unit number of sinex file
 
      integer*4 unit
 
*   line        - Line read from file
 
      character*(*) line
 
* LOCAL VARIABLES
 
*   ierr        - IOSTAT error
 
      integer*4 ierr
 
*   end_block   - Set true at end of block
 
      logical end_block

*   ver         - Version of program 

      real*8 ver

*   st_yr, st_doy, st_sec - Start of data yr, doy and seconds of day
*   en_yr, en_doy, en_sec - end of data yr, doy and seconds of day
*   cr_yr, cr_doy, cr_sec - Creation yr, doy and seconds of day
*   np          - Number of parameters in solution
*   ns          - Short version of solution number
*   cons_type   - Constraint type

      integer*4 st_yr, st_doy, st_sec, en_yr, en_doy, en_sec, np, ns,
     .          cons_type , cr_yr, cr_doy, cr_sec, jerr

*   type        - Type of input (+added; = combined)
*   prog        - File type for input
*   ver         - Version of file type.
*   creator     - Who created file
*   owner       - Who owned the filed
*   sys_type    - System type P=GPS, R=VLBI, S=SLR
*   anal_type   - Types of parameters in solution

      character*4 type, prog, creator, owner, sys_type
      character*8 anal_type

****  Read down each line saving the information
      cnum_soln_recs = 0
      cnum_comb = 0
      end_block = .false.
      ierr = 0
      ns = 0
      do while ( ierr.eq.0 .and. .not.end_block )
 
          read(unit,'(a)', iostat=ierr) line
          if( ierr.eq.0 .and. line(1:1).eq.' ' ) then
              read(line,120,iostat=jerr) type,prog, ver, 
     .                       creator, cr_yr,cr_doy, 
     .                       cr_sec, owner,st_yr, st_doy, st_sec,
     .                       en_yr, en_doy, en_sec, sys_type, np, 
     .                       cons_type, anal_type
 120          format(1x,a1,a3,f5.2, 1x,a3,1x,i2,1x,i3,1x,i5,1x,
     .               a3,1x,i2,1x,i3,1x,i5,1x,i2,1x,i3,1x,i5,1x,
     .               a1,1x,i5,1x,i1,1x,a5)

*             See if this in input solution record or if it was a combined
*             solution.
              cnum_soln_recs = cnum_soln_recs + 1
* MOD TAH 020628: Check that we don't exceed limit
              if( cnum_soln_recs.gt.max_sln_save ) then
*                 Too many solution records, tell user what to do
                  write(*,140)  max_sln_save
 140              format(/,'**DISASTER** Number of saved solutions ',
     .               ' exceeds limit of ',i6,/,
     .               '             Modify the parameter max_sln_save ',
     .               ' in gg/kf/htoglb/htoglb_comm.h to at least ',
     .               ' number of solution records in SNX file',/,
     .               '             Use make to re-make the glbtosnx.')
                  call report_stat('FATAL','htoglb','dec_input_hist',
     .                            ' ','Too many solution records',0)              
              end if

              if( type(1:1).eq.'=' ) qnum_comb = qnum_comb + 1

*             Save the information.  It will later be written out as soln_def
*             records.
              ns = cnum_soln_recs
              qsprog_gen(ns) = type(1:1) // prog
              qsglb_ver(ns) = ver*100
              qscreator(ns)  = creator
              qsowner(ns)    = owner
              call yds_to_jd( cr_yr,cr_doy,cr_sec, qsrun_time(ns))
              call yds_to_jd( st_yr,st_doy,st_sec, qsepoch_start(ns))
              call yds_to_jd( en_yr,en_doy,en_sec, qsepoch_end(ns))
*
*             Set the data type for this solution
              if( sys_type(1:1).eq.'P') qssys_type(ns) = 2
              if( sys_type(1:1).eq.'R') qssys_type(ns) = 1
              if( sys_type(1:1).eq.'L') qssys_type(ns) = 4
              if( sys_type(1:1).eq.'M') qssys_type(ns) = 8
*             This is combined and really don't know so we set all main
*             bits
              if( sys_type(1:1).eq.'C') qssys_type(ns) = 7

              qsnum_parn(ns) = np
              qscons_type(ns) = cons_type
              qsanal_type(ns) = anal_type

          end if
          if( line(1:1).eq.'-' .or. ierr.ne.0 ) then
              end_block = .true.
          end if
      end do
      
      qnum_soln_recs = cnum_soln_recs
      qnum_comb = cnum_comb
      
****  Thats all
      return
      end
