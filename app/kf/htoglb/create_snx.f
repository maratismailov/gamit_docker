 
CTITLE CREATE_SNX
 
      subroutine create_snx( np, cov_parm, sol_parm )
 
      implicit none

*     Routine to create and write SINEX file from the information stored from
*     a globk binary h-file.
*     Current version for Sinex 0.05
 
      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
 
 
* np     - Number of parameters in this solution
 
 
      integer*4 unit, np
 
* cov_parm(np,np) - Covarince matrix
* sol_parm(np)    - Solution vector
 
 
      real*8 cov_parm(np,np), sol_parm(np)
 
 
* LOCAL VARIABLES
 
* i,j      - Loop counters
* num_in_file = number of solutions in file (1 only for snx solns)
* ierr     - IOSTAT error
* trimlen  - Length of string
 
 
      integer*4 i,  ierr, trimlen
 
*  cyr, cdoy, csec - SNX creation start year, doy and ssecs
*  syr, sdoy, ssec - SNX data start year, doy and ssecs
*  eyr, edoy, esec - SNX data end year, doy and secs
 
 
      integer*4 syr, sdoy, ssec, eyr, edoy, esec, 
     .          cyr, cdoy, csec, date(5), unitc

* MOD TAH 130331: Added secondary comment file specific for the
*     snx being written (implemented with + concatenation).
*     Only +SNX to -SNX lines are copied 

      integer*4 units   ! Secondary comment file
      integer*4 pluspos  ! Position of plus sign 
      character*256 sec_comfile   ! File name
     
      real*8 jd, sectag

      character*4 sys_char
      character*128 line, snxfile
 
***** Start by creating the name of the sinex file to be output
      
      snxfile = hfile
      call gen_snxname(cepoch_expt, cepoch_start, cepoch_end,
     .                 cowner, glb_dir, snxfile )
      unitc = 101
      if( trimlen(snx_comfile).eq.0 ) snx_comfile = 'head.snx'

* MOD TAH 130331: See if we need to split comment file name
      pluspos = index(snx_comfile,'+')
      units = 0
      if( pluspos.gt.0 ) then !   Split the names
          sec_comfile = snx_comfile(pluspos+1:)
          if( pluspos.eq.1 ) then  ! Default head.snx + new name
              snx_comfile = 'head.snx'
          else
              snx_comfile = snx_comfile(1:pluspos-1)
          end if
*         Try to open file
          units = 105
          open(units, file=sec_comfile, iostat=ierr, status='old')
          call report_error('IOSTAT',ierr,'open',sec_comfile,0,
     .          'Create_snx/Secondary comment file')
          if( ierr.ne.0 ) units = 0
      end if

*     Try to open main sinex comments file
      open(unitc, file=snx_comfile, iostat=ierr, status='old')
      
      if( ierr.ne.0 ) then

*         Try to open version in HELP_DIR (use line for name)
          call getenv('HELP_DIR',line)
          line(trimlen(line)+1:) = '/' // snx_comfile
          open(unitc, file=line, iostat=ierr, status='old')
          if ( ierr.ne.0 ) unitc = 0
      end if
     
*
***** Tell user what is happening
      write(*,120) snxfile(1:trimlen(snxfile))
 120  format(' Creating SINEX file ',a)
      unit = 200
      open(unit, file=snxfile, iostat=ierr, status='unknown')
      call report_error('IOSTAT',ierr,'creat',snxfile,1,'create_snx')

****  Generate the first line of the file
      do i = 1,5
         date(i) = crun_time(i)
      end do
      sectag = crun_time(6)
      call ymdhms_to_jd( date, sectag, jd)
      call jd_to_yds( jd, cyr, cdoy, csec)
      call jd_to_yds( cepoch_start, syr, sdoy, ssec)
      call jd_to_yds( cepoch_end, eyr, edoy, esec)

****  Get the type of data in solution
      call sub_null(ccreator)
      call sub_null(canal_type)

* MOD TAH 070723: Check for sites with large sigmas and mark these are XPS
*     sites.
* MOD TAH 21013: Moved code to here since number of parameters could be
*     effected by sites being removed.
      call chksnx_site( np, cov_parm, sol_parm )
      
      call snx_sys_type( csys_type, sys_char)
      write(line, 140) cowner, mod(cyr,100), cdoy, csec, ccreator, 
     .                 mod(syr,100), sdoy, ssec, mod(eyr,100), edoy,
     .           esec, sys_char, qnum_parn, ccons_type, canal_type
 140  format('%=SNX 2.01 ',a3,1x,i2.2,':',i3.3,':',i5.5,
     .       1x,a3,1x,i2.2,':',i3.3,':',i5.5,1x,
     .       i2.2,':',i3.3,':',i5.5,1x,a1,1x,i5.5,1x,i1,1x,a8)
      write(*,'(a)') line(1:trimlen(line))
      write(unit,'(a)',iostat=ierr) line(1:trimlen(line))
      call report_error('IOSTAT',ierr,'SNC write',hfile,1,'create_snx')

****  Now read down the comments file echoing the +SNX -> -SNX lines
*     into the sinex file
      if ( units.eq.0 ) then   ! Get comments from main file, else
                               ! take from secondary
          call cp_comments(unit, unitc, 'SNX')
      else
          call cp_comments(unit, units, 'SNX')
      end if

*     Update the DOMES numbers

      call upd_domes(unitc)
      
****  Create the FILE block

      call cresnx_file( unit, unitc)
      call cresnx_input( unit, unitc )

      call cresnx_site( unit, unitc )
      call cresnx_sat ( unit, unitc )
      call cresnx_soln( unit, unitc, np, cov_parm, sol_parm )

      call cp_comments(unit, unitc, 'ENDSNX')
      write(unit,'(a)') '%ENDSNX'
       
***** Thats all 
      return
      end
      
