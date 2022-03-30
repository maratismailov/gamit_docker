CTITLE CRESNX_INPUT
 
      subroutine cresnx_input (unit, unitc )
 
      implicit none

*     Routine to create the input  blocks from SINEX:
 
      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
 
* PASSED VARIABLES
 
*   unit        - Unit number
*   unitc       - Commnts file unit number
 
      integer*4 unit, unitc

 
* LOCAL VARIABLES
 

 
****  Copy over the input/ACKNOWLEDGMENTS lines
      write(unit,'(a)') '+INPUT/ACKNOWLEDGMENTS'
      call cp_comments(unit, unitc, 'INPUT/ACKNOWLEDGMENTS')
      write(unit,'(a)') '-INPUT/ACKNOWLEDGMENTS'

****  Now do the remaining blocks
      call cre_input_hist( unit, unitc )
      call cre_input_file( unit, unitc )
 
 
****  Thats all
      return
      end
 
CTITLE CRE_INPUT_FILE  

      subroutine cre_input_file(unit, unitc )
 
      implicit none

*     Routine to read the input files block
 
      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
 
* PASSED VARIABLES
 
*   unit        - Unit number of sinex file
*   unitc       - Commnts file unit number
 
      integer*4 unit, unitc
  
 
* LOCAL VARIABLES
 
*   cr_yr, cr_doy, cr_sec - Creation yr, doy and seconds of day
*   sf, ef, lf  - Start and end position for part of Kalobs
*                 file names to be printed and the length of
*                 name.
 
      integer*4  cr_yr, cr_doy, cr_sec, i, sf, ef, lf, trimlen
 
****  Start creating

      write(unit,'(a)' ) '+INPUT/FILES'
      call cp_comments(unit, unitc, 'INPUT/FILES')
 
      do i = 1, cnum_soln_recs 
 
          call jd_to_yds( qsrun_time(i), cr_yr,cr_doy,cr_sec)
         
          call sub_null( qskalobs_file(i) ) 
*         Check the length of the file name
          lf = trimlen(qskalobs_file(i))
          if( lf.gt.29 ) then
              ef = lf 
              sf = lf - 28
          else
              ef = 29
              sf = 1
          end if
          call sub_null( qsexpt_title(i) ) 
          call sub_null( qscreator(i) ) 

          write(unit, 120) qscreator(i), mod(cr_yr,100),cr_doy,
     .            cr_sec, qskalobs_file(i)(sf:ef), qsexpt_title(i)
 120      format(1x, a3, 1x,i2.2,':',i3.3,':',i5.5,1x, a29,1x,a32)


      end do

      write(unit,'(a)' ) '-INPUT/FILES'

****  Thats all
      return
      end
 
CTITLE CRE_INPUT_HIST

      subroutine cre_input_hist(unit, unitc )
 
      implicit none

*     Routine to read the input/history block
 
      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
 
* PASSED VARIABLES
 
*   unit        - Unit number of sinex file
*   unitc       - Commnts file unit number
 
      integer*4 unit, unitc
 
 
* LOCAL VARIABLES
 
*   ierr        - IOSTAT error
 
      integer*4 ierr
 

*   st_yr, st_doy, st_sec - Start of data yr, doy and seconds of day
*   en_yr, en_doy, en_sec - end of data yr, doy and seconds of day
*   cr_yr, cr_doy, cr_sec - Creation yr, doy and seconds of day
*   np          - Number of parameters in solution
*   ns          - Short version of solution number
*   cons_type   - Constraint type

      integer*4 st_yr, st_doy, st_sec, en_yr, en_doy, en_sec,  ns,
     .          cr_yr, cr_doy, cr_sec

*   type        - Type of input (+added; = combined)
*   prog        - File type for input
*   ver         - Version of file type.
*   owner       - Who owned the filed
*   sys_type    - System type P=GPS, R=VLBI, S=SLR
*   anal_type   - Types of parameters in solution

      character*4   sys_type


****  Read down each line saving the information
      write(unit,'(a)' ) '+INPUT/HISTORY'
      call cp_comments(unit, unitc, 'INPUT/HISTORY')
      do ns = 1, cnum_soln_recs
          call jd_to_yds(qsrun_time(ns),  cr_yr,cr_doy,cr_sec)
          call jd_to_yds(qsepoch_start(ns),st_yr,st_doy,st_sec )
          call jd_to_yds(qsepoch_end(ns), en_yr,en_doy,en_sec ) 
          call snx_sys_type( qssys_type(ns), sys_type)

          call sub_null( qsprog_gen(ns)  )
          call sub_null( qscreator(ns)   )
          call sub_null( qsowner(ns)     )
          call sub_null( qsanal_type(ns) )

          if (ns.eq.cnum_soln_recs ) then 
              qsprog_gen(ns) = '=SNX'
c             qsglb_ver(ns)  = 100
          end if

          write(unit,120, iostat=ierr) qsprog_gen(ns), 
     .           qsglb_ver(ns)/100.0, qscreator(ns),
     .           mod(cr_yr,100),cr_doy, cr_sec,
     .           qsowner(ns), mod(st_yr,100), st_doy, st_sec,
     .           mod(en_yr,100), en_doy, en_sec, sys_type,
     .           qsnum_parn(ns), qscons_type(ns), qsanal_type(ns) 

 120          format(1x,a4,f5.2, 1x,a3,1x,i2.2,':',i3.3,':',i5.5,1x,
     .               a3,1x,i2.2,':',i3.3,':',i5.5,1x,
     .               i2.2,':',i3.3,':',i5.5,1x,
     .               a1,1x,i5.5,1x,i1,1x,a5)

      end do
      write(unit,'(a)' ) '-INPUT/HISTORY'
****  Thats all
      return
      end
