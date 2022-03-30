CTITLE CRESNX_SAT

      subroutine cresnx_sat( unit, unitc)

      implicit none

*     Routine to create the satellite information blocks

      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
      include '../includes/glb_hdr_def.h'
      include '../includes/svinf_def.h'

* PASSED VARIABLES
 
*   unit        - Unit number
*   unitc       - Commnts file unit number
 
      integer*4 unit, unitc

****  Call the main routines.  Only write records if they exist
*     and prn number of first satellite is not zero.
      if( crec_svinf.gt.0 .and. qsvi_prn(1).ne.0 ) then 
         call cre_sat_id ( unit, unitc)
         call cre_sat_phs( unit, unitc)
      endif

****  Thats all
      return
      end

CTITLE CRE_SAT_ID

      subroutine cre_sat_id ( unit, unitc)

      implicit none

*     Routine to write the satellite ID block.

      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/svinf_def.h'

* PASSED VARIABLES
 
*   unit        - Unit number of sinex file
*   unitc       - Commnts file unit number
 
      integer*4 unit, unitc
 
 
* LOCAL VARIABLES
      real*8 dtl   ! Total Days since launch

      integer*4 dyrl  ! Years since launch
     .,         ddayl ! Days since launch   
     .,         dsl   ! Seconds part of day since launch
     .,         i     ! Loop counter
     .,         prn, svs   ! PRN and SVS numbers

      character*20 sv_ant_name  ! Full numbe of antenna
      character*1 sys    ! Type of system (GREC for GNSS)
      character*1 offcon    ! Function to return constellation letter 
                            ! based on offset PRN (original PRN for 
                            ! mod(prn,100))

  
      write(unit,'(a)') '+SATELLITE/ID'
      call cp_comments(unit, unitc, 'SATELLITE/ID')
      do i = 1, cnum_svs

* MOD TAH 1906127: Only output entry if site actually used.
         if( gtol_svs(i).ge.0 ) then 

* MOD TAH 051026: Change time to launch date
            dtl   = qsvi_launch(i)
            if( dtl.gt. 2401960.50d0 ) then 
               call jd_to_yds( dtl, dyrl, ddayl, dsl)
            else
               dyrl  = 0
	     ddayl = 0
	     dsl   = 0
            end if

* MOD TAH 180401: Output for GNSS
            sys =  offcon( qsvi_prn(i))
            prn =  mod(qsvi_prn(i),100)    

* MOD TAH 190621: Generate sv_ant_name based on qsvi_block
            call name_to_blk(-1, sv_ant_name, qsvi_block(i))
            
            write(unit, 120) sys, qsvi_svn(i), prn, 'P', dyrl,
     .            ddayl, dsl, sv_ant_name 
 120        format(1x,a1,I3.3,1x,i2.2,1x,9('-'),1x,a1,1x,
     .             I2.2,':',i3.3,':',I5.5,1x,'00:000:00000',1x,
     .             a20)
         endif

      end do
      write(unit,'(a)') '-SATELLITE/ID'

****  Thats all
      return
      end

CTITLE CRE_SAT_PHS

      subroutine cre_sat_phs( unit, unitc)

      implicit none

*     Routine to write the satellite Phase model block.

      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/svinf_def.h'

* PASSED VARIABLES
 
*   unit        - Unit number of sinex file
*   unitc       - Commnts file unit number
 
      integer*4 unit, unitc
 
 
* LOCAL VARIABLES
      integer*4 i
      character*1 offcon    ! Function to return constellation letter 
                            ! based on offset PRN (original PRN for 
                            ! mod(prn,100))
      character*1 sys       ! GREC system.
  
      write(unit,'(a)') '+SATELLITE/PHASE_CENTER'
      call cp_comments(unit, unitc, 'SATELLITE/PHASE_CENTER')
      do i = 1, cnum_svs

* MOD TAH 1906127: Only output entry if site actually used.
         if( gtol_svs(i).ge.0 ) then 
 
* MOD TAH 180401: Output for GNSS
            sys =  offcon( qsvi_prn(i))

            write(unit, 120) sys,qsvi_svn(i), qsvi_ocode(i)(1:1),
     .           qsvi_antpos(3,1,i), qsvi_antpos(1,1,i),
     .           qsvi_antpos(2,1,i), qsvi_ocode(i)(2:2),
     .           qsvi_antpos(3,2,i), qsvi_antpos(1,2,i),
     .           qsvi_antpos(2,2,i), qsvi_antmod(i)(1:10),
     .           qsvi_ocode(i)(3:3), qsvi_ocode(i)(4:4)

 120        format(1x,a1,I3.3,2(1x,a1,3(1x,F6.4)),1x,a10,1x,
     .             a1,1x,a1)
         endif

      end do
      write(unit,'(a)') '-SATELLITE/PHASE_CENTER'

****  Thats all
      return
      end


         

