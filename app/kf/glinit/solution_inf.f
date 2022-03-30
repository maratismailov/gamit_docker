CTITLE SOLUTION_INF
 
      subroutine solution_inf

      implicit none 
 
 
*     Saves counts for number of observations, and deletes
 
      include '../includes/kalman_param.h'
      include '../includes/globk_cntl.h'
      include '../includes/globk_markov.h'
      include '../includes/glb_hdr_def.h'
 
*   i       - Loop counter
*   first_call -- Logical to indicate first call
*   save_ep    -- Saved epoch to get the first or last experiment
*                 information
 
      integer*4 i, cand
      logical first_call
      real*8 save_ep
      
      save save_ep  

c***  added by S. Shapiro (rwk 000814) for IRIX system (ok for others?)
      save first_call
      
      data  first_call / .true. /
      
****  Increment number of observations
 
      gnum_obs = gnum_obs + cnum_obs
      gnum_soln_recs = gnum_soln_recs + cnum_soln_recs
      gsys_type = cand(gsys_type, csys_type)
 
      do i = 1, max_edit_types
          gdelete_count(i) = gdelete_count(i) + cdelete_count(i)
      end do

***** Save the frame etc. information
      if( cgamit_mod.eq.0 ) then
          cgamit_mod = 3*65536 + 3
* MOD TAH 150802: Changed test for swap (not really needed anymore)
      else if ( cgamit_mod.gt.2**30 ) then
          call report_error('SWAP',1,'Swapping','cgamit_mod',0)
          ! call swap_bytes(4, cgamit_mod, 1)
      endif
      if( first_call ) then
          save_ep = cepoch_expt
          gsys_type = csys_type
          ggamit_mod = cgamit_mod
          gload_mod = cload_mod
          first_call = .false.
      else
          if( ggamit_mod.ne.cgamit_mod ) then
              call report_gamit_mod(ggamit_mod,cgamit_mod)
          endif 
      end if

      if( (cepoch_expt-save_ep)*sort_direction.ge.0 ) then
          ggtime   = cgtime
          ggframe  = cgframe
          ggprec   = cgprec
          ggsrpmod = cgsrpmod

* MOD TAH 050622: Save the model information.  In future should
*         check the status of the models and bits to see if
*         they match
          gload_mod = cload_mod
          gspeopmod = cspeopmod 
          getidemod = cetidemod  
          gotidemod = cotidemod  
          goatmlmod = coatmlmod  
          gatmtdmod = catmtdmod  
          ghydromod = chydromod 

* MOD TAH 080103: Save nutation and gravity field
* MOD TAH 140327: Added eradmod and antradmod
          if( cgnut(1:1).eq.'I') then
              ggnut  = cgnut
              gggrav = cggrav
              ggeradmod  = ceradmod
              ggantradmod = cantradmod
          else    ! Save default
              ggnut  = 'IAU80'
              gggrav = 'EGM96'
              ggeradmod  = 'NONE'
              ggantradmod = 'NONE'
          endif

* MOD TAH 140403: Added atm modeling options
          call sub_null(cdryzen)
          call sub_null(cwetzen)
          call sub_null(cdrymap)
          call sub_null(cwetmap)
          call sub_null(cionsrc)
          call sub_null(cmagfield)

          ggdryzen   = cdryzen 
          ggwetzen   = cwetzen 
          ggdrymap   = cdrymap 
          ggwetmap   = cwetmap 
          ggionsrc   = cionsrc 
          ggmagfield = cmagfield 
           
      end if
      
           
***** Thats all
      return
      end

CTITLE REPORT_GAMIT_MOD

      subroutine report_gamit_mod( ggamit_mod, cgamit_mod )

      implicit none 

*     Routine to report differences between gamit models from
*     first hfile and subsequent hfiles

* ggamit_mod  -- Global Gamit models set (bit mapped)
* cgamit_mod  -- Gamit models from hfile being read

      integer*4 ggamit_mod, cgamit_mod

* LOCAL Variables

* kbit  -- Checks bits
* trimlen -- Length of line
* lm      -- Current length of message
* gamit_mod_name -- Function to return name of gamit model name

      integer*4 trimlen, lm, i
      logical kbit
      character*8 gamit_mod_name

* mess  -- Message line written to user
      character*256 mess

****  See which bits differ
      write(mess,120) cgamit_mod, ggamit_mod
 120  format('GAMIT Model options differ: Codes ',2o10, ' Names :')
C     mess = 'GAMIT Model options differ:'
      do i = 1, 32
         if( kbit(cgamit_mod,i).neqv.kbit(ggamit_mod,i)) then
             lm = trimlen(mess) + 2
             mess(lm:) = gamit_mod_name(i)
         end if
      enddo

****  Report to user
      write(*,'(a)') mess(1:lm)
      
      call report_stat('WARNING','GLINIT','solution_inf',' ',
     .                 mess,0)
      
****  Thats all
      return
      end

 
