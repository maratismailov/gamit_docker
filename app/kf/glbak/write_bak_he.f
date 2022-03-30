CTITLE WRITE_BAK_HEADER
 
      subroutine write_bak_header( iout, options )
 
      implicit none 
 
*     Routine to write out the header for the Global back solution.
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
      include '../includes/glb_hdr_def.h'
 
*   expt_ref(5) - Reference epoch for the solution
*   i,j         - Loop counters
*   iout        - Output Lu
*   ierr        - IOSTAT error
*   options     - Options for output
*   TrimLen     - HP function for length of string
 
      integer*4 expt_ref(5), i, iout,
     .    ierr, options, TrimLen
 
*   deletes     - Number of deleted data
*   used        - Number of used data
 
      integer*4 deletes, used
 
*   postfit_chi  - Chisquared of the prefit parameter fit
 
      real*4 postfit_chi
 
*   dec_ref     - Decimal years for reference
*   jd          - Main memory JD
*   sec_tag     - Dummy second tag
 
      real*8 dec_ref, jd, sec_tag
      character*12 hsver, full_ver

*   mod_line   - Line with model names
*   mod_ent    - individual model name
*   gamit_mod_name -- Function to return model names
*   lm         - Length of mod_line
*   Load_use   - Set when load model is available but may have been removed
*                using the appload_mod command (USED or NOT USED when removed)

      integer*4 lm
      character*256 mod_line
      character*8 mod_ent, gamit_mod_name, load_use

*   kbit        - Bit checking function (starts at bit 1)
      logical kbit 
 
***** Write the header
      full_ver = hsver(glbak_version) 
      write(iout,100, iostat=ierr) full_ver(1:trimlen(full_ver)), 
     .           glb_inp_file(1:trimlen(glb_inp_file))
  100 format(//,57('-'),/,
     .          ' GLBAK Ver ',a,1x,a,/,
     .          57('-'),/)
 
*     Get start and end time
      jd = gepoch_expt
      call JD_to_YMDHMS( jd, expt_ref, sec_tag)
      call JD_to_Decyrs( jd, dec_ref )
 
      write(iout, 150) expt_ref, dec_ref, sec_tag,
     .                 expt_ref, dec_ref, sec_tag
 
  150 format(' EXPERIMENT date : ',i4,2('/',i2),1x,i2,':',i2,
     .       4x,'(',f9.4,') [Seconds tag ',F7.3,']',/,
     .       ' Solution refers to     : ',i4,2('/',i2),1x,i2,':',i2,
     .       4x,'(',f9.4,') [Seconds tag ',F7.3,']' )

* MOD TAH 930122: Output the IC epoch so that we can generate gfiles
*     directly from GLBAK output.
*     See if we should output the satellite epochs (Only output if
*     number of satellites in this session is greater than zero.
      if( cnum_svs.gt.0 ) then
          call jd_to_ymdhms( csvs_epoch, expt_ref, sec_tag)
          write(iout,160) expt_ref, sec_tag
  160     format(' Satellite IC epoch     : ',i4,2('/',i2),1x,i2,':',i2,
     .                                       1x,f5.2)
*         Now write out the corrdinate system parameters
          if( trimlen(cgtime).eq.0 .or. ichar(cgtime(1:1)).eq.0 ) then
              cgtime = 'UNKN'
          end if
          if( trimlen(cgframe).eq.0 .or. ichar(cgframe(1:1)).eq.0 ) then
              cgframe = 'UNKN '
          end if
          if( trimlen(cgprec).eq.0 .or. ichar(cgprec(1:1)).eq.0 ) then
              cgprec = 'UNKN '
          end if
          if(trimlen(cgsrpmod).eq.0 .or. ichar(cgsrpmod(1:1)).eq.0) then
              cgsrpmod = 'UNKN '
          end if
          cgnut  = ggnut
          cggrav = gggrav
          if(trimlen(ggnut).eq.0 .or. ichar(ggnut(1:1)).eq.0) then
              cgnut = 'UNKN '
          end if
          if(trimlen(gggrav).eq.0 .or. ichar(gggrav(1:1)).eq.0) then
              cggrav = 'UNKN '
          end if

* MOD TAH 140327: Added ggeradmod, ggantradmod
          if(trimlen(ceradmod).eq.0 .or. 
     .       ichar(ceradmod(1:1)).eq.0) then
              ceradmod = 'UNKN '
          end if
          if(trimlen(cantradmod).eq.0 .or. 
     .       ichar(cantradmod(1:1)).eq.0) then
              cantradmod = 'UNKN '
          end if

* MOD TAH 140327: Added ceradmod and cantradmod to output          
          write(iout,165) cgtime, cgframe, cgprec, cgsrpmod, 
     .                    cgnut, cggrav, ceradmod, cantradmod
  165     format(' GPS System Information : Time ',a4,' Frame ',a5,
     .           ' Precession ',a5,' Radiation model ',a5,
     .           ' Nutation ',a5,' Gravity ',a5,
     .           ' EarthRad ',a5,' AntThrust ',a5)

      end if

* MOD TAH 191222: Added more model output to be consistent with
*     GLOUT/GLORG

*     Write out the gamit model information
      mod_line = ' MODELS Used in Analysis: '
      lm = trimlen(mod_line) + 2
      do i = 1,32
         if( kbit(cgamit_mod,i) ) then
             mod_ent = gamit_mod_name(i)
             mod_line(lm:) = mod_ent // '|'
             lm = lm + 10
         end if
      end do

      write(iout, '(a)' ) mod_line(1:lm)

* MOD TAH 130419: Add reporting of the load model status.  See if 
*     atmload name set
      mod_line = ' LOAD Models Used       : '
      lm = trimlen(mod_line) + 2
      if( trimlen(catmtdmod).gt.0 ) then
          mod_line(lm:) = catmtdmod // '|'
          lm = lm + 10
      end if

      if( trimlen(coatmlmod).gt.0 ) then
*         There is an atm loading model 
          mod_line(lm:) = coatmlmod
*         check to see if used
          load_use = 'USED'
          if( kbit(appload_mod,1) .and. .not. kbit(appload_mod,2) )
     .                                        load_use = 'NOT USED'
          lm = lm + 10
          mod_line(lm:) = load_use // '|'
      end if
*     Now check Hydrolgical model
      if( trimlen(chydromod).gt.0 ) then
*         There is an atm loading model 
          lm = trimlen(mod_line) + 2
          mod_line(lm:) = chydromod
*         check to see if used
          load_use = 'USED'
          if( kbit(appload_mod,1) .and. .not. kbit(appload_mod,3) )
     .                                        load_use = 'NOT USED'
          lm = lm + 10
          mod_line(lm:) = load_use // '|'
      end if

      if ( lm.gt.30 ) write(iout, '(a)' ) trim(mod_line)

* MOD TAH 140403: Add line for Atm delay models used.
      write( iout,166) cdryzen, cwetzen, cdrymap,  cwetmap, 
     .                 cionsrc, cmagfield
 166  format(' ATM Delay Models Used  : ',4(a4,'    | '), 
     .       '2nd Order Ion     | ',2(a8,'| '))
! MODELS Used in Analysis: SD-WOB  | SD-UT1  | IERS10  | E-Tide  | K1-Tide | PoleTide| OC-Load | MPT2010 |  
! ATM Delay Models Used  : UFL     | GP25    | VMF1    | VMF1    | 2nd Order Ion     | GMAP    | IGRF11  | 
       

* MOD TAH 080724: Add output of reference frame (passed through
*     the +REFERENCE_FRAME line in the apriori files
      if( trimlen(reference_frame).eq.0 .or. 
     .    ichar(reference_frame(1:1)).eq.0 ) then
          reference_frame = 'Unknown'
      endif
      write(iout,167) reference_frame
 167  format(' Reference Frame        : ',a)

      write(iout, 190) grun_time
  190 format(' Run time        : ',i4,2('/',i2),1x,i2,':',i2,1x,
     .       i2,'.',i2)
 
*     Write number of data and solutions
*     Get used and deleted data
 
      deletes = 0
      do i = 1, max_edit_types
          deletes = deletes + cdelete_count(i)
      end do
 
      used = cnum_obs - deletes
 
      write(iout, 200) used, deletes, cnum_obs,
     .                 cnum_parn
  200 format(/,' There were ',i7,' data used, ',i7,' data not used and',
     .         i7,' data total',/,
     .         ' There were ',i7,' local parameters estimated')
      write(iout,210) cnum_sites, cnum_sources, cnum_svs
  210 format(' There were ',i4,' stations, ',i4,' radio sources, and ',
     .                      i4,' satellites' )

      
      
 
*     Get the prefit chisquared for the solution
      if ( compute_glb_res ) then
          if( cnum_used.gt.0 ) then
              postfit_chi = sum_loc_chi/cnum_used
          else
              postfit_chi = 0
          end if
 
          write(iout, 250) cnum_used, postfit_chi
  250     format(/' The postfit chisquared for ',i3,' local parameters',
     .            ' is ',f10.3)

      else

*         Output the prefit resisual for this day. 
          write(iout,260) cnum_used, dchi_save
  260     format(' The  prefit chi**2 for ',i7,' input parameters is',
     .           f10.3)

      end if
 
***** Thats all
      return
      end
 
