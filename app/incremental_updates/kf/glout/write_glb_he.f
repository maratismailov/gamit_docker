CTITLE WRITE_GLB_HEADER
 
      subroutine write_glb_header( iout, options, cov_parm, sol_parm )

      implicit none 
 
*     Routine to write out the header for the Global solution.
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'

* PASSED VARIABLES
      real*8   cov_parm(num_glb_parn,num_glb_parn), 
     .         sol_parm(num_glb_parn)     ! Solution and covariance function

 
*   expt_end(5) - Ending epoch of solution
*   expt_start(5)   - Staring epoch of solution
*   expt_ref(5) - Reference epoch for the solution
*   i,j         - Loop counters
*   iout        - Output Lu
*   ierr        - IOSTAT error
*   options     - Options for output
*   TrimLen     - HP function for length of string
 
      integer*4 expt_end(5), expt_start(5), expt_ref(5), i, iout,
     .    ierr, options, TrimLen

*   kbit        - Bit checking function (starts at bit 1)
      logical kbit 
 
*   deletes     - Number of deleted data
*   used_data   - Number of used data
*   used_sites  - Number of sites used in solution (guse_site set).
 
      integer*4 deletes, used_data,  used_sites 
 
*   prefit_chi  - Chisquared of the prefit parameter fit
*   postfit_chi - Postfit chisquared if residuals calculated
 
 
      real*4 prefit_chi, postfit_chi
 
*   num_pchi    - Number of values in prefit/postfit chi**2 (used to
*               - get the value as integer
 
      integer*4 num_pchi
 
*   dec_end     - Decimal years for end
*   dec_start   - Decimal years for start
*   dec_ref     - Decimal years for reference
*   jd          - Main memory JD
*   sec_tag     - Dummy second tag
 
      real*8 dec_end, dec_start, dec_ref, jd, sec_tag

      character*12 full_ver, hsver

*   mod_line   - Line with model names
*   mod_ent    - individual model name
*   gamit_mod_name -- Function to return model names
*   lm         - Length of mod_line
*   Load_use   - Set when load model is available but may have been removed
*                using the appload_mod command (USED or NOT USED when removed)

      integer*4 lm
      character*256 mod_line
      character*8 mod_ent, gamit_mod_name, load_use
 
***** Write the header
      write(iout,'(a)') gdescription(1:trimlen(gdescription))
      full_ver = hsver(globk_version)
      write(iout,100, iostat=ierr) full_ver(1:trimlen(full_ver))
  100 format(//,57('-'),/,
     .          ' GLOBK Ver ',a,', Global solution',/,
     .          57('-'),/)
 
*     Get start and end time
      jd = gepoch_start
      call JD_to_YMDHMS( jd, expt_start, sec_tag)
      call JD_to_Decyrs( jd, Dec_start )
      jd = gepoch_end
      call JD_to_YMDHMS( jd, expt_end, sec_tag)
      call JD_to_Decyrs( jd, Dec_end )
 
*     Reference epoch . Set when globk common is read to be the
*     end or beginning depending on the sort_direction.
      jd = gepoch_out
      call JD_to_YMDHMS( jd, expt_ref, sec_tag)
      call JD_to_Decyrs( jd, Dec_ref )
      write(iout, 150) expt_start, dec_start, expt_end, dec_end,
     .                 expt_ref,   dec_ref, sec_tag
  150 format(' Solution commenced with: ',i4,2('/',i2),1x,i2,':',i2,
     .       4x,'(',f9.4,')',/,
     .       ' Solution ended with    : ',i4,2('/',i2),1x,i2,':',i2,
     .       4x,'(',f9.4,')',/,
     .       ' Solution refers to     : ',i4,2('/',i2),1x,i2,':',i2,
     .       4x,'(',f9.4,') [Seconds tag ',F7.3,']' )

*     See if we should output the satellite epochs
      if( gnum_svs.gt.0 ) then
          call jd_to_ymdhms( svs_epoch(1), expt_ref, sec_tag)
          write(iout,160) expt_ref, sec_tag
  160     format(' Satellite IC epoch     : ',i4,2('/',i2),1x,i2,':',i2,
     .                                       1x,f5.2)
     
*         Now write out the corrdinate system parameters
          if( trimlen(ggtime).eq.0 .or. ichar(ggtime(1:1)).eq.0 ) then
              ggtime = 'UNKN'
          end if
          if( trimlen(ggframe).eq.0 .or. ichar(ggframe(1:1)).eq.0 ) then
              ggframe = 'UNKN '
          end if
          if( trimlen(ggprec).eq.0 .or. ichar(ggprec(1:1)).eq.0 ) then
              ggprec = 'UNKN '
          end if
          if(trimlen(ggsrpmod).eq.0 .or. ichar(ggsrpmod(1:1)).eq.0) then
              ggsrpmod = 'UNKN '
          end if
          if(trimlen(ggnut).eq.0 .or. ichar(ggnut(1:1)).eq.0) then
              ggnut = 'UNKN '
          end if
          if(trimlen(gggrav).eq.0 .or. ichar(gggrav(1:1)).eq.0) then
              gggrav = 'UNKN '
          end if

* MOD TAH 140327: Added ggeradmod, ggantradmod
          if(trimlen(ggeradmod).eq.0 .or. 
     .       ichar(ggeradmod(1:1)).eq.0) then
              ggeradmod = 'UNKN '
          end if
          if(trimlen(ggantradmod).eq.0 .or. 
     .       ichar(ggantradmod(1:1)).eq.0) then
              ggantradmod = 'UNKN '
          end if

* MOD TAH 140327: Added ggeradmod and ggantradmod to output          
          write(iout,165) ggtime, ggframe, ggprec, ggsrpmod, 
     .                    ggnut, gggrav, ggeradmod, ggantradmod
  165     format(' GPS System Information : Time ',a4,' Frame ',a5,
     .           ' Precession ',a5,' Radiation model ',a5,
     .           ' Nutation ',a5,' Gravity ',a5,
     .           ' EarthRad ',a5,' AntThrust ',a5)
      
      end if

*     Write out the gamit model information
      mod_line = ' MODELS Used in Analysis: '
      lm = trimlen(mod_line) + 2
      do i = 1,32
         if( kbit(ggamit_mod,i) ) then
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
      if( trimlen(gatmtdmod).gt.0 ) then
          mod_line(lm:) = gatmtdmod // '|'
          lm = lm + 10
      end if

      if( trimlen(goatmlmod).gt.0 ) then
*         There is an atm loading model 
          mod_line(lm:) = goatmlmod
*         check to see if used
          load_use = 'USED'
          if( kbit(appload_mod,1) .and. .not. kbit(appload_mod,2) )
     .                                        load_use = 'NOT USED'
          lm = lm + 10
          mod_line(lm:) = load_use // '|'
      end if
*     Now check Hydrolgical model
      if( trimlen(ghydromod).gt.0 ) then
*         There is an atm loading model 
          lm = trimlen(mod_line) + 2
          mod_line(lm:) = ghydromod
*         check to see if used
          load_use = 'USED'
          if( kbit(appload_mod,1) .and. .not. kbit(appload_mod,3) )
     .                                        load_use = 'NOT USED'
          lm = lm + 10
          mod_line(lm:) = load_use // '|'
      end if

      if ( lm.gt.30 ) write(iout, '(a)' ) trim(mod_line)

* MOD TAH 140403: Add line for Atm delay models used.
      write( iout,166) ggdryzen, ggwetzen, ggdrymap,  ggwetmap, 
     .                 ggionsrc, ggmagfield
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

* MOD TAH 210509: Output Balance status
      if( kbit(org_opts,32) .and. gepoch_end-gepoch_start.lt.2.0 ) then
         write(iout,'(a)') 
     .      ' Site and PCO covariance balanced by number of sessions'
      endif

****  Write runtime
      write(iout, 170) grun_time
  170 format(' Run time               : ',i4,2('/',i2),1x,i2,':',i2,1x,
     .       i2,'.',i2.2)

 
*     Write number of data and solutions
*     Get used and deleted data
 
      deletes = 0
      do i = 1, max_edit_types
          deletes = deletes + gdelete_count(i)
      end do
 
      used_data = gnum_obs - deletes

* MOD TAH 140109: Compute number of used sites
      used_sites = 0
      do i = 1, gnum_sites
          if( kbit(guse_site,i) )  used_sites = used_sites + 1
      end do
 
      write(iout, 200) gnum_soln_recs-gnum_comb, num_glb_sol, 
     .                 used_data, deletes, 
     .                 gnum_obs,  num_glb_parn
  200 format(/,' There were ',i9,' exps from  ',
     .      i9,' global files in the solution',/,
     .    ' There were ',i12,' data used, ',i7,' data not used and',
     .      1x,i12,' data total',/,
     .    ' There were ',i12,' global parameters estimated')
      write(iout,210) gnum_sites, used_sites, gnum_sources, gnum_svs
  210 format(' There were ',i5,' Sites, ',i4,' Used Sites, ',I4,
     .        ' radio sources, and ', i4,' Satellites' )
 
*     Get the prefit chisquared for the solution
      if( sum_chi_num.gt.0 ) then
          prefit_chi = sum_chi/sum_chi_num
      else
          prefit_chi = 0
      end if
 
      num_pchi = sum_chi_num
      write(iout, 250) num_pchi, prefit_chi
  250 format(/' The  prefit chi**2 for ',i7,' input parameters is',
     .        f10.3)
 
*     Now do some more statistics if we ran back solution with residuals
*     calculated
      if( glb_bak_soln .and. compute_glb_res ) then
 
          if( sum_post_chi_num.gt.0 ) then
              postfit_chi = sum_post_chi/sum_post_chi_num
          else
              postfit_chi = 0
          end if
          num_pchi = sum_post_chi_num
          write(iout, 275) num_pchi, postfit_chi
  275     format(' The postfit chi**2 for ',i7,' input parameters is',
     .            f10.3)
 
****      Now do the type statistics
 
          call out_type_stats( iout )
 
      end if
      if( uni_wght ) write(iout,280)
 280  format(' Weights of sites unified by number of times used')
 
 
*     Output file names
      write(iout, 300) list_file(1:max(1,trimlen(list_file))),
     .                 glb_com_file(1:max(1,trimlen(glb_com_file))),
     .                 glb_mar_file(1:max(1,trimlen(glb_mar_file))),
     .                 glr_cmd_file(1:max(1,trimlen(glr_cmd_file)))
      write(iout, 310) ( glb_apr_file(i)(1:max(1,
     .                   trimlen(glb_apr_file(i)))),i=1,num_apr_files)
   
      write(iout, 320) nut_inp_file(1:max(1,trimlen(nut_inp_file))),
     .                 plan_inp_file(1:max(1,trimlen(plan_inp_file))),
     .                 sd_inp_file(1:max(1,trimlen(sd_inp_file))),
     .                 pmu_inp_file(1:max(1,trimlen(pmu_inp_file))),
     .                 glb_bak_file(1:max(1,trimlen(glb_bak_file))),
     .                 glb_out_file(1:max(1,trimlen(glb_out_file))),
     .                 glb_svs_file(1:max(1,trimlen(glb_svs_file))),
     .                 svs_mar_file(1:max(1,trimlen(svs_mar_file)))
      do i = 1, num_eqfiles
         write(iout,410) eq_inp_file(i)(1:trimlen(eq_inp_file(i)))
      end do
 
  300 format(/' LIST file      : ',a,/,
     .        ' COMMON file    : ',a,/,
     .        ' GLOBK CMD file : ',a,/,
     .        ' GLORG CMD file : ',a)
  310 format( ' APRIORI file   : ',a)  
  320 format( ' NUTATION file  : ',a,/,
     .        ' PLANETARY file : ',a,/,
     .        ' SD ORIENT file : ',a,/,
     .        ' PMU file       : ',a,/,
     .        ' BACK SOLN file : ',a,/,
     .        ' OUTGLOBAL file : ',a,/,
     .        ' SVS EPHEM file : ',a,/,
     .        ' SVS MARKOV file: ',a)
  410 format( ' EARTHQUAKE file: ',a)

****  Report earthquakes
      call report_eq( iout, 'UNIQ' )

****  Write out the comparison of eq estmates
      if( kbit(options,20) ) then
         call summ_eq(iout, 0, cov_parm, sol_parm ) 
      endif

****  Report non-secular terms used in the analysis
      call report_nonsec( iout )

****  Report on stations with no-site corrdinate updates

      call report_noapr( iout )

****  Rreport the Markov process noise for sites if requested
      if( kbit(options,27) ) then
         call report_smar( iout)
      endif 

****  Output the list of the markov file if bit 4 is set
      if( kbit(options,4) ) then
          call list_markov( iout, glb_mar_file, comopt, 
     .                      list_file, 'GLOBK')
          if( trimlen(glr_cmd_file).gt.0 ) then
             call list_markov( iout, glr_cmd_file, comopt, 
     .                         list_file, 'GLORG')
          endif
      end if

* MOD TAH 970706: See if the srt_file is to be listed
      if( kbit(options,15) ) then
          call list_sort( iout, sort_file, num_glb_sol ) 
      end if
 
***** Thats all
      return
      end

CTITLE LIST_SORT

      subroutine list_sort( iout, sort_file, num_glb_sol ) 

*     Routine to list the GDL file as actually used in the solution.

      include '../includes/kalman_param.h' 

* PASSED VARIBALES

* iout  - Output unit number
* num_glb_sol - NUmber of global solutions 
* sort_file -- Name of the sort file

      integer*4 iout, num_glb_sol 
      character*(*) sort_file

* LOVAL VARIABLES

* ierr - IOSTAT error
* i    - Record number  
* max_name_len -- Maximum name length 
* trimlen -- Length of string 

      integer*4 ierr, i, max_name_len, trimlen 

* glb_var   -- Over all scaling factor
* glb_diag  -- Diagonal scaling factor

      real*8 glb_var, glb_diag
      
* MOD TAH 980519: Added reading and writing of forward and 
*   back chi**2/f to srt_file.

*   for_chi, bak_chi -- Forward and backwarsd chi**2/f

      real*4 for_chi, bak_chi

* cname     - Name of global file
* used_str  - String to denote if used.
* header    - Headings line for list 

      character*(sort_recl) cname 
      character*12 used_str
      character*256 header 

*     Open the sort file and read the whole thing to get
*     maximum file name length

      open(300, iostat=ierr, file=sort_file,  status='old',
     .         access='direct', recl=sort_recl+8+8)

      if( ierr.ne.0 ) then
         write(iout, 100) sort_file(1:trimlen(sort_file))
 100     format('Cannot list SRT_FILE ',a,' Open error ',i5)
         RETURN
      end if

****  Now scan the names getting the maxiumim length
      max_name_len = 0
      do i = 1, num_glb_sol
         read(300,iostat=ierr,rec=i) cname, glb_var,
     .                               for_chi, bak_chi
         max_name_len = max( max_name_len, trimlen(cname))
      end do
 
****  Write out the header name
      write(iout,200) sort_file(1:trimlen(sort_file))
 200  format(/,' EXPERIMENT LIST from ',a)
      header = '     #  Name '
      header(max_name_len+12:) = 
     .       '    SCALE Diag PPM  Forw Chi2 Back Chi2 Status '
      write(iout,'(a)') header(1:trimlen(header))
 
*     Now write the contents of the file
      do i = 1, num_glb_sol
         read(300,iostat=ierr,rec=i) cname, glb_var,
     .                                for_chi, bak_chi
 
*        See if actually used
         if( ichar(cname(1:1)).gt.128 ) then
             used_str = 'NOT USED'
C -- wrong order of brackets? (simon) cname(1:1) = char(ichar(cname(1:1)-128))
             cname(1:1) = char(ichar(cname(1:1))-128)
             
* MOD TAH 150514: Added for_chi == -1 to denote that no parameters were found
*            for this h-files (usually when use_site is removing many sites so
*            that some sites have no sites at all). (Set glfor when cnum_used = 0 )
             if( for_chi .eq. -1 ) used_str = 'NO DATA'
         else
             used_str = 'USED'
         end if
 
*        Now break up the scaling factors
         if( glb_var.lt.0 ) then
              glb_diag = (abs(glb_var*1.d3) - int(abs(glb_var*1.d3)))*
     .                   1000.d0
              glb_var = abs(glb_var) - glb_diag/1000.d3
              glb_diag = (glb_diag-1.d0)*1.d6 
              if( glb_diag.lt.0.d0 ) glb_diag = 0.d0 
          else  if ( glb_var.eq.0 ) then
              glb_var = 1.d00
              glb_diag = 0.d0
          else
              glb_diag = 0.d0
          end if 

          write(iout,250) i, cname(1:max_name_len), glb_var, glb_diag,
     .                  for_chi, bak_chi, used_str 
 250      format(i6,1x,a,F13.6,F8.3,2(1x,F9.3),3x,a)
      end do
 
***** Thats all
      close(300)
      return
      end
 
