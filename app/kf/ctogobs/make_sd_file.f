CTITLE WRITE_SD_FILE
 
      subroutine write_sd_file( ep, cf, prn, L1r_phs_cse, L2r_phs_cse,
     .    L1_cyc_cse, L2_cyc_cse,  L1r_rng_cse, L2r_rng_cse,
     .    ctol_cse, data_flag_cse, 
     .    params_cse, par_flag_cse  )

      implicit none

*     This routine will write out single difference (station-station)
*     files which can be read by mon_data for processing carrier
*     phase point position determinations.
 
* INCLUDES
 
      include '../includes/const_param.h'
      include '../includes/kalman_param.h'
      include '../includes/cfile_def.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES

*   ep  - Epoch number
*   cf  - Cfile numbers
*   prn - prn number 
*   data_flag_cse(num_chan, num_cfiles, num_ep) - Data flag for each
*                   - measurement (same definition as in Gobs File)
*   ctol_cse(num_chan, num_cfiles, num_ep)  - Conversion from
*                   - channel number to satellite number
*   par_flag_cse(num_param, num_ep)     - Parameter estimate quality
*                   - flags.
 
 
      integer*4 ep, cf, prn, 
     .    data_flag_cse(num_chan, num_cfiles, num_ep),
     .    ctol_cse(num_chan, num_cfiles, num_ep),
     .    par_flag_cse(num_param, num_ep)
 
*   L1_cyc_cse(num_chan, num_cfiles, num_ep)    - Number of cycles
*                   - needed for each L1 phase measurement .  May be
*                   - fracttional for half cycle units.
*   L2_cyc_cse(num_chan, num_cfiles, num_ep)    - number of cycles
*                   - needed for each L2 phase measurement
*   L1r_phs_cse(num_chan, num_cfiles, num_ep)  - L! phase residuals
*                   - cylces at L1
*   L2r_phs_cse(num_chan, num_cfiles, num_ep)  - L2 phase residuals
*                   - cycles at L2
*   params_cse(num_param, num_ep)       - Clock parameter estimates
*                   - by epoch.
 
 
      real*8 L1_cyc_cse(num_chan, num_cfiles, num_ep),
     .    L2_cyc_cse(num_chan, num_cfiles, num_ep),
     .    L1r_phs_cse(num_chan, num_cfiles, num_ep),
     .    L2r_phs_cse(num_chan, num_cfiles, num_ep),
     .    L1r_rng_cse(num_chan, num_cfiles, num_ep),
     .    L2r_rng_cse(num_chan, num_cfiles, num_ep),
     .    params_cse(num_param, num_ep)

* LOCAL VARIABLES
 
*   c1, c2   - Channel numbers
*   ltoc     - Function to convert list number to channel number
*   ls       - list number
*   curr_seq(max_gprn) - Current sequecnce number by PRN (incremented
*              and new file opened for each bias flag).
*   curr_sum(max_gprn) - Number of data written to each file.  Used
*              to see if we should just delete an empty file.
*   curr_cf  - Current cfile number 
*   ierr     - IOstat error
 
      integer*4 c1, c2, ltoc , ls, j, trimlen, curr_seq(max_gprn),
     .          curr_sum(max_gprn), curr_cf, ierr
 
*   data_ok  - Logoical function returns true is data OK
*   kbit     - Checks if bit is set.
 
      logical data_ok, kbit

*   dL1, dL2  - L1 and L2  residuals.
*   dR1,dR2   - Range residual
 
      real*8 dL1, dL2, dR1,dR2
 
*   outfile - Name of output file
 
      character*128 outfile

      data curr_cf /  0   /

****  If this is this first call then open all of the output
*     files.
      if( curr_cf.ne.cf ) then
 
          do j = 1, num_sat
              curr_seq(prn_list(j)) = 0
              curr_sum(prn_list(j)) = 0
          end do
          curr_cf = cf         

*         Open the unit for the position estimates currently
*         in the cfile.
          write(outfile,100) 
     .            sng_diff_root(1:trimlen(sng_diff_root)),
     .            cf_codes(1), cf_codes(cf)
 100      format(a,'POS_',a4,'_',a4)
          open(400,file=outfile,iostat=ierr, status='unknown')
          call report_error('IOSTAT',ierr,'open',outfile,0,
     .                      'WRITE_SD_FILE')
          if( ierr.ne.0 ) RETURN
          write(400,110) cf_codes(cf), cf_codes(1)
 110      format('* Apriori position estimates for ',a,' and ',a,/,
     .           '* Epoch       X Y and Z (meters)')
      end if

*     Now write out the entries for each satellite.
      ls = prntol(prn)
      c1 = ltoc(ctol_cse(1,1,ep),ls,actual_max_chan)
      c2 = ltoc(ctol_cse(1,cf,ep),ls,actual_max_chan)

*     see if we havedata
      if( c1.gt.0 .and. c2.gt.0 ) then

*         See if data OK

          if( data_ok(data_flag_cse(c1,1,ep),0,phs_mask) .and.
     .        data_ok(data_flag_cse(c2,cf,ep),0,phs_mask) ) then

*              See if bias flag set (close ould file and open a new
*              one.
               if( kbit(data_flag_cse(c1,1,ep),31) .or.
     .             kbit(data_flag_cse(c1,1,ep),32) .or.
     .             kbit(data_flag_cse(c2,cf,ep),31) .or.
     .             kbit(data_flag_cse(c2,cf,ep),32) ) then

                   if( curr_sum(prn).eq.0 ) then
C                      close(400+prn, iostat=ierr, status='delete')
                       close(400+prn)
                       write(*,9202) prn, curr_seq(prn)
9202                   format('Closing sd file for PRN ',i2,
     .                        ' current sequence is ',i4,' Delete')
                   else
                       close(400+prn)
                   end if
                   curr_sum(prn) = 0
                   curr_seq(prn) = curr_seq(prn) + 1
                   write(outfile,205) 
     .                    sng_diff_root(1:trimlen(sng_diff_root)),
     .                    cf_codes(1), cf_codes(cf), prn,curr_seq(prn) 
 205               format(a,a4,'_',a4,'_PRN',i2.2,'.',I3.3)
                   open(400+prn, file = outfile,  status='unknown', 
     .                  iostat=ierr)
                   write(*,210) outfile(1:trimlen(outfile))
 210               format(' Creating SD file ',a)
                   write(400+prn,220,iostat=ierr) cf_codes(1),  
     .                        cf_codes(cf),  prn, curr_seq(prn)
 220               format(' Single difference file for ',a4,'-',a4,
     .                    ' PRN ',i2.2,' Bias ',i3,/,
     .                    ' YMID=   0.e0')
               end if

*              compute residuals
               dL1 =(L1r_phs_cse(c2,cf,ep) + L1_cyc_cse(c2,cf,ep)
     .              -  params_cse(cf,ep) ) -
     .              (L1r_phs_cse(c1,1,ep) + L1_cyc_cse(c1,1,ep)
     .              -  params_cse(1,ep) )
               dL2 =(L2r_phs_cse(c2,cf,ep) + L2_cyc_cse(c2,cf,ep)
     .              -  params_cse(cf,ep)*fl2(ls)/fl1(ls) ) -
     .              (L2r_phs_cse(c1,1,ep) + L2_cyc_cse(c1,1,ep)
     .              -  params_cse(1,ep)*fl2(ls)/fl1(ls) )
               dR1 = (L1r_rng_cse(c2,cf,ep)-params_cse(cf,ep)) -
     .               (L1r_rng_cse(c1,1,ep)-params_cse(1,ep))
               dR2 = (L2r_rng_cse(c2,cf,ep)-params_cse(cf,ep)*
     .                  fL2(ls)/fL1(ls)) -
     .               (L2r_rng_cse(c1,1,ep)-params_cse(1,ep)*
     .                  fL2(ls)/fL1(ls))

               if( abs(dR2).gt.1.d6 ) dR2 = 0

*              Now write the file
               curr_sum(prn) = curr_sum(prn) + 1
               write(400+prn,300, iostat=ierr) cf_elev(1)*180/pi,
     .                     cf_azimuth(1)*180/pi,  ep, dL1, dL2, 
     .                     dR1, dR2, data_flag_cse(c2,cf,ep),
     .                     data_flag_cse(c1,1,ep)
 300           format(2f10.4,1x,I5,1x,2(F12.4,1x),2(F10.2,1x),2(1x,o12))
          end if
      else

*         We must be missing the data point at the other station.  If
*         we have a bias flag at this station then create the single 
*         difference file
          if( data_ok(data_flag_cse(c2,cf,ep),0,phs_mask) ) then

*             See if bias flag.  If so close old file and start a new
*             one.
              if( kbit(data_flag_cse(c2,cf,ep),31) .or.
     .            kbit(data_flag_cse(c2,cf,ep),32) ) then

                   if( curr_sum(prn).eq.0 ) then
C                      close(400+prn, iostat=ierr, status='delete')
                       close(400+prn, iostat=ierr)
                       write(*,*) 'Closing PRN',prn,' Status delete?'
                   else
                       close(400+prn)
                   end if

                   curr_seq(prn) = curr_seq(prn) + 1
                   curr_sum(prn) = 0
                   write(outfile,205) 
     .                    sng_diff_root(1:trimlen(sng_diff_root)),
     .                    cf_codes(1), cf_codes(cf), prn,curr_seq(prn) 
                   open(400+prn, file = outfile,  status='unknown')
                   write(*,210) outfile(1:trimlen(outfile))
                   write(400+prn,220,iostat=ierr) 
     .                     cf_codes(1), cf_codes(cf), 
     .                          prn, curr_seq(prn)
               end if
           end if
      end if
 
***** Thats all
      return
      end
 
