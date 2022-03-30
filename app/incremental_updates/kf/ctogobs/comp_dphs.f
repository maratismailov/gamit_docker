
CTITLE GEN_MFILE_NAME
 
      subroutine gen_mfile_name( line, indx, df_name, mf_name)

      implicit none
 
*     This routine will generate the mfile name from the df_name.
*     If no name is passed then form assumed is m<expt>a.<ddd>;
*     where <expt> amd <ddd> come from dfile name
*     if one character passed, then the character subs for a above;
*     if less than or equal seven characters are passed then the
*     characters passed substitute into the dfile name.
*     Otherwize the name as passed is used.
 
*
 
* PASSED VARIABLES
 
*         indx  - Current position in line
 
      integer*4 indx
 
*   line        - Line read from command file
*   df_name     - Name of dfile
*   mf_name     - Name of mfile
 
      character*(*) line, df_name, mf_name
 
* LOCAL VARIABLES
 
*   len_mfr     - Length of read portion of mfile root
*   trimlen     - Length of string.
 
      integer*4 len_mfr, trimlen
 
*   mf_root     - part of mf name read from command
 
      character*128 mf_root
 
****  Get the mf_root from the line
      call GetWord(line, mf_root, indx )
 
****  See what to do based on length of string
      len_mfr = trimlen(mf_root)
      if( len_mfr .le.1 ) then
          if( len_mfr.eq.0 ) mf_root = 'a'
          mf_name = 'm' // df_name(2:5) // mf_root(1:1) //
     .                    df_name(7:)
      else if( len_mfr.lt.7 ) then
          mf_name = mf_root(1:len_mfr) // df_name(len_mfr+1:)
      else
          mf_name = mf_root
      end if
 
***** Thats all
      return
      end
 
CTITLE COMP_DPHS
 
      subroutine comp_dphs( ep, cf, lcf, dphs_L1)

      implicit none
 
*     Routine to compute the phase correction at L1 based on the
*     adjustments to the parameters in the mfile.
 
      include '../includes/kalman_param.h'
      include '../includes/cfile_def.h'
      include '../includes/mfile_def.h'
 
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
 
*         ep    - Epoch number
*   cf      - Cfile number in mfile
*   lcf     - Local cfiles number
 
      integer*4 ep, cf, lcf
 
*   dphs_L1 - Correction to phase (cycles) at L1.
*   ff      - Fraction of time from atm boundary to current
*           - epoch
*   da      - Difference in atmospheric delay between knots
*           - on PWL atmospheric delay estimates.
 
      real*8 dphs_L1, ff, da 
 
* LOCAL VARIABLES
 
*         i,j       - Loop counters
*   svs_indx        - Parameter number for first orbital
*                   - element in current satellite
*   svs_cfix        - parameter of first orbital element by satellite
*                     in the cfile preval list minus 1
*   sof_indx        - Parameter number for satellite antenna offsets
*   sof_cfix        - parameter number in cfile for svs antenna offsets.
*   svs_clk_indx    - Parameter number for clock in current
*                   - current satellite.
*   eop_part        - Points to eop partials in cfile tm_part list
*                     (This numbers depends on if Berne model is used
*                      for orbits).
*   iel             - Element is zenith delay PWL or CON model 
*                     appropriate for current ep.
 
      integer*4 i,j, svs_indx, svs_clk_indx, eop_part, iel, svs_cfix,
     .          sof_indx, sof_cfix
     
*   linout          - Counter on number of lines written

      integer*4 linout
 
*      dt   - TIme difference in seconds from clock epoch time
*           - to current time
*   grad_part(2)    - Gradinent partials in NS and EW directions.
  
      real*8 dt, grad_part(2)

* MOD TAH 971204 Added change to apriori possibility      
*   dp_tot  - Total change in parameter value.  This is the 
*             mf_adj - (cf_apr_save-mf_aprval)

      real*8 dp_tot
      
      data linout / 0 /
 
****  Start, loop over all of the parameters
 
      dphs_L1 = 0.d0
   
*     Site positions
      if( site_np.gt.0 .and. cf.gt.0 ) then
          do i = 1, 3
              dp_tot  = mf_adjust(site_np+i-1) -
     .                  (cf_apr_save(i,lcf) - mf_aprval(site_np+i-1))         
*              dphs_L1 = dphs_L1 + cf_tmpart(i)*mf_adjust(site_np+i-1)
              if ( abs(dp_tot*cf_tmpart(i)).gt.100 .and.
     .             ep.lt.100 ) then
                  write(*,*) 'BAD SITE ADJ ',ep, cf, i, 
     .                       dp_tot, mf_adjust(site_np+i-1),
     .                       cf_apr_save(i,lcf), mf_aprval(site_np+i-1)
              end if 
              dphs_L1 = dphs_L1 + cf_tmpart(i)*dp_tot
          end do
      end if

* MOD TAH 040702: Check to see if average ATM delay needs to be added
      if( acon_np.gt.0 .and. cf.gt.0 ) then
          dphs_L1 = dphs_L1 + cf_tmpart(4)*mf_adjust(acon_np)
          if( ep.eq.-1 )
     .       print *,'acon_np ',cf, acon_np, mf_adjust(acon_np)
      endif
      
*     Atmospheric zenith delay parameters: PWL part
      if( atm_np.gt.0 .and. cf.gt.0 ) then
 
*         See if single zenith delay
* MDO TAH 040702: Test on single zenith delay not needed (and should not be
*         used) because of acon_np now accounts for signle zenith delay
C         if( mf_numzen.eq.1 ) then
C             dphs_L1 = dphs_L1 + cf_tmpart(4)*mf_adjust(atm_np)
C         else if( mf_zenmod.eq.'CON' ) then
          if( mf_zenmod.eq.'CON' ) then  
 
*             Get the correct index for the parameter
              call get_atm_iel(ep, mf_idtzen, mf_numzen, iel, 'CON')
 
              dphs_L1 = dphs_L1 + cf_tmpart(4)*
     .                            mf_adjust(atm_np+iel-1)
          else
              call get_atm_iel(ep, mf_idtzen, mf_numzen, iel, 'PWL')
              ff = float(ep-mf_idtzen(iel))/
     .            float(mf_idtzen(2)-mf_idtzen(1))
              da = mf_adjust(atm_np+iel) - mf_adjust(atm_np+iel-1)
              dphs_L1 = dphs_L1 + cf_tmpart(4)*
     .                        (mf_adjust(atm_np+iel-1)+da*ff)
          end if
     
      end if
 
*     See of gradients estimated
      do j = 1,2
          if( grad_np(j).gt.0 .and. cf.gt.0 ) then
              call comp_grad_part( grad_part, cf_elev(1), cf_azimuth(1))
              if( mf_numgrad.eq.1 ) then
                  dphs_L1 = dphs_L1 + grad_part(j)*
     .                      mf_adjust(grad_np(j))
              else if( mf_zenmod.eq.'CON' ) then
                 call get_atm_iel(ep, mf_idtgrad, mf_numgrad, 
     .                            iel, 'CON')
                     dphs_L1 = dphs_L1 + grad_part(j)*
     .                               mf_adjust(grad_np(j)+iel-1)
c                    write(*,*) 'GRA2', j, iel,grad_np(j),grad_part(j)*
c     .                               mf_adjust(grad_np(j)+iel-1)
              else
                 call get_atm_iel(ep, mf_idtgrad, mf_numgrad, 
     .                            iel, 'PWL')
                 ff = float(ep-mf_idtgrad(iel))/
     .                float(mf_idtgrad(2)-mf_idtgrad(1))
                 da = mf_adjust(grad_np(j)+iel) - 
     .                mf_adjust(grad_np(j)+iel-1)
                     dphs_L1 = dphs_L1 + grad_part(j)*
     .                         (mf_adjust(grad_np(j)+iel-1)+da*ff)
c                    write(*,*) 'GRA3', j, iel,grad_np(j),grad_part(j)*
c     .                               mf_adjust(grad_np(j)+iel-1)
              end if
          end if
      end do 
       

*     See if site clock estimated
      if( site_clk_np.gt.0 .and. cf.gt.0 ) then
          dt = (ep-1)*mf_inter
          if( cf_nversn.lt.980 ) then
             dphs_L1 = dphs_L1 + cf_tmpart(5)*(mf_adjust(site_clk_np) +
     .               dt*mf_adjust(site_clk_np+1) +
     .               dt**2/2*mf_adjust(site_clk_np+2)/86400.d0)
          else
             dphs_L1 = dphs_L1 + cf_tmpart(5)*mf_adjust(site_clk_np) 
          endif 
          if( ep.eq.-1 ) print *,'SITE_CLK_NP ',site_clk_np,
     .        mf_adjust(site_clk_np)               
      end if
 
*     Add orbital adjustments
      if( svs_np.gt.0 ) then
 
*         Do the IC's (always present)
          svs_indx = svs_np + mf_norb*(prntol(cf_iprn)-1)
          if( cf_nversn.lt.980 ) then
             svs_cfix = 8 + mf_norb*(prntol(cf_iprn)-1)
          else          
             svs_cfix = 6 + mf_norb*(prntol(cf_iprn)-1)
          endif
          if( ep.eq.-1 ) then
              print *,'SVS_CFIX ', cf, svs_cfix, svs_indx, mf_norb, 
     .                            cf_iprn, cf_nversn,
     .              (mf_adjust(svs_indx+i-1) -
     .             (cf_apr_save(svs_cfix+i-1,lcf) - 
     .              mf_aprval(svs_indx+i-1)),i=1,mf_norb) 
          end if

          do i = 1, mf_norb 
              dp_tot  = mf_adjust(svs_indx+i-1) -
     .                    (cf_apr_save(svs_cfix+i-1,lcf) - 
     .                     mf_aprval(svs_indx+i-1) )         
              dphs_L1 = dphs_L1 + cf_tmpart(5+i)*dp_tot
     

* MOD TAH 980220: Increased tolerance to 500 cycles (~100 meter
*             orbit adjusts).               
              if ( abs(dp_tot*cf_tmpart(5+i)).gt.500.0 .and.
     .             linout.lt.100  ) then
                  linout = linout + 1
                  write(*,*) '**WARNING** Large Orbit Adj ',
     .                       ep, cf, i, svs_cfix, svs_indx, dp_tot
              end if 
          end do
      end if
      
*     Do the Satellite antenna offset contributions
      if( sof_np.gt.0 ) then
 
*         Do the 3 values for the offsets 
          sof_indx = sof_np + 3*(prntol(cf_iprn)-1)
          sof_cfix = 6 + cf_norb*cf_nsat + 3*(prntol(cf_iprn)-1)
          if( ep.eq.-1 ) then
              print *,'SOF_INDX ',sof_indx, sof_cfix, sof_np, 
     .                 prntol(cf_iprn), cf_norb, cf_nsat,
     .                 (mf_adjust(sof_indx+i-1),
     .                 cf_apr_save(sof_cfix+i-1,lcf), 
     .                     mf_aprval(sof_indx+i-1),i=1,3)
          endif

          do i = 1,3 
              dp_tot  = mf_adjust(sof_indx+i-1) -
     .                    (cf_apr_save(sof_cfix+i-1,lcf) - 
     .                     mf_aprval(sof_indx+i-1) )         
              dphs_L1 = dphs_L1 + cf_tmpart(5+cf_norb+i)*dp_tot
          end do
      endif
      
*     Satellite clocks
      if( svs_clk_np.gt.0 ) then
          dt = (ep-1)*mf_inter
          svs_clk_indx = svs_clk_np + 3*(prntol(cf_iprn)-1)
          dphs_L1 = dphs_L1 + cf_tmpart(5)*(mf_adjust(svs_clk_indx) +
     .            dt*mf_adjust(svs_clk_indx+1) +
     .            dt**2/2*mf_adjust(svs_clk_indx+2)/86400.d0)
          if( ep.eq.-1 ) then
             print *,'svs_clk_indx ',svs_clk_indx, svs_clk_np,
     .            mf_adjust(svs_clk_indx)
          endif
      end if
 
*     EOP Parameters
      if( eop_np.gt.0 ) then
          eop_part = cf_npart - 6
          do i = 1,6
              dphs_L1 = dphs_L1 + cf_tmpart(eop_part+i)*
     .                    mf_adjust(eop_np+i-1)
          end do
          if( ep.eq.-1 ) then
              print *,'EOP_NP ',eop_np, (mf_adjust(eop_np+i-1),i=1,6)
          endif

      end if 
      if( ep.eq.-1 ) then
         print *,'DPHS 6 ',ep, dphs_L1, 6 + cf_norb*cf_nsat, cf_nparam
         do i = 6 + cf_norb*cf_nsat, cf_nparam
             print *,'cf_apr_save ',i,cf_apr_save(i,lcf)
         end do
       endif 

***** Thats all
      return
      end
      
CTITLE GET_MFS

      subroutine get_mfs( code, mfs )

      implicit none
      
*     Routine to get the site number for the cfile being processed
*     from the list of cfiles in the mfile
            
      include '../includes/kalman_param.h'
      include '../includes/mfile_def.h'
      
* PASSED VARIABLES

*  code  - 4-character code of cfile
*  mfs   - Site number

      integer*4 mfs
      character*4 code
      
* LOCAL VARIABLES

*  i  - loop counter

      integer*4 i

***** Scan the list of mfile cfiles to match code
      mfs = -1
      do i = 1, mf_nrfile
         if( code(1:4).eq.mf_rfname(i)(2:5) ) then
            mfs = i
         end if
      end do
      
****  Thats all
      return
      end      
 
CTITLE SET_NP
 
      subroutine set_np( mfs )

      implicit none
 
*     Routine to set up the pointers to the starts of the parameters
*     for the current site.
 
      include '../includes/kalman_param.h'
      include '../includes/cfile_def.h'
      include '../includes/mfile_def.h'
 
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
 
*         mfs   - Mfile site number (could be different if
*               - a list of cfiles are given) 

      integer*4 mfs
      
* Slot numbers for the different parameter types.
* MOD TAH 040702: Added acon_sl: slot number for average part of atm delay
* MOD TAH 210710: Added uupdates for 1072 version with new slot numbers to
*     accommodate upto 45 Beidou satellites.

      integer*4 site_sl, atm_sl, grad_sl(2), site_clk_sl, svs_sl,
     .          eop_sl, svs_clk_sl, i, sof_sl, acon_sl
 
****  Set up the slot numbers that we expect
 
      site_sl = mfs
* MOD TAH 040702: Removed the testing on numzen since new m-files
*     can have both average constant delay and PWL variation
      acon_sl = 300 + mfs
* MOD TAH 210701: New 1072 slot numbers
      if ( mf_nversn.lt.1072 ) then 
         atm_sl  = 11501 + mf_numzen*(mfs-1)
         grad_sl(1)  = 14001 + mf_numgrad*(mfs-1)
         grad_sl(2)  = 16501 + mf_numgrad*(mfs-1)
      else
         atm_sl  = 21501 + mf_numzen*(mfs-1)
         grad_sl(1)  = 24001 + mf_numgrad*(mfs-1)
         grad_sl(2)  = 26501 + mf_numgrad*(mfs-1)
      endif
 
      site_clk_sl = 400 + mfs
      svs_sl      = 501
      if( cf_nversn.lt.980 ) then
         sof_sl   = 0
         svs_clk_sl  = 2201
      else
* MOD TAH 190601: see if this earlier than mf_nversn 10.61
         if( mf_nversn.lt.1061 ) then  ! Pre-10.61 GAMIT version with ECOMC
            sof_sl      = 2001
            svs_clk_sl  = 0
         else
* MOD TAH 190601: Changed satellite PCO slot to 2400 from 2000 in eaarlier
*           gamit versions.
            sof_sl      = 2401
            svs_clk_sl  = 0
         endif
      endif
      eop_sl      = 80001
 
****  Now loop over all the parameters setting the parameter numbers
* MOD 040702: Added constant (average) atm delay
      site_np     = 0
      atm_np      = 0
      acon_np     = 0
      grad_np(1)  = 0
      grad_np(2)  = 0
      site_clk_np = 0
      svs_np      = 0
      sof_np      = 0
      eop_np      = 0
      svs_clk_np  = 0
      
      do i = 1, mf_mtpart
          if( mf_islot_all(i).eq. site_sl     ) site_np     = i
          if( mf_islot_all(i).eq. atm_sl      ) atm_np      = i
          if( mf_islot_all(i).eq. acon_sl     ) acon_np     = i
          if( mf_islot_all(i).eq. grad_sl(1)  ) grad_np(1)  = i
          if( mf_islot_all(i).eq. grad_sl(2)  ) grad_np(2)  = i
          if( mf_islot_all(i).eq. site_clk_sl ) site_clk_np = i
          if( mf_islot_all(i).eq. svs_sl      ) svs_np      = i
          if( mf_islot_all(i).eq. sof_sl      ) sof_np      = i
          if( mf_islot_all(i).eq. eop_sl      ) eop_np      = i
          if( mf_islot_all(i).eq. svs_clk_sl  ) svs_clk_np  = i
      end do

****  Thats all
      return
      end
 
CTITLE COMP_GRAD_PART
 
      subroutine comp_grad_part( grad_part, cf_elev, cf_azimuth)

      implicit none
 
*     Routine to compute the gradient partial derivatives
 
* PASSED VARIABLES
 
*      grad_part(2) - Gradient NS and EW partial.
      real*8 grad_part(2)
 
*   cf_elev         - Elevation (rads)
*   cf_azimuth      - Azimiuth (rads)
 
      real*4 cf_elev, cf_azimuth
 
* LOCAL VARIABLES
 
 
*     Compute the partials
* MOD TAH 191007: Fixed gradient partial (tan(el) missing)
 
      grad_part(1) = cos(cf_azimuth)/(sin(cf_elev)*tan(cf_elev)+0.003d0)
      grad_part(2) = sin(cf_azimuth)/(sin(cf_elev)*tan(cf_elev)+0.003d0)
 
***** Thats all
      return
      end
 
 
 
 
CTITLE GET_ATM_IEL
 
      subroutine get_atm_iel( ep, mf_idtzen,mf_numzen, iel, type)

      implicit none
 
*     Routine to get correct index for atmosphere parameters
 
* PASSED VARIABLES
 
*   ep      - Current epoch
*   mf_numzen       - Number of zenith delays
*   mf_idtzen(mf_numzen)    - Epoch boundaries for atm estimates.
*           - (Taken as ep .ge. correct and .lt. next entry)
*   iel     - Index in mf_idtzen that should be taken as left
*           - boundary
 
      integer*4 ep, mf_numzen, mf_idtzen(mf_numzen), iel
 
*             type  - Type of estimation CON or PWL.
 
      character*(*) type
 
* LOCAL VARIABLES
 
*       found   - Indicates that we have found the correct
*               - index.
 
      logical found
 
***** Loop up the index list until we find the one immediately
*     before our epoch
 
      iel = 0
      found = .false.
      do while ( .not.found )
          iel = iel + 1
          if( ep.lt.mf_idtzen(iel+1).and.
     .        ep.ge. mf_idtzen(iel)  ) then
              found = .true.
          else if( type(1:3).eq.'CON' .and.
     .            iel.eq.mf_numzen-1     ) then
*             Since with CON zenith delays, there is no
*             last entry in the idtzen list, trap the last
*             entry case specially.
              found = .true.
              iel = mf_numzen          
          else if ( iel.eq.mf_numzen ) then
 
*             We have a problem.  Should have found the
*             entry by now.
              call report_stat('warning','autcln','get_atm_iel',
     .            'Unexpected end of zenith delay epoch list',
     .            ' ',ep)
              iel = mf_numzen -1
              found = .true.
 
          end if
      end do
 
***** Thats all
      return
      end
 
 
