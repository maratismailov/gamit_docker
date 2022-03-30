CTITLE READ_CF1
 
      subroutine read_cf1( unit, option, ierr )

      implicit none 
 
*     This routine will read the type 1 record of a C-file.  It is
*     the user responsibility to act on any error which is returned.
*     Here the error will be reported but the program will continue
*     executing.  The option passed indicates how much of the record
*     will be read.

* INCLUDE FILES
 
      include '../includes/cfile_def.h'
 
* PASSED VARIABLES
 
*   unit        - Unit number to which the file is attached.
*   ierr        - IOSTAT error in general but if the record type
*           - is wrong, then -1000 is returned as the error.
 
      integer*4 unit, ierr
 
*   option  - Option for read:  There are
*           - SHORT - Just read read the flag
*           - ALL   - Read complete record (default)
 
      character*(*) option
 
* LOCAL VARIABLES
 
*   flag        - First I*4 value in record and indicates the
*           - type of record.
*   jerr        - Error generated in check c-file
*           - -1000 - Type not correct
*           - -1001 - Dimension too small
*   i       - Loop counter
 
      integer*4 flag, jerr, i
 
*   loc_opt - Local version of option so that we can
*           - casefold it.
 
      character*4 loc_opt
 
***** Check the option and see what we should do:
      loc_opt = option
      call casefold(loc_opt)
 
      if( loc_opt(1:2).eq.'SH' ) then
* MOD TAH 950629: Read the version number.
          read(unit,iostat=ierr) flag, cf_nversn
          cf_ntext = 0
      else
          read(unit,iostat=ierr) flag, cf_nversn, cf_ntext,
     .        (cf_text(i), i=1, cf_ntext )
      end if                           
      
* MOD TAH 950629: Check the version number to see if OK
      if( cf_nversn.lt.900 ) then
          write(*,100) cf_nversn/100.d0
 100      format('**DISASTER** Old version of cfile.  Version',
     .           'number read ',F6.2)
          stop 'READ_CF1: Obsolete cfile version'
      end if
 
****  See what happened and check flags and sizes.
      if( ierr.eq.0 ) then
          jerr = 0
          call check_cflag( flag, 1, jerr )
          call check_csize( cf_ntext, cf_maxtxt, 'TEXT', jerr )
 
*         Save the jerr error in the return error flag
          ierr = jerr
      else
          call report_error('IOSTAT',ierr,'read',
     .        'Cfile Type 1 record', 0, 'READ_CF1')
      end if
 
****  Thats all
      return
      end
 
CTITLE READ_CF2
 
      subroutine read_cf2( unit, option, ierr )
                        
      implicit none 
 
*     This routine will read the type 2 record of a C-file.  It is
*     the user responsibility to act on any error which is returned.
*     Here the error will be reported but the program will continue
*     executing.  The option passed indicates how much of the record
*     will be read.
 
* INCLUDE FILES
 
      include '../includes/cfile_def.h'
 
* PASSED VARIABLES
 
*   unit        - Unit number to which the file is attached.
*   ierr        - IOSTAT error in general but if the record type
*           - is wrong, then -1000 is returned as the error.
 
      integer*4 unit, ierr
 
*   option  - Option for read:  There are
*           - DATA - Just read enough of this record so
*           -        the data can be read
*           - ALL   - Read complete record (default)
 
      character*(*) option
 
* LOCAL VARIABLES
 
*   flag        - First I*4 value in record and indicates the
*           - type of record.
*   jerr        - Error generated in check c-file
*           - -1000 - Type not correct
*           - -1001 - Dimension too small
*   i,j     - Loop counters
 
      integer*4 flag, jerr, i,j
 
*   loc_opt - Local version of option so that we can
*           - casefold it.
 
      character*4 loc_opt
 
***** Check the option and see what we should do:
      loc_opt = option
      call casefold(loc_opt)
            
cd      print *,'cf_nversn loc_opt ',cf_nversn,loc_opt
* MOD TAH 980915: Check for older versions of cfiles
      if( cf_nversn.lt.980 ) then 
         if( loc_opt(1:2).eq.'DA' ) then
             read(unit,iostat=ierr) flag, cf_sitnam, cf_rcvrsw, 
     .           cf_swver, cf_npart, cf_nsat,
     .           (cf_ischan(i),i=1,cf_nsat),
     .           cf_ndat, (cf_dattyp(i),i=1,cf_ndat),
     .           ((cf_lambda(i,j),i=1,cf_nsat),j=1,cf_ndat),
     .           cf_skd, cf_nepoch, cf_inter, cf_ircint, cf_mtime,
     .           cf_isessn, cf_iy, cf_im, cf_id, cf_ihr, cf_min, cf_sec
 
             cf_nslip = 0
             cf_nextra = 0
         else
             read(unit,iostat=ierr) flag, cf_sitnam, cf_rctype,
     .           cf_rcvnum, cf_rcvrsw, cf_swver,
     .           cf_anttyp, cf_antnum, 
     .           cf_npart, cf_nsat, (cf_ischan(i),i=1,cf_nsat),
     .           cf_ndat, (cf_dattyp(i),i=1,cf_ndat),
     .           ((cf_lambda(i,j),i=1,cf_nsat),j=1,cf_ndat),
     .           cf_skd, cf_nepoch, cf_inter, cf_ircint, cf_mtime,
     .           cf_isessn, cf_iy, cf_im, cf_id, cf_ihr, cf_min, cf_sec,
     .           (cf_offarp(i),i=1,3),
     .           (cf_offsL1(i),i=1,3), (cf_offsL2(i),i=1,3),
     .           cf_obfiln, cf_tfiln, cf_jfiln,
     .           cf_frame, cf_precmod, cf_nutmod, cf_gravmod,
     .           cf_srpmod, cf_ietide, cf_isptide, cf_antmod,
     .           cf_elvcut, cf_nclock, (cf_clock(i),i=1,cf_nclock),
     .           cf_jde, cf_te, cf_jdr, cf_tr,
     .           cf_ut1, cf_xp, cf_yp,  cf_psi,cf_eps,
     .           cf_avlmet, cf_serwgt,
     .           cf_nslip, (cf_islip(i),i=1,cf_nslip),
     .           (cf_slpst(i),i=1,cf_nslip),
     .           cf_nextra, (cf_extra(i),i=1,cf_nextra)
          end if
*         Set the value now read from cfile.
          cf_norb = 15
      elseif( cf_nversn.lt. 1020 ) then
*         Read the 980 version of the cfile record.      
          if( loc_opt(1:2).eq.'DA' ) then
             read(unit,iostat=ierr) flag, cf_sitnam, cf_rcvrsw, 
     .           cf_swver, cf_npart, cf_norb, cf_nsat,
     .           (cf_ischan(i),i=1,cf_nsat),
     .           cf_ndat, (cf_dattyp(i),i=1,cf_ndat),
     .           ((cf_lambda(i,j),i=1,cf_nsat),j=1,cf_ndat),
     .           cf_skd, cf_nepoch, cf_inter, cf_ircint, cf_mtime,
     .           cf_isessn, cf_iy, cf_im, cf_id, cf_ihr, cf_min, cf_sec
 
             cf_nslip = 0
             cf_nextra = 0
          else
             read(unit,iostat=ierr) flag, cf_sitnam, cf_rctype,
     .           cf_rcvnum, cf_rcvrsw, cf_swver,
     .           cf_anttyp, cf_antnum, cf_npart, cf_norb,
     .           cf_nsat, (cf_ischan(i),i=1,cf_nsat),
     .           cf_ndat, (cf_dattyp(i),i=1,cf_ndat),
     .           ((cf_lambda(i,j),i=1,cf_nsat),j=1,cf_ndat),
     .           cf_skd, cf_nepoch, cf_inter, cf_ircint, cf_mtime,
     .           cf_isessn, cf_iy, cf_im, cf_id, cf_ihr, cf_min, cf_sec,
     .           (cf_offarp(i),i=1,3),
     .           (cf_offsL1(i),i=1,3), (cf_offsL2(i),i=1,3),
     .           cf_obfiln, cf_tfiln, cf_jfiln,
     .           cf_frame, cf_precmod, cf_nutmod, cf_gravmod,
     .           cf_srpmod, cf_ietide, cf_isptide, cf_antmod,
     .           cf_elvcut, cf_nclock, (cf_clock(i),i=1,cf_nclock),
     .           cf_jde, cf_te, cf_jdr, cf_tr,
     .           cf_ut1, cf_xp, cf_yp,  cf_psi,cf_eps,
     .           cf_avlmet, cf_serwgt,
     .           cf_nslip, (cf_islip(i),i=1,cf_nslip),
     .           (cf_slpst(i),i=1,cf_nslip),
     .           cf_nextra, (cf_extra(i),i=1,cf_nextra)
          end if
* MOD TAH 100902: Test for newer 1040 version of cfiles        
      elseif( cf_nversn.lt. 1040 ) then  
*         Read the 1020 version of the cfile record.      
          if( loc_opt(1:2).eq.'DA' ) then
             read(unit,iostat=ierr) flag, cf_sitnam, cf_rcvrsw, 
     .           cf_swver, cf_npart, cf_norb, cf_nsat,
     .           (cf_ischan(i),i=1,cf_nsat),
     .           cf_ndat, (cf_dattyp(i),i=1,cf_ndat),
     .           ((cf_lambda(i,j),i=1,cf_nsat),j=1,cf_ndat),
     .           cf_skd, cf_nepoch, cf_inter, cf_ircint, cf_mtime,
     .           cf_isessn, cf_iy, cf_im, cf_id, cf_ihr, cf_min, cf_sec
 
             cf_nslip = 0
             cf_nextra = 0
          else
             read(unit,iostat=ierr) flag, cf_sitnam, cf_rctype,
     .           cf_rcvnum, cf_rcvrsw, cf_swver,
     .           cf_anttyp, cf_antnum, cf_npart, cf_norb,
     .           cf_nsat, (cf_ischan(i),i=1,cf_nsat),
     .           cf_ndat, (cf_dattyp(i),i=1,cf_ndat),
     .           ((cf_lambda(i,j),i=1,cf_nsat),j=1,cf_ndat),
     .           cf_skd, cf_nepoch, cf_inter, cf_ircint, cf_mtime,
     .           cf_isessn, cf_iy, cf_im, cf_id, cf_ihr, cf_min, cf_sec,
     .           (cf_offarp(i),i=1,3),
     .           (cf_offsL1(i),i=1,3), (cf_offsL2(i),i=1,3),
* MOD TAH 200126: Read just the L1 satellite PCO values.
     .           (cf_svantdx(:,1,j), j=1,cf_nsat),
     .           cf_obfiln, cf_tfiln, cf_jfiln,
     .           cf_frame, cf_precmod, cf_nutmod, cf_gravmod,
     .           cf_srpmod, cf_ietide, cf_isptide, 
     .           cf_speopmod, cf_etidemod, cf_otidemod, cf_atmtide,
     .           cf_atmlmod, cf_hydrolmod, 
     .           (cf_atmlavg(i),i=1,3), (cf_hydrolavg(i),i=1,3),
     .           cf_drymap, cf_wetmap, cf_antmod,
     .           (cf_iblk(i),i=1,cf_nsat),(cf_svantmod(i),i=1,cf_nsat),
     .           cf_elvcut, cf_nclock, (cf_clock(i),i=1,cf_nclock),
     .           cf_jde, cf_te, cf_jdr, cf_tr,
     .           cf_ut1, cf_xp, cf_yp,  cf_psi,cf_eps,
     .           cf_avlmet, cf_serwgt,
     .           cf_nslip, (cf_islip(i),i=1,cf_nslip),
     .           (cf_slpst(i),i=1,cf_nslip),
     .           cf_niextra, (cf_iextra(i),i=1,cf_niextra),
     .           cf_nextra,  (cf_extra(i), i=1,cf_nextra),
     .           cf_ncextra, (cf_cextra(i),i=1,cf_ncextra)
          end if
* MOD TAH 100902: Read the new default version
      elseif( cf_nversn.lt. 1041 ) then 
*         Read the 1040 version of the cfile record.      
          if( loc_opt(1:2).eq.'DA' ) then
             read(unit,iostat=ierr) flag, cf_sitnam, cf_rcvrsw, 
     .           cf_swver, cf_npart, cf_norb, cf_nsat,
     .           (cf_ischan(i),i=1,cf_nsat),
     .           cf_ndat, (cf_dattyp(i),i=1,cf_ndat),
     .           ((cf_lambda(i,j),i=1,cf_nsat),j=1,cf_ndat),
     .           cf_skd, cf_nepoch, cf_inter, cf_ircint, cf_mtime,
     .           cf_isessn, cf_iy, cf_im, cf_id, cf_ihr, cf_min, cf_sec
 
             cf_nslip = 0
             cf_nextra = 0
         else
             read(unit,iostat=ierr) flag, cf_sitnam, cf_rctype,
     .           cf_rcvnum, cf_rcvrsw, cf_swver,
     .           cf_anttyp, cf_antnum, cf_npart, cf_norb,
     .           cf_nsat, (cf_ischan(i),i=1,cf_nsat),
     .           cf_ndat, (cf_dattyp(i),i=1,cf_ndat),
     .           ((cf_lambda(i,j),i=1,cf_nsat),j=1,cf_ndat),
     .           cf_skd, cf_nepoch, cf_inter, cf_ircint, cf_mtime,
     .           cf_isessn, cf_iy, cf_im, cf_id, cf_ihr, cf_min, cf_sec,
     .           (cf_offarp(i),i=1,3),
     .           (cf_offsL1(i),i=1,3), (cf_offsL2(i),i=1,3),
* MOD TAH 200126: Read just the L1 satellite PCO values.
     .           ((cf_svantdx(i,1,j),i=1,3), j=1,cf_nsat),
     .           cf_obfiln, cf_tfiln, cf_jfiln,
     .           cf_frame, cf_precmod, cf_nutmod, cf_gravmod,
     .           cf_srpmod, cf_ietide, cf_isptide, 
     .           cf_speopmod, cf_etidemod, cf_otidemod, cf_atmtide,
     .           cf_atmlmod, cf_hydrolmod, 
     .           (cf_atmlavg(i),i=1,3), (cf_hydrolavg(i),i=1,3),
     .           cf_drymap,cf_wetmap, cf_antmod_snx, cf_antmod,
     .           (cf_iblk(i),i=1,cf_nsat),
     .           (cf_svantmod_snx(i),i=1,cf_nsat), 
     .           (cf_svantmod(i),i=1,cf_nsat),
     .           cf_elvcut, cf_nclock, (cf_clock(i),i=1,cf_nclock),
     .           cf_jde, cf_te, cf_jdr, cf_tr,
     .           cf_ut1, cf_xp, cf_yp,  cf_psi,cf_eps,
     .           cf_avlmet, cf_serwgt,
     .           cf_nslip, (cf_islip(i),i=1,cf_nslip),
     .           (cf_slpst(i),i=1,cf_nslip),
     .           cf_niextra, (cf_iextra(i),i=1,cf_niextra),
     .           cf_nextra,  (cf_extra(i), i=1,cf_nextra),
     .           cf_ncextra, (cf_cextra(i),i=1,cf_ncextra)

          end if 
* MOD TAH 130119: Added for 10.41 version       
      elseif( cf_nversn.lt. 1042 ) then   
*         Read the 1041 version of the cfile record.(cf_dryzen 
*         cf_wetzen added).    
          if( loc_opt(1:2).eq.'DA' ) then
             read(unit,iostat=ierr) flag, cf_sitnam, cf_rcvrsw, 
     .           cf_swver, cf_npart, cf_norb, cf_nsat,
     .           (cf_ischan(i),i=1,cf_nsat),
     .           cf_ndat, (cf_dattyp(i),i=1,cf_ndat),
     .           ((cf_lambda(i,j),i=1,cf_nsat),j=1,cf_ndat),
     .           cf_skd, cf_nepoch, cf_inter, cf_ircint, cf_mtime,
     .           cf_isessn, cf_iy, cf_im, cf_id, cf_ihr, cf_min, cf_sec
 
             cf_nslip = 0
             cf_nextra = 0
         else
             read(unit,iostat=ierr) flag, cf_sitnam, cf_rctype,
     .           cf_rcvnum, cf_rcvrsw, cf_swver,
     .           cf_anttyp, cf_antnum, cf_npart, cf_norb,
     .           cf_nsat, (cf_ischan(i),i=1,cf_nsat),
     .           cf_ndat, (cf_dattyp(i),i=1,cf_ndat),
     .           ((cf_lambda(i,j),i=1,cf_nsat),j=1,cf_ndat),
     .           cf_skd, cf_nepoch, cf_inter, cf_ircint, cf_mtime,
     .           cf_isessn, cf_iy, cf_im, cf_id, cf_ihr, cf_min, cf_sec,
     .           (cf_offarp(i),i=1,3),
     .           (cf_offsL1(i),i=1,3), (cf_offsL2(i),i=1,3),
* MOD TAH 200126: Read just the L1 satellite PCO values.
     .           ((cf_svantdx(i,1,j),i=1,3), j=1,cf_nsat),
     .           cf_obfiln, cf_tfiln, cf_jfiln,
     .           cf_frame, cf_precmod, cf_nutmod, cf_gravmod,
     .           cf_srpmod, cf_ietide, cf_isptide, 
     .           cf_speopmod, cf_etidemod, cf_otidemod, cf_atmtide,
     .           cf_atmlmod, cf_hydrolmod, 
     .           (cf_atmlavg(i),i=1,3), (cf_hydrolavg(i),i=1,3),
     .           cf_dryzen, cf_wetzen, cf_drymap,cf_wetmap, 
     .           cf_antmod_snx, cf_antmod,
     .           (cf_iblk(i),i=1,cf_nsat),
     .           (cf_svantmod_snx(i),i=1,cf_nsat), 
     .           (cf_svantmod(i),i=1,cf_nsat),
     .           cf_elvcut, cf_nclock, (cf_clock(i),i=1,cf_nclock),
     .           cf_jde, cf_te, cf_jdr, cf_tr,
     .           cf_ut1, cf_xp, cf_yp,  cf_psi,cf_eps,
     .           cf_avlmet, cf_serwgt,
     .           cf_nslip, (cf_islip(i),i=1,cf_nslip),
     .           (cf_slpst(i),i=1,cf_nslip),
     .           cf_niextra, (cf_iextra(i),i=1,cf_niextra),
     .           cf_nextra,  (cf_extra(i), i=1,cf_nextra),
     .           cf_ncextra, (cf_cextra(i),i=1,cf_ncextra)
          endif
* MOD TAH 140327: Added for 10.42 version       
      elseif( cf_nversn.lt. 1060 ) then  
*         Read the 1042 version of the cfile record. Added 
*         cf_eradmod, cf_antradmod
          if( loc_opt(1:2).eq.'DA' ) then
             read(unit,iostat=ierr) flag, cf_sitnam, cf_rcvrsw, 
     .           cf_swver, cf_npart, cf_norb, cf_nsat,
     .           (cf_ischan(i),i=1,cf_nsat),
     .           cf_ndat, (cf_dattyp(i),i=1,cf_ndat),
     .           ((cf_lambda(i,j),i=1,cf_nsat),j=1,cf_ndat),
     .           cf_skd, cf_nepoch, cf_inter, cf_ircint, cf_mtime,
     .           cf_isessn, cf_iy, cf_im, cf_id, cf_ihr, cf_min, cf_sec
 
             cf_nslip = 0
             cf_nextra = 0
          else
             read(unit,iostat=ierr) flag, cf_sitnam, cf_rctype,
     .           cf_rcvnum, cf_rcvrsw, cf_swver,
     .           cf_anttyp, cf_antnum, cf_npart, cf_norb,
     .           cf_nsat, (cf_ischan(i),i=1,cf_nsat),
     .           cf_ndat, (cf_dattyp(i),i=1,cf_ndat),
     .           ((cf_lambda(i,j),i=1,cf_nsat),j=1,cf_ndat),
     .           cf_skd, cf_nepoch, cf_inter, cf_ircint, cf_mtime,
     .           cf_isessn, cf_iy, cf_im, cf_id, cf_ihr, cf_min, cf_sec,
     .           (cf_offarp(i),i=1,3),
     .           (cf_offsL1(i),i=1,3), (cf_offsL2(i),i=1,3),
* MOD TAH 200126: Read just the L1 satellite PCO values.
     .           ((cf_svantdx(i,1,j),i=1,3), j=1,cf_nsat),
     .           cf_obfiln, cf_tfiln, cf_jfiln,
     .           cf_frame, cf_precmod, cf_nutmod, cf_gravmod,
     .           cf_srpmod, cf_eradmod, cf_antradmod, 
     .           cf_ietide, cf_isptide, 
     .           cf_speopmod, cf_etidemod, cf_otidemod, cf_atmtide,
     .           cf_atmlmod, cf_hydrolmod, 
     .           (cf_atmlavg(i),i=1,3), (cf_hydrolavg(i),i=1,3),
     .           cf_dryzen, cf_wetzen, cf_drymap,cf_wetmap, 
     .           cf_ionsrc, cf_magfield, 
     .           cf_antmod_snx, cf_antmod,
     .           (cf_iblk(i),i=1,cf_nsat),
     .           (cf_svantmod_snx(i),i=1,cf_nsat), 
     .           (cf_svantmod(i),i=1,cf_nsat),
     .           cf_elvcut, cf_nclock, (cf_clock(i),i=1,cf_nclock),
     .           cf_jde, cf_te, cf_jdr, cf_tr,
     .           cf_ut1, cf_xp, cf_yp,  cf_psi,cf_eps,
     .           cf_avlmet, cf_serwgt,
     .           cf_nslip, (cf_islip(i),i=1,cf_nslip),
     .           (cf_slpst(i),i=1,cf_nslip),
     .           cf_niextra, (cf_iextra(i),i=1,cf_niextra),
     .           cf_nextra,  (cf_extra(i), i=1,cf_nextra),
     .           cf_ncextra, (cf_cextra(i),i=1,cf_ncextra)
          endif

* MOD TAH 200126: Read 1061 version of c-file. 
      elseif( cf_nversn.lt. 1071 ) then  
*         Read the 1060 version of the cfile record. Added cf_gnss, 
*         cf_fL1, cf_fL2, cv_svantbody; remove cf_iblk, cf_serwgt 
          if( loc_opt(1:2).eq.'DA' ) then
             read(unit,iostat=ierr) flag, cf_sitnam, cf_rctype,
     .           cf_rcvnum, cf_rcvrsw, cf_swver, 
     .           cf_anttyp, cf_antnum, cf_npart, cf_norb,
     .           cf_gnss, cf_nsat, (cf_ischan(i),i=1,cf_nsat), 
     .           (cf_fL1(i),cf_fL2(i),i=1,cf_nsat),
     .           cf_ndat, (cf_dattyp(i),i=1,cf_ndat),
     .           ((cf_lambda(i,j),i=1,cf_nsat),j=1,cf_ndat),
     .           cf_skd, cf_nepoch, cf_inter, cf_ircint, cf_mtime,
     .           cf_isessn, cf_iy, cf_im, cf_id, cf_ihr, cf_min, cf_sec
             cf_nslip = 0
             cf_nextra = 0
          else
             read(unit,iostat=ierr) flag, cf_sitnam, cf_rctype,
     .           cf_rcvnum, cf_rcvrsw, cf_swver,
     .           cf_anttyp, cf_antnum, cf_npart, cf_norb,
     .           cf_gnss, cf_nsat, (cf_ischan(i),i=1,cf_nsat),
     .           (cf_fL1(i),cf_fL2(i),i=1,cf_nsat),
     .           cf_ndat, (cf_dattyp(i),i=1,cf_ndat),
     .           ((cf_lambda(i,j),i=1,cf_nsat),j=1,cf_ndat),
     .           cf_skd, cf_nepoch, cf_inter, cf_ircint, cf_mtime,
     .           cf_isessn, cf_iy, cf_im, cf_id, cf_ihr, cf_min, cf_sec,
     .           (cf_offarp(i),i=1,3),
     .           (cf_offsL1(i),i=1,3), (cf_offsL2(i),i=1,3),
* MOD TAH 200126: Read just the L1 satellite PCO values.
     .           ((cf_svantdx(i,1,j),i=1,3), j=1,cf_nsat),
     .           cf_obfiln, cf_tfiln, cf_jfiln,
     .           cf_frame, cf_precmod, cf_nutmod, cf_gravmod,
     .           cf_srpmod, cf_eradmod, cf_antradmod, 
     .           cf_ietide, cf_isptide, 
     .           cf_speopmod, cf_etidemod, cf_otidemod, cf_atmtide,
     .           cf_atmlmod, cf_hydrolmod,
     .           (cf_atmlavg(i),i=1,3), (cf_hydrolavg(i),i=1,3),
     .           cf_dryzen, cf_wetzen, cf_drymap,cf_wetmap, 
     .           cf_ionsrc, cf_magfield, 
     .           cf_antmod_snx, cf_antmod,
     .           (cf_svantbody(i),i=1,cf_nsat),
     .           (cf_svantmod_snx(i),i=1,cf_nsat), 
     .           (cf_svantmod(i),i=1,cf_nsat),
     .           cf_elvcut, cf_nclock, (cf_clock(i),i=1,cf_nclock),
     .           cf_jde, cf_te, cf_jdr, cf_tr,
     .           cf_ut1, cf_xp, cf_yp,  cf_psi,cf_eps,
     .           cf_avlmet, 
     .           cf_nslip, (cf_islip(i),i=1,cf_nslip),
     .           (cf_slpst(i),i=1,cf_nslip),
     .           cf_niextra, (cf_iextra(i),i=1,cf_niextra),
     .           cf_nextra,  (cf_extra(i), i=1,cf_nextra),
     .           cf_ncextra, (cf_cextra(i),i=1,cf_ncextra)
          end if        
      else  
*         Read the 1071 version of the cfile record. L1/L2 satellte
*         PCO values cf_svantdx(3,2,nsat) 
          if( loc_opt(1:2).eq.'DA' ) then
             read(unit,iostat=ierr) flag, cf_sitnam, cf_rctype,
     .           cf_rcvnum, cf_rcvrsw, cf_swver, 
     .           cf_anttyp, cf_antnum, cf_npart, cf_norb,
     .           cf_gnss, cf_nsat, (cf_ischan(i),i=1,cf_nsat), 
     .           (cf_fL1(i),cf_fL2(i),i=1,cf_nsat),
     .           cf_ndat, (cf_dattyp(i),i=1,cf_ndat),
     .           ((cf_lambda(i,j),i=1,cf_nsat),j=1,cf_ndat),
     .           cf_skd, cf_nepoch, cf_inter, cf_ircint, cf_mtime,
     .           cf_isessn, cf_iy, cf_im, cf_id, cf_ihr, cf_min, cf_sec
             cf_nslip = 0
             cf_nextra = 0
          else
*  MOD TAH 200205: Added cf_antdaz and size of cf_snavtdx
             read(unit,iostat=ierr) flag, cf_sitnam, cf_rctype,
     .           cf_rcvnum, cf_rcvrsw, cf_swver,
     .           cf_anttyp, cf_antnum, cf_npart, cf_norb,
     .           cf_gnss, cf_nsat, (cf_ischan(i),i=1,cf_nsat),
     .           (cf_fL1(i),cf_fL2(i),i=1,cf_nsat),
     .           cf_ndat, (cf_dattyp(i),i=1,cf_ndat),
     .           ((cf_lambda(i,j),i=1,cf_nsat),j=1,cf_ndat),
     .           cf_skd, cf_nepoch, cf_inter, cf_ircint, cf_mtime,
     .           cf_isessn, cf_iy, cf_im, cf_id, cf_ihr, cf_min, cf_sec,
     .           (cf_offarp(i),i=1,3),
     .           (cf_offsL1(i),i=1,3), (cf_offsL2(i),i=1,3), cf_antdaz,
     .           (cf_svantdx(:,:,j), j=1,cf_nsat),
     .           cf_obfiln, cf_tfiln, cf_jfiln,
     .           cf_frame, cf_precmod, cf_nutmod, cf_gravmod,
     .           cf_srpmod, cf_eradmod, cf_antradmod, 
     .           cf_ietide, cf_isptide, 
     .           cf_speopmod, cf_etidemod, cf_otidemod, cf_atmtide,
     .           cf_atmlmod, cf_hydrolmod,
     .           (cf_atmlavg(i),i=1,3), (cf_hydrolavg(i),i=1,3),
     .           cf_dryzen, cf_wetzen, cf_drymap,cf_wetmap, 
     .           cf_ionsrc, cf_magfield, 
     .           cf_antmod_snx, cf_antmod,
     .           (cf_svantbody(i),i=1,cf_nsat),
     .           (cf_svantmod_snx(i),i=1,cf_nsat), 
     .           (cf_svantmod(i),i=1,cf_nsat),
     .           cf_elvcut, cf_nclock, (cf_clock(i),i=1,cf_nclock),
     .           cf_jde, cf_te, cf_jdr, cf_tr,
     .           cf_ut1, cf_xp, cf_yp,  cf_psi,cf_eps,
     .           cf_avlmet, 
     .           cf_nslip, (cf_islip(i),i=1,cf_nslip),
     .           (cf_slpst(i),i=1,cf_nslip),
     .           cf_niextra, (cf_iextra(i),i=1,cf_niextra),
     .           cf_nextra,  (cf_extra(i), i=1,cf_nextra),
     .           cf_ncextra, (cf_cextra(i),i=1,cf_ncextra)
          end if        
      end if

* MOD TAH 190711: Save a copy of the frequencies from the c-file if this 
*     GLONASS FDMA data
      do i = 1, cf_nsat
         sv_fL1(i) = cf_fL1(i)
         sv_fL2(i) = cf_fL2(i)
      end do

* MOD TAH 200126: If less than 1071 version, copy L1 PCO to L2 slots
      if( cf_nversn.lt. 1071 ) then
         do i = 1,cf_nsat
           cf_svantdx(:,2,i) = cf_svantdx(:,1,i)
         enddo
      endif

****  See what happened and check flags and sizes.
      if( ierr.eq.0 ) then
          jerr = 0
          call check_cflag( flag, 2, jerr )
          call check_csize( cf_nsat, cf_maxsat, 'MAX SAT', jerr )
          call check_csize( cf_ndat, cf_maxdat, 'MAX DAT', jerr )
          call check_csize( cf_nslip, cf_maxcsb, 'MAX SLIPS', jerr )
          call check_csize( cf_nextra, cf_maxext, 'MAX EXTRA', jerr )
          call check_csize( cf_nclock, cf_maxclk, 'MAX NCLOCK', jerr )
 
*         Save the jerr error in the return error flag
          ierr = jerr
      else
          call report_error('IOSTAT',ierr,'read',
     .        'Cfile Type 2 record', 0, 'READ_CF2')
      end if
 
****  Thats all
      return
      end
 
CTITLE READ_CF3
 
      subroutine read_cf3( unit, option, ierr )
                        
      implicit none 
 
*     This routine will read the type 3 record of a C-file.  It is
*     the user responsibility to act on any error which is returned.
*     Here the error will be reported but the program will continue
*     executing.  The option passed indicates how much of the record
*     will be read.
 
* INCLUDE FILES
 
      include '../includes/cfile_def.h'
 
* PASSED VARIABLES
 
*   unit        - Unit number to which the file is attached.
*   ierr        - IOSTAT error in general but if the record type
*           - is wrong, then -1000 is returned as the error.
 
      integer*4 unit, ierr
 
*   option  - Option for read:  There are
*           - DATA - Just read enough of this record so
*           -        the data can be read
*           - ALL   - Read complete record (default)
 
      character*(*) option
 
* LOCAL VARIABLES
 
*   flag        - First I*4 value in record and indicates the
*           - type of record.
*   jerr        - Error generated in check c-file
*           - -1000 - Type not correct
*           - -1001 - Dimension too small
*   i,j     - Loop counters
 
      integer*4 flag, jerr, i
 
*   loc_opt - Local version of option so that we can
*           - casefold it.
 
      character*4 loc_opt
 
***** Check the option and see what we should do:
      loc_opt = option
      call casefold(loc_opt)
 
      if( loc_opt(1:2).eq.'DA' ) then
          read(unit,iostat=ierr) flag
          cf_nlabel = 0
          cf_nparam = 0
      else
          read(unit,iostat=ierr) flag,cf_nlabel,cf_nparam,
     .        (cf_islot (i),i=1,cf_nlabel),
     .        (cf_idms  (i),i=1,cf_nlabel),
     .        (cf_rlabel(i),i=1,cf_nlabel),
     .        (cf_preval(i),i=1,cf_nparam)
 
      end if
 
****  See what happened and check flags and sizes.
      if( ierr.eq.0 ) then
          jerr = 0
          call check_cflag( flag, 3, jerr )
          call check_csize( cf_nlabel, cf_maxlab, 'MAX LABEL', jerr )
          call check_csize( cf_nparam, cf_maxprm, 'MAX PARAM', jerr )
 
*         Save the jerr error in the return error flag
          ierr = jerr
      else
          call report_error('IOSTAT',ierr,'read',
     .        'Cfile Type 3 record', 0, 'READ_CF3')
      end if
 
****  Thats all
      return
      end
 
CTITLE READ_CF4
 
      subroutine read_cf4( unit, option, ierr )
                        
      implicit none 
 
*     This routine will read the type 4 record of a C-file.  It is
*     the user responsibility to act on any error which is returned.
*     Here the error will be reported but the program will continue
*     executing.  The option passed indicates how much of the record
*     will be read.
 
* INCLUDE FILES
 
      include '../includes/cfile_def.h'
 
* PASSED VARIABLES
 
*   unit        - Unit number to which the file is attached.
*   ierr        - IOSTAT error in general but if the record type
*           - is wrong, then -1000 is returned as the error.
 
      integer*4 unit, ierr
 
*   option  - Option for read:  There are
*           - DATA - Just read enough of this record so
*           -        the data can be read
*           - ALL   - Read complete record (default)
 
      character*(*) option
 
* LOCAL VARIABLES
 
*   flag        - First I*4 value in record and indicates the
*           - type of record.
*   jerr        - Error generated in check c-file
*           - -1000 - Type not correct
*           - -1001 - Dimension too small
*   i,j     - Loop counters
 
      integer*4 flag, jerr, i
 
*   loc_opt - Local version of option so that we can
*           - casefold it.
 
      character*4 loc_opt
 
***** Check the option and see what we should do:
      loc_opt = option
      call casefold(loc_opt)
* MOD TAH 100902: Add in test for 1040 version where cf_atmlod added 
      if( cf_nversn.lt. 1040 ) then 
         if( loc_opt(1:2).eq.'DA' ) then
             read(unit,iostat=ierr) flag,
     .           cf_msat, (cf_mprn(i),i=1,cf_msat),
     .           cf_iepoch, cf_iyr, cf_idoy, cf_sod, cf_rclock,
     .           cf_okmet, cf_pres, cf_temp, cf_relhum
             cf_nsave = 0
         else
             read(unit,iostat=ierr) flag,
     .           cf_msat, (cf_mprn(i),i=1,cf_msat),
     .           cf_iepoch, cf_iyr, cf_idoy, cf_sod, cf_rclock,
     .           cf_okmet, cf_pres, cf_temp, cf_relhum,
     .           cf_nsave, (cf_save(i),i=1,cf_nsave),
     .           cf_kflag, cf_ksite, cf_klatr, cf_klonr, cf_kradk,
     .           cf_L1Z, cf_L1N, cf_L1E, cf_L2Z, cf_L2N, cf_L2E,
     .           cf_antaz
         endif
      elseif( cf_nversn.lt. 1060 ) then
         if( loc_opt(1:2).eq.'DA' ) then
             read(unit,iostat=ierr) flag,
     .           cf_msat, (cf_mprn(i),i=1,cf_msat),
     .           cf_iepoch, cf_iyr, cf_idoy, cf_sod, cf_rclock,
     .           cf_okmet, cf_pres, cf_temp, cf_relhum
             cf_nsave = 0
         else
             read(unit,iostat=ierr) flag,
     .           cf_msat, (cf_mprn(i),i=1,cf_msat),
     .           cf_iepoch, cf_iyr, cf_idoy, cf_sod, cf_rclock,
     .           cf_okmet, cf_pres, cf_temp, cf_relhum, cf_atmlod,
     .           cf_nsave, (cf_save(i),i=1,cf_nsave),
     .           cf_kflag, cf_ksite, cf_klatr, cf_klonr, cf_kradk,
     .           cf_L1Z, cf_L1N, cf_L1E, cf_L2Z, cf_L2N, cf_L2E,
     .           cf_antaz
          endif    
                                                
      else   
* MOD RWK 141215; Read 1060 version adding cf_zendel, removing kflag, 
*                 renaming ksite, klatr, klonr, krad to sitecd. latr_sph, 
*                 lonr, radius, and reordering the last 2 lines.
         if( loc_opt(1:2).eq.'DA' ) then
             read(unit,iostat=ierr) flag,
     .           cf_msat, (cf_mprn(i),i=1,cf_msat),
     .           cf_iepoch, cf_iyr, cf_idoy, cf_sod, cf_rclock,
     .           cf_okmet, cf_pres, cf_temp, cf_relhum
             cf_nsave = 0
         else
             read(unit,iostat=ierr) flag,
     .           cf_msat, (cf_mprn(i),i=1,cf_msat),
     .           cf_iepoch, cf_iyr, cf_idoy, cf_sod, cf_rclock,
     .           cf_zendel, cf_okmet, cf_pres, cf_temp, cf_relhum,
     .           cf_atmlod, cf_sitecd, cf_latr_sph, cf_lonr, cf_radius,  
     .           cf_L1Z, cf_L1N, cf_L1E, cf_L2Z, cf_L2N, cf_L2E,
     .           cf_antaz, cf_nsave, (cf_save(i),i=1,cf_nsave)
          endif 

      end if
 
****  See what happened and check flags and sizes.
      if( ierr.eq.0 ) then
          jerr = 0
          call check_cflag( flag, 4, jerr )
          call check_csize( cf_msat, cf_maxsat, 'MAX SAT', jerr )
          call check_csize( cf_nsave, cf_maxsav, 'MAX SAVE', jerr )
 
*         Save the jerr error in the return error flag
          ierr = jerr
      else
          call report_error('IOSTAT',ierr,'read',
     .        'Cfile Type 4 record', 0, 'READ_CF4')
      end if
 
****  Thats all
      return
      end
 
 
CTITLE READ_CF5
 
      subroutine read_cf5( unit, option, ierr )
                        
      implicit none 
 
*     This routine will read the type 5 record of a C-file.  It is
*     the user responsibility to act on any error which is returned.
*     Here the error will be reported but the program will continue
*     executing.  The option passed indicates how much of the record
*     will be read.
 
* INCLUDE FILES
 
      include '../includes/cfile_def.h'
      include '../includes/const_param.h'  ! Debug temporary (RF_CFACT)

* PASSED VARIABLES
 
*   unit        - Unit number to which the file is attached.
*   ierr        - IOSTAT error in general but if the record type
*           - is wrong, then -1000 is returned as the error.
 
      integer*4 unit, ierr
 
*   option  - Option for read:  There are
*           - DATA - Just read enough of this record so
*           -        the data can be read
*           - ALL   - Read complete record (default)
 
      character*(*) option
 
* LOCAL VARIABLES
 
*   flag        - First I*4 value in record and indicates the
*           - type of record.
*   jerr        - Error generated in check c-file
*           - -1000 - Type not correct
*           - -1001 - Dimension too small
*   i,j     - Loop counters
*   cf_ndats4   - I*4 version for passing to check routine
*   cf_nspare4  - I*4 version for passing to check routine
*   cf_nparts4  - I*4 version for passing to check routine
 
      integer*4 flag, jerr, i, cf_ndats4, cf_nspare4, cf_nparts4

* cf_prntolv - Function to return list number from PRN number.
*     (prntol array is set in ctogobs_comm so not used here). 
* lv  -- Satellite sequential number 
      integer*4 cf_prntolv, lv

* GLOSNASS reference frencies: (These are svaed in ctoobs_com
*     bu are not included here so save locally 

* fL1_R0 -- GLONASS L1 center Hz
* fL2_R0 -- GLONASS L2 center Hz
      real*8 fL1_R0, fL2_R0

*   loc_opt - Local version of option so that we can
*           - casefold it.
 
      character*4 loc_opt

* Data statements for GLONASS frequencies

      data fL1_R0 / 1602.0d6 /, fL2_R0 / 1246.0d6 / 

      save fL1_R0, fL2_R0
 
***** Check the option and see what we should do:
      loc_opt = option
      call casefold(loc_opt)

      if( cf_nversn.lt.980 ) then 
          if( loc_opt(1:2).eq.'DA' ) then
              read(unit,iostat=ierr) flag, cf_iprn,cf_elev, cf_azimuth,
     .            cf_okwvr, cf_wvrdel, cf_atmdel, cf_tau, cf_drate,
     .            cf_ierfl, cf_data_flag, cf_ndats,
     .            (cf_obsv(i), i=1, cf_ndats),
     .            (cf_omcs(i), i=1, cf_ndats),
     .            (cf_isnr(i), i=1, cf_ndats)
              cf_nspare = 0
              cf_nparts = 0
          else
              read(unit,iostat=ierr) flag, cf_iprn,cf_elev, cf_azimuth,
     .            cf_okwvr, cf_wvrdel, cf_atmdel, cf_tau, cf_drate,
     .            cf_ierfl, cf_data_flag, cf_ndats,
     .            (cf_obsv(i), i=1, cf_ndats),
     .            (cf_omcs(i), i=1, cf_ndats),
     .            (cf_isnr(i), i=1, cf_ndats), cf_ampl1, cf_ampl2,
     .            (cf_obswgt(i), i=1, cf_ndats),
     .            cf_nspare,( cf_spare(i), i=1, cf_nspare),
     .            cf_nparts,( cf_tmpart(i), i=1, cf_nparts)
          end if
          cf_svcepc = 0.d0
          cf_svcL1  = 0.d0
      elseif (cf_nversn.lt.1060 ) then
          if( loc_opt(1:2).eq.'DA' ) then
              read(unit,iostat=ierr) flag, cf_iprn,cf_elev, cf_azimuth,
     .            cf_okwvr, cf_wvrdel, cf_atmdel, 
     .            cf_svcepc, cf_svcL1, cf_tau, cf_drate,
     .            cf_ierfl, cf_data_flag, cf_ndats,
     .            (cf_obsv(i), i=1, cf_ndats),
     .            (cf_omcs(i), i=1, cf_ndats),
     .            (cf_isnr(i), i=1, cf_ndats)
              cf_nspare = 0
              cf_nparts = 0
          else
              read(unit,iostat=ierr) flag, cf_iprn,cf_elev, cf_azimuth,
     .            cf_okwvr, cf_wvrdel, cf_atmdel, 
     .            cf_svcepc, cf_svcL1, cf_tau, cf_drate,
     .            cf_ierfl, cf_data_flag, cf_ndats,
     .            (cf_obsv(i), i=1, cf_ndats),
     .            (cf_omcs(i), i=1, cf_ndats),
     .            (cf_isnr(i), i=1, cf_ndats), cf_ampl1, cf_ampl2,
     .            (cf_obswgt(i), i=1, cf_ndats),
     .            cf_nspare,( cf_spare(i), i=1, cf_nspare),
     .            cf_nparts,( cf_tmpart(i), i=1, cf_nparts)
          end if                      
      else      
* MOD RWK 141215: Read 1060 version,adding cf_nadang, and removing cf_okwvr, cf_wvrdel, cf_obswgt
         if( loc_opt(1:2).eq.'DA' ) then
              read(unit,iostat=ierr) flag, cf_iprn,cf_elev, cf_azimuth,
     .            cf_nadang, cf_atmdel,
     .            cf_svcepc, cf_svcL1, cf_tau, cf_drate,
     .            cf_ierfl, cf_data_flag, cf_ndats,
     .            (cf_obsv(i), i=1, cf_ndats),
     .            (cf_omcs(i), i=1, cf_ndats),
     .            (cf_isnr(i), i=1, cf_ndats)
              cf_nspare = 0
              cf_nparts = 0
          else
              read(unit,iostat=ierr) flag, cf_iprn,cf_elev, cf_azimuth,
     .            cf_nadang, cf_atmdel, 
     .            cf_svcepc, cf_svcL1, cf_tau, cf_drate,
     .            cf_ierfl, cf_data_flag, cf_ndats,
     .            (cf_obsv(i), i=1, cf_ndats),
     .            (cf_omcs(i), i=1, cf_ndats),
     .            (cf_isnr(i), i=1, cf_ndats), cf_ampl1, cf_ampl2,
     .            cf_nspare,( cf_spare(i), i=1, cf_nspare),
     .            cf_nparts,( cf_tmpart(i), i=1, cf_nparts)
          end if  

* MOD TAH 180311: Remap the GLOSNASS frequencies to 1 one frequency.  This
*        will allow the standard estimates of clocks etc.  Ultimately, factors
*        used in compuring LC, EX-WL and MW-WL will be modfified to accoint 
*        for this.
         if( cf_gnss.eq.'R' ) then

***         Loop over cf_omcs and ch_obsv (phase values)
            if( cf_ndats.ne.4 .and. cf_iepoch.eq.1 ) then
               write(*,310) cf_ndats
 310           format('**WARNING** Only mapping 4 for ',i2,
     .             ' data types for GLONASS')
            endif
            lv = cf_prntolv(cf_iprn) 

C MOD TAH 190711: Use the saved versions of the frequencies because
C           cf_fl1 and cf_fl2 are replaced with constant values when
C           cfiles are written out.
            cf_omcs(1) = cf_omcs(1)*fL1_R0/sv_fl1(lv)
            cf_omcs(2) = cf_omcs(2)*fL2_R0/sv_fl2(lv)
            cf_omcs(3) = cf_omcs(3)*fL1_R0/sv_fl1(lv)
            cf_omcs(4) = cf_omcs(4)*fL2_R0/sv_fl2(lv)
*           Now phase observations.
            cf_obsv(1) = cf_obsv(1)*fL1_R0/sv_fl1(lv)
            cf_obsv(2) = cf_obsv(2)*fL2_R0/sv_fl2(lv)
         end if
!        if( unit.le.102 )
!       .write(*,200) cf_iprn, cf_elev(1)*180/pi, cf_azimuth(1)*180/pi,
!    .                cf_omcs(1:cf_ndats)
!200     format('RD_CFACT ',I3,2F8.2,1x,10(F20.4,1x))
      end if

****  See what happened and check flags and sizes.
*     Assign I*2 size variables to I*4 one for processing in checks.
      cf_ndats4  = cf_ndats
      cf_nspare4 = cf_nspare
      cf_nparts4 = cf_nparts
 
      if( ierr.eq.0 ) then
          jerr = 0
          call check_cflag( flag, 5, jerr )
          call check_csize( cf_ndats4, cf_maxdat, 'MAX DAT', jerr )
          call check_csize( cf_nspare4, cf_maxspr, 'MAX SPARE', jerr )
          call check_csize( cf_nparts4, cf_maxlab, 'MAX LABELS', jerr )
 
*         Save the jerr error in the return error flag
          ierr = jerr
      else
          call report_error('IOSTAT',ierr,'read',
     .        'Cfile Type 5 record', 0, 'READ_CF5')
      end if
 
****  Thats all
      return
      end
 
CTITLE CF_PRNTOLV 
 
      integer*4 function cf_prntolv(prn) 

      implicit none

*     Function to return list number for PRN number using cf_ischan
*     array.
 
* INCLUDE FILES
 
      include '../includes/cfile_def.h'
 
* PASSED VARIABLES

      integer*2 prn   ! PRN number of vechicle.

* LOCAL VARUABLES
      integer*4 k     ! Loop couunter

****  Search ovef cf_ischan to find PRN number
      cf_prntolv = 0 
      do k = 1, cf_nsat
         if ( prn.eq.cf_ischan(k)) then
            cf_prntolv = k
            exit
         endif
      end do

*     Make sure found PRN
      if( cf_prntolv.eq.0 ) then
         call report_stat('FATAL','autcln','cf_perntolv', 
     .          cf_sitnam,'PRN Not found in cf_ischan',prn)
      endif

****  That all
      return 
      end


  

 
 
