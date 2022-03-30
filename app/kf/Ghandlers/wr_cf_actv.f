 
CTITLE WRITE_CF1
 
      subroutine write_cf1( unit, option, ierr )
                        
      implicit none 
 
*     This routine will write the type 1 record of a C-file.  It is
*     the user responsibility to act on any error which is returned.
*     Here the error will be reported but the program will continue
*     executing.  The option passed indicates how much of the record
*     will be write.
 
* INCLUDE FILES
 
      include '../includes/cfile_def.h'
 
* PASSED VARIABLES
 
*   unit        - Unit number to which the file is attached.
*   ierr        - IOSTAT error in general but if the record type
*           - is wrong, then -1000 is returned as the error.
 
      integer*4 unit, ierr
 
*   option  - Option for write:  NONE
 
      character*(*) option
 
* LOCAL VARIABLES
 
*   flag        - First I*4 value in record and indicates the
*           - type of record.
*   i       - Loop counter
 
      integer*4 flag, i
 
*   loc_opt - Local version of option so that we can
*           - casefold it.
 
      character*4 loc_opt
 
***** Check the option and see what we should do:
      loc_opt = option
      call casefold(loc_opt)

      flag = 1
 
      write(unit,iostat=ierr) flag, cf_nversn, cf_ntext,
     .        (cf_text(i), i=1, cf_ntext )
 
      call report_error('IOSTAT',ierr,'write',
     .        'Cfile Type 1 record', 0, 'WRITE_CF1')
 
****  Thats all
      return
      end
 
CTITLE WRITE_CF2
 
      subroutine write_cf2( unit, option, ierr )
                        
      implicit none 
 
*     This routine will write the type 2 record of a C-file.  It is
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
 
*   option  - Option for write:  NONE
 
      character*(*) option
 
* LOCAL VARIABLES
 
*   flag        - First I*4 value in record and indicates the
*   i,j     - Loop counters
 
      integer*4 flag, i,j
 
*   loc_opt - Local version of option so that we can
*           - casefold it.
 
      character*4 loc_opt

* MOD TAH 180311: Add GLONASS frequencies (valuaes are saved in
*     ctogobs_com.f and not in cfile_def.f
* fL1_R0 -- GLONASS L1 center Hz
* fL2_R0 -- GLONASS L2 center Hz
      real*8 fL1_R0, fL2_R0

* Data statements for GLONASS frequencies

      data fL1_R0 / 1602.0d6 /, fL2_R0 / 1246.0d6 / 

      save fL1_R0, fL2_R0
 
***** Check the option and see what we should do:
      loc_opt = option
      call casefold(loc_opt)

      flag = 2

*     Check the version of the cfile
      if( cf_nversn.lt.980 ) then  
         write(unit,iostat=ierr) flag, cf_sitnam, cf_rctype,
     .        cf_rcvnum, cf_rcvrsw, cf_swver,
     .        cf_anttyp, cf_antnum,
     .        cf_npart, cf_nsat, (cf_ischan(i),i=1,cf_nsat),
     .        cf_ndat, (cf_dattyp(i),i=1,cf_ndat),
     .        ((cf_lambda(i,j),i=1,cf_nsat),j=1,cf_ndat),
     .        cf_skd, cf_nepoch, cf_inter, cf_ircint, cf_mtime,
     .        cf_isessn, cf_iy, cf_im, cf_id, cf_ihr, cf_min, cf_sec,
     .        (cf_offarp(i),i=1,3),
     .        (cf_offsL1(i),i=1,3), (cf_offsL2(i),i=1,3),
     .        cf_obfiln, cf_tfiln, cf_jfiln,
     .        cf_frame, cf_precmod, cf_nutmod, cf_gravmod,
     .        cf_srpmod, cf_ietide, cf_isptide, cf_antmod,
     .        cf_elvcut, cf_nclock, (cf_clock(i),i=1,cf_nclock),
     .        cf_jde, cf_te, cf_jdr, cf_tr,
     .        cf_ut1, cf_xp, cf_yp,  cf_psi,cf_eps,
     .        cf_avlmet, cf_serwgt,
     .        cf_nslip, (cf_islip(i),i=1,cf_nslip),
     .        (cf_slpst(i),i=1,cf_nslip),
     .        cf_nextra, (cf_extra(i),i=1,cf_nextra)
      elseif( cf_nversn.lt. 1020 ) then
*        Write the 980 version of the cfile.
         write(unit,iostat=ierr) flag, cf_sitnam, cf_rctype,
     .        cf_rcvnum, cf_rcvrsw, cf_swver,
     .        cf_anttyp, cf_antnum,cf_npart, cf_norb, cf_nsat,
     .        (cf_ischan(i),i=1,cf_nsat),
     .        cf_ndat, (cf_dattyp(i),i=1,cf_ndat),
     .        ((cf_lambda(i,j),i=1,cf_nsat),j=1,cf_ndat),
     .        cf_skd, cf_nepoch, cf_inter, cf_ircint, cf_mtime,
     .        cf_isessn, cf_iy, cf_im, cf_id, cf_ihr, cf_min, cf_sec,
     .        (cf_offarp(i),i=1,3),
     .        (cf_offsL1(i),i=1,3), (cf_offsL2(i),i=1,3),
     .        cf_obfiln, cf_tfiln, cf_jfiln,
     .        cf_frame, cf_precmod, cf_nutmod, cf_gravmod,
     .        cf_srpmod, cf_ietide, cf_isptide, cf_antmod,
     .        cf_elvcut, cf_nclock, (cf_clock(i),i=1,cf_nclock),
     .        cf_jde, cf_te, cf_jdr, cf_tr,
     .        cf_ut1, cf_xp, cf_yp,  cf_psi,cf_eps,
     .        cf_avlmet, cf_serwgt,
     .        cf_nslip, (cf_islip(i),i=1,cf_nslip),
     .        (cf_slpst(i),i=1,cf_nslip),
     .        cf_nextra, (cf_extra(i),i=1,cf_nextra)
* MOD TAH 100902: Test for newer 1040 version of cfiles        
      elseif( cf_nversn.lt. 1040 ) then  
*         Write the 1020 version of the cfile record.      
          write(unit,iostat=ierr) flag, cf_sitnam, cf_rctype,
     .        cf_rcvnum, cf_rcvrsw, cf_swver,
     .        cf_anttyp, cf_antnum, cf_npart, cf_norb,
     .        cf_nsat, (cf_ischan(i),i=1,cf_nsat),
     .        cf_ndat, (cf_dattyp(i),i=1,cf_ndat),
     .        ((cf_lambda(i,j),i=1,cf_nsat),j=1,cf_ndat),
     .        cf_skd, cf_nepoch, cf_inter, cf_ircint, cf_mtime,
     .        cf_isessn, cf_iy, cf_im, cf_id, cf_ihr, cf_min, cf_sec,
     .        (cf_offarp(i),i=1,3),
     .        (cf_offsL1(i),i=1,3), (cf_offsL2(i),i=1,3),
* MOD TAH 200126: Write just the L1 PCO values
     .        ((cf_svantdx(i,1,j),i=1,3), j=1,cf_nsat),
     .        cf_obfiln, cf_tfiln, cf_jfiln,
     .        cf_frame, cf_precmod, cf_nutmod, cf_gravmod,
     .        cf_srpmod, cf_ietide, cf_isptide, 
     .        cf_speopmod, cf_etidemod, cf_otidemod, cf_atmtide,
     .        cf_atmlmod, cf_hydrolmod, 
     .        (cf_atmlavg(i),i=1,3), (cf_hydrolavg(i),i=1,3),
     .        cf_drymap,cf_wetmap, cf_antmod,
     .        (cf_iblk(i),i=1,cf_nsat),(cf_svantmod(i),i=1,cf_nsat),
     .        cf_elvcut, cf_nclock, (cf_clock(i),i=1,cf_nclock),
     .        cf_jde, cf_te, cf_jdr, cf_tr,
     .        cf_ut1, cf_xp, cf_yp,  cf_psi,cf_eps,
     .        cf_avlmet, cf_serwgt,
     .        cf_nslip, (cf_islip(i),i=1,cf_nslip),
     .        (cf_slpst(i),i=1,cf_nslip),
     .        cf_niextra, (cf_iextra(i),i=1,cf_niextra),
     .        cf_nextra,  (cf_extra(i), i=1,cf_nextra),
     .        cf_ncextra, (cf_cextra(i),i=1,cf_ncextra)
      elseif( cf_nversn.lt. 1041 ) then 
*         Write the 1040 version of the cfile record.      
          write(unit,iostat=ierr) flag, cf_sitnam, cf_rctype,
     .        cf_rcvnum, cf_rcvrsw, cf_swver,
     .        cf_anttyp, cf_antnum, cf_npart, cf_norb,
     .        cf_nsat, (cf_ischan(i),i=1,cf_nsat),
     .        cf_ndat, (cf_dattyp(i),i=1,cf_ndat),
     .        ((cf_lambda(i,j),i=1,cf_nsat),j=1,cf_ndat),
     .        cf_skd, cf_nepoch, cf_inter, cf_ircint, cf_mtime,
     .        cf_isessn, cf_iy, cf_im, cf_id, cf_ihr, cf_min, cf_sec,
     .        (cf_offarp(i),i=1,3),
     .        (cf_offsL1(i),i=1,3), (cf_offsL2(i),i=1,3),
* MOD TAH 200126: Write just the L1 PCO values
     .        ((cf_svantdx(i,1,j),i=1,3), j=1,cf_nsat),
     .        cf_obfiln, cf_tfiln, cf_jfiln,
     .        cf_frame, cf_precmod, cf_nutmod, cf_gravmod,
     .        cf_srpmod, cf_ietide, cf_isptide, 
     .        cf_speopmod, cf_etidemod, cf_otidemod, cf_atmtide,
     .        cf_atmlmod, cf_hydrolmod, 
     .        (cf_atmlavg(i),i=1,3), (cf_hydrolavg(i),i=1,3),
     .        cf_drymap,cf_wetmap, cf_antmod_snx, cf_antmod,
     .        (cf_iblk(i),i=1,cf_nsat),
     .        (cf_svantmod_snx(i),i=1,cf_nsat), 
     .        (cf_svantmod(i),i=1,cf_nsat),
     .        cf_elvcut, cf_nclock, (cf_clock(i),i=1,cf_nclock),
     .        cf_jde, cf_te, cf_jdr, cf_tr,
     .        cf_ut1, cf_xp, cf_yp,  cf_psi,cf_eps,
     .        cf_avlmet, cf_serwgt,
     .        cf_nslip, (cf_islip(i),i=1,cf_nslip),
     .        (cf_slpst(i),i=1,cf_nslip),
     .        cf_niextra, (cf_iextra(i),i=1,cf_niextra),
     .        cf_nextra,  (cf_extra(i), i=1,cf_nextra),
     .        cf_ncextra, (cf_cextra(i),i=1,cf_ncextra)
      elseif( cf_nversn.lt. 1042 ) then 
*         Write the 1041 version of the cfile record.      
          write(unit,iostat=ierr) flag, cf_sitnam, cf_rctype,
     .        cf_rcvnum, cf_rcvrsw, cf_swver,
     .        cf_anttyp, cf_antnum, cf_npart, cf_norb,
     .        cf_nsat, (cf_ischan(i),i=1,cf_nsat),
     .        cf_ndat, (cf_dattyp(i),i=1,cf_ndat),
     .        ((cf_lambda(i,j),i=1,cf_nsat),j=1,cf_ndat),
     .        cf_skd, cf_nepoch, cf_inter, cf_ircint, cf_mtime,
     .        cf_isessn, cf_iy, cf_im, cf_id, cf_ihr, cf_min, cf_sec,
     .        (cf_offarp(i),i=1,3),
     .        (cf_offsL1(i),i=1,3), (cf_offsL2(i),i=1,3),
* MOD TAH 200126: Write just the L1 PCO values
     .        ((cf_svantdx(i,1,j),i=1,3), j=1,cf_nsat),
     .        cf_obfiln, cf_tfiln, cf_jfiln,
     .        cf_frame, cf_precmod, cf_nutmod, cf_gravmod,
     .        cf_srpmod, cf_ietide, cf_isptide, 
     .        cf_speopmod, cf_etidemod, cf_otidemod, cf_atmtide,
     .        cf_atmlmod, cf_hydrolmod, 
     .        (cf_atmlavg(i),i=1,3), (cf_hydrolavg(i),i=1,3),
     .        cf_dryzen,cf_wetzen, cf_drymap,cf_wetmap, 
     .        cf_antmod_snx, cf_antmod,
     .        (cf_iblk(i),i=1,cf_nsat),
     .        (cf_svantmod_snx(i),i=1,cf_nsat), 
     .        (cf_svantmod(i),i=1,cf_nsat),
     .        cf_elvcut, cf_nclock, (cf_clock(i),i=1,cf_nclock),
     .        cf_jde, cf_te, cf_jdr, cf_tr,
     .        cf_ut1, cf_xp, cf_yp,  cf_psi,cf_eps,
     .        cf_avlmet, cf_serwgt,
     .        cf_nslip, (cf_islip(i),i=1,cf_nslip),
     .        (cf_slpst(i),i=1,cf_nslip),
     .        cf_niextra, (cf_iextra(i),i=1,cf_niextra),
     .        cf_nextra,  (cf_extra(i), i=1,cf_nextra),
     .        cf_ncextra, (cf_cextra(i),i=1,cf_ncextra)  
      elseif( cf_nversn.lt. 1060 ) then   
* MOD TAH 140327: Write the 1042 version of the cfile record.  
         write(unit,iostat=ierr) flag, cf_sitnam, cf_rctype,
     .        cf_rcvnum, cf_rcvrsw, cf_swver,
     .        cf_anttyp, cf_antnum, cf_npart, cf_norb,
     .        cf_nsat, (cf_ischan(i),i=1,cf_nsat),
     .        cf_ndat, (cf_dattyp(i),i=1,cf_ndat),
     .        ((cf_lambda(i,j),i=1,cf_nsat),j=1,cf_ndat),
     .        cf_skd, cf_nepoch, cf_inter, cf_ircint, cf_mtime,
     .        cf_isessn, cf_iy, cf_im, cf_id, cf_ihr, cf_min, cf_sec,
     .        (cf_offarp(i),i=1,3),
     .        (cf_offsL1(i),i=1,3), (cf_offsL2(i),i=1,3),
* MOD TAH 200126: Write just the L1 PCO values
     .        ((cf_svantdx(i,1,j),i=1,3), j=1,cf_nsat),
     .        cf_obfiln, cf_tfiln, cf_jfiln,
     .        cf_frame, cf_precmod, cf_nutmod, cf_gravmod,
     .        cf_srpmod, cf_eradmod, cf_antradmod, 
     .        cf_ietide, cf_isptide, 
     .        cf_speopmod, cf_etidemod, cf_otidemod, cf_atmtide,
     .        cf_atmlmod, cf_hydrolmod, 
     .        (cf_atmlavg(i),i=1,3), (cf_hydrolavg(i),i=1,3),
     .        cf_dryzen,cf_wetzen, cf_drymap,cf_wetmap, 
     .        cf_ionsrc, cf_magfield, 
     .        cf_antmod_snx, cf_antmod,
     .        (cf_iblk(i),i=1,cf_nsat),
     .        (cf_svantmod_snx(i),i=1,cf_nsat), 
     .        (cf_svantmod(i),i=1,cf_nsat),
     .        cf_elvcut, cf_nclock, (cf_clock(i),i=1,cf_nclock),
     .        cf_jde, cf_te, cf_jdr, cf_tr,
     .        cf_ut1, cf_xp, cf_yp,  cf_psi,cf_eps,
     .        cf_avlmet, cf_serwgt,
     .        cf_nslip, (cf_islip(i),i=1,cf_nslip),
     .        (cf_slpst(i),i=1,cf_nslip),
     .        cf_niextra, (cf_iextra(i),i=1,cf_niextra),
     .        cf_nextra,  (cf_extra(i), i=1,cf_nextra),
     .        cf_ncextra, (cf_cextra(i),i=1,cf_ncextra)
      elseif( cf_nversn.lt.1071 ) then    
* MOD RWK 141215: Write the 1060 version of the cfile record. 
* MOD TAH 180311: For GLONASS set the frequencies back to constant
*         values so the cview will not re-map
          if( cf_gnss.eq.'R' ) then
             do i = 1, cf_nsat
                cf_fL1(i) = fL1_R0
                cf_fl2(i) = fL2_R0
             end do
          end if
          
          write(unit,iostat=ierr) flag, cf_sitnam, cf_rctype,
     .        cf_rcvnum, cf_rcvrsw, cf_swver,
     .        cf_anttyp, cf_antnum, cf_npart, cf_norb,
     .        cf_gnss, cf_nsat, (cf_ischan(i),i=1,cf_nsat),
     .        (cf_fL1(i),cf_fl2(i),i=1,cf_nsat),
     .        cf_ndat, (cf_dattyp(i),i=1,cf_ndat),
     .        ((cf_lambda(i,j),i=1,cf_nsat),j=1,cf_ndat),
     .        cf_skd, cf_nepoch, cf_inter, cf_ircint, cf_mtime,
     .        cf_isessn, cf_iy, cf_im, cf_id, cf_ihr, cf_min, cf_sec,
     .        (cf_offarp(i),i=1,3),
     .        (cf_offsL1(i),i=1,3), (cf_offsL2(i),i=1,3),
* MOD TAH 200126: Write just the L1 PCO values
     .        ((cf_svantdx(i,1,j),i=1,3), j=1,cf_nsat),
     .        cf_obfiln, cf_tfiln, cf_jfiln,
     .        cf_frame, cf_precmod, cf_nutmod, cf_gravmod,
     .        cf_srpmod, cf_eradmod, cf_antradmod, 
     .        cf_ietide, cf_isptide, 
     .        cf_speopmod, cf_etidemod, cf_otidemod, cf_atmtide,
     .        cf_atmlmod, cf_hydrolmod, 
     .        (cf_atmlavg(i),i=1,3), (cf_hydrolavg(i),i=1,3),
     .        cf_dryzen,cf_wetzen, cf_drymap,cf_wetmap, 
     .        cf_ionsrc, cf_magfield, 
     .        cf_antmod_snx, cf_antmod,
     .        (cf_svantbody(i),i=1,cf_nsat),
     .        (cf_svantmod_snx(i),i=1,cf_nsat), 
     .        (cf_svantmod(i),i=1,cf_nsat),
     .        cf_elvcut, cf_nclock, (cf_clock(i),i=1,cf_nclock),
     .        cf_jde, cf_te, cf_jdr, cf_tr,
     .        cf_ut1, cf_xp, cf_yp,  cf_psi,cf_eps,
     .        cf_avlmet, 
     .        cf_nslip, (cf_islip(i),i=1,cf_nslip),
     .        (cf_slpst(i),i=1,cf_nslip),
     .        cf_niextra, (cf_iextra(i),i=1,cf_niextra),
     .        cf_nextra,  (cf_extra(i), i=1,cf_nextra),
     .        cf_ncextra, (cf_cextra(i),i=1,cf_ncextra) 
      else
*         Write 10.71 version with 2 frequency SVS PCO values
*         and cf_antdaz (antenna orientation).  
*         values so the cview will not re-map
          if( cf_gnss.eq.'R' ) then
             do i = 1, cf_nsat
                cf_fL1(i) = fL1_R0
                cf_fl2(i) = fL2_R0
             end do
          end if
* MOD TAH 200205: Added cf_antdaz and size of cf_svantdx       
          write(unit,iostat=ierr) flag, cf_sitnam, cf_rctype,
     .        cf_rcvnum, cf_rcvrsw, cf_swver,
     .        cf_anttyp, cf_antnum, cf_npart, cf_norb,
     .        cf_gnss, cf_nsat, (cf_ischan(i),i=1,cf_nsat),
     .        (cf_fL1(i),cf_fl2(i),i=1,cf_nsat),
     .        cf_ndat, (cf_dattyp(i),i=1,cf_ndat),
     .        ((cf_lambda(i,j),i=1,cf_nsat),j=1,cf_ndat),
     .        cf_skd, cf_nepoch, cf_inter, cf_ircint, cf_mtime,
     .        cf_isessn, cf_iy, cf_im, cf_id, cf_ihr, cf_min, cf_sec,
     .        (cf_offarp(i),i=1,3),
     .        (cf_offsL1(i),i=1,3), (cf_offsL2(i),i=1,3), cf_antdaz,
     .        (cf_svantdx(:,:,j), j=1,cf_nsat),
     .        cf_obfiln, cf_tfiln, cf_jfiln,
     .        cf_frame, cf_precmod, cf_nutmod, cf_gravmod,
     .        cf_srpmod, cf_eradmod, cf_antradmod, 
     .        cf_ietide, cf_isptide, 
     .        cf_speopmod, cf_etidemod, cf_otidemod, cf_atmtide,
     .        cf_atmlmod, cf_hydrolmod, 
     .        (cf_atmlavg(i),i=1,3), (cf_hydrolavg(i),i=1,3),
     .        cf_dryzen,cf_wetzen, cf_drymap,cf_wetmap, 
     .        cf_ionsrc, cf_magfield, 
     .        cf_antmod_snx, cf_antmod,
     .        (cf_svantbody(i),i=1,cf_nsat),
     .        (cf_svantmod_snx(i),i=1,cf_nsat), 
     .        (cf_svantmod(i),i=1,cf_nsat),
     .        cf_elvcut, cf_nclock, (cf_clock(i),i=1,cf_nclock),
     .        cf_jde, cf_te, cf_jdr, cf_tr,
     .        cf_ut1, cf_xp, cf_yp,  cf_psi,cf_eps,
     .        cf_avlmet, 
     .        cf_nslip, (cf_islip(i),i=1,cf_nslip),
     .        (cf_slpst(i),i=1,cf_nslip),
     .        cf_niextra, (cf_iextra(i),i=1,cf_niextra),
     .        cf_nextra,  (cf_extra(i), i=1,cf_nextra),
     .        cf_ncextra, (cf_cextra(i),i=1,cf_ncextra) 
      
      end if
 
 
****  See what happened and check flags and sizes.
      call report_error('IOSTAT',ierr,'write',
     .        'Cfile Type 2 record', 0, 'WRITE_CF2')
 
****  Thats all
      return
      end
 
CTITLE WRITE_CF3
 
      subroutine write_cf3( unit, option, ierr )
                        
      implicit none 
 
*     This routine will write the type 3 record of a C-file.  It is
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
 
*   option  - Option for write:  NONE
 
      character*(*) option
 
* LOCAL VARIABLES
 
*   flag        - First I*4 value in record and indicates the
*           - type of record.
*   i,j     - Loop counters
 
      integer*4 flag, i
 
*   loc_opt - Local version of option so that we can
*           - casefold it.
 
      character*4 loc_opt
 
***** Check the option and see what we should do:
      loc_opt = option
      call casefold(loc_opt)

      flag = 3
 
      write(unit,iostat=ierr) flag,cf_nlabel,cf_nparam,
     .        (cf_islot (i),i=1,cf_nlabel),
     .        (cf_idms  (i),i=1,cf_nlabel),
     .        (cf_rlabel(i),i=1,cf_nlabel),
     .        (cf_preval(i),i=1,cf_nparam)
 
 
****  See what happened and check flags and sizes.
      call report_error('IOSTAT',ierr,'writ',
     .        'Cfile Type 3 record', 0, 'WRITE_CF3')
 
****  Thats all
      return
      end
 
CTITLE WRITE_CF4
 
      subroutine write_cf4( unit, option, ierr )
                        
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
 
*   option  - Option for write:  NONE
 
      character*(*) option
 
* LOCAL VARIABLES
 
*   flag        - First I*4 value in record and indicates the
*           - type of record.
*   i,j     - Loop counters
 
      integer*4 flag,  i
 
*   loc_opt - Local version of option so that we can
*           - casefold it.
 
      character*4 loc_opt
 
***** Check the option and see what we should do:
      loc_opt = option
      call casefold(loc_opt)

      flag = 4
* MOD TAH 100902: Add in test for 1040 version where cf_atmlod added 
      if( cf_nversn.lt. 1040 ) then 
          write(unit,iostat=ierr) flag,
     .           cf_msat, (cf_mprn(i),i=1,cf_msat),
     .           cf_iepoch, cf_iyr, cf_idoy, cf_sod, cf_rclock,
     .           cf_okmet, cf_pres, cf_temp, cf_relhum,
     .           cf_nsave, (cf_save(i),i=1,cf_nsave),
     .           cf_kflag, cf_ksite, cf_klatr, cf_klonr, cf_kradk,
     .           cf_L1Z, cf_L1N, cf_L1E, cf_L2Z, cf_L2N, cf_L2E,
     .           cf_antaz
      elseif( cf_nversn.lt.1060 ) then
* MOD TAH 100902: Add in test for 1040 version where cf_atmlod added 
          write(unit,iostat=ierr) flag,
     .           cf_msat, (cf_mprn(i),i=1,cf_msat),
     .           cf_iepoch, cf_iyr, cf_idoy, cf_sod, cf_rclock,
     .           cf_okmet, cf_pres, cf_temp, cf_relhum, cf_atmlod,
     .           cf_nsave, (cf_save(i),i=1,cf_nsave),
     .           cf_kflag, cf_ksite, cf_klatr, cf_klonr, cf_kradk,
     .           cf_L1Z, cf_L1N, cf_L1E, cf_L2Z, cf_L2N, cf_L2E,
     .           cf_antaz  
      else
* MOD RWK 141215; Write 1060 version adding cf_zendel, removing kflag, 
*                 renaming ksite, klatr, klonr, krad to sitecd. latr_sph, 
*                 lonr, radius, and reordering the last 2 lines.
          write(unit,iostat=ierr) flag,   
     .           cf_msat, (cf_mprn(i),i=1,cf_msat),
     .           cf_iepoch, cf_iyr, cf_idoy, cf_sod, cf_rclock,
     .           cf_zendel, cf_okmet, cf_pres, cf_temp, cf_relhum,
     .           cf_atmlod, cf_sitecd, cf_latr_sph, cf_lonr, cf_radius,  
     .           cf_L1Z, cf_L1N, cf_L1E, cf_L2Z, cf_L2N, cf_L2E,
     .           cf_antaz, cf_nsave, (cf_save(i),i=1,cf_nsave)  

      end if
 
 
****  See what happened and check flags and sizes.
      call report_error('IOSTAT',ierr,'write',
     .        'Cfile Type 4 record', 0, 'WRITE_CF4')
 
****  Thats all
      return
      end
 
 
CTITLE WRITE_CF5
 
      subroutine write_cf5( unit, option, ierr )
                        
      implicit none 
 
*     This routine will write the type 5 record of a C-file.  It is
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
 
*   option  - Option for write: 
*             'SH' Mean write no partials 
*             'AL' Means write all the data (defaults)
 
      character*(*) option
 
* LOCAL VARIABLES
 
*   flag        - First I*4 value in record and indicates the
*           - type of record.
*   i,j     - Loop counters
 
      integer*4 flag, i
 
*   loc_opt - Local version of option so that we can
*           - casefold it.
 
      character*4 loc_opt
 
***** Check the option and see what we should do:
      loc_opt = option
      call casefold(loc_opt)
 
      flag = 5

      if( loc_opt(1:2).eq.'SH' ) then
         cf_nparts = 0
      end if 

      if( cf_nversn.lt.980 ) then  
          write(unit,iostat=ierr) flag, cf_iprn,cf_elev, cf_azimuth,
     .            cf_okwvr, cf_wvrdel, cf_atmdel, cf_tau, cf_drate,
     .            cf_ierfl, cf_data_flag, cf_ndats,
     .            (cf_obsv(i), i=1, cf_ndats),
     .            (cf_omcs(i), i=1, cf_ndats),
     .            (cf_isnr(i), i=1, cf_ndats), cf_ampl1, cf_ampl2,
     .            (cf_obswgt(i), i=1, cf_ndats),
     .            cf_nspare,( cf_spare(i), i=1, cf_nspare),
     .            cf_nparts,( cf_tmpart(i), i=1, cf_nparts)
      elseif( cf_nversn.lt. 1060 ) then
          write(unit,iostat=ierr) flag, cf_iprn,cf_elev, cf_azimuth,
     .            cf_okwvr, cf_wvrdel, cf_atmdel,
     .            cf_svcepc, cf_svcL1, cf_tau, cf_drate,
     .            cf_ierfl, cf_data_flag, cf_ndats,
     .            (cf_obsv(i), i=1, cf_ndats),
     .            (cf_omcs(i), i=1, cf_ndats),
     .            (cf_isnr(i), i=1, cf_ndats), cf_ampl1, cf_ampl2,
     .            (cf_obswgt(i), i=1, cf_ndats),
     .            cf_nspare,( cf_spare(i), i=1, cf_nspare),
     .            cf_nparts,( cf_tmpart(i), i=1, cf_nparts) 
      else
* MOD RWK 141215: Read 1060 version,adding cf_nadang, and removing cf_okwvr, cf_wvrdel, cf_obswgt
            write(unit,iostat=ierr) flag, cf_iprn,cf_elev, cf_azimuth,
     .            cf_nadang, cf_atmdel, 
     .            cf_svcepc, cf_svcL1, cf_tau, cf_drate,
     .            cf_ierfl, cf_data_flag, cf_ndats,
     .            (cf_obsv(i), i=1, cf_ndats),
     .            (cf_omcs(i), i=1, cf_ndats),
     .            (cf_isnr(i), i=1, cf_ndats), cf_ampl1, cf_ampl2,
     .            cf_nspare,( cf_spare(i), i=1, cf_nspare),
     .            cf_nparts,( cf_tmpart(i), i=1, cf_nparts)
      end if

 
      call report_error('IOSTAT',ierr,'write',
     .        'Cfile Type 5 record', 0, 'WRITE_CF5')
 
****  Thats all
      return
      end
 
 
