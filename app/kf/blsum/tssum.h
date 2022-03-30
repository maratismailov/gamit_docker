*     DECLARATIONS of tssum
      integer*4 max_ent   ! Maximum number of total entries
      integer*4 max_site  ! Maximum number of sites
      integer*4 max_ts    ! Maximum number of entries per site

      character*(*) ts_ver   ! Version of time series files

C     parameter ( max_ent  = 10000000 )
C MOD TAH 180306 Increased to support UNR stations.
C MOD TAH 180615: Put back to 10000000 due to strange behavior with larger value
C MOD TAH 180615: Seems OK with larger value.  Problem may be with build.
      parameter ( max_ent  = 40000000 )   ! Increased for UNR files (7Gb to run).
C MOD TAH 210203: Increased max_site to 9000 from 8192 to handle GAGE number of stations.
      parameter ( max_site =    9000 )
      parameter ( max_ts   =   20000 )
*
* MOD TAH 051129: Added version 1.0.1: Version number in file and type 
*     designation (rapid/final)
* MOD TAH 121214: Version 1.1.0: Added headers and unified header files.
      parameter ( ts_ver   = '1.1.0' )

* Declaratons
      integer*4 num_ent, num_site, num_code, num_ts

* Data records from .org files
      integer*4 in_ns(max_ent), in_cs(max_ent)
      integer*4 ln_tsdir, ln_prod_id

      integer*4 date_rel(5)

      integer*4 ts_edt(3,max_ts)  ! Edit status: 
                                  !    0 -- OK; 
                                  ! 1  1 -- max_sigma 
                                  ! 2  2 -- rename edit
                                  ! 3  4 -- outlier
                                  ! 4  8 -- Outside time range
                                  ! 5 16 -- Data on day of earthquake

      real*8    sec_rel

      real*8 in_mjd(max_ent), in_xyz(3,max_ent), in_xyz_std(6,max_ent),
     .       in_llu(3,max_ent), 
     .       in_neu(3,max_ent),in_neu_std(6,max_ent) 

      real*8 ts_mjd(max_ts), ts_xyz(3,max_ts), ts_xyz_std(6,max_ts),
     .       ts_llu(3,max_ts), 
     .       ts_neu(3,max_ts),ts_neu_std(6,max_ts), 
     .       ts_neu_res(max_ts,3), ts_neu_sig(max_ts,3),
     .       ts_neu_ures(max_ts,3)

      real*8 sv_mjd, sv_xyz(3), sv_xyz_std(6),
     .       sv_llu(3), 
     .       sv_neu(3),sv_neu_std(6) 

      real*8 ts_first, ts_last

      real*8 ref_xyz(3), ref_llu(3), ref_neu(3)
      real*8 save_ref_xyz(3,max_site) ! Saved values of the reference coordinates

      character*8  in_site(max_site)
      character*4  in_code(max_site), ts_code
      character*128 tsdir, prod_id
      character*5  ts_type(max_ts)     ! Based on prod_id, rapid/final and any new
                                ! type passed in the prod_id.  Normally updated from
                                ! the value read in.
     .,            sv_type
     .,            ts_ref_type  ! New reference type based on prod_id
          
      character*5 in_type(max_ent)    ! Solution type: Mostly set to the
                 ! ts_ref_type but when .pos are read, these values are
                 ! perserved.

     
      character*256 in_file, ts_file
      character*16 in_full(max_site), ts_full

      character*16 ts_ver_read
      character*16 reference_frame  ! Reference frame read for org files.

      character*16 jn_full

      character*32  solnstr  ! Name of solution and station

      character*8  tsprog    ! Name of program (tssum/tscon/tsfit)

      logical rf_out  ! Set true when reference frame line written

! ts_recomp -- Force re-computation of the reference position (tsssum -C)
! ts_keep   -- Keeps bad entries (old behavior) but now removes them

      logical ts_refresh, ts_recomp, ts_keep, new_ts


      common / ts_i4 / num_ent, num_site, num_code, num_ts,
     .      in_ns, in_cs, date_rel, ln_tsdir, ln_prod_id,
     .      ts_refresh, ts_recomp, new_ts, rf_out, ts_edt,
     .      ts_keep

      common / ts_r8 /sec_rel,  in_mjd, in_xyz, in_xyz_std,
     .       in_llu, in_neu,in_neu_std, ts_mjd, ts_xyz, ts_xyz_std,
     .       ts_llu, ts_neu,ts_neu_std, ts_first, ts_last, 
     .       ref_xyz, ref_llu, ref_neu, ts_neu_res, ts_neu_sig,
     .       ts_neu_ures, save_ref_xyz,
     .       sv_mjd, sv_xyz, sv_xyz_std, sv_llu, sv_neu,sv_neu_std 

      common / ts_ch /in_site,in_code, ts_code, tsdir, 
     .       prod_id, in_file, ts_file, in_full, ts_full, ts_ver_read,
     .       ts_type, sv_type, ts_ref_type, reference_frame,
     .       jn_full, solnstr, in_type, tsprog
