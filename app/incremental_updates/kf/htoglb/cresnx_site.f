CTITLE CRESNX_SITE
 
      subroutine cresnx_site(unit, unitc )
 
      implicit none 

*     Routine to creode the site blocks from SINEX:
 
      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
 
* PASSED VARIABLES
 
*   unit        - Unit number
*   unitc       - Commnts file unit number
 
      integer*4 unit, unitc
   

****  Start creoding the types of site blocks
*     Generate the unique site occ and pt values
      call qname_to_codeoccpt
      call cre_site_id ( unit, unitc )
      call cre_site_rcv( unit, unitc )
      call cre_site_ant( unit, unitc )
      call cre_site_phs( unit, unitc )
      call cre_site_ecc( unit, unitc )
 
****  Thats all
      return
      end
 
CTITLE CRE_SITE_ID
 
      subroutine cre_site_id(unit, unitc  )

      implicit none 
 
*     Routine to create the site ID block
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
 
* PASSED VARIABLES
 
*   unit        - Unit number of sinex file
*   unitc       - Commnts file unit number
 
      integer*4 unit, unitc
 
 
* LOCAL VARIABLES
 
*   ns      - number of stations found
*   lat_deg, lat_min, lng_deg, lng_min  -- Lat and long degrees and minutes
*   sign_lat  - Sign of the latitude.
 
      integer*4 ns, lat_deg, lat_min, lng_deg, lng_min, sign_lat, occ

*   lat_sec, lng_sec, loc_coord(3), rot_mat(3,3) -- Lat and long seconds,
*       and local coordinates in radians.

      real*8 lat, lng, lat_sec, lng_sec, loc_coord(3), RM(3,3)

 
*   domes       - DOMES number for site (9 characters needed)
*   type        - P=GPS, R=VLBI, S=SLR
*   pt          - Point number.  This is added to the end of the four 
*                 character code.  If it is not A we will make _GP[Pt]

      character*12 domes
      character*4 type, pt, code

****  Start creating the records.

      write(unit,'(a)') '+SITE/ID'
      call cp_comments(unit, unitc, 'SITE/ID')

      do ns = 1, cnum_sites

* MOD TAH 1906127: Only output entry if site actually used.
         if( gtol_sites(ns).ge.0 ) then 
 
*           Slip out the point and domes number from the site name
            call name_to_id(qsite_names, ns, code, pt, occ)
            domes = qfull_names(ns)(23:31)
            call sub_null(domes)
            type  = qfull_names(ns)(32:32)

*           Check to see if whole entry is null.  This can
*           happen when there are not site records in SINEX
*           files.
            call sub_null(type)
            call sub_null(qfull_names(ns))
            qfull_names(ns)(32:32) = type

****        Now get the approximate lat and long
            if( (site_pos(1,ns)**2+site_pos(2,ns)**2+ 
     .           site_pos(3,ns)**2).gt.6000000.d0**2  ) then
                 call xyz_to_geod( RM, site_pos(1,ns), loc_coord)
            else 
                 loc_coord(1) = 0
                 loc_coord(2) = 0
                 loc_coord(3) = 0
            end if
                
            lat = (pi/2 - loc_coord(1))*180.d0/pi
            if( lat.lt.0 ) then
                sign_lat = -1
                lat = abs(lat)
            else
                sign_lat = +1
            end if
            lat_deg = int(lat)
            lat_min = (lat - lat_deg)*60.d0
            lat_sec = (lat - lat_deg - lat_min/60.d0)*3600.d0
            if( lat_deg.gt.0 ) then
                lat_deg = sign_lat*lat_deg
            else if ( lat_min.gt.0 ) then
                lat_min = sign_lat*lat_min
            else
                lat_sec = sign_lat*lat_sec
            end if

            lng = loc_coord(2)*180.d0/pi
            lng_deg = int(lng)
            lng_min = (lng - lng_deg)*60.d0
            lng_sec = (lng - lng_deg - lng_min/60.d0)*3600.d0

            write(unit, 120) qsite_names(ns), pt, domes, type, 
     .                       qfull_names(ns), lng_deg, lng_min, lng_sec,
     .                       lat_deg, lat_min, lat_sec, loc_coord(3)
 120        format(1x,a4,2x,a1,1x,a9,1x,a1,1x,a22, 1x,i3.3,1x,i2.2,1x,
     .             F4.1,1x,i3,i3.2,F5.1,1x,f7.1)
         endif
      end do
      write(unit,'(a)') '-SITE/ID'
****  Thats all
      return
      end
 
CTITLE CRE_SITE_RCV
 
      subroutine cre_site_rcv(unit, unitc  )

      implicit none 
 
*     Routine to read the site RECEIVER block
 
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
*   ns      - number of stations found
*   os      - Original station number of Occ changes
*   refs_yr, refs_doy, refs_sec - Starting yr, day-of-year and
*             seconds
*   refe_yr, refe_doy, refe_sec - Ending yr, day-of-year and
*             seconds
*   occ     - Occupation number 
*   indx    - Position in string
 
      integer*4 ns,  refs_yr, refs_doy, refs_sec ,
     .          refe_yr, refe_doy, refe_sec, occ
 
*   code        - 4 character code from which we generate name
*   type        - P=GPS, R=VLBI, S=SLR
*   pt          - Point number.  This is added to the end of the four character
*                 code.  Later if it is not A we will make _GP[Pt]

      character*4 code, type, pt

*   recv_ty     - Reveiver type
*   recv_sn     - Receiver serial number
*   recv_fw     - Receiver firmware
*   recv_ty_full - Full name of receiver type written to file.
*   recv_fw_full - Full firrm ware written to file.
* MOD TAH 101016: Made recv_ty full length 

      character*16 recv_fw
      character*20 recv_ty, recv_ty_full
      character*11 recv_fw_full 
      character*8  recv_sn

****  Start creoding
      write(unit,'(a)') '+SITE/RECEIVER'
      call cp_comments(unit, unitc, 'SITE/RECEIVER')

      do ns = 1, cnum_sites

* MOD TAH 1906127: Only output entry if site actually used.
         if( gtol_sites(ns).ge.0 ) then 

****        Start creating the output information
            call name_to_id( qsite_names, ns, code, pt, occ)
            type  = qfull_names(ns)(32:32)

****        OK Save the information about the receiver
            call jd_to_yds( qrecv_st(ns), refs_yr,refs_doy,refs_sec )
            call jd_to_yds( qrecv_en(ns), refe_yr,refe_doy,refe_sec )
            recv_ty = qrecv_ty(ns)
            if( recv_ty(1:4).eq.'----' .and. 
     .          qsite_names(ns+1)(1:4).eq.qsite_names(ns)(1:4) ) 
     .          recv_ty = qrecv_ty(ns+1)
            recv_sn = qrecv_sn(ns)
            if( recv_sn(1:4).eq.'----' .and. 
     .          qsite_names(ns+1)(1:4).eq.qsite_names(ns)(1:4) ) 
     .          recv_sn = qrecv_sn(ns+1)
            recv_fw = qrecv_fw(ns)
            if( recv_fw(1:4).eq.'----' .and. 
     .          qsite_names(ns+1)(1:4).eq.qsite_names(ns)(1:4) ) 
     .          recv_fw = qrecv_fw(ns+1)
            call sub_null(recv_ty)
            call sub_null(recv_sn)
            call sub_null(recv_fw(1:11))
            if( ichar(recv_fw(12:12)).ne.0 .and.  
     .          recv_fw(12:12).ne.'-'             ) then 
                recv_ty_full = recv_ty // recv_fw(12:15)
            else
                recv_ty_full = recv_ty 
            end if 
            recv_fw_full = recv_fw(1:11)
            write(unit, 120) code, pt, occ, type, mod(refs_yr,100),
     .                       refs_doy, refs_sec, mod(refe_yr,100),
     .                       refe_doy, refe_sec, recv_ty_full,
     .                       recv_sn,  recv_fw_full 
 120        format(1x,a4,2x,a1,2x,i3,1x,a1,1x,i2.2,':',i3.3,':',i5.5,
     .                  1x,i2.2,':',i3.3,':',i5.5,1x,a20,1x,a5,1x,a11)
         endif
      end do

      write(unit,'(a)') '-SITE/RECEIVER'

****  Thats all
      return
      end
      
CTITLE CRE_SITE_ANT
 
      subroutine cre_site_ant(unit, unitc )
 
      implicit none 

*     Routine to create the site ANTENNA block
 
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
*   ns      - number of stations found
*   refs_yr, refs_doy, refs_sec - Starting yr, day-of-year and
*             seconds
*   refe_yr, refe_doy, refe_sec - Ending yr, day-of-year and
*             seconds
*   occ     - Occupation number 
*   indx    - Position in string
 
      integer*4 ns,  refs_yr, refs_doy, refs_sec ,
     .          refe_yr, refe_doy, refe_sec, occ 
 
*   code        - 4 character code from which we generate name
*   type        - P=GPS, R=VLBI, S=SLR
*   pt          - Point number.  This is added to the end of the four character
*                 code.  Later if it is not A we will make _GP[Pt]

      character*4 code, type, pt

*   ante_ty     - Reveiver type
*   ante_sn     - Receiver serial number
*   ante_fw     - Receiver firmware

      character*16 ante_ty
      character*20 ante_ty_full 
      character*8  ante_sn
      character*5  radome

****  Start creoding
      write(unit,'(a)') '+SITE/ANTENNA'
      call cp_comments(unit, unitc, 'SITE/ANTENNA')

      do ns = 1, cnum_sites
 
* MOD TAH 1906127: Only output entry if site actually used.
         if( gtol_sites(ns).ge.0 ) then 
 
****        Create the line
            call name_to_id( qsite_names, ns, code, pt, occ)
            type  = qfull_names(ns)(32:32)

****        OK Save the information about the receiver
            call jd_to_yds( qante_st(ns), refs_yr,refs_doy,refs_sec )
            call jd_to_yds( qante_en(ns), refe_yr,refe_doy,refe_sec )

            ante_ty = qante_ty(ns) 
            if( ante_ty(1:4).eq.'----' .and. 
     .          qsite_names(ns+1)(1:4).eq.qsite_names(ns)(1:4) ) 
     .          ante_ty = qante_ty(ns+1)
            ante_sn = qante_sn(ns)
            if( ante_sn(1:4).eq.'----' .and. 
     .          qsite_names(ns+1)(1:4).eq.qsite_names(ns)(1:4) ) 
     .          ante_sn = qante_sn(ns+1)
            call sub_null(ante_ty)
            call sub_null(ante_sn(1:5))

            if( ichar(ante_sn(6:6)).ne.0 .and. 
     .          ante_sn(6:6).ne.'-' )   then
                ante_ty_full = ante_ty // ante_sn(6:8) 
            else 
                ante_ty_full = ante_ty 
            end if 
* MOD TAH 041229: Move characters at end of antenna name accross one
*           spot
            radome = ante_ty_full(16:19)
* MOD TAH 200805: Fixed issue with radome being truncated with determined
*           from extended part of serial number. (Normally saved in
*           qradome_ty but does not happen if original SINEX has spaces
*           for radome.
            if( radome(2:4).eq.'NON' ) radome(2:) = 'NONE'
            if( qradome_ty(ns).ne.'----' ) then
                ante_ty_full = ante_ty // qradome_ty(ns)
            else
                ante_ty_full(16:16) = ' '
* MOD TAH 200805: Fixed 2:5 of radome (char 1 is +/-/ )
                ante_ty_full(17:20) =  radome(2:5)
            endif
* MOD TAH 051027: Clear the +- sign that says that the radome type was
*           found.
            if( ante_ty_full(16:16).eq. '+' .or.
     .          ante_ty_full(16:16).eq. '-'  ) ante_ty_full(16:16) = ' '
            
            write(unit, 120) code, pt, occ, type, mod(refs_yr,100),
     .                           refs_doy, refs_sec, mod(refe_yr,100),
     .                           refe_doy, refe_sec, ante_ty_full, 
     .                           ante_sn,int(qantdaz(ns))
 120        format(1x,a4,2x,a1,2x,i3,1x,a1,1x,i2.2,':',i3.3,':',i5.5,
     .                  1x,i2.2,':',i3.3,':',i5.5,1x,a20,1x,a5,1x,I4) 
         endif
      end do
!+SITE/ANTENNA
!*CODE PT SOLN T _DATA START_ __DATA_END__ ____ANTENNA_TYPE____ _S/N_ _DAZ
! MKEA  A ---- P 96:221:00000 15:023:00000 AOAD/M_T        NONE 482      0
! MKEA  A ---- P 15:023:00000 19:351:79200 JAVRINGANT_DM   SCIS 00983   90
!-SITE/ANTENNA

      write(unit,'(a)') '-SITE/ANTENNA'
****  Thats all
      return
      end

CTITLE CRE_SITE_PHS
 
      subroutine cre_site_phs(unit, unitc  )
 
      implicit none 

*     Routine to read the site PHASE CENTER MODEL block
 
      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
 
* PASSED VARIABLES
 
*   unit        - Unit number of sinex file
*   unitc       - Commnts file unit number
 
      integer*4 unit, unitc
 
 
* LOCAL VARIABLES
 
*   i, j        - Loop counters
 
      integer*4  i, j
      
*   L1a_ecc, L2a_ecc  -- L1 and L2 phase center offsets (UNE form)

      real*4 L1a_ecc(3), L2a_ecc(3)
 
*   not_same   - Set true when a different phase center model is found
 
      logical not_same
 
*   ante_ty     - Reveiver type
*   ante_sn     - Receiver serial number
*   ant_mod     - Type of annetta phase model used

      character*16 ante_ty
      character*20 ante_ty_full 
      character*8  ante_sn
      character*10  ant_mod

****  Start creating phase center records.  Only echo values when they
*     differ.
      write(unit,'(a)') '+SITE/GPS_PHASE_CENTER'
      call cp_comments(unit, unitc, 'SITE/GPS_PHASE_CENTER')

      do i = 1, cnum_sites

* MOD TAH 1906127: Only output entry if site actually used.
         if( gtol_sites(i).ge.0 ) then 
 
            not_same = .true.
            do j = 1, i-1
* MOD TAH 200108: Added check that j site was also used.
               if( qante_ty(j).eq.qante_ty(i) .and.
     .             qante_sn(j).eq.qante_sn(i) .and.
     .             qant_mod(j).eq.qant_mod(i) .and.
     .             qradome_ty(j).eq.qradome_ty(i)  .and.
     .             qL1a_ecc(1,j).eq.qL1a_ecc(1,i)  .and.
     .             qL1a_ecc(2,j).eq.qL1a_ecc(2,i)  .and.
     .             qL1a_ecc(3,j).eq.qL1a_ecc(3,i)  .and.
     .             qL2a_ecc(1,j).eq.qL2a_ecc(1,i)  .and.
     .             qL2a_ecc(2,j).eq.qL2a_ecc(2,i)  .and.
     .             qL2a_ecc(3,j).eq.qL2a_ecc(3,i)  .and.
     .             gtol_sites(j).ge.0  ) then
                   not_same = .false.
                end if
            end do

            if( not_same ) then
                ante_ty = qante_ty(i)
                ante_sn = qante_sn(i)
*               MOD TAH 960926: convert code to full model name
*               MOD TAH 101015: See if only 4-Char code
                if( qant_mod(i)(5:5).eq.' ' ) then  ! Old code            
                    call phsmod_to_code( ant_mod,qant_mod(i),'DECODE') 
                else
                    ant_mod = qant_mod(i)
                end if
               
                call sub_null( ante_ty )
                call sub_null( ante_sn )
                call sub_null( ant_mod )

                if( ichar(ante_sn(6:6)).ne.0 .and. 
     .              ante_sn(6:6).ne.'-' )   then
                    ante_ty_full = ante_ty // ante_sn(6:8) 
                else 
                    ante_ty_full = ante_ty 
                end if
                if( qradome_ty(i).ne.'----' ) then
                     ante_ty_full = ante_ty // qradome_ty(i)
                endif
*               See if the radome type was found; denoted by +- sign
*               in character 1
C               print *,'I ',qsite_names(i),'|',ante_ty_full,'|',
C    .                  ante_sn,'|',qante_ty(i),'|',qradome_ty(i),'|',
C    .                      qante_sn(i),'|'
                if( ante_ty_full(16:16).eq.'-' ) then
                    ante_ty_full(17:) = 'NONE'
                    ante_ty_full(16:16) = ' '
                elseif ( ante_ty_full(16:16).eq.'+' ) then
                    ante_ty_full(16:16) = ' '
                endif

                L1a_ecc(2) = qL1a_ecc(1,i) 
                L1a_ecc(3) = qL1a_ecc(2,i) 
                L1a_ecc(1) = qL1a_ecc(3,i) 
                L2a_ecc(2) = qL2a_ecc(1,i) 
                L2a_ecc(3) = qL2a_ecc(2,i) 
                L2a_ecc(1) = qL2a_ecc(3,i) 

****            Creode the line
                write(unit,120) ante_ty_full, ante_sn, L1a_ecc, L2a_ecc,
     .                          ant_mod
 120            format(1x,a20,1x,a5,6(1x,f6.4),1x,a10) 
            end if
         end if
      end do

      write(unit,'(a)') '-SITE/GPS_PHASE_CENTER'

****  Thats all
      return
      end

CTITLE CRE_SITE_ECC

      subroutine cre_site_ecc(unit, unitc )
 
      implicit none 

*     Routine to create the site ECCENTRICITY block
 
      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
 
* PASSED VARIABLES
 
*   unit        - Unit number of sinex file
*   unitc       - Commnts file unit number
 
      integer*4 unit, unitc
 
* LOCAL VARIABLES

*   arp_ecc(3)  - Antenna monument to ARP eccentricity 

       real*4 arp_ecc(3) 
 
*   ierr        - IOSTAT error
*   ns      - number of stations found
*   os      - Original station number of Occ changes
*   refs_yr, refs_doy, refs_sec - Starting yr, day-of-year and
*             seconds
*   refe_yr, refe_doy, refe_sec - Ending yr, day-of-year and
*             seconds
*   occ     - Occupation number 
*   indx    - Position in string
 
      integer*4 ns, refs_yr, refs_doy, refs_sec ,
     .          refe_yr, refe_doy, refe_sec, occ
 
 
*   code        - 4 character code from which we generate name
*   type        - P=GPS, R=VLBI, S=SLR
*   pt          - Point number.  This is added to the end of the four character
*                 code.  Later if it is not A we will make _GP[Pt]
*   ecc_type    - Method for specifying the ECC (onlu UNE accepted at
*                 the moment)

      character*4 code, type, pt, ecc_type


****  Start creoding
      write(unit,'(a)') '+SITE/ECCENTRICITY'
      call cp_comments(unit, unitc, 'SITE/ECCENTRICITY')

      do ns = 1, cnum_sites

* MOD TAH 1906127: Only output entry if site actually used.
         if( gtol_sites(ns).ge.0 ) then 
 
            call name_to_id( qsite_names, ns, code, pt, occ)
            type  = qfull_names(ns)(32:32)

****        OK Save the information about the receiver.  Assume that
*           these are the same as the antenna
            call jd_to_yds( qante_st(ns), refs_yr,refs_doy,refs_sec )
            call jd_to_yds( qante_en(ns), refe_yr,refe_doy,refe_sec )

            arp_ecc(2) = qarp_ecc(1,ns)
            arp_ecc(3) = qarp_ecc(2,ns)
            arp_ecc(1) = qarp_ecc(3,ns)
            ecc_type   = 'UNE'
            if( qante_ty(ns)(1:4).eq.'----' .and. 
     .          qsite_names(ns+1)(1:4).eq.qsite_names(ns)(1:4) ) then
               arp_ecc(2) = qarp_ecc(1,ns+1)
               arp_ecc(3) = qarp_ecc(2,ns+1)
               arp_ecc(1) = qarp_ecc(3,ns+1)
            endif


            write(unit, 120) code, pt, occ, type, mod(refs_yr,100),
     .                       refs_doy, refs_sec, mod(refe_yr,100),
     .                       refe_doy, refe_sec, ecc_type,
     .                       arp_ecc
 120        format(1x,a4,2x,a1,2x,i3,1x,a1,1x,i2.2,':',i3.3,':',i5.5,
     .                  1x,i2.2,':',i3.3,':',i5.5,1x,a3,3f9.4)
         endif 
      end do
      write(unit,'(a)') '-SITE/ECCENTRICITY'

****  Thats all
      return
      end

CTITLE CHKSNX_SITE
 
      subroutine chksnx_site(np, cov_parm, sol_parm )
 
      implicit none 

*     Routine to check the sigmas on the site positions and rename the
*     site to ...._X.. if the sigma is too large (>sqrt(20)=4.5 meters)
 
      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
 
* PASSED VARIABLES
 
*   np      - number of parameters expected
 
      integer*4  np
 
* cov_parm(np,np) - Covarince matrix
* sol_parm(np)    - Solution vector
 
 
      real*8 cov_parm(np,np), sol_parm(np)
 
 
* LOCAL VARIABLES
      integer*4 i, j, k, ni, no
* MOD TAH 210113: Updated to allow removal of low accuracy site
      logical reset   ! Set true if we are eliminating a site



****  Loop over the sites and check the size of the sites
      reset = .false.
      do i = 1, cnum_sites

* MOD TAH 1906127: Only output entry if site actually used.
         if( gtol_sites(i).ge.0 ) then 
            do j = 1,3 
               ni = qparn_sites(j,i)
               if( ni.gt.0 ) then
* MOD TAH 210112: Added Var_Limit for value of large sigma
*                 (use -V=<Var> to set)
                   if( cov_parm(ni,ni).gt.Var_Limit ) then
                       if( qsite_names(i)(5:6) .ne. '_X' ) then
                           write(*,120) i,qsite_names(i), j,ni, 
     .                                  sqrt(cov_parm(ni,ni)),
     .                                  sqrt(Var_Limit)
 120                       format('HighSigma Site ',i4,1x,a8,1x,
     .                            'Comp ',i1,' NP ',i4,1x,' Sigma ',
     .                             F6.2,' m; Limit ',F6.2,' m')
                           qsite_names(i)(5:6) = '_X'
* MOD TAH 210113: Set to show reset of parameter mapping needed.
                           reset = .true.
                       endif
                   endif
               endif
            enddo
* MOD TAH 210113: If the site name has been marked with _X remove
*           from gtol_sites (set value to -1) 
            if( qsite_names(i)(5:6).eq.'_X' ) then
               gtol_sites(i) = -1
               qparn_sites(:,i) = 0  ! Reset all pointers to 0
            endif
         endif
      enddo

****  See if we need to re-map parameters (code from SR_CODES)
      if( .not. reset ) RETURN

*     Originally created in SR_CODES but if we eliminate sites we
*     need to redo here
      write(*,220) qnum_parn
 220  format('RESETING Parameter pointers: qnum_parn ',I8)
****  Now get the mapping of the output parameters to the input ones.
      no = 0
      do i = 1, cnum_sites
         do j = 1, 3
            if( qparn_sites(j,i).ne.0 ) then
                no = no + 1
                otoi(no) = qparn_sites(j,i)
                itoo(qparn_sites(j,i)) = no
                atos(qparn_sites(j,i)) = i 
            end if
         end do
         do j = 1, 3
            if( qparn_vel(j,i).ne.0 ) then
                no = no + 1
                otoi(no) = qparn_vel(j,i)
                itoo(qparn_vel(j,i)) = no
                atos(qparn_vel(j,i)) = i
            end if
         end do
      end do

* MOD TAH 0506022: Check satelliet offsets
      do i = 1, cnum_svs
         do j = 1,3
            if( qparn_svs(max_svs_elem-3+j,i).ne.0 ) then
                no = no + 1
                otoi(no) = qparn_svs(max_svs_elem-3+j,i)
                itoo(qparn_svs(max_svs_elem-3+j,i)) = no
                atos(qparn_svs(max_svs_elem-3+j,i)) = i
            end if
          end do
      end do   
                
****  Do translation and rate
      do j = 1,2 
         do i = 1,3
            if( qparn_tran(i,j).ne.0 ) then
                no = no + 1
                otoi(no) = qparn_tran(i,j)
                itoo(qparn_tran(i,j)) = no
            end if
         end do
      end do

*     Now do the scale and rate of change    
      do i = 1,2
         if( qparn_scale(i).ne.0 ) then
             no = no + 1
             otoi(no) = qparn_scale(i)
             itoo(qparn_scale(i)) = no
         end if
      end do

****  Now do the PMU parameters
      do i = 1, 3
         do j = 1, 2
            if( qparn_pmu(j,i).ne.0 ) then
                no = no + 1
                otoi(no) = qparn_pmu(j,i)
                itoo(qparn_pmu(j,i)) = no
            end if
         end do
      end do
      
****  Now do the multi-PMU parameters
      do i = 1, 3
         do j = 1, 2
            do k = 1, qnum_mul_pmu(j,i) 
               if( qparn_mul_pmu(j,i,k).ne.0 ) then
                   no = no + 1
                   otoi(no) = qparn_mul_pmu(j,i,k)
                   itoo(qparn_mul_pmu(j,i,k)) = no
               end if
            end do
         end do
      end do

***** Save the number of parameters to be output to the SINEX file
      qnum_parn = no

      write(*,240) qnum_parn
 240  format('DONE resetting pointers    : qnum_parn ',I8)

***** Thats all
      return
      end

                     

            

