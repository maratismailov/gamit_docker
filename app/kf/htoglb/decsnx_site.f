CTITLE DECSNX_SITE
 
      subroutine decsnx_site(unit, line )
 
      implicit none 

*     Routine to decode the site blocks from SINEX:
 
      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
 
* PASSED VARIABLES
 
*   unit        - Unit number
 
      integer*4 unit
 
*   line        - Line read from input
 
      character*(*) line
 
* LOCAL VARIABLES
 
*   indx        - Pointer in string
*   trimlen     - Length of string
 
      integer*4 indx, trimlen
 
*   block_found - True is block found
 
      logical block_found

****  Start decoding the types of site blocks
      block_found = .false.
 
      indx = index(line,'SITE/ID')
      if( indx.gt.0 ) then
          block_found = .true.
          call dec_site_id( unit, line  )
      end if
 
      indx = index(line,'SITE/RECEIVER')
      if( indx.gt.0 ) then
          block_found = .true.
          call dec_site_rcv( unit, line  )
      end if
      indx = index(line,'SITE/ANTENNA')
      if( indx.gt.0 ) then
          block_found = .true.
          call dec_site_ant( unit, line  )
      end if
 
      indx = index(line,'SITE/GPS_PHASE_CENTER')
      if( indx.gt.0 ) then
          block_found = .true.
          call dec_site_phs( unit, line )
      end if
 
      indx = index(line,'SITE/ECCENTRICITY')
      if( indx.gt.0 ) then
          block_found = .true.
          call dec_site_ecc( unit, line )
      end if
 
 
****  Thats all the blocks so far
      if( .not.block_found ) then
          write(*,500) line(1:trimlen(line))
500       format('DECSNX_SITE: Unknown block type',/a)
          call snx_finbl(unit)
      end if
 
****  Thats all
      return
      end
 
CTITLE DEC_SITE_ID
 
      subroutine dec_site_id(unit, line )
 
      implicit none 

*     Routine to read the site ID block
 
      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
 
* PASSED VARIABLES
 
*   unit        - Unit number of sinex file
 
      integer*4 unit
 
*   line        - Line read from file
 
      character*(*) line
 
* LOCAL VARIABLES
 
*   ierr        - IOSTAT error
*   ns      - number of stations found
 
      integer*4 ierr
 
*   occ_found   - Occ found in case of duplicate numbers
*   end_block   - Set true at end of block
 
      logical end_block, occ_found
 
*   domes       - DOMES number for site (9 characters needed)
*   type        - P=GPS, R=VLBI, S=SLR
*   pt          - Point number.  This is added to the end of the four character
*                 code.  Later if it is not A we will make _GP[Pt]
*   test_name   - Test name to see if we already name.  If so then program
*                 assumes that occupation number should be incremented.
*   test_occ    - Test occ number
*   os, indx    - SIte number (if gt 0 then we already have so try the
*                 next occupation number; indx - position in string
      integer*4 test_occ, os, indx, i

      character*12 domes
      character*4 type, pt
      character*8 test_name, test_read
      character*32 test_full

****  Start decoding
 
      end_block = .false.
      ierr = 0
      do while ( ierr.eq.0 .and. .not.end_block )
 
          read(unit,'(a)', iostat=ierr) line
          if( ierr.eq.0 .and. line(1:1).eq.' ' ) then
 
****          Decode the line
              read(line, 120) test_read, pt, domes, type, 
     .                        test_full
 120          format(1x,a4,1x,a2,1x,a9,1x,a1,1x,a22)

* MOD TAH 970716: Check for blank pt names, replace with ' A'
              if( pt(2:2).eq.' ' .or. pt(2:2).eq.'-' ) pt = ' A'

****          Generate test name
              occ_found = .false.
              test_occ  = 1
              
              write(test_name,130) test_read, test_occ,
     .              pt(2:2)
 130         format(a4,i3.3,a1)

****         See if we have this name
             call casefold(test_name)
             indx = 1
             call get_cmd(test_name, qsite_names, qnum_sites,
     .                         os, indx ) 

             call casefold(type)
             if( os.gt.0 ) then
                 call trimlead(test_full)
                 qfull_names(os) = test_full
                 qfull_names(os)(23:32) = domes(1:9) // type(1:1)
             else 
*                Try to see if unique with out the point code
                 os = 0
                 do i = 1, qnum_sites
                    if( test_name(1:4).eq.qsite_names(i)(1:4) .and.
     .                  test_name(8:8).eq.qsite_names(i)(8:8) ) then
 
*                       This matches, if os is still zero save
                        os = i
                        qfull_names(os) = test_full  
                        qfull_names(os)(23:32) = domes(1:9) // type(1:1)
                   end if
                 end do
                 if( os.eq.0 ) then
                     write(*,210) test_name
 210                 format('Site ',a8,' Not used in ESTIMATES,',
     .                             ' or in Occ')
                 end if
            end if
          end if
          if( line(1:1).eq.'-' .or. ierr.ne.0 ) then
              end_block = .true.
          end if
      end do
 
****  Thats all
      return
      end
 
CTITLE DEC_SITE_RCV
 
      subroutine dec_site_rcv(unit, line )

      implicit none 
 
*     Routine to read the site RECEIVER block
 
      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
 
* PASSED VARIABLES
 
*   unit        - Unit number of sinex file
 
      integer*4 unit
 
*   line        - Line read from file
 
      character*(*) line
 
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
 
      integer*4 ierr, ns,  refs_yr, refs_doy, refs_sec ,
     .          refe_yr, refe_doy, refe_sec, occ, indx
 
*   end_block   - Set true at end of block
*   time_OK     - Set true is values in time range of experiment.
 
      logical end_block, time_OK

*   st_ep, en_ep - Start and stop epoch (checked for range)

      real*8 st_ep, en_ep
 
*   code        - 4 character code from which we generate name
*   type        - P=GPS, R=VLBI, S=SLR
*   pt          - Point number.  This is added to the end of the four character
*                 code.  Later if it is not A we will make _GP[Pt]
*   chocc       - Character version of occupation number.

      character*4 code, type, pt, chocc

*   recv_ty     - Reveiver type
*   recv_sn     - Receiver serial number
*   recv_fw     - Receiver firmware
*   full_name   - Full name of site (with occ+pt)
*   recv_ty_full - Full name of receiver type read from file.
*   recv_fw_full - Full firrm ware read from file. 

      character*8 full_name
* MOD TAH 101015: Made recv_ty full C*20; last 4-characters will be
*     copied to recv_ty_end later
      character*16 recv_fw
      character*20 recv_ty, recv_ty_full
      character*11 recv_fw_full 
      character*8  recv_sn

****  Start decoding
 
      end_block = .false.
      ierr = 0
      ns = 0
      do while ( ierr.eq.0 .and. .not.end_block )
 
          read(unit,'(a)', iostat=ierr) line
          if( ierr.eq.0 .and. line(1:1).eq.' ' ) then
 
****          Decode the line
              if( cglb_vers.eq.5 ) then
                  read(line, 120) code, pt, chocc, type, refs_yr,
     .                            refs_doy, refs_sec, refe_yr,
     .                            refe_doy, refe_sec, recv_ty,
     .                            recv_sn,  recv_fw
 120              format(1x,a4,1x,a2,1x,a4,1x,a1,1x,i2,1x,i3,1x,i5,
     .                   1x,i2,1x,i3,1x,i5,1x,a16,1x,a5,1x,a15)
              else
                  read(line, 125) code, pt, chocc, type, refs_yr,
     .                            refs_doy, refs_sec, refe_yr,
     .                            refe_doy, refe_sec, recv_ty_full,
     .                            recv_sn,  recv_fw_full 
 125              format(1x,a4,1x,a2,1x,a4,1x,a1,1x,i2,1x,i3,1x,i5,
     .                   1x,i2,1x,i3,1x,i5,1x,a20,1x,a5,1x,a11)
                  call trimlead(recv_sn)
                  call trimlead(recv_ty_full)
                  recv_ty = recv_ty_full
                  recv_fw(1:11) = recv_fw_full
                  recv_fw(12:15) = recv_ty_full(17:20)
              end if

              call check_full_name( code, chocc, pt, full_name,
     .                             'RECEIVER' ) 
              indx = 1
              call get_cmd(full_name, qsite_names, qnum_sites, 
     .                         ns,indx ) 

* MOD TAH 100825: Make sure we have the name.  When if occ code
*             is different from 1 we may need to add entry
              if( ns.le.0 ) then
                   write(*,*) 'REC ',trim(line), full_name
                   call Add_site(full_name, ns)
              end if


              call yds_to_jd( refs_yr,refs_doy,refs_sec, st_ep)
              call yds_to_jd( refe_yr,refe_doy,refe_sec, en_ep)
              call check_rng( st_ep, en_ep, time_OK )
              if( ns.gt.0 .and. time_OK ) then

****              OK Save the information about the receiver
                  qrecv_st(ns) = st_ep
                  qrecv_en(ns) = en_ep
                  qrecv_ty(ns) = recv_ty
                  qrecv_sn(ns) = recv_sn
                  qrecv_fw(ns) = recv_fw
              end if

          end if
          if( line(1:1).eq.'-' .or. ierr.ne.0 ) then
              end_block = .true.
          end if
      end do
 
****  Thats all
      return
      end
      
CTITLE DEC_SITE_ANT
 
      subroutine dec_site_ant(unit, line )
 
      implicit none 

*     Routine to read the site ANTENNA block
 
      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
 
* PASSED VARIABLES
 
*   unit        - Unit number of sinex file
 
      integer*4 unit
 
*   line        - Line read from file
 
      character*(*) line
 
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
 
      integer*4 ierr, ns,  refs_yr, refs_doy, refs_sec ,
     .          refe_yr, refe_doy, refe_sec, occ, indx
 
*   end_block   - Set true at end of block
*   time_OK     - Set true is values in time range of experiment.
 
      logical end_block, time_OK

*   st_ep, en_ep - Start and stop epoch (checked for range)

      real*8 st_ep, en_ep

* MOD TAH 200205: Added Antenna aligment to True north
      integer*4 iantdaz   ! Antenna aligment to True north _DAZ
 
*   code        - 4 character code from which we generate name
*   type        - P=GPS, R=VLBI, S=SLR
*   pt          - Point number.  This is added to the end of the four character
*                 code.  Later if it is not A we will make _GP[Pt]

      character*4 code, type, pt

*   ante_ty     - Reveiver type
*   ante_sn     - Receiver serial number
*   ante_fw     - Receiver firmware
*   full_name   - Full name of site (with occ+pt)

      character*8 full_name

      character*16 ante_ty
      character*20 ante_ty_full 
      character*8  ante_sn
      character*4  chocc, chocc_copy


****  Start decoding
 
      end_block = .false.
      ierr = 0
      ns = 0
      do while ( ierr.eq.0 .and. .not.end_block )
 
          read(unit,'(a)', iostat=ierr) line
          if( ierr.eq.0 .and. line(1:1).eq.' ' ) then
 
****          Decode the line
              if( cglb_vers.eq.5 ) then
                  read(line, 120) code, pt, chocc, type, refs_yr,
     .                            refs_doy, refs_sec, refe_yr,
     .                            refe_doy, refe_sec, ante_ty,
     .                            ante_sn
 120              format(1x,a4,1x,a2,1x,a4,1x,a1,1x,i2,1x,i3,1x,i5,
     .                   1x,i2,1x,i3,1x,i5,1x,a16,1x,a5)
              else
* MOD TAH 200205: Added antdaz for azimuth alignment.
                  read(line, 125) code, pt, chocc, type, refs_yr,
     .                            refs_doy, refs_sec, refe_yr,
     .                            refe_doy, refe_sec, ante_ty_full,
     .                            ante_sn, iantdaz
 125              format(1x,a4,1x,a2,1x,a4,1x,a1,1x,i2,1x,i3,1x,i5,
     .                   1x,i2,1x,i3,1x,i5,1x,a20,1x,a5,1x,I4)    
                  call trimlead(ante_ty_full)
                  ante_ty = ante_ty_full(1:16)
                  call trimlead(ante_sn)
                  ante_sn(6:8) = ante_ty_full(17:19) 
              end if

              call check_full_name( code, chocc, pt, full_name,
     .                             'ANTENNA' ) 
              indx = 1
              call get_cmd(full_name, qsite_names, qnum_sites, 
     .                         ns,indx ) 
* MOD TAH 100825: Make sure we have the name.  When if occ code
*             is different from 1 we may need to add entry
              if( ns.le.0 ) then
                   write(*,*) 'ANT ',trim(line), full_name
                  call Add_site(full_name, ns)
              end if

              call yds_to_jd( refs_yr,refs_doy,refs_sec, st_ep)
              call yds_to_jd( refe_yr,refe_doy,refe_sec, en_ep)
              call check_rng( st_ep, en_ep, time_OK )
              if( ns.gt.0 .and. time_OK ) then

***               OK Save the information about the receiver
                  qante_st(ns) = st_ep
                  qante_en(ns) = en_ep
                  qante_ty(ns) = ante_ty
* MOD TAH 050709: Save the radome type
                  qradome_ty(ns) = ante_ty_full(17:20)
                  qante_sn(ns) = ante_sn
* MOD TAH 200205: save the antenna alignent.
                  qantdaz(ns) = iantdaz
               end if

          end if
          if( line(1:1).eq.'-' .or. ierr.ne.0 ) then
              end_block = .true.
          end if
      end do
 
****  Thats all
      return
      end

CTITLE DEC_SITE_PHS
 
      subroutine dec_site_phs(unit, line )
 
      implicit none 

*     Routine to read the site PHASE CENTER MODEL block
 
      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
 
* PASSED VARIABLES
 
*   unit        - Unit number of sinex file
 
      integer*4 unit
 
*   line        - Line read from file
 
      character*(*) line
 
* LOCAL VARIABLES
 
*   ierr        - IOSTAT error
*   i, j        - Loop counters
 
      integer*4 ierr, i, jerr, trimlen

*   L1a_ecc, L2a_ecc  -- L1 and L2 phase center offsets (UNE form)

      real*4 L1a_ecc(3), L2a_ecc(3)
 
*   end_block   - Set true at end of block
*   Ant_set(max_glb_sites)  - Set true once a full antenna model has been set
 
      logical end_block,Ant_set(max_glb_sites) 

*   ante_ty     - Reveiver type
*   ante_sn     - Receiver serial number
*   ant_mod     - Type of annetta phase model used
*   ant_cod     - CODE for the model.  phsmod_to_cod converts back and
*                 forth.

      character*16 ante_ty
      character*20 ante_ty_full 
      character*8  ante_sn
      character*4  ant_cod
      character*16 ant_mod

****  Initialize the Ant_set logicals
      do i = 1, qnum_sites
         Ant_set(i) = .false.
      end do

****  Start decoding
 
      end_block = .false.
      ierr = 0
      do while ( ierr.eq.0 .and. .not.end_block )
 
          read(unit,'(a)', iostat=ierr) line
          if( ierr.eq.0 .and. line(1:1).eq.' ' ) then
 
****          Decode the line
              if( cglb_vers.eq.5 ) then
                  read(line,120) ante_ty, ante_sn, L1a_ecc, L2a_ecc,
     .                           ant_mod
 120              format(1x,a16,1x, a5,6f7.4,1x,a10)
              else
                  read(line,125,iostat=jerr) ante_ty_full, ante_sn, 
     .                           L1a_ecc,  L2a_ecc, ant_mod
C125              format(1x,a20,1x,a5,6(1x,f6.4),1x,a8)
* MOD TAH 970806: Changed to allow for mis-placed minus sign in
*                 SIO generated sinex files.
 125              format(1x,a20,1x,a5,6(f7.4),1x,a10)
                  if( jerr.ne.0 ) then
                     read(line,125,iostat=jerr) ante_ty_full, ante_sn, 
     .                           L1a_ecc,  L2a_ecc, ant_mod
                  end if

                  ante_ty = ante_ty_full(1:16)
                  call trimlead(ante_sn)
                  ante_sn(6:8) = ante_ty_full(17:19)
              endif

              if ( trimlen(ante_sn(1:5)).eq.0 ) ante_sn(1:5) = '-----' 

              call trimlead( ant_mod )        
              call phsmod_to_code( ant_mod, ant_cod, 'ENCODE')

****          Now match the antenna to the sites.
              do i = 1, qnum_sites

*                If we have not yet set this site; set the value for no-radome
*                case
                 if( .not. Ant_set(i) ) then
*                   First match antenna with no-radome case
                    if( qante_ty(i)(1:15).eq.ante_ty(1:15) .and.
     .                  (ante_ty_full(17:20).eq.'    '  .or.
     .                   ante_ty_full(17:20).eq.'NONE')  .and. 
     .                  (qante_sn(i)(1:5).eq.ante_sn(1:5) .or. 
     .                   ante_sn(1:5).eq.'-----') ) then

*                      Save the values re-arranging from UNE to
*                      NEU
                       qL1a_ecc(1,i) = L1a_ecc(2)
                       qL1a_ecc(2,i) = L1a_ecc(3)
                       qL1a_ecc(3,i) = L1a_ecc(1)
                       qL2a_ecc(1,i) = L2a_ecc(2)
                       qL2a_ecc(2,i) = L2a_ecc(3)
                       qL2a_ecc(3,i) = L2a_ecc(1)

                       qant_mod(i) = ant_mod
                       qante_ty(i)(16:16) = '-'
                    end if
*****               Now try matching with radome type as well.  Only
*                   in the latter case do we set the antenna as known         

                    if( qante_ty(i)(1:15).eq.ante_ty(1:15) .and.
     .                  qradome_ty(i).eq.ante_ty_full(17:20) .and. 
     .                  (qante_sn(i)(1:5).eq.ante_sn(1:5) .or. 
     .                   ante_sn(1:5).eq.'-----') ) then

*                      Save the values re-arranging from UNE to
*                      NEU
                       qL1a_ecc(1,i) = L1a_ecc(2)
                       qL1a_ecc(2,i) = L1a_ecc(3)
                       qL1a_ecc(3,i) = L1a_ecc(1)
                       qL2a_ecc(1,i) = L2a_ecc(2)
                       qL2a_ecc(2,i) = L2a_ecc(3)
                       qL2a_ecc(3,i) = L2a_ecc(1)

                       qant_mod(i) = ant_mod
                       Ant_set(i) = .true.
                       qante_ty(i)(16:16) = '+'
                    end if
                 end if
              end do

          end if
          if( line(1:1).eq.'-' .or. ierr.ne.0 ) then
              end_block = .true.
          end if
      end do
 
****  Thats all
      return
      end

CTITLE DEC_SITE_ECC

      subroutine dec_site_ecc(unit, line )
 
      implicit none 

*     Routine to read the site ECCENTRICITY block
 
      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
 
* PASSED VARIABLES
 
*   unit        - Unit number of sinex file
 
      integer*4 unit
 
*   line        - Line read from file
 
      character*(*) line
 
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
 
      integer*4 ierr, ns,  refs_yr, refs_doy, refs_sec ,
     .          refe_yr, refe_doy, refe_sec, occ, indx
 
*   end_block   - Set true at end of block
*   time_OK     - Set true is values in time range of experiment.
 
      logical end_block, time_OK

*   st_ep, en_ep - Start and stop epoch (checked for range)

      real*8 st_ep, en_ep
 
*   code        - 4 character code from which we generate name
*   type        - P=GPS, R=VLBI, S=SLR
*   pt          - Point number.  This is added to the end of the four character
*                 code.  Later if it is not A we will make _GP[Pt]
*   ecc_type    - Method for specifying the ECC (onlu UNE accepted at
*                 the moment)
*   chocc       - Character version of occupation number.

      character*4 code, type, pt, ecc_type, chocc


      character*8 full_name

****  Start decoding
 
      end_block = .false.
      ierr = 0
      ns = 0
      do while ( ierr.eq.0 .and. .not.end_block )
 
          read(unit,'(a)', iostat=ierr) line
          if( ierr.eq.0 .and. line(1:1).eq.' ' ) then
 
****          Decode the line
              read(line, 120) code, pt, chocc, type, refs_yr,
     .                        refs_doy, refs_sec, refe_yr,
     .                        refe_doy, refe_sec, ecc_type,
     .                        arp_ecc
 120          format(1x,a4,1x,a2,1x,a4,1x,a1,1x,i2,1x,i3,1x,i5,
     .               1x,i2,1x,i3,1x,i5,1x,a3,3f9.4)

              call check_full_name( code, chocc, pt, full_name,
     .                             'ECCENTRICITY' ) 
              indx = 1
              call get_cmd(full_name, qsite_names, qnum_sites, 
     .                         ns,indx ) 
* MOD TAH 100825: Make sure we have the name.  When if occ code
*             is different from 1 we may need to add entry
              if( ns.le.0 ) then
                  write(*,*) 'ECC ',trim(line), full_name
                  call Add_site(full_name, ns)
              end if

              call yds_to_jd( refs_yr,refs_doy,refs_sec, st_ep)
              call yds_to_jd( refe_yr,refe_doy,refe_sec, en_ep)
              call check_rng( st_ep, en_ep, time_OK )

              if( ecc_type(1:3).eq. 'UNE' .and. ns.gt.0 .and.
     .            time_OK ) then
                  qarp_ecc(1,ns) = arp_ecc(2)
                  qarp_ecc(2,ns) = arp_ecc(3)
                  qarp_ecc(3,ns) = arp_ecc(1)
              else if( ecc_type(1:3).eq. 'UEN' .and. ns.gt.0 .and.
     .            time_OK ) then
                  qarp_ecc(1,ns) = arp_ecc(3)
                  qarp_ecc(2,ns) = arp_ecc(2)
                  qarp_ecc(3,ns) = arp_ecc(1)
              else if( ns.gt.0 .and. time_OK) then
                   write(*,180) ecc_type, full_name
 180               format('**WARNING** Unknown ECC type in SITE/ECC.',
     .                    ' Type is ',a4,' Site ',a)
               end if
          end if
          if( line(1:1).eq.'-' .or. ierr.ne.0 ) then
              end_block = .true.
          end if
      end do
 
****  Thats all
      return
      end

      
CTITLE CHECK_RNG

      subroutine check_rng( st, en, OK )

      implicit none 

*     Routine to check that the start and stop ranges
*     full within the sinex file being processed.

      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'

* PASSED VARIABLES

*  st, en  - Strart and End JD for the site paramter

      real*8 st, en

*  OK      - Set true if the values are with range

      logical OK

*  LOCAL VARIABLES

*  zero_jd  -  JD when zero's are passed as date.  Indicates
*      value is unknown.

      real*8 zero_jd

      data zero_jd / 2415019.5d0 /

****  Start checking ranges
      OK = .true.
      if ( en.lt. qstart_epoch .and. en.gt. zero_jd ) OK = .false.
      if ( st.gt. qend_epoch   .and. st.gt. zero_jd ) OK = .false.

****  Thats all
      return
      end

CTITLE Add_site

      subroutine Add_site(full_name, ns)

      implicit none 

*     Routine to add a new site name when the occ code changes and
*     no site name exists for that code.  full_name has 
*     <4-character code><3-digit Occ><pt code>

      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'

      character*(*) full_name
      integer*4 ns      ! Site number

      integer*4 i 

****  First see if we can find a site with the name 4-char and PT
      ns = -1
      do i = 1, qnum_sites
         if( full_name(1:4).eq.qsite_names(i)(1:4) .and.
     .       full_name(8:8).eq.qsite_names(i)(8:8) ) then
             ns = i
             exit
         end if
      end do

****  See if we found a match
      if( ns.gt.0 ) then
         write(*,120) full_name, ns, qsite_names(ns)
 120     format('Found match for ',a,' for site ',i4,1x,a)

****     Add site into list
         qnum_sites = qnum_sites + 1
         qsite_names(qnum_sites) = full_name
         qfull_names(qnum_sites) =  qfull_names(ns)
         ns =   qnum_sites
      else
         write(*,140) full_name
 140     format('Unable to find match for ',a,' Adding')
         qnum_sites = qnum_sites + 1
         qsite_names(qnum_sites) = full_name
         qfull_names(qnum_sites) = full_name(1:4) // 
     .                             '_GPS             -----M---P'
         ns =   qnum_sites
 
      end if

***** Thats all
      return 
      end

CTITLE CHECK_FULL_NAME

      subroutine check_full_name( code, chocc, pt, full_name, block )

      implicit none

*     Routine to generate full name 

      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'

* PASSED IN:
      character*(*) code   ! 4-char code
      character*(*) chocc  ! Occupany or soln number
      character*(*) pt     ! PT-code (A-Z)
      character*(*) block  ! BLOCK Calling this routine

* PASSED OUT
      character*(*) full_name !  Full name of site

* LOCAL
      character*4 chocc_copy  ! Copy of chocc to see if '-----'
      integer*4   occ         ! Numeric value of occupancy/soln code
      integer*4   indx        ! Counter in string
      integer*4   ns          ! Site number (<= 0 if problem)

***** Start creating full name

      chocc_copy = chocc
      if( chocc.eq.'----' ) chocc = '   1'
      read(chocc,'(i4)') occ    

* MOD TAH 970716: Check for blank pt names, replace with ' A'
      if( pt(2:2).eq.' ' .or. pt(2:2).eq.'-' ) pt = ' A'
      
****  OK, construct the full name for the site and
*     see if we find it.
      write(full_name, 140) code, occ, pt(2:2)
 140  format(a4,i3.3,a1)
      indx = 1
      call get_cmd(full_name, qsite_names, qnum_sites, 
     .                 ns,indx ) 
* MOD TAH 100825: Make sure we have the name.  When if occ code
*      is different from 1 we may need to add entry
       if( ns.le.0 ) then
* MOD TAH 171002: See if just a "default" soln code (chocc_copy)
           if ( chocc_copy.eq.'----' ) then
*             Try to match just the 4-char code and warn user.
              indx = 1
              call get_cmd(code, qsite_names, qnum_sites, 
     .                  ns,indx ) 
              if( ns.gt.0 ) then 
                 write(*,210) chocc_copy, code, qsite_names(ns),
     .                        block
 210             format('Default SOLN code ',a,' Site ',a4,
     .                  ' Matched to ',a,' SINEX Block ',a)
                 full_name = qsite_names(ns)
              else
                 write(*,220) code, chocc_copy, pt, block
 220             format('Non-unique named decoded from ',a,1x,a,1x,a,/,
     .                  'Check ',a,' BLOCK')
                 stop 'HTOGLB: Non-unique name in REC/ANT/ECC BLOCK'
              endif
           endif
       endif

****   Thats all
       return
       end

