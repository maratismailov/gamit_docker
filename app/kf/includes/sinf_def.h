*-----------------------------------------------------------------
* sinf_def.h 
*     Include file to define the station information to be included
*     in version 1.0 global files.

*     ssdata_st  - Start JD for the data in this solution for this
*                  site.  This nee
*     ssdata_en  - End JD for the data in this solution.

*     ssrecv_st  - Start JD for the receiver characteristics (not
*                  the data spanned in this solution)
*     ssrecv_en  - End JD for  receiver characteristics
*     ssante_st  - Start JD for the antenna characteristics (not
*                  the data spanned in this solution) 
*     ssante_en  - End JD for antenna characteristics

      real*8 ssdata_st, ssdata_en, ssrecv_st, ssrecv_en, 
     .       ssante_st, ssante_en
 
*      ssarp_ecc(3) - Single station eccentricity of the
*           - antenna with respect ground mark
*           - (NEU m)
*     ssL1a_ecc(3)  - L1 phase center to ARP (NEU m).
*     ssL2a_ecc(3)  - L2 phase center to ARP (NEU m).
*     sselev_cut    - Elevation cutoff angle (rads)
*     sscons_size   - Ratio of final to apriori variance (used to
*                     get the constraint type).
 
      real*4 ssarp_ecc(3), ssL1a_ecc(3), ssL2a_ecc(3), sselev_cut,
     .       sscons_size
      real*4 ssantdaz ! Antenna aligment from True N (deg). Value
                      ! removes one slot from ssres (now a dummy array). 

* Log estimate information
      real*8 slog_ep   ! JD for for log for site (JD)
      real*4 slog_tau  ! Time constant for log (days)

* MOD TAH 050622: Added load contributions for sites
*   satmload(3) -- Average daily load applied at this site (NEU mm)
*   shydload(3) -- Average hydrographic load (NEU mm)

       real*4 satmload(3),  shydload(3)
 
 
*         ssnum_zen - Number of zenith delay parameters
 
      integer*4 ssnum_zen
 
*           ssant_mod   - Antenna model used for this station
 
      character*4 ssant_mod
 
*           sscode  - 8 character code name for station.  This
*           - needs to match the gsite_names array.
*   ssante_sn   - Antenna serial number (----- if not known)
*   ssrecv_sn   - Serial number for the receiver (----- if
*           - not known)
 
      character*8 sscode, ssante_sn, ssrecv_sn
 
*            ssrecv_sty  - Type of receiver
*   ssante_ty   - Antenna type
*   ssrecv_fw   - Firm ware version for the receiver.
 
      character*16 ssrecv_ty, ssante_ty, ssrecv_fw

* MOD TAH 050622: Add radome type to model
      character*4 ssradome_ty

* MOD TAH 101005: Added sti_antmod with name of antenna model
*     SINEX needs C*10, we add ELEV, AZEL for C*16. Old ssant_mod
*     set to FULL when full name is available.
*     ssrecv_ty_end is last 4-characters of receiver type to allow
*     character*20 definition 
      character*16 sti_antmod
      character*4  ssrecv_ty_end 
* MOD TAH 190627: Add ss_lnum for local site number when BACK type
*     GLX file is created in GLBAK.  Used to see if values to be
*     output.
      integer*4 ss_lnum
 
*         ssres(21) - We have 21 I*4 words free for future
*           - expansion.
* Remove 3 I*4 words TAH 030615: For slog_ep (R*8) and slog_tau (R*4) 
* Remove 6 R*4 words, 1 C*4 TAH 050622: For load contributions and radome
C     integer*4 ssres(21)
C     integer*4 ssres(18)
C     integer*4 ssres(11)
C Remove 5 I*4 words for sti_antmod, and ssrecv_ty_end
C      integer*4 ssres(6)
C Remove 1 I*4 word for ss_lnum.
C      integer*4 ssres(5)
C Remove 1 R*4 for ssantdaz (alignment of antenna)
       integer*4 ssres(4)

*-------------------------------------------------------------------------
 
      common / ssdef_rec / ssdata_st, ssdata_en, ssrecv_st, ssrecv_en, 
     .       ssante_st, ssante_en,
     .       ssarp_ecc, ssL1a_ecc, ssL2a_ecc, sselev_cut, sscons_size,
     .       ssnum_zen , ssant_mod, sscode, ssante_sn, ssrecv_sn,
     .       ssrecv_ty, ssante_ty, ssrecv_fw, slog_tau, slog_ep, 
     .       satmload,  shydload, ssradome_ty, ss_lnum, ssantdaz, 
     .       ssres, sti_antmod, ssrecv_ty_end
 
