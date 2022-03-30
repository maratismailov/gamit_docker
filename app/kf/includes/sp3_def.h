*    Include file for definitions of Multi-GNSS SP3 files.
*    (Merger from track and svsp3 include files).

* PRNs are mapped with the additions of constants
* GPS     G    +0
* GLONASS R  +100
* GALIELO E  +200
* BEIDOU  C  +300
* QZSS    J  +400
* SBAS    S  +500
* IRNSS   I  +600
 
*  max_sat       - Maximum number of satellites
*  max_eph_eps   - Maximum number of ephemeris epochs in the sp3 files.
 
      integer*4 max_sat, max_eph_eps

* MOD TAH 200303: Changed from 100 to 150 for maximum number of satellites 
      parameter ( max_sat = 150 )
      parameter ( max_eph_eps = 1000 )

      integer*4 max_gnss  ! Maximum number of GNSS types (make constinent
                          ! with function conoff which gives the constelation
                          ! PRN offset
      parameter ( max_gnss = 7 )


      integer*4 max_dna         ! Max SVS dnadir values
* MOD TAH 180705: Inceased max nadir values to 18.  Added check to read
* MOD TAH 180712: Inceased max nadir values to 41 to accomadate Galileo chamber calibration. 
      parameter ( max_dna = 41  )

* Multi-GNSS Frequencies

      real*8 fL1(max_sat), fL2(max_sat)   !  Frequnecies for each
                  ! satellite (GLONASS determined by SVNAV.dat).
      real*8 fL5(max_sat)   ! Frequency for possible 3rd frequency
                  ! (Frequencies are set when the SP3 file is read).

*  num_sat       - Number of satelites found
*  num_sp3       - Number of epoch entries in the sp3 file
*  prn_sp3(max_sat)  - PRN numbers of sateliites (0 initially)

      integer*4 num_sat, num_sp3, prn_sp3(max_sat)

      logical out_sp3freqs  ! Set .true. to output SP3 satellite frequencies
      logical out_sp3PCV    ! Set .true. to output antenna offsets

*  sp3_xyz(3,max_sat, max_eph_eps) - Earth fixed coordinates of satellites
*                          (XYZ, m)
*  sp3_time(maxmax_ephs_eps) - JD's for the ephemeris entries
*  sp3_clk(max_sat,max_eph_eps)       - Clock error by time (us)
*  sp3_sec(max_ephs_eps)  - Seconds of day part of the SP3 time tags running
*       from the first day of the file.
*  sp3_refmjd  - Reference MJD for SP3 file

      real*8  sp3_xyz(3,max_sat, max_eph_eps), 
     .    sp3_time(max_eph_eps), sp3_clk(max_sat,max_eph_eps),
     .    sp3_sec(max_eph_eps) , sp3_refmjd

* num_svs_dph(2,max_sat) -- Number of phase entries in dnadir and azimuth directions
*      for satellites
      integer*4 num_svs_dph(2,max_sat)
* svs_dphs(max_dna,1, 2,max_sat) -- Phase center adjustments at function of Nadir, azimuth,
*      L1/L2 by satellite (mm)
* svs_dna(3,2,max_sat) -- zen1, zen2, dzen, az1, az2, daz for each satellite (deg)
* svs_L12(3,2,max_sat) -- "North, East, Up" offsets at L1/L2 for each satellite (mm)

      real*4 svs_dphs(max_dna,1, 2,max_sat)
      real*8 svs_dna(3,2,max_sat), svs_L12(3,2,max_sat)


****  COMMON DECLARATIONS

      common / sp3_R8 / sp3_xyz, sp3_time, sp3_clk, sp3_sec,  
     .                  sp3_refmjd, fL1, fL2, svs_dna, svs_L12

      common / sp3_I4 / num_sat, num_sp3, prn_sp3,  num_svs_dph, 
     .       out_sp3freqs, out_sp3PCV

      common / sp3_R4 / svs_dphs
