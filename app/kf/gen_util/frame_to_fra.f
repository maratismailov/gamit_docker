CTITLE FRAME_TO_FRAME
 
      subroutine frame_to_frame(sys_frame, out_frame, rot_vec )

      implicit none 
 
 
*     Routine to get the rotation rate bewteen two frames.  The
*     frames are defined in the frame_names array and all rotation
*     rates are given with respect to NUVEL_NNR field.
* MOD TAH 031031: Added automatic checking of frames.dat or 
*     ~/gg/tables/frames.dat to see if sys_frame or out_frame can
*     be found, if they are not found in standard list.  First
*     ocrrurence of frame definition in frames.dat used.
* MOD TAH 120120: Added LIST option to sys_frame which will list
*     Available frames
* MOD TAH 180117: Updated poles for ITRF2014.  Added 11 plates plus
*     3 for ITRF14, IGS14 and NAM14.
*     Citation: Altamimi, Zuheir, Laurent Métivier, Paul Rebischung, Hélène Rouby, and 
*               Xavier Collilieux. 2017. "ITRF2014 Plate Motion Model." 
*               Geophysical Journal International 209 (3). Oxford University Press: 
*               1906-12. doi:10.1093/gji/ggx136.


      include '../includes/const_param.h'
 
 
*   max_frames  - Maximum number of frames supported
 
      integer*4 max_frames
 
*     parameter ( max_frames = 63 )
* MOD TAH 180117: Added 14 ITRF2014 plate names and poles
      parameter ( max_frames = 78 )
 
* PASSED VARIABLES
 
*   rot_vec(3)  - Rotation vector to use
 
      real*8 rot_vec(3)
 
*   sys_frame   - Frame of current system
*   out_frame   - Desired frame
 
      character*(*) sys_frame, out_frame
 
* LOCAL VARIABLES
 
*   i       - Loop counter
*   nfs, nfo    - Numbers of the two frames passed
*   lenf    - Length of frame name
*   trimlen - Length of string
 
      integer*4 i, nfs, nfo, lenf, trimlen, ierr, jerr, lenh,
     .          indx, ifnd
 
*   frame_data(3,max_frames)    - Data for the frames based
*             on IERS standards (deg/Myr)
*   scs, sco   - Scale factors to get NUVEL_1A for the system
*             and output frames
 
      real*8 frame_data(3,max_frames), scs, sco, frame_ev(3),
     .       frame_es(3), frame_eo(3)
 
*   frame_names(max_frames) - Names of frames supported
 
      character*8 frame_names(max_frames)

*   newframe -- Framce read from frames.dat file
*   frame_un -- Units for frame vectors (default is deg/Myr but
*               using rad will interpret values as rad/Myr.
*   cdum     -- Dummy string for multiread

      character*8 newframe, frame_un, cdum

      character*128 line, frame_file, home_dir
 
      data frame_names / 'PCFC    ', 'COCO    ', 'NAZC    ', 'CARB    '
     .                ,  'SAFD    ', 'ANTA    ', 'INDI    ', 'AUST    '
     .                ,  'AFRC    ', 'ARAB    ', 'EURA    ', 'NAFD    '
     .                ,  'JUAN    ', 'PHIL    ', 'RIVERA  ', 'SCOTIA  '
     .                ,  'NUV-NNR ', 'AM-02   ', 'ITRF93  ', 'ITRF94  '
     .                ,  'GG_PCFC ', 'ANTA_I00', 'AUST_I00', 'EURA_I00'
     .                ,  'NOAM_I00', 'PCFC_I00', 'SOAM_I00', 'ITRF00  '
     .                ,  'ARAB_MCC', 'ARAB_M06', 'AMUR_I05', 'ANTA_I05'
     .                ,  'ARAB_I05', 'AUST_I05', 'CARB_I05', 'EURA_I05'
     .                ,  'INDI_I05', 'NAZC_I05', 'NOAM_I05', 'NUBI_I05'
     .                ,  'OKHT_I05', 'PCFC_I05', 'SOAM_I05', 'SOMA_I05'
     .                ,  'YANG_I05', 'ITRF05  ', 'ITRF08  ', 'NUBI_I08'
     .                ,  'ANTA_I08', 'ARAB_I08', 'AUST_I08', 'CARB_I08'
     .                ,  'EURA_I08', 'INDI_I08', 'NAZC_I08', 'NOAM_I08' 
     .                ,  'PCFC_I08', 'SOAM_I08', 'SUND_I08', 'SOMA_I08'
     .                ,  'AMUR_I08', 'ANTA_I14', 'ARAB_I14', 'AUST_I14'
     .                ,  'EURA_I14', 'INDI_I14', 'NAZC_I14', 'NOAM_I14'
     .                ,  'NUBI_I14', 'PCFC_I14', 'SOAM_I14', 'SOMA_I14'
     .                ,  'NAM08   ', 'IGS08   ', 'ITRF14  ', 'IGS14   '
     .                ,  'NAM14   ', 'ANT14   ' /
     

*     Nuvel-1A Rates in degrees per million year.  Does not include the
*     revised pacific rate from Demets, 3545-3548, 1995.  From this paper
*     it is not clear what the rate should be.  Tables from IERS convertions
*     convert to deg/My  (*57.2958)
      data frame_data /-0.08652d0,   0.27731d0,  -0.57124d0, 
     .                 -0.59731d0,  -1.23788d0,   0.62596d0, 
     .                 -0.08778d0,  -0.49143d0,   0.55056d0, 
     .                 -0.01020d0,  -0.19395d0,   0.09058d0, 
     .                 -0.05947d0,  -0.08680d0,  -0.04985d0, 
     .                 -0.04704d0,  -0.09746d0,   0.21234d0, 
     .                  0.38216d0,   0.00229d0,   0.38904d0, 
     .                  0.44914d0,   0.29358d0,   0.35993d0, 
     .                  0.05105d0,  -0.17756d0,   0.22471d0, 
     .                  0.38302d0,  -0.02985d0,   0.38732d0, 
     .                 -0.05621d0,  -0.13722d0,   0.18065d0, 
     .                  0.01478d0,  -0.20621d0,  -0.00877d0, 
     .                  0.29794d0,   0.49332d0,  -0.33346d0, 
     .                  0.57811d0,  -0.41024d0,  -0.55405d0, 
     .                 -0.53801d0,  -1.77388d0,   0.69041d0, 
     .                 -0.02349d0,  -0.15241d0,  -0.07277d0, 
     .                    0.0000d0,  0.0000d0,  0.0000d0,
     .                   -0.0178d0,  0.0128d0,  0.0049d0,
     .                    0.0333d0,  0.0806d0, -0.0056d0,
     .                    0.0000d0,  0.0000d0,  0.0000d0,
     .                 -0.106303d0,  0.27194d0,  -0.60755d0,
* Following value from ITRF 2000 paper [Altamimi et al. JGR 107, 2214,doi:10.1029/2001/JB000561, 2002] 
* Plate        Phi    Long    Rate (deg/myr)
* ANTA_I00   61.830 -125.574 0.231  2.143 3.689 0.015
* AUST_I00   32.327   39.437 0.614  0.652 0.816 0.006
* EURA_I00   57.965  -99.374 0.260  1.211 2.710 0.005
* NOAM_I00   -5.036  -83.144 0.194  1.142 1.945 0.003
* PCFC_I00  -64.176  110.194 0.666  0.404 1.345 0.005
* SOAM_I00  -21.457 -134.631 0.113  2.806 4.762 0.005
     .                  -0.06344d0,    -0.08870d0,     0.20364d0,
     .                   0.40071d0,     0.32958d0,     0.32834d0,
     .                  -0.02246d0,    -0.13607d0,     0.22041d0,
     .                   0.02307d0,    -0.19187d0,    -0.01703d0,
     .                  -0.10015d0,     0.27228d0,    -0.59949d0,
     .                  -0.07388d0,    -0.07484d0,    -0.04134d0, 
* ITRF00
     .                   0.0000d0,      0.0000d0,      0.0000d0,  
* Following value for ARAB MCC from GPS CONSTRAINTS ON AFRICA (NUBIA) AND ARABIA PLATE MOTIONS, 
*   S. McClusky, R. Reilinger, S. Mahmoud, D. Ben Sari, A. Tealeb. GJI, 2003 
*   www-gpsg.mit.edu/edocs/MCCLUSKY_2003_GJI/africa_arabia_GJI-format+figs.pdf
*   lat 27.4, Long 18.4, Rate 0.4 deg/Myrs (AR-EU: Converted to NNR frame with 
*   ITRF00 Eurasia pole.
     .                   0.31451d0,    -0.02398d0,     0.40449d0,    
* Following value for ARAB M06 from From Robert Reilinger, Simon McClusky, Philippe Vernant, 
* Shawn Lawrence, * Semih Ergintav, Rahsan Cakmak, Haluk Ozener, Fakhraddin Kadirov, Ibrahim Guliev, 
* Ruben Stepanyan. (2006) GPS constraints on continental deformation in the 
* Africa-Arabia-Eurasia continental collision zone and implications for the 
* dynamics of plate interactions. Journal of Geophysical Research 111:b5, B05411
     .                   0.33478d0,    -0.01723d0,     0.42398d0,  
* Following values from ITRF2005 paper (Altamimi et al. JGR 112, B09401, doi:10.1029/2007/JB004949, 2007).
* Plate        Phi    Long    Rate (deg/myr)
* AMUR_I05    56.263  -102.789  0.269
* ANTA_I05    59.813  -125.315  0.223
* ARAB_I05    49.642     5.061  0.579   
* AUST_I05    32.407    37.367  0.628
* CARB_I05    39.318  -104.279  0.241
* EURA_I05    56.330   -95.979  0.261
* INDI_I05    49.823    21.841  0.614
* NAZC_I05    45.101  -101.441  0.642
* NOAM_I05    -4.291   -87.385  0.192
* NUBI_I05    49.955   -82.501  0.269
* OKHT_I05   -32.041  -132.910  0.083
* PCFC_I05   -62.569   112.873  0.682
* SOAM_I05   -16.800  -129.631  0.121
* SOMA_I05    53.661   -89.542  0.309 
* YANG_I05    59.425  -109.737  0.310  
     .               -0.0331d0, -0.1457d0,  0.2237d0,
     .               -0.0648d0, -0.0915d0,  0.1928d0,
     .                0.3735d0,  0.0331d0,  0.4412d0,
     .                0.4214d0,  0.3218d0,  0.3366d0,
     .               -0.0460d0, -0.1807d0,  0.1527d0,
     .               -0.0151d0, -0.1439d0,  0.2172d0,
     .                0.3677d0,  0.1474d0,  0.4691d0,
     .               -0.0899d0, -0.4442d0,  0.4548d0,
     .                0.0087d0, -0.1913d0, -0.0144d0,
     .                0.0226d0, -0.1716d0,  0.2059d0,
     .               -0.0479d0, -0.0515d0, -0.0440d0,
     .               -0.1221d0,  0.2895d0, -0.6053d0,
     .               -0.0739d0, -0.0892d0, -0.0350d0,
     .                0.0015d0, -0.1831d0,  0.2489d0,
     .               -0.0533d0, -0.1484d0,  0.2669d0,
* ITRF05
     .                    0.d0,  0.d0    ,  0.d0,
* ITRF08  
     .                    0.d0,  0.d0    ,  0.d0,
* Following values from ITRF2008 analysis, sent to RWK by Z. Altamimi,
* 6 Jan 2012, and distributed to GAMIT users.  However, there was one
* further revision of the paper itself (Altamimi et al. J. Geophys. Res.
* 117, B07402, doi: 10.1029/2011JB008930, 2012), which we've now used to 
* supersede the first set.  Note that Altamimi recommends using translation 
* rates (Tx Ty Tz) = (0.41 0.22 0.41)  mm/yr with these rotation vectors.  
* The translation rates are not yet included in this  other routines.
*       Wx       Wy      Wz (mas/yr) W (deg/My)  Long       Lat     
*  Values and sigmas 120112
* NUBI   0.0973  -0.6175   0.7306      0.2671    -81.0437     49.4493
* Sigma  0.0076   0.0043   0.0043      0.0012      0.7058      0.2448
* 
* ANTA  -0.2446  -0.3145   0.6354      0.2083   -127.8712     57.9091
*        0.0030   0.0032   0.0087      0.0023      0.4124      0.3182
* 
* ARAB   1.1915  -0.0805   1.4650      0.5250     -3.8664     50.8141
*        0.0624   0.0805   0.0474      0.0201      4.0523      0.4666
* 
* AUST   1.5156   1.1770   1.2164      0.6311     37.8322     32.3697
*        0.0048   0.0045   0.0040      0.0010      0.1839      0.0786
* 
* CARB   0.0524  -1.1227   0.6898      0.3663    -87.3300     31.5387
*        0.0877   0.1773   0.0617      0.0518      4.0458      1.8847
* 
* EURA  -0.0788  -0.5310   0.7673      0.2601    -98.4423     55.0253
*        0.0072   0.0016   0.0081      0.0016      0.7519      0.3829
* 
* INDI   1.2464   0.2972   1.5344      0.5553     13.4131     50.1350
*        0.0169   0.0653   0.0147      0.0084      2.6768      0.4575
* 
* NAZC  -0.3408  -1.5548   1.6370      0.6342   -102.3625     45.8045
*        0.0073   0.0218   0.0080      0.0035      0.2446      0.4926
* 
* NOAM   0.0275  -0.6752  -0.0729      0.1888    -87.6655     -6.1542
*        0.0021   0.0066   0.0058      0.0017      0.1688      0.5394
* 
* PCFC  -0.4088   1.0495  -2.1659      0.6781    111.2833    -62.5232
*        0.0054   0.0030   0.0026      0.0007      0.2847      0.0620
* 
* SOAM  -0.2430  -0.3250  -0.1524      0.1204   -126.7906    -20.5772
*        0.0067   0.0068   0.0045      0.0012      1.2383      0.5389
* 
* SUND   0.0109  -0.7900   0.9649      0.3464    -89.2070     50.6898
*        0.1502   0.6195   0.0192      0.1059     10.2694     22.5729
* 
* SOMA  -0.0578  -0.7451   0.8869      0.3222    -94.4351     49.8809
*        0.0184   0.0216   0.0069      0.0051      1.2874      0.7480
* 
* AMUR  -0.1825  -0.4334   0.8970      0.2813   -112.8392     62.3344
*        0.0255   0.0312   0.0346      0.0062      4.2836      1.8613
*
* Values and sigmas from Altamimi et al. (2012)   
* Plate	NSa	omegax (mas/a)	omegay (mas/a)	omegaz (mas/a)	omega (degree/Ma)	WRMS E	N
* AMUR	3	-0.190 +/- 0.040	-0.442 +/- 0.051	0.915 +/- 0.049	0.287 +/- 0.008	0.14	0.24
* ANTA	9	-0.252 +/- 0.008	-0.302 +/- 0.006	0.643 +/- 0.009	0.209 +/- 0.003	0.40	0.29
* ARAB	4	1.202 +/- 0.082	    -0.054 +/- 0.100	1.485 +/- 0.063	0.531 +/- 0.027	0.23	0.15
* AUST	19	1.504 +/- 0.007	     1.172 +/- 0.007	1.228 +/- 0.007	0.630 +/- 0.002	0.29	0.25
* CARB	2	0.049 +/- 0.201	    -1.088 +/- 0.417	0.664 +/- 0.146	0.354 +/- 0.122	0.06	0.04
* EURA	69	-0.083 +/- 0.008	-0.534 +/- 0.007	0.750 +/- 0.008	0.257 +/- 0.002	0.34	0.28
* INDI	4	1.232 +/- 0.031	     0.303 +/- 0.128	1.540 +/- 0.030	0.554 +/- 0.017	0.55	0.55
* NAZC	3	-0.330 +/- 0.011	-1.551 +/- 0.029	1.625 +/- 0.013	0.631 +/- 0.005	0.09	0.08
* NOAM	44	0.035 +/- 0.008	    -0.662 +/- 0.009   -0.100 +/- 0.008	0.186 +/- 0.002	0.27	0.32
* NUBI	11	0.095 +/- 0.009	    -0.598 +/- 0.007	0.723 +/- 0.009	0.262 +/- 0.003	0.26	0.35
* PCFC	23	-0.411 +/- 0.007	 1.036 +/- 0.007   -2.166 +/- 0.009	0.677 +/- 0.002	0.42	0.44
* SOAM	10	-0.243 +/- 0.009	-0.311 +/- 0.010   -0.154 +/- 0.009	0.118 +/- 0.002	0.44	0.34
* SOMA	3	-0.080 +/- 0.028	-0.745 +/- 0.030	0.897 +/- 0.012	0.325 +/- 0.007	0.28	0.21
* SUND	2	0.047 +/- 0.381	    -1.000 +/- 1.570	0.975 +/- 0.045	0.388 +/- 0.308	0.08	0.05
* ITRF2008-PMM						0.33	0.31
*                                                     
*  These are the Jan 2012 values, now commented out:
*              Wx           Wx          Wx   in deg/My
*    .     0.2703D-01, -0.1715D+00,  0.2029D+00,
*    .    -0.6794D-01, -0.8736D-01,  0.1765D+00,
*    .     0.3310D+00, -0.2236D-01,  0.4069D+00,
*    .     0.4210D+00,  0.3269D+00,  0.3379D+00,
*    .     0.1456D-01, -0.3119D+00,  0.1916D+00,
*    .    -0.2189D-01, -0.1475D+00,  0.2131D+00,
*    .     0.3462D+00,  0.8256D-01,  0.4262D+00,
*    .    -0.9467D-01, -0.4319D+00,  0.4547D+00,
*    .     0.7639D-02, -0.1876D+00, -0.2025D-01,
*    .    -0.1136D+00,  0.2915D+00, -0.6016D+00,
*    .    -0.6750D-01, -0.9028D-01, -0.4233D-01,
*    .     0.3028D-02, -0.2194D+00,  0.2680D+00,
*    .    -0.1606D-01, -0.2070D+00,  0.2464D+00,
*    .    -0.5069D-01, -0.1204D+00,  0.2492D+00 /
* These from the paper:
     .     0.2639D-01, -0.1661D+00,  0.2008D+00,
     .    -0.7000D-01, -0.8389D-01,  0.1786D+00,
     .     0.3339D+00, -0.1500D-01,  0.4125D+00,
     .     0.4178D+00,  0.3256D+00,  0.3411D+00,
     .     0.1361D-01, -0.3022D+00,  0.1844D+00,
     .    -0.2306D-01, -0.1483D+00,  0.2083D+00,
     .     0.3422D+00,  0.8417D-01,  0.4278D+00,
     .    -0.9167D-01, -0.4308D+00,  0.4514D+00,
     .     0.9722D-02, -0.1839D+00, -0.2778D-01,
     .    -0.1142D+00,  0.2878D+00, -0.6017D+00,
     .    -0.6750D-01, -0.8639D-01, -0.4278D-01,
     .     0.1306D-01, -0.2778D+00,  0.2708D+00,
     .    -0.2222D-01, -0.2069D+00,  0.2492D+00,
     .    -0.5278D-01, -0.1228D+00,  0.2542D+00,
* Values from Table-S1.txt from Altimini et al., 2017 (see reference
* above) Supplemental file.
* Altamimi, Z., L. Metivier, P. Rebischung, H. Rouby, X. Collilieux; 
*     ITRF2014 plate motion model, Geophysical Journal International, 
*     Volume 209, Issue 3, 1 June 2017, Pages 1906-1912, 
*     https://doi.org/10.1093/gji/ggx136
* Table S1: ITRF2014-PMM rotation pole cartesian components in DEG/My.
*--------------------------------------
*   plate  omega_x   omega_y   omega_z
*          -----------DEG/My----------
*--------------------------------------
*   ANTA   -0.0688   -0.0900    0.1874
*   ARAB    0.3205   -0.0378    0.4011
*   AUST    0.4194    0.3284    0.3375
*   EURA   -0.0235   -0.1476    0.2140
*   INDI    0.3205   -0.0014    0.4038
*   NAZC   -0.0925   -0.4290    0.4508
*   NOAM    0.0066   -0.1928   -0.0176
*   NUBI    0.0274   -0.1704    0.2037
*   PCFC   -0.1135    0.2907   -0.6025
*   SOAM   -0.0751   -0.0835   -0.0389
*   SOMA   -0.0336   -0.2206    0.2454
     .      -0.0688D00,  -0.0900D00,   0.1874D00,
     .       0.3205D00,  -0.0378D00,   0.4011D00,
     .       0.4194D00,   0.3284D00,   0.3375D00,
     .      -0.0235D00,  -0.1476D00,   0.2140D00,
     .       0.3205D00,  -0.0014D00,   0.4038D00,
     .      -0.0925D00,  -0.4290D00,   0.4508D00,
     .       0.0066D00,  -0.1928D00,  -0.0176D00,
     .       0.0274D00,  -0.1704D00,   0.2037D00,
     .      -0.1135D00,   0.2907D00,  -0.6025D00,
     .      -0.0751D00,  -0.0835D00,  -0.0389D00,
     .      -0.0336D00,  -0.2206D00,   0.2454D00,
*
     .     0.9722D-02, -0.1839D+00, -0.2778D-01,   ! NAM08 == NOAM_I08
     .           0.d0,        0.d0,        0.d0,   ! IGS08 == ITRF2008
     .           0.d0,        0.d0,        0.d0,   ! ITRF2014
     .           0.d0,        0.d0,        0.d0,   ! IGS14 == ITRF2014
     .       0.0066D00,  -0.1928D00,  -0.0176D00,  ! NAM14 = NOAM_I14
     .      -0.0688D00,  -0.0900D00,   0.1874D00 / ! ANT14 = ANTA_I14


* Not valid values: ITRF94 is nominally NNR-NUV so itrf94 set to zero.     
*    .                   -0.0383d0, -0.0617d0,  0.0089d0 /
* MOD TAH 970713: Added ITRF94 values from IERS TN 20; values
*     converted from -0.138, -0.222, 0.032 mas/yr to deg/Myr.
* MOD TAH 030816: Added ARAB_MCC from McClusky et al [2003] below
* MOD ??? ??????: Added ARAB_M06 from Reilinger et al. [2006]
* MOD RWK 071227: Added ITRF05 values from Altamimi et al. [2007]
* MOD RWK 120106: Added ITRF08 values from Altamimi et al. [2012]
* MOD TAH 120120: Added LIST option
      if( sys_frame(1:4).eq.'LIST' ) then
          write(*,120) max_frames
 120      format('* Frame data is available for ',i3,' frames',/,
     .           '*  #   Name           wx         wy',
     .           '         wz  (deg/Myr)')
          do i = 1, max_frames
             write(*,140) i, frame_names(i), frame_data(:,i)
 140         format('* ',i3,1x,a8,1x,3(F10.6,1x))
          end do
          call getenv('HOME',home_dir)
          frame_file = home_dir(1:lenh) // 
     .               '/gg/tables/frames.dat'
          write(*,160) trim(frame_file)
 160      format('* Also check ./frames.dat and ',a,
     .           ' for additional frames')

          RETURN
      end if


****  See if we are not using a frame (set if no file name passed
*     for this routine)
      do i = 1,3
          rot_vec(i) = 0.d0
      end do
 
      if( sys_frame(1:4).eq.'NONE' ) RETURN
 
****  Scan do the list of frames available and get number for
*     sys frame
 
      nfs = 0
      nfo = 0
 
      do i = 1, max_frames
*          lenf = trimlen(frame_names(i))
           lenf = 8
          if( sys_frame(1:lenf).eq.frame_names(i)(1:lenf)) then
              nfs = i
          end if
          if( out_frame(1:lenf).eq.frame_names(i)(1:lenf)) then
              nfo = i
          end if
      end do

***** See if we found the frames.  If we have not, check to see if
*     frames.dat exists either locally on in ~/gg/tables
      if( nfs.eq.0 .or. nfo.eq.0 ) then
         frame_file = 'frames.dat'
         open(107,file=frame_file,iostat=ierr,status='old')
         if( ierr.ne.0 ) then
*            Try a version in $HOME/gg/tables
             call getenv('HOME',home_dir)
             lenh = trimlen(home_dir)
             if( lenh.gt.0 ) then
                frame_file = home_dir(1:lenh) // 
     .               '/gg/tables/frames.dat'
                open(107,file=frame_file,iostat=ierr,status='old')
             end if
         end if
         if( ierr.eq.0 ) then
            write(*,180) frame_file(1:trimlen(frame_file))
 180        format(' Reading additional frame definitions from ',a)
         end if
*        If file opened OK, see if we can find the frames we are
*        looking for
***** SCM 04/6/16 Added ifnd index to indicate when correct plate entry is found. 
         ifnd = 0
         do while ( ierr.eq.0 )
            read(107,'(a)',iostat=ierr) line
            if( ierr.eq.0 .and. line(1:1).eq.' ' .and.
     .          trimlen(line).gt.0  .and. ifnd.eq.0) then
               indx = 0
               call GetWord(line,newframe, indx)  
               call casefold(newframe)
*              See if this is one we are looking for
               if( sys_frame(1:lenf).eq.newframe(1:lenf)) nfs = -1
               if( out_frame(1:lenf).eq.newframe(1:lenf)) nfo = -1
               if( nfs.lt.0 .or. nfo.lt.0 ) then
*                  OK: Found a match to frame name.  Read the rest
*                  of the line
                   call multiread(line, indx, 'R8',jerr, frame_ev,
     .                   cdum, 3) 
*                  See if units passed
                   call GetWord(line,frame_un, indx)
                   call casefold(frame_un)
*                  Default units are deg/Myr.  Rad/Myrs may be
*                  used in file and results are converted to 
*                  deh/Myr (actually gets converted back below)
                   if( index(frame_un,'RAD').eq.1 ) then
                       do i = 1,3
                          frame_ev(i) = frame_ev(i)*180.d0/pi
                       end do
                   endif
*                  Assign the values.  Only the first occurrence will
*                  be used.
                   if( nfs.eq.-1 ) then
                       do i = 1,3
                          frame_es(i) = frame_ev(i)
                       end do
                       nfs = -2
                   end if
                   if( nfo.eq.-1 ) then
                       do i = 1,3
                          frame_eo(i) = frame_ev(i)
                       end do 
                       nfo = -2                      
                   end if
                   if( nfs.lt.0 .and. nfo.lt.0 ) ifnd = 1
               end if
            end if
         end do
         close(107)
      end if

****  Check to see if :A is in the name.  If so set the scale
*     for NUVEL_1A model.  
* MOD TAH 970713: Make the default system be NUVEL-1A. Adding
*     :O to name will yield original NUVEL-1.
*     Do not modify the ITRF values becuase these are already in
*     NUVEL-1A.
      scs = 1.d0 
      sco = 1.d0 
      if( index(sys_frame,':O').gt.0 .and. 
     .    sys_frame(1:4).ne.'ITRF' ) scs = 1.d0/0.9562d0   
      if( index(out_frame,':O').gt.0 .and.
     .    out_frame(1:4).ne.'ITRF' ) sco = 1.d0/0.9562d0    
 
****  If we can't find frame, tell user about continue (rot_vec
*     will be zero)
 
      if( nfs.eq.0 .or. nfo.eq.0 ) then
          write(*,200) sys_frame, out_frame, nfs, nfo
  200     format(' *** WARNING *** Could not find frame.  Frames',
     .           ' were ',a10, 1x,a10,'. Numbers were ',2i3,/,
     .           '                 Continuing with no frame change')
          RETURN
      end if
 
****  Get the rotation between the frames
      if( nfs.gt.0 ) then
         do i = 1,3
            frame_es(i) = frame_data(i,nfs)
         end do
      endif
      if( nfo.gt.0 ) then
         do i = 1,3
            frame_eo(i) = frame_data(i,nfo)
         end do
      endif


      do i = 1,3
          rot_vec(i) = frame_es(i)*scs*pi/180.d6 - 
     .                 frame_eo(i)*sco*pi/180.d6
      end do
 
****  Thats all
      return
      end
 
