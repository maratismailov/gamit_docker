      Subroutine ant_alias( antcod_in,antcod_out )

c     Convert the input antenna code into a unique internal standard
c     code for GAMIT.  See also table rcvant.dat 

c     R. King  11 January 1997; last modified by R. King  22 August 2013 

c     Antenna                                 Std code (output)  Aliases allowed (input)
c     ----------------------                   ----------------  -----------------------

c     Altus Integrated GPS + GLONASS          ALTPS3             TIA83G  (TI Asahi)
                                                                   
c     Ashtech Marine L1 only                  ASHMAR             none
c     Ashtech Marine L1/L2 ver A              ATMR2A             ASHMRA
c     Ashtech Marine L1/L2 ver B              ATMR2B             ASHMRB
c     Ashtech Marine L1/L2 ver C              ATMR2C             ASHMRC
c     Ashtech Geodetic L                      ASHL12             none
c     Ashtech Geodetic L Extended GP          ASHLEX             none
c     Ashtech Geodetic P                      ASHP12             ASHXII  ATGEOE TOPP12 (Topcon)
c     Ashtech Geodetic III                    ASHGD3             ASHGDA  ATGE3B TOPGD3  (Topcon)
c     Ashtech Geodetic III USCG (radome)      ASHG3R             none
c     Ashtech C-R 700936A  (no radome)        ASHDMG             DMGASH
c     Ashtech C-R 700936A  (radome)           ASHDMR             DMRASH
c     Ashtech C-R 700936.02 Rev E (no radome) ASHDMC             none
c     Ashtech C-R 700936.02 Rev E  (radome)   ASHDMD             none    
c     Ashtech C-R GPS/GLONASS 701073.1        ASHGG1             none
c     Ashtech C-R GPS/GLONASS 701073.3        ASHGG2             none
c     Ashtech C-R GPS/GLONASS 701941.1        ATGGD1             none
c     Ashtech C-R GPS/GLONASS 701941.2        ATGGD2             none
c     Ashtech C-R GPS/GLONASS 701941.A        ATGGDA             none
c     Ashtech C-R GPS/GLONASS 701941.B        ATGGDB             none
c     Ashtech C-R 701945-NN Rev A             ATDMRA             none
c     Ashtech C-R 701945-NN Rev B             ATDMRB             ATDMR2
c     Ashtech C-R 701945-01 Rev C             ATDM1C             none    
c     Ashtech C-R 701945-01 Rev D             ATDM1D             none
c     Ashtech C-R 701945-01 Rev E             ATDM1E             none
c     Ashtech C-R GPS/GLONASS 701946-01 Rev 2 ATDMG2             ATDMGG
c     Ashtech C-R GPS/GLONASS 701946-01 Rev32 ATDMG3             none
c     Ashtech Geodetic IV Rev A               ATGD4A             none
c     Ashtech Geodetic IV Rev A w/ GP         ATGP4A             none    
c     Ashtech ProMark 800                     AT147A             ATP800

c     JPL Dorne-Margolin R                    ROGSNR             DMRCHR  ROGDMR JPLDMR   
c     JPL DM-R retrofit with Ashtech LNA      ROGMRA             RDGMRA 
c     AOA Dorne-Margoiin B                    ROGAOA             DMBCHR  ROGDMB AOADMB
c     AOA Dorne-Margolin T                    TRBROG             DMTCHR  ROGDMT AOADMT  
c     AOAD/M_TA_NGS                           TRBMTA             TRBMRA 

c     FRPA-2 (w/ TI400)                       FRPA-2             none

c     JAVTRIANT_A                             JAVTRA             JAVTRI

c     Leica internal                          LC299I             SR299I SR399I LEIGRT LC_299
c     Leica AT201                             LC_201             AT201E
c     Leica AT202/302                         LC202N             AT202E SR299E
c     Leica AT202/302 GP                      LC202G             AT202G SR299G    
c     Leica AT302 (same as AT202)             LC302N             AT302E SR399E
c     Leica AT302 (same as AT302 GP)          LC302G             AT302G SR399G   
c     Leica AT303 (micropule choke-ring)      LC_303             AT303G AT503G LEIDMG
c     Leica AT501                             LC_501             AT501E 
c     Leica AT502                             LC_502             AT502E 
c     Leica AT503 (same as AT303)             LC_503             AT503E

c     Leica AT504G (DM-T choke-ring)          LC_504             AT504G  

c     Macrometer X-dipole                     MINXDP             MIN6AT     
c     Macrometer patch                        MINPCH but not yet coded in hisub

c     Sercel TR5S                             SRTR5S             none
c     Sercal NR52                             SRN52              none

c     Sokia 700                               SOKGD3             ASHGD3

c     TI 4100 100/2000 series                 TI_100             TI4100 
c     TI 4100 2000 series (same as 100s?)     TI2000             TI4100
c     TI 4100 4000 series                     TI4000             none
 
c     Topcon HIPER GD                         TPSHIP             TOPHIP                                              
c     Topcon Choke Ring GD                    TPSC3D             TOPC3D JNSC3D JNSZCT
c     Topcon CHoke Ring w/ SNOW Radome        TPSC3R             TOPC3R JNSC3R JNSZCR
c        These two sold by JNS as 'Zero-centered' not 'C3R'

c     Trimble 4000SX Micro                    TRMSXD             none
c     Trimble 4000SLD L1/L2                   TRMSLD             none
c     Trimble 4000SST L1/L2 Geodetic          TRMSST             none
c     Trimble TR Geodetic L1/L2 (SSE/SSi)     TRMSSE             none
c     Trimble TR Geodetic L1/L2 w/o GP        TRMKIN but not yet coded in hisub
c     Trimble choke-ring                      TRMDMG             DMGTRM 
c     Trimble Zephyr - no ground plane        TRMZEP             TRZEPH
c     Trimble Zephyr - with ground plane      TRMZGP             TRZEPG

c     Wild-Magnovox WM-102                    WM-102             none

c     Geotracer 2000                          GEO200             none
c     Geotracer 2200                          GEO220             none

c     Micropulse MP-1372W REGP (UNAVCO CMC)   MPREGP             none


      implicit none

      character*6 antcod_in,antcod_out
                          
      call uppers(antcod_in)
                         
      if( antcod_in .eq. 'TI4100' ) then
         antcod_out = 'TI_100' 

      elseif( antcod_in .eq. 'ASHXII' .or. antcod_in.eq.'ATGEOE' .or.
     .    antcod_in .eq. 'TOPP12' ) then
         antcod_out = 'ASHP12'

      elseif( antcod_in .eq. 'ASHGDA' .or. antcod_in.eq.'ATGE3B' .or. 
     .        antcod_in .eq. 'TOPGD3' .or. antcod_in.eq.'ASHGD3' ) then
         antcod_out = 'ASHGD3' 

      elseif( antcod_in .eq. 'ASHMRA' ) then
         antcod_out = 'ATMR2A'
      elseif( antcod_in .eq. 'ASHMRB' ) then
         antcod_out = 'ATMR2B'
      elseif( antcod_in .eq. 'ASHMRC' ) then
         antcod_out = 'ATMR2C' 

      elseif( antcod_in .eq. 'ATDMR2' ) then
         antcod_out = 'ATDMRB'
 
      elseif( antcod_in .eq. 'DMGASH' ) then
         antcod_out=  'ASHDMG'

      elseif( antcod_in .eq. 'DMRASH' )  then
         antcod_out=  'ASHDMR'
                                          
      elseif( antcod_in .eq. 'ATDMGG' ) then 
         antcod_out = 'ATDMG2'

      elseif( antcod_in .eq.'ATP800' ) then
         antcod_out = 'AT147A'

      elseif( antcod_in .eq. 'DMRCHR' .or. antcod_in.eq.'ROGDMR' .or.
     .        antcod_in .eq. 'JPLDMR' )  then
         antcod_out = 'ROGSNR'                                

      elseif( antcod_in.eq. 'RDGMRA' ) then
         antcod_out = 'ROGMRA'

      elseif( antcod_in .eq. 'DMBCHR' .or. antcod_in.eq.'ROGDMB' .or.
     .        antcod_in .eq. 'AOADMB' )  then
         antcod_out = 'ROGAOA'

      elseif( antcod_in .eq. 'DMTCHR' .or. antcod_in.eq.'ROGDMT' .or.
     .        antcod_in .eq. 'AOADMT' )  then
         antcod_out = 'TRBROG'     

      elseif( antcod_in .eq. 'TRBMRA' )  then
         antcod_out = 'TRBMTA'

      elseif( antcod_in .eq. 'MIN6AT' )  then
         antcod_out = 'MINXDP'

      elseif( antcod_in .eq.'TOPC3D' .or. antcod_in.eq.'JNSC3D' .or.
     .        antcod_in .eq.'JNSZCT' ) then
         antcod_out = 'TPSC3D'

      elseif( antcod_in .eq.'TOPC3R' .or. antcod_in.eq.'JNSC3R' .or.
     .        antcod_in .eq.'JNSZ3R' ) then
         antcod_out = 'TPSC3R'

      elseif( antcod_in .eq. 'DMGTRM' )  then
         antcod_out = 'TRMDMG' 

      elseif( antcod_in .eq. 'TRZEPH' ) then
         antcod_out = 'TRMZEP'

      elseif( antcod_in .eq. 'TRZEPG' ) then
         antcod_out = 'TRMZGP'     
                        
      elseif( antcod_in .eq. 'SR299I' .or. antcod_in.eq. 'LC_299' .or.
     .        antcod_in .eq. 'LEIGRT' ) then
         antcod_out = 'LC299I'

      elseif( antcod_in .eq. 'SR399I' .or. antcod_in.eq. 'LC_399' ) then
         antcod_out = 'LC399I'

      elseif( antcod_in .eq. 'SR299E' .or. antcod_in.eq. 'AT202E' ) then
         antcod_out = 'LC202N'

      elseif( antcod_in .eq. 'SR399E' .or. antcod_in.eq. 'AT302E' ) then
         antcod_out = 'LC302N'

      elseif( antcod_in .eq. 'SR299G' .or. antcod_in.eq. 'AT202G' ) then
         antcod_out = 'LC202G'   

      elseif( antcod_in .eq. 'SR399G' .or. antcod_in.eq. 'AT302G' ) then
         antcod_out = 'LC302G'

      elseif( antcod_in .eq. 'LEIDMG' .or. antcod_in.eq. 'AT503G' .or.
     .        antcod_in .eq. 'AT303G' ) then
         antcod_out = 'LC_303'
 
      elseif( antcod_in .eq. 'AT201E' ) then
              antcod_out = 'LC_201'

      elseif( antcod_in .eq. 'AT501E' ) then
         antcod_out = 'LC_501'

      elseif( antcod_in .eq. 'AT502E' ) then
         antcod_out = 'LC_502'

      elseif( antcod_in .eq. 'AT503E' ) then
         antcod_out = 'LC_503'

      elseif( antcod_in .eq. 'AT504G' ) then
         antcod_out = 'LC_504'

      elseif( antcod_in.eq.'SOKGD3' ) then
         antcod_out = 'ASHGD3'

      elseif( antcod_in.eq.'TIA83G' ) then 
         antcod_out= 'ALTPS3'
	 
      elseif( antcod_in.eq.'JAVTRI' ) then
         antcod_out= 'JAVTRA'

      else
         antcod_out = antcod_in

      endif
                   
      return
      end
