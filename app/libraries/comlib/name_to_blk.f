CTITLE NAME_TO_BLK

      subroutine name_to_blk(dir, antbody, code)

      implicit none

*     Convert to satellite antenna/boby string to numerical code
*     Names are given below
*     dir set direction: +1 from name to code; -1 from code to name

      character*(*) antbody
      integer*4 dir         ! Direction: +ve name to (numeric) code
      integer*4 code        ! Numeric code
 
      integer*4 num_abtypes
* MOD TAH 180331: Extendded for GNSS, added BLOCK III
      parameter ( num_abtypes = 45 )
      character*22 antbody_types(num_abtypes)

* MOD TAH 191219: Mod to handle the GLONASS-M+ designation which
*     only appears in someplaces.  Solution is to make the M+ be M.
* MOD TAH 200602: M+, K1A and K1B body types now in svnav.dat igs_metadata.snx
*        so this mod not longer need but we need to add add new 
*        block types so that code can be mapped back to body type.  
*        Impact in gamit/arc/earthradTUM.f which is coded with 
*        numeric body types.
C     character*22 mod_antbody  ! antbody with M+ -> M

* MOD TAH 190602: Changes BLOCK III to BLOCK IIIA

* MOD TAH 200302: Added BEIDOU-3M-CAST and BEIDOU-3M-SECM types 
* MOD TAH 200603: GLONASS changes for M->M+, K1->K1A/B:
*     Old     New
*     M   12  12
*     M+  --  13
*     K1  13  14
*     K1A --  15
*     K1B --  16
* MOD TAH 200704: BEIDOU-3M-SECM-A and BEIDOU-3M-SECM-B body types
*     replaces the BEIDOU-3M-SECM designation. (Not codes in 
*     earthradTUM.f yet) 
* MOD MAF 20210113: BEIDOU-3I body type added (not coded in earthradTUM.f yet),
*                   replacing dummy placeholder "BEI NEXT 07" (ID 37)
*     
      integer*4 j

      data antbody_types / 'BLOCK I               '   !  1
     .,                    'BLOCK II              '   !  2
     .,                    'BLOCK IIA             '   !  3
     .,                    'BLOCK IIR-A           '   !  4
     .,                    'BLOCK IIR-B           '   !  5
     .,                    'BLOCK IIR-M           '   !  6
     .,                    'BLOCK IIF             '   !  7
     .,                    'BLOCK IIIA            '   !  8
     .,                    'GPS NEXT 09           '   !  9
     .,                    'GPS NEXT 10           '   ! 10
     .,                    'GLONASS               '   ! 11
     .,                    'GLONASS-M             '   ! 12
     .,                    'GLONASS-M+            '   ! 13  --
     .,                    'GLONASS-K1            '   ! 14  13
     .,                    'GLONASS-K1A           '   ! 15
     .,                    'GLONASS-K1B           '   ! 16
     .,                    'GLO NEXT 07           '   ! 17
     .,                    'GLO NEXT 08           '   ! 18
     .,                    'GLO NEXT 09           '   ! 19
     .,                    'GLO NEXT 10           '   ! 20
     .,                    'GALILEO-0A            '   ! 21
     .,                    'GALILEO-0B            '   ! 22
     .,                    'GALILEO-1             '   ! 23
     .,                    'GALILEO-2             '   ! 24
     .,                    'GAL NEXT 05           '   ! 25
     .,                    'GAL NEXT 06           '   ! 26
     .,                    'GAL NEXT 07           '   ! 27
     .,                    'GAL NEXT 08           '   ! 28
     .,                    'GAL NEXT 09           '   ! 29
     .,                    'GAL NEXT 10           '   ! 30
     .,                    'BEIDOU-2G             '   ! 31
     .,                    'BEIDOU-2I             '   ! 32
     .,                    'BEIDOU-2M             '   ! 33
     .,                    'BEIDOU-3M-CAST        '   ! 34 
     .,                    'BEIDOU-3M-SECM-A      '   ! 35 
     .,                    'BEIDOU-3M-SECM-B      '   ! 36 
     .,                    'BEIDOU-3I             '   ! 37 
     .,                    'BEI NEXT 08           '   ! 38 
     .,                    'BEI NEXT 09           '   ! 39 
     .,                    'BEI NEXT 10           '   ! 40 
     .,                    'EXTRA-41              '   ! 41
     .,                    'EXTRA-42              '   ! 42
     .,                    'EXTRA-43              '   ! 43
     .,                    'EXTRA-44              '   ! 44
     .,                    'EXTRA-45              ' / ! 45

* 
*     See what type we have
      if( dir.gt.0 ) then
* MOD TAH 191219: Make copy of antbody and replace M+ with M
* MOD TAH 200602: M+, K1A and K1B body types now in svnav.dat igs_metadata.snx
*        so this mod not longer need/
C        mod_antbody = antbody
C        if( antbody(1:10).eq.'GLONASS-M+'  )  antbody = 'GLONASS-M'
C        if( antbody(1:11).eq.'GLONASS-K1A' )  antbody = 'GLONASS-K1'
C        if( antbody(1:11).eq.'GLONASS-K1B' )  antbody = 'GLONASS-K1'
         code = 0
         do j = 1, num_abtypes
            if( antbody.eq.antbody_types(j) ) then
                code = j
                exit
            endif
         end do
       else
         antbody = antbody_types(code)
       endif

****  Thats all
      return
      end

 
