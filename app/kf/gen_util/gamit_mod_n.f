CTITLE GAMIT_MOD_NAME

      character*8 function gamit_mod_name(mn)

      implicit none 

*     Function to return the name of gamit model number.
* MOD TAH 110512: Added IERS2010 Mean Pole tide model (MPT2010)
*     Bit 23 of gamit_mod
* MOD TAH 150213: Replaced OceanPT with OPTD96 (ocean pole tide applied to
*     IERS96 mean pole).
*     Bit 24 of gamit_mod: OPTD96
*     Bit 25 of gamit_mod: OceanPT ocean pole tide applied to IERS2010 mean pole.
* MOD TAH 180205: Started to use the end of the GAMIT Model bits to show
*     htoglb covaraince matrix changes.  Bit 14/15 used for RCOV and TCOV to
*     show rotation and translation covariance added.  (Can be used to remove
*     later if desired--not implemented as of 180205.
* MOD TAH 200220: Updated codes to reflect the new IERS2020 secular pole
*     location. MPT2020 -- Mean pole IERS2020; OPTD20 ocean pole tide applied
*     to the same mean pole.
* MOD TAH 200505: Added the entries for UT1-Libration and the Desai and Sibios 
*     diurnal/semidiurnal ocean tide model.

* mn  -- Model number

      integer*4 mn

* LOCAL Variables
* mod_names(32)  -- Names of the models 

      character*8 mod_names(32)
                              

      data mod_names  / 'SD-WOB  ','SD-UT1  ','RAY-MOD ',
     .                  'IERS10  ','UT1-LIBR','GIPSON17',
     .                  'D&S2016 ','UNKNOWN ','UNKNOWN ',
     .                  'UNKNOWN ','UNKNOWN ','UNKNOWN ',
     .                  'UNKNOWN ','UNKNOWN ','RotCov  ',
     .                  'TranCov ','E-Tide  ','K1-Tide ',
     .                  'PTide   ','OC-Load ','IERS96  ',
     .                  'AtmTide ','IERS10  ','OPTD96  ',
     .                  'OPTD10  ','IERS20  ','OPTD20  ',
     .                  'UNKNOWN ','UNKNOWN ','UNKNOWN ',
     .                  'UNKNOWN ','UNKNOWN '             /

****  Assign the result
      gamit_mod_name = mod_names(mn)

****  Thats all
      return
      end


