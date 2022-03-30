*-----------------------------------------------------------------------
* svinf_def.h
*     Satellite information records
      integer*4 svi_prn   ! PRN number
     .,         svi_svn   ! Satellite number (unique)
     .,         svi_block ! Satellite block
* MOD TAH 190627: Used former dummy variable (svi_dum) to save local SV number
*               In BACK glsave runs this value will be zero if satellite not used.
     .,         svi_lnum  ! Local satellite number (if -1 not used
                          ! in this solution.  Only used for GLBAK created GLX files.

      real*8 svi_antpos(3,2)  ! XYZ position of phase center relative
                            ! to center of mass (two sets of values)
      real*8 svi_launch     ! JD for laucch of satellite

      character*16 svi_antmod  ! Model used of phase center
      character*4 svi_ocode    ! Codes for offsets and types of models used
                               ! 1 -- 1,2,5 for L1, L2, LC first entry
                               ! 2 -- 1,2,5 for L1, L2, LC second entry
                               ! 3 -- PCV Type R or A
                               ! 4 -- PCV application E (elevation) F full

*     Total I*4 word used 14+1+6+2
      integer*4 svi_fill(9)   ! Available extra space.  This makes record
* 32 I*4 words long and therefore 4 of these records in each binary file 
* record.

      common / svidef_rec / svi_prn,  svi_svn,  svi_block,  svi_lnum  
     .,         svi_antpos, svi_launch, svi_antmod, svi_ocode  
     .,         svi_fill

*-----------------------------------------------------------------------


