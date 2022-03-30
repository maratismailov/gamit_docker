*-----------------------------------------------------------------------
*     Added block for satellite information: Entries taken from htoglb_comm.h
*     and so can't be mixed.
*     Created TAH 131111 to allow antenna offsets to be read into globk
*
      integer*4 qsvi_prn(max_glb_svs)    ! PRNs of saved entries
     .,         qsvi_svn(max_glb_svs)    ! SV number
     .,         qsvi_block(max_glb_svs)  ! Block number

      real*8 qsvi_antpos(3,2,max_glb_svs)  ! Antenna positions read from input files
      real*8 qsvi_launch(max_glb_svs)      ! Launch date for satellite
      character*16 qsvi_antmod(max_glb_svs)  ! Antenna model for each satellite (SINEX
                                         ! only need 10 charcaters in name)
      character*4 qsvi_ocode(max_glb_svs)  ! Codes for satellites.
                                         ! 1/2 Frequencues; A/R and E/F
 
      common / qsvinf_com / qsvi_prn, qsvi_svn, qsvi_block
     .,      qsvi_antpos, qsvi_launch, qsvi_antmod, qsvi_ocode
*-----------------------------------------------------------------------
