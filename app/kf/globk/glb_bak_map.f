CTITLE GLB_BAK_MAP
 
      subroutine glb_bak_map( cnp, vma_i8wrd )
 
      implicit none 

 
*     This routine will set up the mapping the back global solution.
*     Some parts of the mapping depend specifically on the size
*     of the solution which is about to be read.
 
      include '../includes/kalman_param.h'
      include '../includes/globk_cntl.h'
      include '../includes/glb_hdr_def.h'

* cnp    - Number of local parameters that will be used in the
*          in the biggest experiment read (found in the glinit
*          run)
* vma_i8wrd  - Number of i*4 words of memory need.
* MOD TAH 150524: Introduced integer*8 version for large 
*     solutions >20000 parameters

      integer*4 cnp
      integer*8 vma_i8wrd

 
*   blocks      - Number of blocks needed for the solution
*               - matrix to be read from disk
*   cov_obs_space   - Space needed for a copy of cov_obs if the are
*                   - computing postfit residuals
*   gain_size       - Amount of space we should allow for kgain if
*                   - we need cov_obs copy as well
*   gain_space      - Space neeed for gain matrices
 
 
      integer*4 blocks

* MOD TAH 190511: Introduced when number of parameters > 32767.
      integer*8 cov_obs_space, gain_size, gain_space

* MOD TAH 190511: Introduced when number of parameters > 32767.
      integer*8 i8  !  i8 value of 1 to force I8 calculations
      data I8 / 1 /

***** Add up the space for the matrices
      icov_sav_parm = istart_vma
 
      isol_sav_parm = icov_sav_parm + 2*(i8*num_glb_parn)*num_glb_parn
 
      blocks     = (isol_sav_parm + 2*num_glb_parn 
     .                        + 1 - istart_vma)/128 + 1
 
      ijmat      = 128*blocks + istart_vma
 
      icov_parm  = ijmat      + 2*(i8*num_glb_parn)*num_glb_parn
 
      isol_parm  = icov_parm  + 2*(i8*num_glb_parn)*num_glb_parn
 
      irow_copy  = isol_parm  + 2*num_glb_parn
 
      ipart_pnt  = irow_copy  + 2*num_glb_parn
      ia_part    = ipart_pnt  + 2*cnp*max_glb_deriv
 
      itemp_gain = ia_part    + 2*cnp*max_glb_deriv
      ikgain     = itemp_gain + 2*(I8*cnp)*num_glb_parn
 
*     If we are computing postfut residuals then allow space for a copy
*     of cov_obs
 
      gain_size = 2*(I8*cnp)*num_glb_parn
 
      if( compute_glb_res ) then
          cov_obs_space = 2*(i8*cnp)*cnp
*                                                   ! Space for Kgain and
          gain_space    = 4*(I8*cnp)*num_glb_parn
*                                                   ! temp_gain
*                                                   ! Increase gain size to
          if( cov_obs_space.gt.gain_space ) then
*                                                   ! allow for cov_obs
              gain_size = cov_obs_space - 2*(I8*cnp)*num_glb_parn
          end if
*                     ! We need to allow space for cov_obs
      end if
 
      icov_obs   = ikgain     + 2*gain_size
      isol_obs   = icov_obs   + 2*(I8*cnp)*cnp

*     Get size of VMA.  Add 128 to allow for maxiumum
*     amount of record at the end of sol_obs when it is read
      vma_i8wrd = (isol_obs+2*cnp-istart_vma+128)

C     write(*,100) istart_vma, num_glb_parn, max_glb_deriv, cnp,
C    .             vma_i8wrd

C100  format('BAK_MAP: istart_vma ',i12, ' num_glb_parn ',i6,
C    .       ' max_glb_deriv ',i6,' Most cparn ',i6,/,
C    .       '         i*4 words needed ',i10 )              
 
***** Thats all
      return
      end
 
