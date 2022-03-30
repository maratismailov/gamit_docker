CTITLE GLB_FOR_MAP
 
      subroutine glb_for_map( cnp, vma_i8wrd )
 
      implicit none 

 
*     This routine will set up the mapping the foward global solution.
*     Some parts of the mapping depend specifically on the size
*     of the solution which is about to be read.

* MOD TAH: 940615 to allow for true dynamic mapping.  The full size
*     is only set by the most local parameters and the routine
*     is called once at the beginning of the run to get the full
*     size and then subsequently to set for each day. 
 
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

* MOD TAH 190520: Mod to allow more than 32767x32767 matrices
      integer*8 I8   ! Needed for large numbers of parameters

      data I8 / 1 /

***** Add up the space for the matrices
 
      icov_parm  = istart_vma
      isol_parm  = icov_parm  + 2*(I8*num_glb_parn)*num_glb_parn
 
      irow_copy  = isol_parm  + 2*num_glb_parn
 
      ipart_pnt  = irow_copy  + 2*num_glb_parn
      ia_part    = ipart_pnt  + 2*cnp*max_glb_deriv
 
      itemp_gain = ia_part    + 2*cnp*max_glb_deriv
      ikgain     = itemp_gain + 2*(I8*cnp)*num_glb_parn
 
      icov_obs   = ikgain     + 2*(I8*cnp)*num_glb_parn
      isol_obs   = icov_obs   + 2*(I8*cnp)*cnp
 
*     Get size of VMA.  Add 128 to allow for maxiumum
*     amount of record at the end of sol_obs when it is read
      vma_i8wrd = (isol_obs+2*cnp-istart_vma+128) 

C     write(*,100) istart_vma, num_glb_parn, max_glb_deriv, cnp,
C    .             vma_i8wrd

 100  format('FOR_MAP: istart_vma ',i18, ' num_glb_parn ',i6,
     .       ' max_glb_deriv ',i6,' Most cparn ',i6,/,
     .       '         I*4 words needed ',i18 )
 
***** Thats all
      return
      end

 
