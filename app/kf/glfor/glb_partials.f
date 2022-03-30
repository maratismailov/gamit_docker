CTITLE GLB_PARTIALS
 
      subroutine glb_partials( part_pnt, a_part )

      implicit none
 
 
*     Routine to compute the partials derivatives of the current
*     solution parameters with respect to the the global parameters
*
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
 
      include '../includes/glb_hdr_def.h'
 
 
*   indx        - Index from parameter code (eg site/source number)
*   i           - Loop counter
*   part_pnt(2,max_glb_deriv,cnum_parn) - Pointer to which partials
*               - are to be used.  The values are arranged in pairs
*               - the first giving parameter number, the next giving
*               - the number of contiguous parameters to be used.
*   type        - Type of parameter
 
 
      integer*4 indx, i, part_pnt(2,max_glb_deriv,cnum_parn), type
      integer*4 act_max_deriv, j, numd
 
*   a_part(max_glb_deriv, cnum_parn) - Non-zero partials to be used.
*   amp_deriv(2,4,max_coeff_periods)   - Amplitude partials for each of the
*                periods at which we can estimate.
*   fcn_period  - Real*8 version of nut_period
 
 
      real*8 a_part(max_glb_deriv, cnum_parn), 
     .       amp_deriv(2,4,max_coeff_periods), fcn_period
 
***** loop over all the parameters and compute the partials
 
      call clear_pnt( cnum_parn, part_pnt )

***** Cpmpute the amplitude derivatives before we start 

      fcn_period = nut_period

      call gamp_parts( gepoch_expt, fcn_period, amp_deriv)

* MOD TAH 051202: Debug
C     write(*,998) guse_site
C998  format('GUSE_SITE ',(5o12))
C     write(*,999) (i,ltog_sites(i),i=1,cnum_sites)
C999  format('LTOG: ',(20i4))
      act_max_deriv = 0

      do i = 1, cnum_parn
 
          call decode_code( gpar_codes(i), type, indx )

          call compute_partials( i, type, indx, part_pnt(1,1,i),
     .                           a_part(1,i), amp_deriv )
          numd = 0
          do j = 1, indx_pnt(i)
             numd = numd + part_pnt(2,j,i)
          end do
          
          act_max_deriv = max(act_max_deriv,numd)
 
      end do
* MOD TAH 140106: Only  report if error.
      if( act_max_deriv.gt. max_glb_deriv ) then
         write(*,120) act_max_deriv, max_glb_deriv
         if( log_unit.ne.6 )
     .   write(log_unit,120) act_max_deriv, max_glb_deriv 
 120     format('** ERROR *** Used ',i3,' derivatives; max allowed ',
     .           i3)

*        Detailed debug to see what is happening:
         act_max_deriv = 0
         call clear_pnt( cnum_parn, part_pnt )

         write(*,220) cnum_parn
 220     format('* PARTIAL count for ',i5,' local parameters')

         do i = 1, cnum_parn
 
             call decode_code( gpar_codes(i), type, indx )

             call compute_partials( i, type, indx, part_pnt(1,1,i),
     .                           a_part(1,i), amp_deriv )
             numd = 0
             do j = 1, indx_pnt(i)
                numd = numd + part_pnt(2,j,i)
             end do
             act_max_deriv = max(act_max_deriv,numd)
* MOD TAH 160324: Only output for offending parameters.
             if( numd.gt.max_glb_deriv ) then
                write(*,240) i, act_max_deriv, numd, indx_pnt(i), 
     .               part_pnt(:,1:indx_pnt(i),i)
 240            format('I ',i4,1x,i3,1x,i3,1x,i4,' : ',20(1x,i6,1x,i2))
                write(*, 260) gpar_codes(i), type, indx, 
     .               ltog_sites(indx), gsite_names(ltog_sites(indx))
 260            format('PAR_CODE ',i12,1x,' TYPE ', I3,' INDX ', i5,
     .                 ' GLB Site ',i5,' Name ',a)
             endif

         end do

         stop 'Error in derivative bounds'
*        
      endif 

****  Thats all
      return
      end
 
