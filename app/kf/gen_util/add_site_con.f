CTITLE ADD_SITE_CONT
 
      subroutine add_site_cont( mask, avail, type, theoretical,
     .                          cont, in, cg,cp,cs )

      implicit none 
 
 
*     Routine to add the site dependent contributions to the
*     theoretical delays.  This routine assumes that the contribution
*     (CONT) is passed as a 2*2 matrix.  The first index runs over
*     site and the second index over delay and rate.  The contribution
*     is assumed to have the sign of the contribution for the site.
*     Hence, the contribution is added for the second site and
*     substracted for the first site.  (For the phase calibration
*     contribution this sign convension is opposite to that adopted
*     by GSFC.  The signs are changed at READIN stage.)
*
*     The user also passes three integer values which are used to
*     stop some contributions being applied to specific delays.
*     e,g, Type 2 contributions (Feed rotation corr.) are applied only
*     to the phase delay.
*
*     All contributions are assumed to be in ps for delays and fs/s
*     for rates
*                                     10:17 PM TUE., 10 Feb., 1987
*
*   avail       - Indicates if the contribution is available at
*               - this contribution is available
*   cg          - Multiplicitive factor for group delay cont.
*               - By passing zero, the contribution will not be
*               - applied to the group delay.
*   cp          - Factor for phase delays (see CG above)
*   cs          - Factor for SB delay (see CG above)
*   in          - Baseline index (either 1 or 2) used to get
*               - contribution from cont array
*   mask        - Indicates if this contribution is to be applied
*               - for this site.
*   type        - Indicates which contribution is to be applied.
*               - (Used in bit masking avial and mask)
 
      integer*4 avail, cg, cp, cs, in, mask, type
 
*   kbit        - Bit checking function
 
      logical kbit
 
*   cont(2,2)   - The contribution to be applied, arranged by
*               - site (index 1) and delay and rate(index 2)
 
      real*4 cont(2,2)
 
*   sign        - Sign of the contribution (plus for site 2
*               - (in=2), minus for site 1, (in=1))
*   theoretical(4)  - The group, phase and SB delays and rate
*               - (ps and fs/s)
 
 
      real*8 sign, theoretical(4)
 
***** See if we want to apply this contribution
 
*                                ! yes, we do want to apply
      if( kbit(mask,type) ) then
 
*         See if contribution is available
*                                     ! Yes, it is, so apply
          if( kbit(avail,type) ) then
 
*             Get sign
              sign = 2*in - 3
 
              theoretical(1) = theoretical(1) + sign*cg*cont(in,1)
              theoretical(2) = theoretical(2) + sign*cp*cont(in,1)
              theoretical(3) = theoretical(3) + sign*cs*cont(in,1)
              theoretical(4) = theoretical(4) + sign*cont(in,2)
 
*                                     ! Contribution not available.
          else
*                                     ! Set mask so it wont be tried
*                                     ! again.
              call sbit( mask,type,0 )
 
          end if
      end if
 
***** Thats all
      return
      end
 
