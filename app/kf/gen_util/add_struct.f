CTITLE ADD_STRUCT
 
      subroutine add_struct( cont_structure,avail_structure,
     .           source, theoretical,source_struc_cont)

      implicit none 
 
 
*     Routine to add source structure contribution if it is
*     avaialable.  If it is not available then the cont_struture
*     bit is turned off.
 
*   avail_structure - Bit mapped word saying which sources
*           - actually have source structure corrections
*   cont_structure  - Bit mapped word saying which sources
*           - are to have structure corrections applied
*   source  - source number for this observation
 
      integer*4 avail_structure, cont_structure, source
 
*   kbit    - Bit checking function
 
      logical kbit
 
*   source_struc_cont(3)    - structure contribution for
*           - grup delay, phase delay and phase delay rate.
 
      real*4 source_struc_cont(3)
 
*   theoretical(4)  - the four theoretical values (grp,phs,
*           - Sb delays and rates
 
      real*8 theoretical(4)
 
***** See if we are to apply
 
*                                             ! Yes, we are to apply
      if( kbit(cont_structure,source) ) then
 
*         See if it is available
*                                                  ! yes it is
          if( kbit(avail_structure,source) ) then
              theoretical(1) = theoretical(1) + source_struc_cont(1)
              theoretical(2) = theoretical(2) + source_struc_cont(2)
              theoretical(3) = theoretical(3) + source_struc_cont(1)
              theoretical(4) = theoretical(4) + source_struc_cont(3)
*                                                  ! No it is not
          else
*                                                  ! reset so we wont
              call sbit(cont_structure,source,0)
*                                                  ! try again
          end if
      end if
 
***** Thats all
      return
      end
 
