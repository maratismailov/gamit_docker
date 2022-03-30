CTITLE SVCODE_to_CODE

      integer*4 function svcode_to_code( old_code )

      implicit none

*     Function to convert satellite radiation parameter
*     codes from GLX file 106 version 107 (cglb_vers)
*     This routine should only be called if cglb_vers
*     is <107 and code type 51 (satellite orbit parameters).
*     Assumption is all <107 codes are ECOM only so we 
*     clear bits 9-20 and then set bits 9-15.
*     First 8-bits are type and bits 9-31 are mapped 
*     to 1-23 below (mutliplying by 256 will shift the
*     bits to the correct position).
* 
      include '../includes/kalman_param.h'
      include '../includes/glb_hdr_def.h'

* Mapping routine
*     Element codes ares:
*           NEW (ECOMC)         (OLD)                    Map
*      1  - X POS             - X POS                    1
*      2  - Y POS             - Y POS                    2
*      3  - Z POS             - Z POS                    3
*      4  - X VEL             - X VEL                    4
*      5  - Y VEL             - Y VEL                    5
*      6  - Z VEL             - Z VEL                    6
*      7  - RAD PRES DIRECT   - RAD PRES DIRECT          7
*      8  - Y AXIS BIAS       - Y AXIS BIAS              8
*      9  - B AXIS BIAS       - Z AXIS BIAS (never used) 0
*     10  - COS DIRECT        - B AXIS BIAS              9
*     11  - SIN DIRECT        - X AXIS BIAS (never used) 0
*     12  - COS Y BIAS        - COS DIRECT               10
*     13  - SIN Y BIAS        - SIN DIRECT               11
*     14  - COS B BIAS        - COS Y BIAS               12
*     15  - SIN B BIAS        - SIN Y BIAS               13
*     16  - COS 2U DIRECT     - COS B BIAS               14
*     17  - SIN 2U DIRECT     - SIN B BIAS               15
*     18  - COS 4U DIRECT     - SIN X1 BIAS (never used) 0
*     19  - SIN 4U DIRECT     - SIN X3 BIAS (never used) 0
*     20  - <OPEN>            - SIN Z1 BIAS (never used) 0
*     21  - SVANT X OFF       - SVANT X OFF              21
*     22  - SVANT Y OFF       - SVANT Y OFF              22
*     23  - SVANT Z OFF       - SVANT Z OFF              23
*
* PASSED VARIABLE
      integer*4 old_code   ! parameter or apriori code from
                           ! binary hfile. 

* LOCAL VARIABLES
      integer*4 code_map(max_svs_elem)  ! Mapps pre-107 GLX version
                           ! satellite codes to 107 and later
      integer*4 type,   ! Parameter type (51)
     .          indx,   ! decode with orbital element and SVnum
     .          old_orb_el,   ! Old orbital element
     .          orb_el,       ! New orbital element
     .          sv_num        ! Satellite number
 
      character*128 message ! Error message if code appears wrong.

      data code_map /  1,  2,  3,  4,  5,  6,  7,  8,  0,  9,  0, 
     .            10, 11, 12, 13, 14, 15,  0,  0,  0, 21, 22, 23 /



****  Initially just copy the code in case there is an error later
      svcode_to_code = old_code

*     Check the cglb_vers
      if( cglb_vers .ge. 107 ) then 
*         Skip printing warnings because this code used to check all
*         codes.
!         print *,'**WARNING** svcode_to_code GLX version ',cglb_vers
          RETURN
      endif 

*     Slit the code apart: Divides into type and then orbital element
*     and satellite number
      call decode_code( old_code, type, indx )
      call decode_code( indx,     old_orb_el, sv_num )

*     Check correct parameter type (51) and that orbital element
*     number mapps successfully
      if( old_orb_el.le.max_svs_elem .and. old_orb_el.gt.0 .and.
     .    type.eq.51  ) then
          orb_el = code_map(old_orb_el)
          if( orb_el.eq.0 ) then
*             Mapping is off (orbital element should not have been used)
!             write(*,120) old_code, cglb_vers, type,old_orb_el, sv_num 
 120          format('**WARNING** Orbital element mapping problem: ',
     .               ' CODE (Oct) ',O10,' CGLB_VERS ',I5,' TYPE ',i3, 
     .               ' ORB_EL ',i3,' SV_NUM ',I4)
* MOD TAH 191201: Encode the new value (orb_el == 0 will be detected)
              svcode_to_code = type + orb_el*256 + sv_num*65536
           else
              svcode_to_code = type + orb_el*256 + sv_num*65536
           endif
      else
*         Something is not correct! Report what we know
!         write(*,120) old_code, cglb_vers, type,old_orb_el, sv_num 
      endif
          
*     Thats all
      return
      end




