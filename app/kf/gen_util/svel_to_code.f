      subroutine svel_to_code( svel, gkel, arel, code, direct)

      implicit none 

*     Routine to map and unmap the Satellite radiation pressure
*     svel names to code number and visa-versa.  Direct
*     set which way the convertion is done by giving NTOC 
*     (name to code), CTON (code to name), GTOC (Globk name to
*     code).  ATOC gives the arc to code conversion.
*     The other character strings are returned when
*     ever a code number is given or returned.

* MOD TAH 190610: Re-order the radition parameters to be fit
*     the ECOMC model (replaced terms from ECOM2 that we never
*     used).

      integer*4 num_svel   ! Sets the number of SV type codes.
* MOD TAH 981020: Increased number to 17 to account new Berne
*     model parameters and satellite antenna offsets. 
*     NOTE:  The antenna offsets are not used in arc.

      parameter ( num_svel = 17 )

* PASSED VARIABLES

      integer*4 code     ! Numeric code for radiation pressure
                         ! term (value offset by 6 for other 
                         ! orbital svels)

* svel,  - Name of svel from hfile
* gkel   - Name of element from globk solution
* arel   - Name of element for gfile
* direct    - direction for conversion.  Options are:
*             NTOC -- H-file name to code
*             GTOC -- Globk solution name to code
*             ATOC -- Arc t-file name to code
*             CTON -- Code number to each of the names.

      character*(*) svel, gkel, arel, direct  

* LOCAL VARIABLES

* i,        - Loop counter
* len_el    - Length of passed svel.
      integer*4 i,  len_el  

      character*4 dir_loc  ! Local name that can be casefolded

* svel_names  - Radiation parameter names in ascii h-files
* gkel_names  - Names used in globk output
* arel_names  - Names used in arc gfiles

* MOD TAH 190610: Added 6 names of solve strings to allow for
*     rename of once-per-rev terms (mapped to new names here).
      character*16 svel_names(num_svel+6)
      character*14 gkel_names(num_svel)
      character*4  arel_names(num_svel)

C     data  svel_names / 'RAD PRES DIRECT',
C    .                   'Y AXIS BIAS    ',
C    .                   'Z AXIS BIAS    ',
C    .                   'B AXIS BIAS    ',
C    .                   'X AXIS BIAS    ',
C    .                   'COS DIRECT     ',
C    .                   'SIN DIRECT     ',
C    .                   'COS Y BIAS     ',
C    .                   'SIN Y BIAS     ',
C    .                   'COS B BIAS     ',
C    .                   'SIN B BIAS     ',
C    .                   'SIN  U X AXIS  ',
C    .                   'SIN 3U X AXIS  ',
C    .                   'SIN  U Z AXIS  ',
C    .                   'SVANT X AXIS   ',
C    .                   'SVANT Y AXIS   ',
C    .                   'SVANT Z AXIS   '  /
C
* MOD TAH 190610: New names from solve      
      data  svel_names / 'RAD PRES DIRECT',
     .                   'Y AXIS BIAS    ',
     .                   'B AXIS BIAS    ',
     .                   'COS U DIRECT   ',
     .                   'SIN U DIRECT   ',
     .                   'COS U Y AXIS   ',
     .                   'SIN U Y AXIS   ',
     .                   'COS U B AXIS   ',
     .                   'SIN U B AXIS   ',
     .                   'COS 2U DIRECT  ',
     .                   'SIN 2U DIRECT  ',
     .                   'COS 4U DIRECT  ',
     .                   'SIN 4U DIRECT  ',
     .                   'UNUSED         ',
     .                   'SVANT X AXIS   ',
     .                   'SVANT Y AXIS   ',
     .                   'SVANT Z AXIS   ',  
* 6 old names
     .                   'COS DIRECT     ',
     .                   'SIN DIRECT     ',
     .                   'COS Y BIAS     ',
     .                   'SIN Y BIAS     ',
     .                   'COS B BIAS     ',
     .                   'SIN B BIAS     ' /

* OLD GLOBK names (non-used values are removed). 
* Strictly requires new globk common version;
* implemememt with warning because only orbits will be
* affected
C     data gkel_names /  'Direct Rad    '
C    .,                  'Y Axis Bias   '
C    .,                  'Z Axis Bias   '
C    .,                  'B Axis Bias   '
C    .,                  'X Axis Bias   '
C    .,                  'Cos Direct    '
C    .,                  'Sin Direct    '
C    .,                  'Cos Y Bias    '
C    .,                  'Sin Y Bias    '
C    .,                  'Cos B Bias    '
C    .,                  'Sin B Bias    '
C    .,                  'Sin U X Bias  '
C    .,                  'Sin 3UX Bias  '
C    .,                  'Sin U Z Bias  '
C    .,                  'SVAnt X Off   '
C    .,                  'SVAnt Y Off   '
C    .,                  'SVAnt Z Off   ' /
* MOD TAH 190610: Add 2U and 4U names and removed non-used values.
*     These names need to copied to write_glb/loc_parn.f routines as
*     well.
      data gkel_names /  'Direct Rad    '
     .,                  'Y Axis Bias   '
     .,                  'B Axis Bias   '
     .,                  'Cos Direct    '
     .,                  'Sin Direct    '
     .,                  'Cos Y Bias    '
     .,                  'Sin Y Bias    '
     .,                  'Cos B Bias    '
     .,                  'Sin B Bias    '
     .,                  'Cos 2U Direct '
     .,                  'Sin 2U Direct '
     .,                  'Cos 4U Direct '
     .,                  'Sin 4U Direct '
     .,                  'Not Used      '
     .,                  'SVAnt X Off   '
     .,                  'SVAnt Y Off   '
     .,                  'SVAnt Z Off   ' /

* OLD ARC codes
C     data arel_names /  'DRAD'
C    .,                  'YRAD'
C    .,                  'ZRAD'
C    .,                  'BRAD'
C    .,                  'XRAD'
C    .,                  'DCOS'
C    .,                  'DSIN'
C    .,                  'YCOS'
C    .,                  'YSIN'
C    .,                  'BCOS'
C    .,                  'BSIN' 
C    .,                  'X1SN'
C    .,                  'X3SN'
C    .,                  'Z1SN'
C    .,                  'XOFF'
C    .,                  'YOFF'
C    .,                  'ZOFF' /
* MOD TAH 190610: Added new arc codes and removed non-used ones.
      data arel_names /  'DRAD'
     .,                  'YRAD'
     .,                  'BRAD'
     .,                  'DCOS'
     .,                  'DSIN'
     .,                  'YCOS'
     .,                  'YSIN'
     .,                  'BCOS'
     .,                  'BSIN' 
     .,                  'D2CS'
     .,                  'D2SN'
     .,                  'D4CS'
     .,                  'D4SN'
     .,                  '????'
     .,                  'XOFF'
     .,                  'YOFF'
     .,                  'ZOFF' /

****  Check which direction we are going
      dir_loc = direct
      call casefold(dir_loc)

      if( dir_loc.eq.'NTOC' ) then

*         Default to name to code mode
          len_el = min(len(svel),len(svel_names(1)))
          code = 0
* MOD TAH 190610: Also search for old names and adjust code
*         as needed (6 old names in solve h-files).
          do i = 1, num_svel+6
             if( svel(1:len_el).eq.svel_names(i)(1:len_el) ) then
                code = i + 6  ! Add in offset set 6 IC elements.
             end if
*            See if old code found and reduce the number accordingly
             if( code-6.gt.num_svel ) then ! Map back i.e.,
                     ! code-6 == 18 -> code-6 = 4
                     ! code-6 == 19 -> code-6 = 5 ...
                     ! ....
                     ! code-6 == 23 -> code-6 = 9  
                 code = code - 14
             end if
          end do
      else if( dir_loc.eq.'GTOC' ) then

*         Default to name to code mode
          len_el = min(len(gkel),len(gkel_names(1)))
          code = 0
          do i = 1, num_svel
             if( gkel(1:len_el).eq.gkel_names(i)(1:len_el) ) then
                code = i + 6
             end if
          end do
      else if ( dir_loc.eq.'ATOC' ) then
*         Arc name to code number      
          len_el = min(len(arel),len(arel_names(1)))
          code = 0
          do i = 1, num_svel
             if( arel(1:len_el).eq.arel_names(i)(1:len_el) ) then
                code = i + 6
             end if
          end do
      end if

****  Assign the names
      if( code.ge.7 .and. code.le.num_svel+6 ) then
          svel = svel_names(code-6)
          gkel = gkel_names(code-6)
          arel = arel_names(code-6)
      else if ( direct.eq.'CTON' ) then
          write(*,220) code, svel, gkel, arel, num_svel+6
 220      format('** ERROR ** svel_to_code: Code number ',i4,
     .           ' out of range. Elements ',a6,1x,a14,1x,a4,
     .           ' Valid Values 7 to ',i2)
      end if

****  Thats all
      return
      end
