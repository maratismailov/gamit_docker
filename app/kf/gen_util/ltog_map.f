CTITLE LTOG_MAP
 
 
      integer*4 function ltog_map ( ltog, in )

      implicit none 
 
 
*     Routine to convert index IN in ltog to site or source number.
*     If the first element of ltog is -1 then, ltog_map returns IN
*     Otherwize it returns ltog(in)
 
*   in      - Index in ltog to be checked (normally a local site
*           - or source number)
*   ltog(1) - Array which maps local to global site or source number
 
 
      integer*4 in, ltog(1)
 
****  Check to see if we should map
 
*                                     ! No mapping needed
      if( ltog(1).eq.-1 ) then
          ltog_map = in
*                                     ! use mapping from array
      else
          ltog_map = ltog(in)
      end if
 
****  Thats all
      return
      end
 
 
