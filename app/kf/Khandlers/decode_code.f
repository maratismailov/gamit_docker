CTITLE DECODE_CODE
 
      subroutine decode_code ( code, type, indx )
 
      implicit none
 
*     Routine to decode the code into type and indx.  The type is
*     stored in the lower 8 bits, and the indx in the upper 8 bits.
 
*   code        - Apriori or parameter code (see
*               - GLB_HEADER_DEF.FTNI)
*   indx        - Pointer to station/source/coefficent
*   type        - Lower 8 bits with type
 
      integer*4 code, indx, type
 
*   full_code   - Full code.  Full code is made positive even when
*               - Bit 16 of code is set (as for coefficient codes)
 
      integer*4 full_code
 
***** Convert code to a positive number
      if( code.gt.0 ) then
          full_code = code
      else
          full_code = 65536 + code
      end if
 
***** Split code into the two parts
 
      indx = full_code/256
      type = full_code - indx*256
 
***** Thats all
      return
      end
 
 
