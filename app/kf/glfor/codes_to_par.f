CTITLE CODES_TO_PARN
 
      subroutine codes_to_parn( type, codes, numc, parn, dimp )

      implicit none  
 
*     This routine will check the parmeter codes given in CODES
*     and save the parameter number they correspond to in the
*     PARN array (in its index postion), if the parameter type matches TYPE.
*
*     The user should ensure that PARN is cleared before this routine
*     is called.
*
 
*   codes(1)    - Parameter codes to be checked
*   dimp        - First dimension of the Parn array
*   i           - Loop counter
*   numc        - number of codes to check
*   type        - TYPE of parameter to match
*   type_test   - TYPE for the code being tested
*   indx        - Index from the code being tested
*   parn(dimp,1)    - Parameter number array to save results in
 
 
      integer*4 codes(1), dimp, i, numc, type, type_test, indx,
     .    parn(dimp,1)
 
***** Loop over the codes, and try to match TYPE
 
      do i = 1, numc
          call decode_code( codes(i), type_test, indx )
 
*                                          ! Save parameter number
          if( type_test.eq.type ) then
              parn(1,indx) = i
          end if
      end do
 
***** Thats all
      return
      end
 
 
