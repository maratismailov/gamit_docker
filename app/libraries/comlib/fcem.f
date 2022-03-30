 
      integer*2 function fcih ( ieee, ihp )
 
 
      implicit none

*     Routine to convert IEEE to HP floating point.
****  NOT IMPLEMENTED ****  -1 error returned
 
 
      real*4 ieee, ihp
 
      fcih = -1
 
      return
      end
 
 
      integer*2 function efcidh ( ieee, ihp )
 
 
*     Routine to convert IEEE to HP floating point
****  NOT IMPLEMENTED ****  -1 error returned
 
 
      real*4 ieee, ihp
 
      efcidh = -1
 
      return
      end
 
      integer*2 function dfcih ( ieee, ihp )
 
 
*     Routine to convert IEEE to HP floating point
****  NOT IMPLEMENTED ****  -1 error returned
 
 
      real*4 ieee, ihp
 
      dfcih  = -1
 
      return
      end
 
      integer*2 function fchi ( ihp, ieee )
 
 
*     Conversion from HP floating real*4 to IEEE real*4 floating
*     point.
 
*         ihp    - HP floating point (a integer for manipuation)
*         mantissa   - Mantissa of floating point
*         exponent   - Exponent of floating point
*         sign_man   - Sign of mantissa
*         sign_exp   - sign of exponent
*         mask_man   - Mask for mantissa
*         mask_exp   - Mask for exponent
*         mask_exps  - Mask for exponent sign
*         mask_mans  - Mask for mantissa sign
 
      integer*4 ihp, mantissa, exponent, sign_man, sign_exp, mask_man,
     .    mask_exp, mask_exps, mask_mans, cand
 
*         ieee   - Should the IEEE floating point value corresponding
*                - to ihp
 
      real*4 ieee
 
*                       12345678901
      data mask_man  / O'17777777400' /
     .,    mask_exp  / O'00000000376' /
     .,    mask_exps / O'00000000001' /
     .,    mask_mans / O'20000000000' /
 
***** Get the mantissa
      mantissa = cand(ihp, mask_man )/128
      exponent = cand(ihp, mask_exp )/2
      sign_exp = cand(ihp, mask_exps )
      if( sign_exp.eq.0 ) then
          sign_exp = 1
      else
          sign_exp = -1
          exponent = 128 - exponent
      end if
 
      sign_man = cand(ihp, mask_mans )
      if( sign_man.eq.0 ) then
          sign_man = 1
      else
          sign_man = -1
          mantissa = 16777216 - mantissa
      end if
 
***** Generate the floating point value
      ieee = sign_man*(mantissa/16777216.d0)*(2.d0)**(sign_exp*exponent)
      fchi = 0
 
***** Thats all
      return
      end
 
 
 
      integer*2 function dfchi ( ihp, ieee )
 
 
*     Conversion from HP floating real*8 to IEEE real*8 floating
*     point.
 
*         ihp(2) - HP floating point (a integer for manipuation)
*         mantissa(2)- Mantissa of floating point in two parts
*         exponent   - Exponent of floating point
*         sign_man   - Sign of mantissa
*         sign_exp   - sign of exponent
*         mask_man(2)- Mask for mantissa in two parts
*         mask_exp   - Mask for exponent
*         mask_exps  - Mask for exponent sign
*         mask_mans  - Mask for mantissa sign
*         man2       - Tempory storage
 
      integer*4 ihp(2), mantissa(2), exponent, sign_man, sign_exp,
     .    mask_man(2), mask_exp, mask_exps, mask_mans, man2, cand
 
*         ieee   - Should the IEEE floating point value corresponding
*                - to ihp
*         manf   - Floating point version of the mantissa
*         bias   - 2**56
 
      real*8 ieee, manf, bias
 
*                       12345678901
      data mask_man  / O'17777777777' , O'37777777400' /
     .,    mask_exp  / O'00000000376' /
     .,    mask_exps / O'00000000001' /
     .,    mask_mans / O'20000000000' /
 
***** Get the mantissa
*     bias = 2.d0**56
      bias = 0.72057594037927940d17/2
 
*     write(*,100) bias - 2.d0**55
* 100 format(' Bias = ',d24.17)
      mantissa(1) = cand(ihp(1), mask_man(1) )
      mantissa(2) = cand(ihp(2), mask_man(2) )/256
      man2 = mantissa(2)
      if( mantissa(2).lt.0 ) then
          mantissa(2) = 16777216 + mantissa(2)
      end if
*     manf = mantissa(1)*4294967296.d0 + mantissa(2)
c     manf = mantissa(1)*33554432.d0   + mantissa(2)
      manf = mantissa(1)*33554432.d0/2   + mantissa(2)
*     write(*,110) mantissa, man2
* 110 format(3o16)
*     write(*,120) mantissa, man2
  120 format(3d25.16)
      exponent = cand(ihp(2), mask_exp )/2
      sign_exp = cand(ihp(2), mask_exps )
      if( sign_exp.eq.0 ) then
          sign_exp = 1
      else
          sign_exp = -1
          exponent = 128 - exponent
      end if
 
      sign_man = cand(ihp(1), mask_mans )
      if( sign_man.eq.0 ) then
          sign_man = 1
      else
          sign_man = -1
          manf = bias - manf
      end if
 
***** Generate the floating point value
      ieee = sign_man*(manf/bias)*(2.d0)**(sign_exp*exponent)
      dfchi = 0
 
***** Thats all
      return
      end
 
 
      integer*2 function efchi ( ihp, ieee )
 
 
*     Conversion from HP floating real*6 to IEEE real*4 floating
*     point.
 
*         ihp(2) - HP floating point (a integer for manipuation)
*                - (**WARNING** the lower bytes of ihp(2) do not
*                -    belong to this routine (ie. only 6 bytes were
*                -    passed), therefore never modify ihp(2))
*         mantissa(2)- Mantissa of floating point (most significant
*                    - part in mantissa(1), mantissa(2) is effectively
*                    - an I*2 variable)
*         exponent   - Exponent of floating point
*         sign_man   - Sign of mantissa
*         sign_exp   - sign of exponent
*         mask_man(2)- Mask for mantissa in two parts.  The second
*                    - part of the mask maskes the lower two bytes
*                    - of ihp(2) as well (effectively converting it
*                    - to a I*2 value (expect for two's complement
*                    - sign change)
*         mask_exp   - Mask for exponent
*         mask_exps  - Mask for exponent sign
*         mask_mans  - Mask for mantissa sign
 
      integer*4 ihp(2), mantissa(2), exponent, sign_man, sign_exp,
     .    mask_man(2), mask_exp, mask_exps, mask_mans, cand
 
*         ieee   - Should the IEEE floating point value corresponding
*                - to ihp
 
      real*4 ieee
 
*         manf   - Floating point version of the mantissa
*         bias   - 2**39.  Number of bits in the mantissa
 
 
      real*8 manf, bias
 
*                       12345678901
      data mask_man  / O'17777777777' , O'37777600000' /
     .,    mask_exp  / O'00077400000' /
     .,    mask_exps / O'00000200000' /
     .,    mask_mans / O'20000000000' /
 
***** Get the mantissa
      bias = 2.d0**39
*     bias = 0.72057594037927940d17/2
 
*     write(*,100)  2.d0**39, 2.d0**24
  100 format(' Bias = ',2d24.17)
      mantissa(1) = cand(ihp(1), mask_man(1) )
      mantissa(2) = cand(ihp(2), mask_man(2) )/2.d0**24
      if( mantissa(2).lt.0 ) then
          mantissa(2) = 256 + mantissa(2)
      end if
      manf = mantissa(1)*2.d0**8    + mantissa(2)
*     write(*,110) mantissa
  110 format(3o16)
*     write(*,120) mantissa
  120 format(3d25.16)
*                                                 ! = 2**17
      exponent = cand(ihp(2), mask_exp )/131072
      sign_exp = cand(ihp(2), mask_exps )
      if( sign_exp.eq.0 ) then
          sign_exp = 1
      else
          sign_exp = -1
          exponent = 128 - exponent
      end if
 
      sign_man = cand(ihp(1), mask_mans )
      if( sign_man.eq.0 ) then
          sign_man = 1
      else
          sign_man = -1
          manf = bias - manf
      end if
 
***** Generate the floating point value
      ieee = sign_man*(manf/bias)*(2.d0)**(sign_exp*exponent)
      efchi = 0
 
***** Thats all
      return
      end
 
 
      integer*2 function efchid( ihp, ieee )
 
 
*     Conversion from HP floating real*6 to IEEE real*4 floating
*     point.
 
*         ihp(2) - HP floating point (a integer for manipuation)
*                - (**WARNING** the lower bytes of ihp(2) do not
*                -    belong to this routine (ie. only 6 bytes were
*                -    passed), therefore never modify ihp(2))
*         mantissa(2)- Mantissa of floating point (most significant
*                    - part in mantissa(1), mantissa(2) is effectively
*                    - an I*2 variable)
*         exponent   - Exponent of floating point
*         sign_man   - Sign of mantissa
*         sign_exp   - sign of exponent
*         mask_man(2)- Mask for mantissa in two parts.  The second
*                    - part of the mask maskes the lower two bytes
*                    - of ihp(2) as well (effectively converting it
*                    - to a I*2 value (expect for two's complement
*                    - sign change)
*         mask_exp   - Mask for exponent
*         mask_exps  - Mask for exponent sign
*         mask_mans  - Mask for mantissa sign
 
      integer*4 ihp(2), mantissa(2), exponent, sign_man, sign_exp,
     .    mask_man(2), mask_exp, mask_exps, mask_mans, cand
 
*         ieee   - Should the IEEE floating point value corresponding
*                - to ihp
 
      real*8 ieee
 
*         manf   - Floating point version of the mantissa
*         bias   - 2**39.  Number of bits in the mantissa
 
 
      real*8 manf, bias
 
*                       12345678901
      data mask_man  / O'17777777777' , O'37700000000' /
     .,    mask_exp  / O'00077400000' /
     .,    mask_exps / O'00000200000' /
     .,    mask_mans / O'20000000000' /
 
***** Get the mantissa
      bias = 2.d0**39
*     bias = 0.72057594037927940d17/2
 
*     write(*,100)  2.d0**39, 2.d0**24
  100 format(' Bias = ',2d24.17)
      mantissa(1) = cand(ihp(1), mask_man(1) )
      mantissa(2) = cand(ihp(2), mask_man(2) )/2.d0**24
      if( mantissa(2).lt.0 ) then
          mantissa(2) = 256 + mantissa(2)
      end if
      manf = mantissa(1)*2.d0**8    + mantissa(2)
*     write(*,110) mantissa
  110 format(3o16)
*     write(*,120) mantissa
  120 format(3d25.16)
*                                                 ! = 2**17
      exponent = cand(ihp(2), mask_exp )/131072
      sign_exp = cand(ihp(2), mask_exps )
      if( sign_exp.eq.0 ) then
          sign_exp = 1
      else
          sign_exp = -1
          exponent = 128 - exponent
      end if
 
      sign_man = cand(ihp(1), mask_mans )
      if( sign_man.eq.0 ) then
          sign_man = 1
      else
          sign_man = -1
          manf = bias - manf
      end if
 
***** Generate the floating point value
      ieee = sign_man*(manf/bias)*(2.d0)**(sign_exp*exponent)
      efchid = 0
 
***** Thats all
      return
      end
 
 
