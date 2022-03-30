
*------------------------------------------------------------------------
*  This include file has the hard-wired GPS frequency quantities once
*  held in const_param.h.   For GNSS it has been supplanted by 
*  libraries/freq_def.h and the arrays in kf/ctogobs_com.h but for the
*  two programs that have not yet been modified to use arrays--kf/svpos
*  and kf/track, this file in temporarily needed.   RWK/TAH 150514

*   fL1, fL2    - GPS frequencies in Hz at L1 and L2

      real*8 fl1, fl2, dfsf, sfdf, lcf1, lcf2, lgf1, lgf2, pcf1, pcf2
     .     , exf1, exf2


      parameter ( fL1 = 154*10.23d6 )    
      parameter ( fL2 = 120*10.23d6 )    
 
      parameter ( dfsf = (fL1-fL2)/(fL1+fL2) )
      parameter ( sfdf = (fL1+fL2)/(fL1-fL2) )

      parameter ( lcf1 = 1.d0/(1.d0 - (fL2/fL1)**2) )
      parameter ( lcf2 = -(fL2/fL1)/(1.d0 - (fL2/fL1)**2) )

      parameter ( lgf1 = -fL2/fL1)
      parameter ( lgf2 = 1.d0 )

      parameter ( exf1 = 1.d0     )
      parameter ( exf2 = -fL1/fL2 )

      parameter ( pcf1 =  fL1**2/(fL1**2-fL2**2) )
      parameter ( pcf2 = -fL2**2/(fL1**2-fL2**2) )

* ION conversions
*  l1tecu - Change in range (m) per TECU at L1
*  l2tecu - Change in range (m) per TECU at L2
      real*8 l1tecu, l2tecu 

      parameter ( l1tecu = 40.3d0/fL1**2*1.d16 )
      parameter ( l2tecu = 40.3d0/fL2**2*1.d16 )


