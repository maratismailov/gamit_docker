*   Fundamental freqencies for GNSS in Hz       
 
      real*8 gps_f1,gps_f2,gps_f5, glonass_f1,glonass_f2,glonass_f3
     .     , glonass_df1,glonass_df2,glonass_df3
     .     , beidou_f2,beidou_f5,beidou_f7,beidou_f6
     .     , galileo_f1,galileo_f5,galileo_f6,galileo_f7,galileo_f8
     .     , irnss_f9,irnss_f5
      parameter ( gps_f1 = 154*10.23d6 )   ! GPS L1 1575.42 MHz  
      parameter ( gps_f2 = 120*10.23d6 )   ! GPS L2 1227.60 MHz
      parameter ( gps_f5 = 115*10.23d6 )   ! GPS L5 1176.45 MHz

      parameter ( glonass_f1  = 1602.d6 )  ! GLONASS L1 FDMA center
      parameter ( glonass_df1 = 562.5d3 )  ! GLONASS L1 FDMA separation
      parameter ( glonass_f2  = 1246.d6 )  ! GLONASS L2 FDMA center
      parameter ( glonass_df2 = 437.5d3 )  ! GLONASS L2 FDMA separation
      parameter ( glonass_f3  = 1201.d6 )   ! GLONASS L3 FDMA center
      parameter ( glonass_df3 = 475.5d3 )  ! GLONASS L3 FDMA separation (check this!)

      parameter ( beidou_f2 = 1561.098d6 ) ! Beidou B1I (http://en.beidou.gov.cn/SYSTEMS/Officialdocument/201902/P020190227601370045731.pdf)
      parameter ( beidou_f5 = 1176.45d6 )  ! Beidou B2a (http://en.beidou.gov.cn/SYSTEMS/Officialdocument/201806/P020180608525870555377.pdf)
      parameter ( beidou_f7 = 1207.14d6 )  ! Beidou B2b (http://en.beidou.gov.cn/SYSTEMS/Officialdocument/202008/P020200803544811195696.pdf)
      parameter ( beidou_f6 = 1268.520d6 ) ! Beidou B3I (http://en.beidou.gov.cn/SYSTEMS/Officialdocument/201806/P020180608525869304359.pdf)

      parameter ( galileo_f1 = 1575.42d6 ) ! Galileo E1 
      parameter ( galileo_f5 = 1176.45d6 ) ! Galileo E5a
      parameter ( galileo_f6 = 1278.75d6 ) ! Galileo E6 
      parameter ( galileo_f7 = 1207.140d6) ! Galileo E5b
      parameter ( galileo_f8 = 1191.795d6) ! Galileo E5(E5a+E5b)

      parameter ( irnss_f9 = 2492.028d6 )  ! IRNSS S-band  (L9)
      parameter ( irnss_f5 = 1176.45d6 )   ! IRNSS I5 

*     Multipliers now defined in model, autcln, solve, and track
 
c      parameter ( dfsf = (f1-fL2)/(f1+f2) )
c      parameter ( sfdf = (f1+fL2)/(f1-f2) )

c      parameter ( lcf1 = 1.d0/(1.d0 - (f2/f1)**2) )
c      parameter ( lcf2 = -(f2/f1)/(1.d0 - (f2/f1)**2) )

c      parameter ( lgf1 = -f2/f1)
c      parameter ( lgf2 = 1.d0 )

c      parameter ( exf1 = 1.d0     )
c      parameter ( exf2 = -f1/f2 )

c      parameter ( pcf1 =  f1**2/(f1**2-f2**2) )
c      parameter ( pcf2 = -f2**2/(f1**2-f2**2) )


