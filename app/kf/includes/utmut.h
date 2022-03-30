*     Include file containg constants that need to be used 
*     with UTM conversions

*     Constant for input system
      real*8 F_ut   ! Flattening
      real*8 RF_ut  ! Recipocal flattening
      real*8 ER_ut  ! Semimajor axis
      real*8 ESQ_ut 
      real*8 EPS_ut 
      real*8 PR_ut
      real*8 EN_ut, EN2_ut, EN3_ut, EN4_ut
      real*8 A_ut
      real*8 B_ut
      real*8 C_ut
      real*8 U_ut, U0_ut, U2_ut, U4_ut, U6_ut
      real*8 V_ut, V0_ut, V2_ut, V4_ut, V6_ut
      real*8 W_ut
      real*8 R_ut
      real*8 OMO_ut
      real*8 SO_ut
      real*8 SF_ut
      real*8 OR_ut 
      real*8 FE_ut, FN_ut
      real*8 dXYZ_ut(3)   ! Offset of center of system 

      character*16 datum_ut

*     Common declaration
      common / datum_defs /  F_ut, RF_ut, ER_ut, ESQ_ut, 
     .         EPS_ut, PR_ut, EN_ut, EN2_ut, EN3_ut, EN4_ut, 
     .         A_ut, B_ut, C_ut,
     .         U_ut, U0_ut, U2_ut, U4_ut, U6_ut,
     .         V_ut, V0_ut, V2_ut, V4_ut, V6_ut,
     .         W_ut, R_ut, OMO_ut, SO_ut, 
     .         SF_ut, OR_ut, FE_ut, FN_ut, dXYZ_ut

      common / datum_names / datum_ut



