 
*     This include concatinates the two parts of the globk common
*     blocks.  The two parts are similar to the SOLVK common in
*     that the control part contains the information needed by the
*     filter segments (this part is as small as possible).  The other
*     part contains the remainder of the infomation needed by the
*     GLOBK Package.
*
*     NOTE: This common uses EMA/VMA and therefore a $ema / glb_ema_block /
*     directive is needed whenever the common block is included.
* MOD TAH 090930: Changed globk_cntl.h  
*
*                                  4:12 PM  WED., 29  JULY, 1987
 
      include '../includes/globk_cntl.h'
 
      include '../includes/globk_markov.h'
 
      include '../includes/globk_ema.h'
 
