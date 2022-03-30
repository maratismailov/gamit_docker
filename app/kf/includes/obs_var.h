*---------------------------------------------------------------------*
*                                                       OBS_VAR.FTNI  *
*                                                                     *
*     Include file with number of variabes in observation record and  *
*     resulting number of 16-bit words for flagging which of these to *
*     write.                                                          *
*                                                                     *
*     J.L. Davis 870311                                               *
*                                                                     *
* MOD JLD 870417 Added room for DB_OBSERVATION and DB_SIGMA           *
*                                                                     *
*---------------------------------------------------------------------*
 
*           dim_var             - The dimension of the flags
*   ,       num_var             - Number of variables in observation
 
      integer*4 dim_var, num_var
 
      parameter (num_var = 61)
 
      parameter (dim_var = num_var / 16 + 1)
 
*           list_name(num_var)  - List of variable names
 
      character*20 list_name(num_var)
 
