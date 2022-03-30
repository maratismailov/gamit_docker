*@BAK_PARAM.FTNI
************************************************************************
*                                                                      *
*     J.L. Davis  870413                                               *
*                                                                      *
*     Parameter file for the back-file handlers                        *
*                                                                      *
*     NOTE:  ONLY parameters which are specific to handling the        *
*            back files should be include here                         *
*                                                                      *
************************************************************************
 
*       bak_array_size              - Size of array for reading
*                                   -   BAKFIL
*   ,   bak_dcb_size                - The size of the DCB for
*                                   - handling the back-files
*   ,   bak_record_size             - The maximum size of a record
*                                   -   in BAKFIL in words
*   ,   max_mar_types               - Max types of Markov parameters
*                                   -   in BAKFIL
*   ,   max_types                   - The max total number of types
*                                   -   which have a name in BAKFIL
*                                   -   (see MAR_TYP_DB.FTN)
 
      integer*4 bak_array_size, bak_dcb_size, bak_record_size,
     .    max_mar_types, max_types
 
      parameter (bak_dcb_size = 16)
 
*                                                     ! Epoch
      parameter (bak_record_size = 2
*                                                     ! Markov values
     .                           + 2 * max_parn
*                                                     ! Residuals
     .                           + 2 * max_baselines
*                                                     ! Elevations
     .                           + max_sites
*                                                     ! Flags
     .                           + max_baselines
*                                                     ! Source number
     .                           + 1
*                                                     ! KALFIL rec #'s
     .                           + max_baselines    )
 
      parameter (bak_array_size =
     .     (bak_record_size / max_headr_size) * bak_record_size
     .   + (max_headr_size / bak_record_size) * max_headr_size  )
*                                 ! BAKFIL array size is the greater of the
*                                 ! maximum header length and the maximum
*                                 ! record size.  Thus, this array can
*                                 ! accomodate both.
 
 
      parameter (max_mar_types = 64)
 
      parameter (max_types = max_mar_types + 4)
 
*---------------------------------------------------------------------------
