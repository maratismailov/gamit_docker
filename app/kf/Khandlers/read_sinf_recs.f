CTITLE READ_SINF_RECS

      subroutine read_sinf_recs

      implicit none

*     Routine to read in the load values from the site information
*     records and save them in the gatm/hydload variables for removal
*     of application later in glfor.


      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
      include '../includes/glb_hdr_def.h'

* LOCAL VARIABLES 
      integer*4 i   ! counter over records
      integer*4 rec ! Record number to read
     .,         len_read  ! Length read
     .,         array(128)   ! Array into which record is read.
      integer*4 ierr

****  Loop over the site recoords
      do i = 1, cnum_sites, 2
          rec = crec_sinf + (i-1)/2 
          call readd(cglb_dcb, ierr, array, 128, len_read, rec)
*         The 'array' record contains two sites so process each
*         half of the record
          call getload_sinf(array(1),ltog_sites(i))
          if( i+1.le.cnum_sites ) then
              call getload_sinf(array(65),ltog_sites(i+1))
          endif
      end do

****  Thats all
      return
      end

CTITLE GETLOAD_SINF

      subroutine getload_sinf ( array, is)

      implicit none

*     Routine to extract load values from site information records
*     and save in the global arrays

* INCLUDE FILES

      include '../includes/kalman_param.h'
      include '../includes/globk_cntl.h'
      include '../includes/glb_hdr_def.h'
      include '../includes/sinf_def.h'

* PASSED VARIABLES 
      integer*4 is   ! global site number
     .,         array(64)   ! Array with file record

* LOCAL
      integer*4 i   ! Counter

****  Move the array into the declaratio common so that we can
*     extract values
      call wmov( array, 1, ssdata_st, 1, 64 )

*     Save the load values
      do i = 1,3
         gatmload(i,is) = satmload(i)
         ghydload(i,is) = shydload(i)
      end do

***** Thats all
      return 
      end






