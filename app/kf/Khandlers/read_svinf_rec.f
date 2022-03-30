CTITLE READ_SVINF_REC
 
      subroutine read_svinf_rec

      implicit none
 
 
*     Routine to read the PRN to SVN number records so that the 
*     svn number can be appended to PRN_pnsv where pn is the prn
*     and sv is the sv number.
*     Results are returned through the  ../glinit/qsvi_rec.h header
*     common 

*     Needs to be called before the rw_name_block subroutine

      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
 
      include '../includes/glb_hdr_def.h'
      include '../glinit/qsvi_rec.h'
 
*   ierr            - File reading error flag
 
      integer*4 ierr
  
      integer*4 rec_copy(128)   ! Buffer for reading svs information
                                ! records
      integer*4 in_rec          ! SVS record number to read
      integer*4 i, len_read 
 
* MOD TAH 150825: Now check the satellite antenna offset values.  Here we 
*     read the records and update if this hfile has a run-time greater than
*     any previously processed ones. 
*     Read over the satellite information records
      if( crec_svinf.gt.0 ) then   ! We have satellite records
          in_rec = crec_svinf
          do i = 1, cnum_svs, 4
             call readd(cglb_dcb, ierr, rec_copy, 128, len_read, in_rec)
             call gr_svinf_rec(i,rec_copy(1))
             call gr_svinf_rec(i+1, rec_copy(33))
             call gr_svinf_rec(i+2, rec_copy(65))
             call gr_svinf_rec(i+3, rec_copy(97))
             in_rec = in_rec + 1
          end do
!         write(*,120) (qsvi_prn(i), qsvi_svn(i), i=1,cnum_svs) 
!120      format('SVINF ',32(I2.2,1x,I2.2,2x))
      endif

      return
      end

CTITLE GR_SVINF_REC

      subroutine gr_svinf_rec(i, record )
      
      implicit none

*     Routine to decode the satellite information records:
*     This GR_SVINF version used to read satellite antenna information in glinit/glist

      include '../includes/kalman_param.h'
      include '../includes/glb_hdr_def.h'
      include '../includes/svinf_def.h'
      include '../glinit/qsvi_rec.h'

* PASSED VALUES

* i  -  Station number
* record(32)  - Record for the Station information

      integer*4 i, j, record(32)


*     Start saving the information
      call wmov( record, 1, svi_prn, 1, 32 )

*     Save the values
      qsvi_prn(i)    = svi_prn
      qsvi_svn(i)    = svi_svn
      qsvi_block(i)  = svi_block

      qsvi_antmod(i) = svi_antmod
      qsvi_ocode(i)  = svi_ocode
      do j = 1,3
         qsvi_antpos(j,1,i) = svi_antpos(j,1)
         qsvi_antpos(j,2,i) = svi_antpos(j,2)
      end do
      qsvi_launch(i) = svi_launch
                             
      return
      end
      
 

