CTITLE EXTENDED_HEADER
 
      subroutine extended_header(lu)

      implicit none
 
*     J.L. Davis                   1:09 PM  FRI.,  5  JUNE, 1987
*
*     Routine to output "extended" header information.
 
 
      include '../includes/kalman_param.h'
      include '../includes/obs_header.h'
 
*       i               - Loop counter
*   ,   lu              - LU for output
*   ,   TrimLen         - HP function
 
      integer*4 i, lu, TrimLen
 
***** KalObs record length
      write(lu,100) obs_rec_len
  100 format(/,' Record length for this KalObs is ',I4,' blocks.')
 
***** Report the edit conditions
      call report_edit(lu)

***** Report clock breaks

      call report_break( lu )
 
***** Write apriori Markov statistics for all sites
      do i = 1, num_sites
 
          call wr_site_markv_apr(lu,i)
 
      end do
 
***** Write PMU descriptor
      if (TrimLen(user_pmu_dsc) .gt. 0)
     .    write(lu,200) user_pmu_dsc(1:TrimLen(user_pmu_dsc))
 
  200 format(A)
 
      end
 
