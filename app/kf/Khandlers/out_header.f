CTITLE OUT_HEADER
 
      subroutine out_header(lu_out, data_file)

      implicit none
 
 
*     J.L. Davis 870311
*
*     Writes the header to the output LU
*
* MOD JLD 870419 To output available met & contribution flags and
*                available baseline flag
 
      include '../includes/kalman_param.h'
      include '../includes/obs_header.h'
 
*       i, j                    - Loop counters
*   ,   len_title               - Length of experiment title
*   ,   len_data_file           - Length of the data file name
*   ,   lu_out                  - LU for writing header
*   ,   date(5)                 - epoch in Y M D H M S format
*   ,   TrimLen                 - HP function
 
      integer*4 i, j, len_title, len_data_file, lu_out, date(5),
     .    TrimLen
 
*   seconds     - Seconds part of epoch
 
      real*8 seconds
 
*   data_file   - Name of the KalObs file
 
      character*(*) data_file
 
***** Write out info
      call jd_to_ymdhms( start_epoch, date, seconds )
 
      len_data_file = max(1, trimlen(data_file))
      write(lu_out,300) data_file(1:len_data_file), kal_ver/100.d0,
     .                  (date(i),i=1,5), nint(seconds)
      write(lu_out,400) data_base, version, CALC_ver
 
      len_title = max(1,TrimLen(expt_title))
      write(lu_out,500) expt_title(1:len_title)
 
      write(lu_out,600)
      do i = 1, num_sites
        write(lu_out,625) site_names(i),mon_names(i),
     .                   (ecc_change(j,i),j=1,3),ellip_hgt(i),
     .                    axis_offsets(i), axis_types(i),
     .                    avail_site(i),avail_met(i)
      end do
 
      write(lu_out,650) (source_names(i),i=1,num_sources)
      write(lu_out,700) num_obs
      write(lu_out,750) avail_baseline
 
*     Read epoch
      call jd_to_ymdhms( read_epoch,    date, seconds)
      write(lu_out,800) 'Read epoch  :', date, nint(seconds)
 
*     Start epoch
      call jd_to_ymdhms( start_epoch,    date, seconds)
      write(lu_out,800) 'Start epoch :', date, nint(seconds)
*     Mid epoch
      call jd_to_ymdhms( mid_epoch,      date, seconds)
      write(lu_out,800) 'Mid epoch   :', date, nint(seconds)
*     End epoch
      call jd_to_ymdhms( end_epoch,      date, seconds)
      write(lu_out,900) 'End epoch   :', date, nint(seconds)
 
  300 format(//,' Data File ',a,t50,'V ',f5.2,
     .          ' Start date ',I4,2('/',I2.2),1X,I2, 2(':',I2.2))
  400 format(' Database ',A,' ver ',I2,' (CALC version ',F6.2,')')
  500 format(1X,A)
  600 format(/,24X,'Eccentricity',8X,'Ellip',7X,'AXIS',
     .       /,3X,'Site',3X,'Monument',4X,'X',7X,'Y',7X,'Z',5X,'Height',
     .         3X,'Offset',1X,'Type',1X,'Avail',1X,'Avail',
     .       /,21X,'(m)',5X,'(m)',5X,'(m)',5X,'(m)',6X,'(m)',8X,'site',
     .         3X,'Met',
     .       /,2(1X,'--------'),3(1X,'-------'),1X,'--------',2X
     .       ,'------',    1X,'----',1X,'-----',1X,'-----')
  625 format(2(1X,A8),3(1X,F7.3),1X,F8.3,1X,F7.3,1X,I4,1X,o5,1X,o5)
  650 format(/,' Sources:',5(:,/,1X,8(:,A,1X)))
  700 format(/,' Number of Observations = ',I5)
  750 format(/,' AVAIL_BASELINE = ',o8)
  800 format(1x,a,t20,I4,2('/',I2.2),1X,I2,2(':',I2.2))
  900 format(1x,a,t20,I4,2('/',I2.2),1X,I2,2(':',I2.2),/)
 
      end
 
