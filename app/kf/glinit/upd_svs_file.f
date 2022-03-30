CTITLE upd_svs_file
 
      subroutine upd_svs_file ( unit, num, ephem_ep, 
     .                          ephem_frame, ephem_prec, ephem_nut,
     .                          svs_frame,   svs_prec,   svs_nut )

      implicit none 
 
*     Routine to writeout current epehemeric elements to an svs_file
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
 
      include '../includes/glb_hdr_def.h'
*   unit  - Unit number for output file
*   num   - Number of peochs that need to written out

      integer*4 unit, num

*   ephem_ep(num)   - Epochs of the experiments that share this ephemeris

      real*8 ephem_ep(*)

* The following values are to allow B1950 and J2000 values to be mixed.
*   ephem_frame(num) - Frames associated with each of the ephemeris
*                 epochs to be output
*   ephem_prec(num)  - Precession modes assocated with each ephemeris
*   ephem_nut(num)   - Nutation mode (IAU80/IAU00)
*   svs_frame(max_svs) - Frames for each IC value saved
*   svs_prec(max_svs)  - Precession modes associated with each IC value saved
*   svs_nut(max_glb_svs) - Nutation mode IAU80/IAU00

      character*8 ephem_frame(*), ephem_prec(*), ephem_nut(*),
     .            svs_frame(max_glb_svs), svs_prec(max_glb_svs), 
     .            svs_nut(max_glb_svs)

 
* LOCAL VARIABLES
 
*   ierr    - IOSTAT error
*   i,j     - Loop counters
*   date(5) - Date of the experiment. Lines output with this
*           - time so that GLOBK can match the ephemeris entry
*           - global file
*   eph_date(4) - ymdh for the ephemeris entry
*   trimlen - Length of string
*   len_out - Length of the output line before and after adding
*             frame and precession modes.
 
      integer*4 ierr, i,j, n, date(5), eph_date(4), trimlen,
     .          len_out
 
*   sectag  - Seconds tag for time
*   Dec_yrs - Deciminal years
 
      real*8 sectag

*   svs_out(6)  - Satellite position and velocity in the output
*                 frame

      real*8 svs_out(6)
      
*    out_line   - output line.  The orbital elements are written
*                 first and then the frame and ephem information
*                 are appended to the end of the line

      character*512 out_line
 
*     Write some header information

      if ( num.eq.0 ) RETURN

      write(unit,100, iostat=ierr) 
 100  format('* EPHEMERIS INFORMATION')
      call report_error('IOSTAT',ierr,'head writ','empheris file',0,
     .    'upd_svs_file/glinit')
 
*     Get the output time to use
      do n = 1, num 
          call jd_to_ymdhms( ephem_ep(n), date, sectag )
          do i = 1,4
              eph_date(i) = date(i)
          end do

          write(unit,110) ephem_frame(n), ephem_prec(n), ephem_nut(n)
 110      format('* Frame: ',a5,' Precession Mode ',a5,' Nutation  ',a5)

*         Now write the satellite ephemeris
          do i = 1, gnum_svs

*             Only write out elements which have been set
*             Check the size of the elements
              do j = 8, max_svs_elem - 3
                 if( abs(apr_val_svs(j,i)).gt.1.0 ) then
                     apr_val_svs(j,i) = 0.d0
                 end if
              end do

* MOD TAH 970317: Check to see if the radiation parameters should be
*             zero'd rather than using the apriori in the hfile
              if( glb_svs_file(trimlen(glb_svs_file)-1:).eq.'_Z' ) then
                  apr_val_svs(7,i) = 1.d0
                  do j = 8,max_svs_elem - 3
                     apr_val_svs(j,i) = 0.d0
                  end do
              end if

*             Write values
              if( apr_val_svs(1,i).ne.0 ) then

*                 Now see if we need to change the frame and prec mode
*                 Do Position
                  call prec_svs(ephem_ep(n), apr_val_svs(1,i), 
     .                          svs_out(1),
     .                          svs_frame(i), ephem_frame(n),  
     .                          svs_prec(i) , ephem_prec(n),
     .                          svs_nut(i),   ephem_nut(n)  )
                  
*                 Do Velocity
                  call prec_svs(ephem_ep(n), apr_val_svs(4,i), 
     .                          svs_out(4),
     .                          svs_frame(i), ephem_frame(n),  
     .                          svs_prec(i) , ephem_prec(n),
     .                          svs_nut(i),   ephem_nut(n)  )

                  write(out_line,300, iostat=ierr) eph_date, 
     .                    gsvs_names(i),
     .                    (svs_out(j), j= 1,6),
     .                    (apr_val_svs(j,i),j=7,max_svs_elem)
 300              format(i5,3i3,1x,a8,6(1x,f15.5),1x,17(1x,f8.5))

*                 Get the length of the output line
                  len_out = trimlen(out_line)
                  write(out_line(len_out+2:),310,iostat=ierr)
     .                    ephem_frame(n),  ephem_prec(n), ephem_nut(n)
 310              format(a5,1x,a5,1x,a5)
                  len_out = trimlen(out_line) 
*                 Now write the string to the file
                  write(unit,'(a)',iostat=ierr) out_line(1:len_out) 
                                             
              end if
          end do
          call report_error('IOSTAT',ierr,'ephemeris writ',
     .            'Ephemeris file',0,'upd_svs_file/glinit')
          if( ierr.ne.0 ) then
             call report_stat('warning','glinit','upd_svs_file',
     .            glb_svs_file,'IOSTAT error writing file',ierr)
          end if
      end do
 
****  Thats all
      return
      end
 
 
 
