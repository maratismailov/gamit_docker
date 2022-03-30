
 
CTITLE fix_apr
 
      program fix_apr
 
*     J.L. Davis                   2:20 PM  WED., 27  MAY , 1987
*
C                                                           in
*     Program to update earth-orientation values (PM & UT1)
*     a KalObs file.  Delay and delay rate theoreticals as well
*     as the table values are changed.  Note that only four points
*     are used for the new tabular values.
 
* MOD TAH 880304  Added code to remove a linear trend (weighted by
*     the time difference from mid_epoch) from the file values of
*     polar motion and UT1.  This is particularly important for
*     UT1 which has a high slope. Changes made in GET_NEW_TABLES.
 
* MOD TAH 880317 Introduced new UT1 type.  Bit 2 of type is set if
*     series is not regularized. ie. if short period UT1 terms have
*     not been removed. (GSFC and NGS are usually like this)
*
* MOD TAH 881228 Added output of corrections to a print device so that
*     we can keep a record of the changes.  Also fixed leap second problem
*     when we are crossing a boundary and using UT1-UTC tables.

* MOD TAH 910115: Added check on date at which fix_apr was run so that
*     runs which did not save the pmu_calc values could be detected.  These
*     runs all occurred before 90/11/25

      include '../includes/kalman_param.h'
      include '../includes/obs_header.h'
 
 
*       num_header          - Number of header lines in EOR file
*                           -   (Default = 10)
*   ,   ut1_def             - Definition of values in UT1 column.
*                           -   1=UT1-UTC (default);2=UT1-TAI
*   ,   printer             - LU for print device
 
      integer*4  printer
 
*  date_bug - Date of not saving pmu_calc bug fixed.

      integer*4 date_bug(5)

*  jd_bug  - Julian date at which the bug of not saving pmu_calc was
*            fixed
*  sectag  - Seconds tag

      real*8 jd_bug, sectag

 
*       pmu_file            - PM-UT1 file name
*   ,   KalObs_file         - KalObs file name
 
      character*256 pmu_file, KalObs_file

      data date_bug / 1990, 11, 25, 0, 0 /
 
***** Get the runstring parameters
      call get_runstring(KalObs_file)

*     Set the bug fix jd

      sectag = 0.d0
      call ymdhms_to_jd( date_bug, sectag, jd_bug )
 
***** Open the KalObs file, get the header info, and read the "a priori"
*     block
      call get_initial_info(KalObs_file)

* MOD TAH 900920: Check to see the number of data points in the tables.
*     If it is greater than 5 then reset the value.
      if( fut1_inf(3).gt.5 .or. fwob_inf(3).gt.5 ) then
          write(printer,100) fut1_inf, fwob_inf
 100      format('*** WARNING *** UT1/Wobble tables are too large.',
     .           ' Current INF array ares:',/,
     .           ' UT1:    ',f12.2,1x,f5.2,1x,f2.0,1x,f5.0,
     .           ' Wobble: ',f12.2,1x,f5.2,1x,f2.0,/,
     .           ' Resetting values to aviod problems.')      
          fut1_inf(3) = 5
          fwob_inf(3) = 5
      end if
 
 
***** Update values for PM and UT1 at mid_epoch
      call update_mid_pmu(Fut1_inf,Fut1_pts,
     .    Fwob_inf,Fwob_pts,printer, KalObs_file)
 
***** Update apriori info in KalObs file
      call update_apriori_info(pmu_file)
 
***** Close all files
      call finish
 
      end
 
CTITLE GET_RUNSTRING
 
      subroutine get_runstring(KalObs_file)
 
*     J.L. Davis                   4:42 PM  WED., 27  MAY , 1987
*
*     Get the runstring parameters and call proper_runstring if
*     necessary
 
 
*   ,   len_KalObs_file     - Length of string descriptor
*   ,   printer             - Print device
*   ,   rcpar               - HP runstring utility
 
      integer*4 len_KalObs_file,
     .    num_header, printer, rcpar, ut1_def
 
*       pmu_file            - File with PM-UT table
*   ,   KalObs_file         - KalObs file name
 
      character*(*) KalObs_file
 
*   ,   help_file*16        - Help file descriptor  -- not used
 
c      character  help_file*20  -- not used
 
c      data help_file / 'fix_apr.hlp' /  -- not used
 
***** Initialize input strings
      KalObs_file    = ' '
 
***** Read in the runstrings
      len_KalObs_file = rcpar(1,KalObs_file   )
 
***** Initialize the numeric quantities with default values
      num_header =  3
      ut1_def    =  1
      printer    =  6

      if( len_KalObs_file.eq.0 ) then
          stop ' Fix_apr <KalObs> -- Fix wob_apr, ut1_apr values'
      end if
 
      end
 
CTITLE GET_INITIAL_INFO
 
      subroutine get_initial_info(KalObs_file)
 
*     J.L. Davis                   4:58 PM  WED., 27  MAY , 1987
*
*     Do some initial opening and reading of the KalObs file
 
 
      include '../includes/kalman_param.h'
      include '../includes/obs_header.h'
 
*       error                   - Error flag
 
      integer*4 error
 
*       KalObs_file             - KalObs file descriptor
 
      character*(*) KalObs_file
 
 
      character*16 prog
 
      data prog / 'GET_INITIAL_INFO' /
 
***** Open KalObs file
      call open_KalObs(KalObs_file,error)
 
***** Error?
      call report_error('FMGR',error,'open',KalObs_file,1,prog)
 
***** Read the KalObs header
      call rw_KalObs_header('R',error)
 
***** Error?
      call report_error('FMGR',error,'read','KalObs header',1,prog)
 
***** Output the header
      call out_header(6,KalObs_file)
 
      end
 
CTITLE FINISH
 
      subroutine finish
 
*     J.L. Davis                   4:06 PM  MON.,  1  JUNE, 1987
*
*     Routine to close files
*
 
 
      integer*4 error
 
      call close_KalObs(error)
 
      end
 
CTITLE UPDATE_MID_PMU
 
      subroutine update_mid_pmu(new_ut1_inf,new_ut1_pts,
     .    new_wob_inf,new_wob_pts, printer, KalObs_file)
 
*     J.L. Davis                   5:25 PM  TUE.,  2  JUNE, 1987
*
*     This routine does the following:
*
*     1. Updates the values of UT1_APR and WOB_APR to reflect
*        the values in the new tables.
*
*     2. Copies the new tables into the location for the old tables
*
*     3. Writes the a priori block to KalObs file
*
*     4. Closes the KalObs file
 
 
      include '../includes/kalman_param.h'
      include '../includes/obs_header.h'
 
*       date(5)                 - Date opf first UT1 point
*       dev                     - Current Output device
*       error                   - Error flag
*       i,j,k                   - Loop counters
*       printer                 - Print device
*       trimlen                 - Return length of string
 
      integer*4 date(5), dev, i,j,k, printer, trimlen
 
*       new_ut1_inf(4)          - New UT1 info array
*   ,   new_ut1_pts(4)          - New UT1 array
*   ,   new_wob_inf(3)          - New wobble info array
*   ,   new_wob_pts(2,4)        - New wobble array
*       Sectag                  - Seconds tag
 
      real*8 new_ut1_inf(4), new_ut1_pts(4), new_wob_inf(3),
     .    new_wob_pts(2,4), Sectag
 
 
c      character*14 prog  -- not used
 
*       KalObs_file            - Name of the KalObs file
 
      character*(*) KalObs_file
 
c      data prog / 'UPDATE_MID_PMU' /  -- not used
 
***** Use Lagrangian interpolation and the new tables to get the
*     "a priori" values for PM-UT1 at the epoch PMU_EPOCH
      call Lagrange_intp(new_wob_inf,new_wob_pts,pmu_epoch,wob_apr,2)
      call Lagrange_intp(new_ut1_inf,new_ut1_pts,pmu_epoch,ut1_apr,1)
 
***** Convert the units of the apriori UT1 value from seconds to mas
      ut1_apr(1) = 15.0D+03 * ut1_apr(1)
      ut1_apr(2) = 15.0D+03 * ut1_apr(2)
 
***** Tell the user the new values being used
      call jd_to_ymdhms( fut1_inf(1), date, sectag )
 
***** Output the correction,  Do LU1 and the printer (if different)
      do k = 1,2
          dev = 0
          if( k.eq.1 ) dev = 6
          if( k.eq.2 .and. printer.ne.6 ) dev = printer
*                                ! Output the values
          if( dev.ne.0 ) then
              write(dev,50) KalObs_file(1:trimlen(KalObs_file))
  50          format(/' For KalObs File ',a)
 
              write(dev,100) (fut1_inf(i),i=1,3), date,
     .                       (new_ut1_inf(i),i=1,3)
  100         format(' Old and New UT1 information array',
     .       /' Start JD ',f12.3,' Step and # ',2(f5.2,1x),' Date ',
     .         I4,'/',i2,'/',i2,1x,i2,':',i2,
     .       /' Start JD ',f12.3,' Step and # ',2(f5.2,1x) )
 
              write(dev,120) (fut1_pts(i), new_ut1_pts(i),
     .             (fut1_pts(i)-new_ut1_pts(i))*15.d3,i=1,4)
  120         format(' UT1 points [old,new (tsec), delta (mas)]',/,
     .        4(2f12.7,2x,f12.3,:/) )
 
 
              write(dev,160) (fwob_inf(i),i=1,3),
     .                       (new_wob_inf(i),i=1,3)
  160         format( ' Old and New wobble information array',
     .       /' Start JD ',f12.3,' Step and # ',2(f5.2,1x),
     .       /' Start JD ',f12.3,' Step and # ',2(f5.2,1x) )
 
              write(dev,180) ((fwob_pts(j,i), new_wob_pts(j,i),
     .              (fwob_pts(j,i)-new_wob_pts(j,i)),j=1,2),i=1,4)
  180         format(' Wobble points (old,new,del) x and y (mas)',/
     .        4( 3f11.5,4x,3f11.5,/) )
*                   ! Need to output
          end if
*                   ! Looping over two devices
      end do

* MOD TAH Modded lengths of moves to match integer*4 variables.
 
***** Zero out the old PM-UT1 tables
c     call ClearBuffer(fwob_pts,20)
c     call ClearBuffer(fut1_pts,10)
 
***** Transfer the info arrays
c     call MoveWords(new_wob_inf,fwob_inf, 6)
c     call MoveWords(new_ut1_inf,fut1_inf, 8)
 
***** Transfer table values
c     call MoveWords(new_wob_pts,fwob_pts,16)
c     call MoveWords(new_ut1_pts,fut1_pts, 8)

* MOD TAH: Don't be fancy just copy arrays
      do i = 1, 3
         fwob_inf(i) = new_wob_inf(i)
      end do
      do i = 1,4
         fut1_inf(i) = new_ut1_inf(i)
      end do
      do i = 1,5
         fut1_pts(i) = new_ut1_pts(i)
         do j = 1,2
            fwob_pts(j,i) = new_wob_pts(j,i)
         end do
      end do

****  Thats all
      return
      end
 
CTITLE UPDATE_APRIORI_INFO
 
      subroutine update_apriori_info(pmu_file)
 
*     J.L. Davis                   1:47 PM  SAT.,  6  JUNE, 1987
*
*     Routine to make the final changes to the apriori block, and
*     to write the block to the KalObs file.
 
 
      include '../includes/kalman_param.h'
      include '../includes/obs_header.h'
 
*       error               - Error flag
*   ,   i                   - Loop counter
 
      integer*4 error, date(5), i, len_file, TrimLen

      real*8 sectag
 
 
      character*20 prog
 
*       pmu_file            - PMU file name
 
      character*(*) pmu_file
 
      data prog / 'UPDATE_APRIORI_INFO' /
 
***** Get the current date
      call systime(date,sectag)
 
***** How many characters does this file name have?
      len_file = TrimLen(pmu_file)
 
***** Create the polar motion descriptor
      write(user_pmu_dsc,100) (date(i),i=1,5)
  100 format('FIX_APR    run on ',I4,'/',I2.2,'/',I2.2,1X,I2.2,':',
     .    I2.2,' using to fix wob_apr and ut1_apr ')
 
***** Write this out
      write(*,150) user_pmu_dsc(1:trimlen(user_pmu_dsc))
  150 format(A)

***** Write A priori block to KalObs
      call increment_KalVer
      call rw_kalobs_header('W', error )
c     call rw_KalObs_block('W','apr',aprioris,error,0)
 
***** Check for error
      call report_error('FMP',error,'writ','a prioris',0,prog)
 
***** Close KalObs file
      call close_KalObs(error)
 
***** Check for error
      call report_error('FMP',error,'clos','KalObs',0,prog)
 
      end
 
CTITLE get_last_pmu

      subroutine get_last_pmu( user_pmu_dsc, jd_last )

*     This routine will read the user_pmu_dsc and get the last time
*     that update pmu was run.  -1 is returned if no last run

* jd_last   - Julian date of last run

      real*8 jd_last

* user_pmu_dsc - Update pmu description line

      character*(*) user_pmu_dsc

* LOCAL

* trimlen  - Length of string
* indx     - Position in string
* date(5)  - Date run from line

      integer*4 trimlen, indx, date(5), i, ierr

* sectag   - Seconds tag
      real*8 sectag

***** First see if we have anything in string

      if( trimlen(user_pmu_dsc).eq.0 .or.
     .    ichar(user_pmu_dsc(1:1)).eq.0 ) then
          jd_last = -1.d0
          return
      end if

***** Get position of time
      indx = index( user_pmu_dsc,'run on') + 6

*     Kill off the / characters
      call sub_char( user_pmu_dsc, '/',' ')

      read(user_pmu_dsc(indx:),*,iostat=ierr) (date(i),i=1,3)
      call report_error('IOSTAT',ierr,'read',user_pmu_dsc,1,
     .                  'fix_apr')

      date(4) = 0
      date(5) = 0
      sectag  = 0
      call ymdhms_to_jd( date, sectag, jd_last )

****  Thats all
      return
      end

