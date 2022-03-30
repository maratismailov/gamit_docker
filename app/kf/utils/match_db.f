      Program match_db
 
*     This program will scan a KO catlg and compare with CfA Tape
*     library to see which data bases currently have no KalObs files
*     associated with them.
 
*   max_exper       - Maxiumum # of experiments
 
      integer*4 max_exper
 
      parameter ( max_exper = 5000 )
 
*
*   dbk_ver(max_exper)  - Data base version associated with
*               - KalObs file
*   dbt_ver(max_exper)  - Data base version associate with
*               - each data base.
*   dbt_tape(max_exper) - Tape number for each data base
*   ct                  - Current tape number
*   dtk(max_exper)  - Data base name entry for each KalObs
*   ktd(max_exper)  - KalObs file for each data base.
*   nk, nd          - Number of KalObs file and data bases
*   i,j,k           - Loop counters
*   ierr            - IOSTAT error
*   rcpar           - Runstring reading
*   len_run         - Length of runstring
*   trimlen         - Length of string
 
      integer*4 dbk_ver(max_exper), dbt_ver(max_exper),
     .    dbt_tape(max_exper), ct, dtk(max_exper), ktd(max_exper),
     .    nk, nd, i,j,k, ierr, rcpar, len_run, trimlen
 
*   kal_match(max_exper)    - KalObs file has matching data base
*   db_match(max_exper) - Data base has a KalObs file
 
      logical kal_match(max_exper), db_match(max_exper)
 
*   dbt_db(max_exper)       - Data base name on each tape
*   kal_db(max_exper)       - Data base name for each KalObs file
 
      character*10 dbt_db(max_exper), kal_db(max_exper)
 
*   ko_name(max_exper)  - KalObs file names
*   kalin               - Name of file with KalObs file
*   dbtin               - Name of file with Data base info
 
      character*128 ko_name(max_exper), kalin, dbtin
 
*   line                - Line used in reading files
 
      character*256 line
 
****  Start, get the file names and open
 
      write(*,'(/,'' Match_db Running'')')
      len_run = rcpar(1, kalin)
      if( len_run.gt.0 ) then
          open(100, file=kalin, iostat=ierr, status='old')
          call report_error('IOSTAT',ierr,'open',kalin,1,'match_db')
      else
          call proper_runstring('match_db.hlp',6,1)
      end if
 
      len_run = rcpar(2, dbtin)
      if( len_run.gt.0 ) then
          open(101, file=dbtin, iostat=ierr, status='old')
          call report_error('IOSTAT',ierr,'open',dbtin,1,'match_db')
      else
          call proper_runstring('match_db.hlp',6,1)
      end if
 
****  Now read the KalObs file list
      i = 0
      do while( ierr.eq.0 )
          read(100,'(a)', iostat=ierr) line
*                                                        ! OK, continue
          if( ierr.eq.0 .and. line(2:4).ne.'KDL' ) then
              i = i + 1
              read(line,'(a)',iostat=ierr) ko_name(i)
              read(100,'(a)', iostat=ierr) line
              j = index(line,'DB')
              read(line(j+3:),100,iostat=ierr) kal_db(i), dbk_ver(i)
 100          format(a10,2x,i3)
          end if
      end do
 
      nk = i
      write(*,150) nk, kalin(1:trimlen(kalin))
 150  format(' There are ',i4,' KalObs files in ',a)
 
      close(100)
 
****  Now read the tape information
      ierr = 0
      i = 0
      do while ( ierr.eq.0 )
 
          read(101,'(a)', iostat=ierr ) line
*                                             ! Get tape number
          if( line(1:4).eq.'Tape' .and. ierr.eq.0 ) then
              j = index(line,'?')
              read(line(j+1:),*, iostat=ierr) ct
*                     ! See if we should db info
          else
              if( trimlen(line).gt.0 .and. ierr.eq.0 ) then
                  i = i + 1
                  read(line,200,iostat=ierr) dbt_db(i), dbt_ver(i)
  200             format(6x,a10,3x,i4)
                  dbt_tape(i) = ct
              end if
          end if
      end do
 
      nd = i
      write(*,220) nd, dbtin(1:trimlen(dbtin))
 220  format('There are ',i4,' data bases in ',a,/,
     .       'Scaning now to find data bases with out KalObs files')
 
 
****  Initialize the pointers
      do i = 1, nk
          dtk(i) = -1
          kal_match(i) = .false.
      end do
 
      do i = 1, nd
          ktd(i) = -1
          db_match(i) = .false.
      end do
 
      do i = 1, nd
 
*         Scan KalObs list
          j = 0
          do while (j.lt.nk )
              j = j + 1
*                                                 ! Found name match
              if( dbt_db(i).eq.kal_db(j) ) then
                  db_match(i) = .true.
                  kal_match(j) = .true.
 
*****             Check to see if we have previous match.  If we do
*                 then use the latest version
                  if( ktd(j).gt.0 ) then
                      if( dbk_ver(j).eq.dbt_ver(i) ) then
                          ktd(j) = i
                      end if
                  end if
                  dtk(i) = j
              end if
          end do
      end do
 
****  Now Output the data bases with out matches
      write(*,300)
 300  format(/' X-Band Data Bases without Matching KalObs files',/,
     .        '   Name    Ver     Tape')
 
      do i = 1, nd
          if( .not. db_match(i) .and. dbt_db(i)(9:9).eq.'X' ) then
              write(*,320) dbt_db(i), dbt_ver(i), dbt_tape(i)
 320          format(1x,a10,1x,i3,1x,i3)
          end if
      end do
 
****  Thats all
      end
 
 
 
