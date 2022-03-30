CTITLE GET_SVSRUN
 
      subroutine get_svsrun

      implicit none
 
*     Routine to get the cunstring.
 
      include 'track_com.h'
 
* LOCAL VARIABLES
 
*  len_run     - Length of the runstring element
*  trimlen     - Length of string
*  i           - Loop counter
*  rcpar       - Gets runstring entry
 
 
      integer*4 len_run, ii, rcpar, trimlen, num_runs, i

       
*  runstring   - Elements of runstring
 
      character*256 runstring

      do i = 1, max_runstr
         runstr(i) = ' '
      enddo
      num_runs = 0
      
      len_run = rcpar(1,runstring)
      if( len_run.le.0.or.runstring(1:1).ne.'-' ) then
          call proper_runstring('track.hlp','track',-1)
      end if
     
      ii=1
      do while(len_run.ne.0)
         len_run = rcpar(ii,runstring)
         if(len_run.ne.0) then
            if(runstring(1:2).eq.'-f') then

*               See if command file names passed
                len_run = rcpar(ii+1,bat_file)
                if( len_run.le.0 ) then
                    call err_msg('-f')
                    bat_file = 'makexk.bat'
                    ii=ii+1
                else
                    ii=ii+2
                endif

            else if ( runstring(1:2).eq.'-a' ) then
                len_run = rcpar(ii+1,ambin_file)
                if( len_run.le.0 ) then
                    call err_msg('-a')
                    ambin_file = ' '
                    ii=ii+1
                else
                    ii=ii+2
                endif

            else if ( runstring(1:2).eq.'-d' ) then
                len_run = rcpar(ii+1,runday)
                if( len_run.le.0 ) then
                    call err_msg('-d')
                    runday = ' '
                    ii=ii+1
                else
                    ii=ii+2
                endif

            else if ( runstring(1:2).eq.'-w' ) then
                len_run = rcpar(ii+1,runweek)
                if( len_run.le.0 ) then
                    call err_msg('-w')
                    runweek = ' '
                    ii=ii+1
                else
                    ii=ii+2
                endif

            elseif ( runstring(1:2).eq.'-s' ) then 
*               Loop over strings to max of 10
                len_run = rcpar(ii+1,runstring)
                do while ( len_run.gt.0 .and. runstring(1:1).ne.'-')
                   num_runs = num_runs + 1
                   runstr(num_runs) = runstring
                   ii = ii + 1
                   if( num_runs.gt.max_runstr ) then
                      write(*,160) num_runs
 160                  format('Too many strings entered. ',i4)
                      stop
                   endif
                   len_run = rcpar(ii+1,runstring)
                enddo
                ! print *,'num_runs ',num_runs
                ! print *,'List ', runstr(1:num_runs)

****        See if help requested.
            else  if(runstring(1:2).eq.'-h') then
               ii=ii+1
               call proper_runstring('track.hlp','track',-1)

            else 
               write(*,200) runstring(1:len_run) 
 200           format('TRACK: Unknown runstring entry: ',a)
               ii = ii + 1   ! Go to next entry

            endif      ! batch file
         
         endif
      end do

      write(*,220) track_version, bat_file(1:trimlen(bat_file))
 220  format(/,' TRACK Version ',a,' GPS Kinematic trajectory program',/,
     .         ' Running with batch file: ',a)
      write(*,240) runday(1:max(1,trimlen(runday))), 
     .             runweek(1:max(1,trimlen(runweek)))
 240  format(' DAY and WEEK Options ',2(a,1x))
      write(*,280) (runstr(i)(1:max(1,trimlen(runstr(i)))),i=1,num_runs)
 280  format(' STRING       Options ',99(a,1x))

****  Thats all
      return
      end

      subroutine err_msg(err)

      implicit none

      character*(*) err
 
      write(*,*)" Command Error: Missing augument after ", err
      write(*,*) "               Use makexk -h for help"

      stop 'TRACK: Incomplete runstring'
      end
