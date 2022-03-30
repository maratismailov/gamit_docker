      program sdrms

      implicit none 
 
*     This program will read output files from the standard
*     ut_sd.gxt and xy_sd.gxt extract scripts and compute
*     rms scatter.
 
*   ierr    - IOSTAT eror
*   i       - Loop cpunter
*   args(6) - Argummets (not saved)
 
      integer*4 ierr, i, args(6)
 
*   ct(11), ca(11), cs(11)  - Cos total, cos adj, cos sig
*   st(11), sa(11), ss(11)  - sin total, sin adj, sin sig
*   stat_all(2)     - Stats for total and adjust for all terms
*   stat_not(2)     - Stats for non main tide lines.
*   stat_std(2)     - Stats for the 8 terms with correcttions.
 
      real*8 ct(11), ca(11), cs(11), st(11), sa(11), ss(11),
     .    stat_all(2), stat_not(2), stat_std(2)
 
*   line            - Line read from input.
 
      character*256 line
 
*     Start to read the file in groups of 11 elememts
 
      ierr = 0
      do while ( ierr.eq.0 )
          do i = 1,11
              read(*,'(a)',iostat=ierr) line
              if( ierr.eq.0 ) then
                  read(line,*) args, ct(i), ca(i),cs(i),
     .                        st(i), sa(i), ss(i)
              end if
          end do
 
***** If no error then accumulate statistics
          if ( ierr.eq.0 ) then
              call all_stat(ct,ca,cs,st,sa,ss, stat_all)
 
*             Now do stats for small tidal lines
              if( index(line, 'Prograde').gt.0 ) then
                  call pro_stat(ct,ca,cs,st,sa,ss, stat_not,stat_std)
              else
                  call ret_stat(ct,ca,cs,st,sa,ss, stat_not,stat_std)
              end if
 
*             Now write out results
              call out_stats(line(72:), stat_all, stat_not,stat_std,
     .                       cs(1))
          end if
      end do
 
      end
 
CTITLE ALL_STAT
 
      subroutine all_stat( ct,ca,cs, st,sa,ss, stat_all)

      implicit none 
 
 
*     Compute rms scatter using all data
 
*   ct(11), ca(11), cs(11)  - Cos total, cos adj, cos sig
*   st(11), sa(11), ss(11)  - sin total, sin adj, sin sig
*   stat_all(2)     - Stats for total and adjust for all terms
 
      real*8 ct(11), ca(11), cs(11), st(11), sa(11), ss(11),
     .    stat_all(2)
 
*   i       - Loop counter
 
      integer*4 i
 
*   sum_t, sum_a    - Sum of totals and adjustments squared
 
      real*8 sum_t, sum_a
 
      sum_t  = 0
      sum_a  = 0
      do i = 1,11
          sum_t = sum_t + ct(i)**2 + st(i)**2
          sum_a = sum_a + ca(i)**2 + sa(i)**2
      end do
 
      stat_all(1) = sqrt(sum_t/22)
      stat_all(2) = sqrt(sum_a/22)
 
      return
      end
 
CTITLE PRO_STAT
 
      subroutine pro_stat( ct,ca,cs, st,sa,ss, stat_not,stat_std)

      implicit none  
 
*     Compute rms scatter using non-tidal data (prograde) data
 
*   ct(11), ca(11), cs(11)  - Cos total, cos adj, cos sig
*   st(11), sa(11), ss(11)  - sin total, sin adj, sin sig
*   stat_not(2)     - Stats for total and adjust for all terms
 
      real*8 ct(11), ca(11), cs(11), st(11), sa(11), ss(11),
     .    stat_not(2), stat_std(2)
 
*   i       - Loop counter
 
      integer*4 i
 
*   sum_t, sum_a    - Sum of totals and adjustments squared
 
      real*8 sum_t, sum_a, sus_t, sus_a
 
      sum_t  = 0
      sum_a  = 0
      sus_t  = 0
      sus_a  = 0
      do i = 1,11
          if( i.ne.1 .and. i.ne.2 .and. i.ne.4 .and. i.ne.6 ) then
              sum_t = sum_t + ct(i)**2 + st(i)**2
              sum_a = sum_a + ca(i)**2 + sa(i)**2
          else
              sus_t = sus_t + ct(i)**2 + st(i)**2
              sus_a = sus_a + ca(i)**2 + sa(i)**2
          end if
      end do
 
      stat_not(1) = sqrt(sum_t/14)
      stat_not(2) = sqrt(sum_a/14)
      stat_std(1) = sqrt(sus_t/ 8)
      stat_std(2) = sqrt(sus_a/ 8)
 
      return
      end
 
CTITLE RET_STAT
 
      subroutine ret_stat( ct,ca,cs, st,sa,ss, stat_not,stat_std)

      implicit none  
 
*     Compute rms scatter using non-tidal data (retrograde) data
 
*   ct(11), ca(11), cs(11)  - Cos total, cos adj, cos sig
*   st(11), sa(11), ss(11)  - sin total, sin adj, sin sig
*   stat_not(2)     - Stats for total and adjust for all terms
 
      real*8 ct(11), ca(11), cs(11), st(11), sa(11), ss(11),
     .    stat_not(2), stat_std(2)
 
*   i       - Loop counter
 
      integer*4 i
 
*   sum_t, sum_a    - Sum of totals and adjustments squared
 
      real*8 sum_t, sum_a, sus_t, sus_a
 
      sum_t  = 0
      sum_a  = 0
      sus_t  = 0
      sus_a  = 0
      do i = 1,11
          if( i.ne.6 .and. i.ne.8 .and. i.ne.10 .and. i.ne.11 ) then
              sum_t = sum_t + ct(i)**2 + st(i)**2
              sum_a = sum_a + ca(i)**2 + sa(i)**2
          else
              sus_t = sus_t + ct(i)**2 + st(i)**2
              sus_a = sus_a + ca(i)**2 + sa(i)**2
          end if
      end do
 
      stat_not(1) = sqrt(sum_t/14)
      stat_not(2) = sqrt(sum_a/14)
 
      stat_std(1) = sqrt(sus_t/ 8)
      stat_std(2) = sqrt(sus_a/ 8)
      return
      end
 
CTITLE OUT_STATS
 
      subroutine out_stats(type, stat_all, stat_not,stat_std, cs)

      implicit none 

*   stat_all(2)     - Stats for total and adjust for all terms
*   stat_not(2)     - Stats for total and adjust for all terms
 
      real*8 stat_all(2), stat_not(2), stat_std(2), cs
 
*   type            - Type of stat
 
      character*(*) type
 
      write(*,100) type, stat_all, cs, stat_not, cs, stat_std, cs
100   format(' For ',a20,' Tot rms    Adj rms  Sig' ,/,
     .,      ' All lines (22) ', 5x,3F10.2,/,
     .,      ' Non-tide  (14) ', 5x,3F10.2,/, 
     .,      ' Corrected ( 8) ', 5x,3F10.2   )
 
      return
      end
