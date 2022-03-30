CTITLE UNIFY_SVS
      program unify_svs
 
      implicit none 

*     The program will read an svs file which has first been sorted by
*     time and then by prn (i.e., sort svs | sort +4 -5) and will make
*     the entries consistent in the sense that all orbital elements at
*     same IC time will be made the same.
*
*     Runstring:
*   % unify_svs in.svs out.svs
*   where in.svs is an input ephemeris file and ou
*         out.svs is the output file.  This output should then be sorted
*                 by time agiain.
 
* VARIABLES
 
*   date(5)     - ymdhm date format
*   ierr        - IOSTAT error
*   trimlen     - Returns the length of a string
*   rcpar       - reads the runstring and returns the entries length
*   len_run     - Length of runstring
*   i           - Loop counter
 
      integer*4 date(5), ierr, trimlen, rcpar, len_run, i
 
*   sectag      - Second tage
*   orb_el(9)   - Orbital elements read from file
*   curr_orb_el(9)  - Current orbital elemenys to be used
*   jd, curr_jd     - Experiment epoch and current experiment epoch
 
      real*8 sectag, orb_el(9), curr_orb_el(9), jd, curr_jd
 
*   prn, curr_prn   - PRN and current prn
 
      character*8 prn, curr_prn
 
*   infile      - Inut file
*   outfile     - Output file
*   line        - Line read from input
 
      character*256 infile, outfile, line
 
****  Read tthe runstring and print the help if not complete
      len_run = rcpar(1,infile)
      if( len_run.le.0 ) then
          call proper_runstring('unify_svs.hlp','unify_svs',1)
      else
          open(100, file=infile, iostat=ierr, status='old')
          call report_error('IOSTAT',ierr,'open',infile,1,'unify_svs')
      end if
 
      len_run = rcpar(2,outfile)
      if( len_run.le.0 ) then
          call proper_runstring('unify_svs','unify_svs.hlp',1)
      else
          open(200, file=outfile, iostat=ierr, status='new')
          call report_error('IOSTAT',ierr,'open',outfile,1,'unify_svs')
      end if
 
****  Start readng the input
      date(5) = 0
      sectag = 0
      curr_jd = 0
      curr_prn = ' '
      do while ( ierr.eq.0 )
          read(100,'(a)', iostat=ierr) line
          if( ierr.eq.0 ) then
              if( line(1:1).eq.' ' .and.trimlen(line).gt.0 ) then
                  read(line,*) (date(i),i=1,4), prn, orb_el
 
*                 See if things have changed
                  call ymdhms_to_jd(date, sectag, jd)
                  if( prn.eq.curr_prn .and.
     .                abs(jd-curr_jd).lt.1.0d0 .and.
     .                abs(orb_el(1)-curr_orb_el(1)).lt.100.d0 .and.
     .                abs(orb_el(2)-curr_orb_el(2)).lt.100.d0 .and.
     .                abs(orb_el(3)-curr_orb_el(3)).lt.100.d0 ) then
 
*                     Write out the current orbital elements, unless this
*                     is a duplicate
                      if( abs(jd-curr_jd).gt.0.001 ) then
                          write(200,200) (date(i),i=1,4), prn,
     .                                    curr_orb_el
                          curr_jd = jd
                          curr_prn = prn
                          do i = 1,9
                              curr_orb_el(i) = orb_el(i)
                          end do
                      end if
 200                  format(1x,i4,3I3,1x,a,6(1x,F13.3),1x,
     .                       f8.4,2(1x,F7.4))
                  else
*                     Write new values and update current values
                      write(200,200) (date(i),i=1,4), prn, orb_el
                      curr_jd = jd
                      curr_prn = prn
                      do i = 1,9
                          curr_orb_el(i) = orb_el(i)
                      end do
                  end if
              else
                  write(200,'(a)') line(1:max(1,trimlen(line)))
              end if
          end if
      end do
 
***** Thats all
      close(100)
      close(200)
      end
 
