CTTITLE DECODE_COMOPT

      subroutine decode_comopt(line, comopt, updated) 

      implicit none

*     Routine to remove from the start of a command line, line,
*     any comopt strings passed by the user.  Line is returned
*     with the comopt strings removed.  Updated is returned true
*     if line modified.
*
* MOD TAH 200602: Added new LAB<str> decode option in which the
*     <LAB> string in a line will be replaced by str e.g.
*     With option LABWkRD the command line
*      bak_file <LAB>_<WEEK>.bak is replaced with
*      bak_file WkRD_<WEEK>.bak

      character*(*) line    ! Input line
      character*(*) comopt  ! Options line in form of
                            ! OPT1+OPT2+OPT3 any of which will
                            ! be processed
      logical updated       ! Set true if line updated

* PARAMETERS
      integer*4 max_opt     ! Maximum number of options
      parameter ( max_opt = 16 ) 

* LOCAL VARIABLES
      character*512 newline ! Updated line
      character*64  opts(max_opt)  ! Each option in line

      integer*4 num_opt   ! Number of options in comopt
     .,         trimlen   ! Length of string
     .,         indx      ! Position in string.
     .,         jndx      ! Second position
     .,         k         ! loop counter
     .,         lopt      ! Length of options components

      logical   done      ! When split string complere
      logical   LAB_rpd   ! See true when missing <LAB> option found

      data LAB_rpd / .false. /

      save  LAB_rpd

****  See what options we have
      updated = .false.

      if( trimlen(comopt).eq.0 ) RETURN   ! Nothing to do

*     Start splitting comopt up
      num_opt = 0
      jndx = 1     ! Start of string
      indx = jndx
      done = .false.
      do while ( .not.done  )
         k = index(comopt(jndx+1:),'+')
         if( k.gt.0 ) then
             jndx = jndx + k
         else
             jndx = trimlen(comopt)+1
             done = .true.
         end if
         num_opt = num_opt + 1
         if( num_opt.gt.max_opt) then
             call report_stat('FATAL','GLOBK','decode_comopt',
     .           comopt, 'Too many options. Max ',max_opt)
         endif

         opts(num_opt) = comopt(indx:jndx-1)
         indx = jndx + 1
       end do

****   Now see if we can find 
       do k = 1, num_opt
          lopt = trimlen(opts(k))
* MOD TAH 200602: See of the start of lopt is LAB.  If so
*         do string replacement instead of removing lopt
          if( opts(k)(1:3).eq.'LAB' ) then
             call sub_char(line,'<LAB>',opts(k)(4:lopt))
          elseif( line(1:lopt).eq.opts(k)(1:lopt) ) then
              updated = .true.
              newline = line
              line = newline(lopt+1:)
          end if
       end do

* MOD TAH 200602 Add check to report <LAB> string still in string
       if( index(line,'<LAB>').gt.0 .and. .not. lab_rpd ) then
          call report_stat('WARNING','GLOBK','decode_compopt',line,
     .                   '<LAB> label but no LABstr option',0)
          lab_rpd = .true.
       endif


*****  Thats all
       return
       end




