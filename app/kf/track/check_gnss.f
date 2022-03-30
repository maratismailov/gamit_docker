CTITLE CHECK_GNSS 

      subroutine check_gnss( prn, tr_gnss, OK )

      implicit none

*     Routine to set OK false the prn is from a system nt being
*     processed.

* PASSED
      integer*2 prn   ! PRN with multiples of 100 added depending
                      ! on constellation
      logical OK      ! Changed to false if not beinf procesed

      character*(*) tr_gnss  ! List of systems being processed in track

* LOCAL

      integer*4 sys    ! System number 

****  Get system
      sys = int(prn/100) + 1
      if( sys.eq.1 .and. index(tr_gnss,'G').eq.0 ) OK = .false. 
      if( sys.eq.2 .and. index(tr_gnss,'R').eq.0 ) OK = .false.
      if( sys.eq.3 .and. index(tr_gnss,'E').eq.0 ) OK = .false.
      if( sys.eq.4 .and. index(tr_gnss,'C').eq.0 ) OK = .false.
      if( sys.ge.5 ) OK = .false.

      RETURN
      end



