CTITLE JD_TO_DECYRS
 
      subroutine JD_to_Decyrs( jdi, decyrs )
 
 
*     Routine to convert JD to decimal years (i.e. to 1986.45)
 
*   date(5) - date array needed by YMDHMS_to_JD
 
      integer*4 date(5)
 
*   decyrs  - Deciminal years values
*   jd      - Julian data
*   jdi     - Input julian date
*   secs    - Seconds tag needed for YMDHMS_to_JD
 
 
      real*8 decyrs, jd, jdi, secs, jde, num_days
 
***** Get years part of decyrs and JD at start of year
 
      call JD_to_YMDHMS ( jdi, date, secs )
 
*                             ! Jan
      date(2) = 1
*                             ! first
      date(3) = 1
*                             ! 0 hrs
      date(4) = 0
*                             ! 0 minutes
      date(5) = 0
 
      secs    = 0
*                                             ! JD at start of year
      call ymdhms_to_jd( date, secs, jd )

*     Get number of days in year
      date(1) = date(1) + 1
      call ymdhms_to_jd( date, secs, jde)
      num_days = jde - jd
      if( num_days.ne.365.d0 .and. num_days.ne.366.d0 ) then
          num_days = 365.d0
      end if
      date(1) = date(1) - 1
 
      decyrs = date(1) + (jdi-jd)/num_days
 
***** Thats all
      return
      end
 
