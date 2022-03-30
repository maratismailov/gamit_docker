CTITLE PTD_STATUS

      subroutine ptd_status( status, mpole, gamit_mod )

      implicit none

*     Routine to return the pole tide status (3 enties) and
*     mean pole model based on the gamit_mod flag (see globk_cntl.h)

* PASSED OUT
      logical status(3)  ! Three entries 
                         ! (1) status known or not (SINEX vs h-files)
                         ! (2) Solid Earth pole-fide applied
                         ! (3) Ocean pole-tide applied
      character*(*) mpole  ! Mean pole model (string to match names 
                         ! used in mean_pole.f

* PASSED IN
      integer*4 gamit_mod  ! GLOBK bit mapped word for GAMIT models applied
                           ! (see ggamit_mod in globk_cntl.h)

* LOCAL
      logical kbit      ! Test if bit is on.

****  Start
      status = .false.   ! Set all entries false first
      mpole = 'UNKN'     ! Initial entry

****  See if pole-tide is on (i.e., applied)
      if( kbit(gamit_mod,19) ) then  ! Model is applied and known
         status(1) = .true.    ! We know status
*        Work backwards to get latest model used for Solid-Earth
*        pole tide
         if( kbit(gamit_mod,26) ) then   ! IERS20 mean pole
            status(2) = .true.
            mpole = 'IERS20'
         elseif( kbit(gamit_mod,23) ) then   ! IERS10 mean pole
            status(2) = .true.
            mpole = 'IERS10'
         elseif( kbit(gamit_mod,21) ) then   ! IERS96 mean pole
            status(2) = .true.
            mpole = 'IERS96'
         else
            status(2) = .true.
            mpole = 'ZERO'     ! Very old h-files. Should not available
         endif

****     See of ocean pole tide applied
         if( kbit(gamit_mod,27) ) then   ! IERS20 mean pole
            status(3) = .true.    ! Should be same pole as SE-P-tide
         elseif( kbit(gamit_mod,25) ) then   ! IERS10 mean pole
            status(3) = .true.    ! Should be same pole as SE-P-tide
         elseif( kbit(gamit_mod,24) ) then   ! IERS96 mean pole
            status(3) = .true.    ! Should be same pole as SE-P-tide
         endif

***** Pole tide not applied.  See if we know status (check model bits)
      else
*       If model bits are on then we know status
        if( kbit(gamit_mod,26).or.kbit(gamit_mod,23).or.
     .      kbit(gamit_mod,21) ) then
           status(1) = .true.
        endif
      endif

****  Thats all 
      return
      end





