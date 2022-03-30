      subroutine get_iut1pol( iut1pol )

      implicit none

*     Routine to return the value of iut1pol bit mapped
*     word from sestbl. or to return a default value which
*     can be set in one location (here).  Herring 200504.

* PASSED (Value returned from call)
      integer*4  iut1pol  ! Bit mapped short period EOP value
                          ! Bit     Meaning
                          !   1   1 Apply to X/Y pole
                          !   2   2 Apply to UT1
                          !   3   4 Old Ray model
                          !   4   8 IERS 2010 standard (May 2020 default)
                          !   5  16 Gipson model (coded by Lei Wang)
                          !   6  32 Desai and Sibiuos (REPRO3 standard)

* LOCAL Variables
      integer*4 ill    ! Error return from rdsest
      integer*4 unit   ! Unit for reading sestbl.  31 seems available
      integer*4 jerr   ! Error return from reading iut1pol.
      integer*4 len    ! Length of runstring
     
      logical reqd   ! Set to default so that sestbl. is not needed
      logical kbit   ! Test bit status

      character*4  aerp ! String returned from sestbl.  

      character*80 prog_name  ! Extracted with rcpar(0,prog_name)

      logical known / .false. /  ! Set true once we know value so won't
                                 ! keep reading sestbl.
      integer*4 kut1pol          ! Known value (restore is known)
      integer*4 rcpar            ! Read runstring function.

      save known, kut1pol

****  See if we know the value from an earilier call
      if( known ) then
         iut1pol = kut1pol
         RETURN
      endif

****  Set default value in case value not set/
      iut1pol = 11 + 16  ! IERS2010 + UT LIBRATION terms
 
****  Open and read sestbl. for the Earth Rotation entey
      reqd = .false.
      unit = 31 
      open( unit, file='sestbl.',status='old',iostat=jerr)
      if( jerr.ne.0 ) then
         len = rcpar(0,prog_name)
         call report_stat('WARNING',prog_name,'orbits/trot','sestbl.',
     .       'Errors in opening session table,set to default',jerr)
      else
         call rdsest( 14,'Earth Rotation',4,aerp,unit,reqd,ill )
         if ( ill.ne.0 ) then
              len = rcpar(0,prog_name)
              call report_stat('WARNING',prog_name,'Earth Rotation',
     .          'sestbl.','Errors in sestbl entries',ill)
         endif
         if(aerp(1:1).ne.' ') then
            read(aerp,*,iostat=jerr) iut1pol
            if ((kbit(iut1pol,1)).or.(kbit(iut1pol,2))) then
* MOD TAH 200504: Added Desai and Sibious model (bit 7).
*              Adopt the latest model and set all lower bits to zero.
               if( kbit(iut1pol,7)) then   !Desai and Sibous
                   call sbit(iut1pol,3,0)
                   call sbit(iut1pol,4,0)
                   call sbit(iut1pol,6,0)
               elseif( kbit(iut1pol,6) ) then
                   call sbit(iut1pol,3,0)
                   call sbit(iut1pol,4,0)
               elseif( kbit(iut1pol,4) ) then
                   call sbit(iut1pol,3,0)
               else   ! Old Ray model; force to IERS2010
                   call report_stat('WARNING','ORBITS','trot',' ',
     .             'old diurnal EOP model--override with IERS10 model',
     .             0)
* MOD    TAH 200504: Use sbit to set IERS model as default.  This might
*                    change later.
                     call sbit(iut1pol,4,1)
               endif
            else
               iut1pol = 0   ! No model to apply
            endif
         endif
      endif

*     Save if there is another call.
      known = .true.
      kut1pol = iut1pol
      close(unit)

***** Thats all
      return 
      end



