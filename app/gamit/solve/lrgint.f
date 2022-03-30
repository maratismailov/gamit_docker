Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
C
      Subroutine LRGINT( PHI,IPHI,FIRST,IERTAG )

      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'

C     re-written TAH 970206.  Previous version seemed to miss
C     removing some large biases for reasons that are not clear.
      
      integer*4 ISATC(MAXSAT),ISTATC(MAXSIT),IPHI(MAXSIT,MAXSAT),
     .          IERTAG(MAXSIT,MAXSAT)
     
      real*8 PHI(MAXSIT,MAXSAT), FIRST(MAXSIT,MAXSAT), 
     .       small  

      
*  ref_sv - Satellite number of the satellite to be used at the
*     reference for computing the large offset in the phase.

      integer*4 ref_sv, i,j,k 

*  Found - Logical to indicate we have found a appropriate reference
*     satellite.
      logical found

      small = 1.0d-12
C
C     COUNT NUMBER STATIONS, SATELLITES AND OBSERVATIONS
      DO 34 I = 1,nsite
   34 ISTATC(I) = nsat
      DO 35 I = 1,nsat
   35 ISATC(I) = nsite

*     Loop over stations and satellites
      DO 22 I = 1,nsat
         DO 22 J = 1,nsite

*           if we see a bias flag set first to zero so that
*           new value will be comuted later in this routine.

            IF(IERTAG(J,I).EQ.10) FIRST(J,I) = 0.0D0
*           Decrement number of satellites visible with this
*           station
            IF(IPHI(J,I).EQ.0) ISTATC(J) = ISTATC(J)-1
 22   continue

*     Decrement number of stations that can be seen from
*     satellite. 
      DO 23 I = 1,nsite
      DO 23 J = 1,nsat
         IF(IPHI(I,J).EQ.0) ISATC(J) = ISATC(J)-1
 23   continue

*     If we have only one satellite visible, ignore the
*     one-way data since it can never be used in a double
*     difference.  (With phase clocks applied in autcln, the
*     LC phase omc is likely to be zero anyway. L1 and L2
*     phases separately will not be zero.)
      DO 3 I = 1,nsat
         IF(ISATC(I).LE.1) THEN
            DO 4 J = 1,nsite
               IPHI(J,I) = 0
               PHI(J,I) = 0.D0
 4          CONTINUE
         ENDIF
 3    CONTINUE
C
C
C  TAKE OFF LARGE INTEGER
      DO I = 1,nsite

*        Make sure that we can see more than one satellite.
         IF(ISTATC(I).gt.1) then     
            DO J = 1,nsat
*              If we have a measurement at this epoch then
*              continue.  Either remove the current offset
*              or estimate a new one based on current data.
               if( iphi(i,j).eq.1 ) then

*                  We have a measurement, see if we should
*                  apply old offset or estimate a new one.
                   if( first(i,j).ne.0 .and. iertag(i,j).ne.10 ) then
* MOD TAH 031204: Donot remove any offset so that autcln widelanes remain correct
* MOD RWK 031209: Make this contingent on LC_AUTCLN                             
                     if( l2flag.ne.4 ) then
                        phi(i,j) = phi(i,j) - first(i,j)
                     endif
                   else

*                      We need to estimate a new offset, see if
*                      we can find another satellite to use.
*                      Ref_sv will be the reference satellite
                       k = 0
                       ref_sv = 0
                       found = .false.
                       do while ( .not.found )
                           k = k + 1
*                          if we are at the last satellite, so
*                          set found true for exit (still check
*                          to see if last satellite is acceptable).
                           if( k.ge.nsat ) then
                               found = .true.
                           end if
                           if( k.ne.j .and. first(i,k).ne.0 .and.
     .                         iphi(i,k).eq.1 ) then
                               ref_sv = k
                               found = .true.
                           end if
                       end do

*                      Now see what we have:
*                      if ref_sv is zero, then use current phi as offset
                       if( ref_sv.eq.0 ) then
                           first(i,j) = nint(phi(i,j))

                       else if( ref_sv.lt.j ) then
*                          if ref_sv is before current satellite then we
*                          have already applied the offset to the phase
*                          at the reference site and the difference we
*                          want is the difference in the phase.
                           first(i,j) = nint(phi(i,j) - phi(i,ref_sv))

                       else 
*                          The only remaining case is ref_sv after j, so
*                          do as above but remove the offset from the
*                          reference site first.
                           first(i,j) = nint(phi(i,j) - 
     .                            (phi(i,ref_sv) - first(i,ref_sv)) )
                       end if

*                      For autcln with phase clocks applied first should
*                      small, so set value to 1 cycle if less than 100
*                      (These are small enough not affect RMS calcualtions
*                      and since all sites and satellites get the same correction
*                      will not affect the bias parameter estimates).
                       if( abs(first(i,j)).lt.100.d0 ) first(i,j) = 1.d0  
* MOD TAH 031204: Again do not remove offset
* MOD RWK 031209: Make this contingent on LC_AUTCLN 
                      if( l2flag.ne.4 ) then
                          phi(i,j) = phi(i,j) - first(i,j)
                       endif
                    end if
                end if
*                       !  Looping over satellites at this site
            end do
*                       !  If we have have more than one satllite observed
         end if
*                       !  Looping over the satellites.
      end do

C     Thats all
      RETURN
      END
