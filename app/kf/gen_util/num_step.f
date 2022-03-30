ctitle
 
      subroutine num_step( data_type, nstep )

      implicit none 
 
c
c     routine to compute the number of kalfile records for each observation
c     There is one record for each data type used.
c
c     Only one delay type can be used at one time (currently) even
c     though there are three to chose from.
c
c Variables
c ---------
c data_type -- word to indicate data types being used. Each bit is set
c     for each type
c nstep -- the number of kalfile records per observation
c kbit  -- solve logical function
c
      integer*4 data_type, nstep
 
c
      logical kbit
 
c
c.... set up nstep by checking data types used.  Increment once for
c     the delay being used and once if rates are used. (see BASE_DATA)
      nstep = 0
*                                             ! Group delay
      if( kbit(data_type,1) ) then
          nstep = nstep + 1
*                                         ! Phase delay
      else if( kbit(data_type,2) ) then
          nstep = nstep + 1
*                                         ! SB delay
      else if( kbit(data_type,3) ) then
          nstep = nstep + 1
      end if
 
*     See if rate is being used.
      if( kbit(data_type,4) ) then
          nstep = nstep + 1
      end if
c
      return
      end
 
c......................................................................
