 
      LOGICAL FUNCTION ifmar(mar, num)
 

      implicit none
c
c     Routine to check if markov process in effect
c
c Variables
c ---------
c mar -- the markov statistics variable
c num -- number of values to be checked.
c ifmar -- logical set if markov process
c
      real*4 mar(1)
 
c
      integer*4 num, i
 
c
c.... Loop over the markov statitics
      ifmar = .false.
      do i = 1, num
         if( mar(i).ne.0 ) ifmar = .true.
      end do
c
      return
      end
 
c......................................................................
