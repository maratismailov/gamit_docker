
CTITLE COMPRESS_site

      subroutine compress_site( num, names, full, obs )

*     This routine will remove from the list of sites those
*     which have no observations.

*   num     - NUmber in list
*   obs(num)    - Number of obsevation per sites

      integer*4 num, obs(num)

*   names(num)  - NAme of the sites
*   full(num)   - Full names of sites

      character*(*) names(num), full(num)

* LOCAL VARIBALES

*   i,j     - Loop counters

      integer*4 i,j

***** First get the total number of observations.  If it zero do
*     nothing since this is an old format h-file.

      j = 0

      do i = 1, num
         j = j + obs(i)
      end do

*     If there is no observations, then just return.

      if( j.eq.0 ) RETURN

****  Scan over list of satelites finding those with no observations

      i = 0
      do while (i.lt.num )
          i = i + 1
*                                 ! Remove from list
        if( obs(i).eq.0 ) then
              do j = i, num-1
                  obs(j) = obs(j+1)
                  names(j) = names(j+1)
                  full(j) = full(j+1)
              end do

*             Decrement number and pointer
              num = num - 1
              i = i - 1
          end if
      end do

****  Thats all
      return
      end



