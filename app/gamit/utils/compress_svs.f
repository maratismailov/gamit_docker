CTITLE COMPRESS_SVS

      subroutine compress_svs( num, names, obs )

*     This routine will remove from the list of satellites those
*     which have no observations.

*   num     - NUmber in list
*   obs(num)    - Number of obsevation per satellites

      integer*4 num, obs(num)

*   names(num)  - NAme of the satellites

      character*(*) names(num)

* LOCAL VARIBALES

*   i,j     - Loop counters

      integer*4 i,j

****  Scan over list of satelites finding those with no observations

      i = 0
      do while (i.lt.num )
          i = i + 1
*                                 ! Remove from list
        if( obs(i).eq.0 ) then
              do j = i, num-1
                  obs(j) = obs(j+1)
                  names(j) = names(j+1)
              end do

*             Decrement number and pointer
              num = num - 1
              i = i - 1
          end if
      end do

****  Thats all
      return
      end



