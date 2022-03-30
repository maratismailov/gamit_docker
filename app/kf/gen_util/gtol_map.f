CTITLE GTOL_MAP

      integer*4 function gtol_map( ig, ltog, gnum )

      implicit none

*     Function to return the local number of site/satellite/source
*     given it global number.  If the site/satellite/source is not
*     in local solution, -1 is returned.

* PASSED 
      integer*4 ig     ! Global entry number
     .,         gnum   ! Number of global entries
     .,         ltog(gnum)  ! Mapping from local number to global value
                       ! (dimension must be at least gnum in size)

* LOCAL 
      integer*4 j      ! Loop counter

*     Scan the ltog array to see if ig is present, if it is return value
      gtol_map = -1    ! Not-found indicator

      do j = 1, gnum
         if( ltog(j).eq.ig ) then   ! We found entry to return
            gtol_map = j
            RETURN
         endif
      enddo

*     Thats all.  Here is value not found
      return
      end

  
