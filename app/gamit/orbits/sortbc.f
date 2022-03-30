      subroutine sortbc( nbrd,nprn,iwk )

c     Sort the elements in the Broadcast Ephemeris arrays
c     according to the pointers in indx.
c     R. W. King    19 November 1991

      include '../includes/dimpar.h'
      include '../includes/orbits.h'

      integer*4 nbrd,indx(maxbrd),nprn(maxbrd),iwk(maxbrd)
     .        , iwksp(maxbrd),i,j

      real*8 ephem(16,maxbrd),wksp(maxbrd)

      common/ephemcom/ephem

c        Convert GPS week, seconds of week into a single argument
c        (precision is not important here)

      do i=1,nbrd
        wksp(i) = iwk(i)*604800.d0 + ephem(1,i)
      enddo

c        Make an index table sorting by time

      call indexx(nbrd,wksp,indx)

c        Reorder all of the arrays using the index table

      do j=1,nbrd
         iwksp(j)= iwk(j)
      enddo
      do j=1,nbrd
         iwk(j) = iwksp(indx(j))
      enddo

      do j=1,nbrd
         iwksp(j)= nprn(j)
      enddo
      do j=1,nbrd
         nprn(j) = iwksp(indx(j))
      enddo

      do i=1,16
         do j=1,nbrd
            wksp(j)= ephem(i,j)
         enddo
         do j=1,nbrd
            ephem(i,j) = wksp(indx(j))
         enddo
      enddo

      return
      end
