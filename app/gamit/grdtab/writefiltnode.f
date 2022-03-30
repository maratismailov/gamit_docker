      subroutine writefiltnode(maxepoch,nepoch,uout,displfilt,leap_yr2
     .                         ,irecbox,nglat,nglon)

c  subroutine to write out the filtered data for a single node for a whole year.
c  We can do this node by node because the grid files are direct-access, so by
c  knowing the starting record for the node we can then just increment by the
c  number of records/epoch for the whole grid and thus write out each epoch for
c  each node, one node at a time.
c
c  P. Tregoning
c  8 October 2008

      implicit none
 
      integer maxepoch,nepoch,uout,i,j,irecbox
     .       ,start_epoch,end_epoch,recl,nglat,nglon
      real*8 displfilt(maxepoch,3)
      real*4  values(3)
      logical leap_yr2

c  regardless of whether it is a leap year or not, the start epoch is
c     line 51*4 + 1 in the filtered array - because we extracted 50 days
c     from the previous year.
      start_epoch = 205
      if (leap_yr2) then
        end_epoch = 366 * 4 + start_epoch 
      else
        end_epoch = 365 * 4 + start_epoch 
      endif

c PT090303: we need at least one more epoch so that the last day of the
c           year can be processed (grdtab reads through to the 0600 epoch
c           of the next day). I've increased the header to 002 of the next
c           year, so do the same here.
      end_epoch = end_epoch + 4

c
c  PT081008: just output to the screen the raw and the filtered time series
c            so that I can plot them and see whether the filter is working
c      print*,'writefiltnode: filtered time series'
c      do i=start_epoch,end_epoch
c        write(*,'(7f10.4)')float(i)/4.0 - 50.25,(displfilt(i,j),j=1,3)
c      enddo

c DEBUG, see if we can read from the output file
c      read(uout,rec=2)values(1)
c      print*,'interval from output file header is',values(1)

c  now we simply loop through the whole year and write out the values into
c  the appropriate records for this node. The first record is going to
c  be irecbox + 3 header lines 
c      recl = irecbox + 3   
c PT/RWK 100902: irecbox is now zero for 90N, so need to add one more
      recl = irecbox + 4 
      do i = start_epoch,end_epoch
        values(1) = displfilt(i,1)
        values(2) = displfilt(i,2)
        values(3) = displfilt(i,3)
c       print*,recl,values
c        recl = 4
        write(uout,rec=recl)values
        recl=recl+nglat*nglon
      enddo

c      stop ' stop in writefiltnode'

      return
      end
