c
      subroutine tapdrv(iunit,istart,iend,ndata,array,mode)
c
c     1-d real array load or unload
c     mode = 1: load, 1-d dumpping
c     mode = 2: load, lower-triangle recording
c     mode = 3: load, rectanguler recording(istart=row,iend=column)
c     mode = 4: unload, 1-d reading
c     mode = 5: unload, lower-triangle reading
c     mode = 6: unload, rectanguler reading(istart=row,iend=column)
c
      implicit none

      integer iunit,istart,iend,ndata,mode,ij,i,k

      real*8  array(ndata)
c
      if (mode.eq.1) then
         write (iunit) (array(i),i = istart,iend)
      endif
c
      if (mode.eq.2) then
         ij = istart*(istart-1)/2
         do 10 i = istart,iend
            write (iunit) (array(ij+k),k = 1,i)
            ij = ij+i
 10      continue
      endif
c
      if (mode.eq.3) then
         ij = 0
         do 30 i = 1,istart
            write (iunit) (array(ij+k),k = 1,iend)
            ij = ij+iend
 30      continue
      endif
c
      if (mode.eq.4) then
         read (iunit) (array(i),i = istart,iend)
      endif
c
      if (mode.eq.5) then
         ij = istart*(istart-1)/2
         do 20 i = istart,iend
            read (iunit) (array(ij+k),k = 1,i)
            ij = ij+i
 20      continue
      endif
c
      if (mode.eq.6) then
         ij = 0
         do 50 i = 1,istart
            read (iunit) (array(ij+k),k = 1,iend)
            ij = ij+iend
 50      continue
      endif
c
      return
      end


