      subroutine wfica(lun,itype,fa,ia,ca,nf,ni,nc)
c
c     write a file in fica (float, int, char (ascii)) format
c     modified from judah levine's wtfica by kf mar 5, 87
c     modified by k feigl, jan 17, 1988 to deal with
c     variable logical unit number (lun)
c
c     this subroutine writes an fica record.  itype is the
c     type, nf,ni and nc are the number of floating variables,
c     integer variables and character variables, respectively.
c
c     fa, ia and ca are the arrays of floating variables, integer
c     variables and character variables, respectively.
c
      implicit none

      integer*4 lun,i
      integer*4 itype,nf,ni,nc
      integer*4 ia(ni)
      real*8 fa(nf)
      character*8 ca(nc)
c
c
c
      write(lun,1) itype,nf,ni,nc
    1 format('BLK  ',4i5)
      if(nf .gt. 0) write(lun,2)(fa(i),i=1,nf)
    2 format(4(1pd20.13))
      if(ni .gt. 0) write(lun,3)(ia(i),i=1,ni)
    3 format(6i12)
      if(nc .gt. 0) write(lun,4)(ca(i),i=1,nc)
    4 format(10a8)

      return
      end
