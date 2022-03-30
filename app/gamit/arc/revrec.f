Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
      subroutine revrec

c     reverse the record order in scratch file on unit 3 so as to
c     reverse the output of a backward integration
c     Rick Abbot - March 1985

      implicit none

      include '../includes/dimpar.h'
      include '../includes/arc.h'

      character*1 upperc

      integer*4 mdum,ndum,i,k,l

      real*8 tdum
      real*8 outvec(maxyt2)

      do  i=1,nrec

        backspace 3
        if (kount.eq.0) go to 600   
        if( upperc(apar).eq."V" ) then    
          read (3) (outvec(k),k=1,neq)
          write(isunit)  (outvec(k),k=1,neq)
        else
          read (3) (outvec(k),k=1,neq/2)
          write (isunit) (outvec(k),k=1,neq/2)
        endif
        goto 1500

c     old way !!!
c      read (3)
c     $          (y(l,4),l=1,3),  (y(l,4),l=7,9),  (y(l,4),l=13,15),
c     $          (y(l,4),l=19,21),(y(l,4),l=25,27),(y(l,4),l=31,33),
c     $          (y(l,4),l=37,39),(y(l,4),l=43,45),(y(l,4),l=49,51),
c     $          (y(l,4),l=55,57)
c      write (isunit)
c     $          (y(l,4),l=1,3),  (y(l,4),l=7,9),  (y(l,4),l=13,15),
c     $          (y(l,4),l=19,21),(y(l,4),l=25,27),(y(l,4),l=31,33),
c     $          (y(l,4),l=37,39),(y(l,4),l=43,45),(y(l,4),l=49,51),
c     $          (y(l,4),l=55,57)
c      go to 1500

  600   read (3)       (y(l,4),l=1,3)
        write (isunit) (y(l,4),l=1,3)
 1500   backspace 3

      enddo
c
      return
      end
