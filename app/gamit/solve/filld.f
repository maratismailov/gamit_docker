      Subroutine FILLD ( isit1,isit2,isat1,isat2,mrowd,mm,ilen)  

c     Danan Dong 1997-88   

c     Fill the D-operator arrays for mapping double-difference biases
              
c     Input:
c        isit1 isit2 isat1 isat2 : site, sat indices 
c        ilen                    : baseline index ordered by lengths

c     Output:
c        ipntd(maxdd) in common /point1/ : indices of mapping operator
c        irowdt(maxnd)   common /point4/ : indices of mapping operator rows transposed
c        d(maxdd)        common /dmat/   : mapping operator elements
c        mrowd                           : row index of d-operator
c        mm                              : sequential index of d-operator elements

      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'

      integer*4 isit1,isit2,isat1,isat2,mrowd,mm
      integer*4 ilen,id1,id2,iau

      logical debug/.false./

      id1 = (isit1-1)*nsat
      id2 = (isit2-1)*nsat
      iau = isigma(id1+isat1)
      ipntd(mm+1) = iau
      d(mm+1) = -1.0d0
      irowdt(iau+2) = irowdt(iau+2)+1   
c      print *,'isit1 isit2 isat1 isat2 mm iau irowdt '
c     .       , isit1,isit2,isat1,isat2,mm+1,iau,irowdt(iau+2)
      iau = isigma(id2+isat1)
      ipntd(mm+2) = iau
      d(mm+2) = 1.0d0
      irowdt(iau+2) = irowdt(iau+2)+1 
c      print *,'isit1 isit2 isat1 isat2 mm,iau irowdt '
c     .       , isit1,isit2,isat1,isat2,mm+2,iau,irowdt(iau+2)
      iau = isigma(id1+isat2)
      ipntd(mm+3) = iau
      d(mm+3) = 1.0d0
      irowdt(iau+2) = irowdt(iau+2)+1   
c      print *,'isit1 isit2 isat1 isat2 mm,iau irowdt '
c     .       , isit1,isit2,isat1,isat2,mm+3,iau,irowdt(iau+2)
      iau = isigma(id2+isat2)
      ipntd(mm+4) = iau
      d(mm+4) = -1.0d0
      irowdt(iau+2) = irowdt(iau+2)+1  
c      print *,'isit1 isit2 isat1 isat2 mm,iau irowdt '
c     .       , isit1,isit2,isat1,isat2,mm+4,iau,irowdt(iau+2)
      mm = mm+4
      mrowd = mrowd+1
      ipntdt(mrowd-2) = ilen
      irowd(mrowd) = mm

      if(debug) print *
     .  ,'FILLD mrowd mm ipntd mrowd ipntdt irowd '
     .         ,mrowd,mm,ipntd(mrowd+4),ipntdt(mrowd-2),irowd(mrowd)

      return
      end
