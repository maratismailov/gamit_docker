c
      Subroutine FILOBS(gearf,iphi,phi,phi2,iones,al1,al2)
c
c     fill observation matrix (by station)
c     form combined observables

      implicit none

      include '../includes/dimpar.h'                             
      include 'solve.h'

      real*8 al1(maxobs),al2(maxobs)
      real*8 phi(maxsit,maxsat),phi2(maxsit,maxsat) 
      integer iphi(maxsit,maxsat)    
      real*8 gearf,g1,g2,phitm1,phitm2
      integer*4 iones,i,j
c     for debug
      integer*4 ii 
 
      logical debug/.false./

      iones = 0

      g1 = gearf/(1.d0-gearf*gearf)
      g2 = 1.d0/(2.d0*gearf)

      do 70 i = 1,nsite
         do 60 j = 1,nsat
            if(iphi(i,j).eq.0) go to 60
            iones = iones+1
            phitm1 = phi(i,j)
            phitm2 = phi2(i,j)
            al1(iones) = phitm1-g1*(phitm2-gearf*phitm1)
            if (l2flag.gt.0)
     .         al2(iones) = phitm1+g2*(phitm2-gearf*phitm1)
         if(debug) print *,'FILOBS sit sat iones al1  '
     .         ,i,j,iones,(al2(ii),ii=1,iones)
 60      continue
 70   continue
c
      return
      end

