      subroutine dcheck(iphi,iphibk,iertag,ipflg,ieflg)
c
c  determine if operator has changed from previous epoch
c  iphi   : current phase matrix (0:no obs, 1: good obs)
c  iphibk : previous phase matrix
c  istat  : total number of stations
c  isat   : total number of satellites
c  ipflg = 0 : good observations same as before
c  ipflg = 1 : need a new operator matrix
c  ieflg = 0 : without bias flag
c  ieflg = 1 : with bias flag
c              
      implicit none

      include '../includes/dimpar.h'     
      include 'solve.h' 

      integer iphi(maxsit,maxsat),iphibk(maxsit,maxsat)
     .      , iertag(maxsit,maxsat),ipflg,ieflg,i,j
c
      ipflg=0
      ieflg=0
      do 10 j=1,nsat
         do 20 i=1,nsite
            if(iphi(i,j).ne.iphibk(i,j)) ipflg=1
            if(iertag(i,j).eq.10.and.iphi(i,j).eq.1) ieflg=1
            if (ipflg.eq.1.and.ieflg.eq.1) goto 40
   20    continue
   10 continue

   40 continue
c
      return
      end
