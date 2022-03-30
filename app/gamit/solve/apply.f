c
      Subroutine APPLY
c
c apply atmospheric constraints to normal matrix
c
      implicit none

      include '../includes/dimpar.h'
      include 'solve.h' 
      include 'parameters.h'

      integer*4 condim,jstat,jparm,ki,kl,iel1,ielem1
     .        , iel2,inorm1,ielem2,inorm2,jel,inorm,iat
     .        , ka(maxsit)

      parameter (condim = maxsit * (maxsit + 1) / 2)
      real*8 outcon(condim)
c
c compute atmospheric constraint matrix
c
      call delcon( coords,nsite,outcon )
c
c determine atmospheric parameter slots
      do 5 jstat=1,nsite
        do 6 jparm=1,npartm(jstat)
          ki=islot2(jparm,jstat)
          kl=islot1(ki)
c check for atmosphere parameter
c  old single zenith-delay code
c          if(kl.ge.301.and.kl.le.400) then
           if(kl.ge.11501.and.kl.le.14000) then
           ka(jstat)=ki
           go to 5
          end if
    6   continue
    5 continue
c
c add constraints to normal matrix
      iel1=0
      do 10 ielem1=1,nsite
        iel2=0
        if (iusest(ielem1).eq.0) goto 10
        iel1=iel1+1
        inorm1=ka(ielem1)
        do 20 ielem2=1,ielem1
          if (iusest(ielem2).eq.0) goto 20
          iel2=iel2+1
          iat=jel(iel1,iel2)
          inorm2=ka(ielem2)
          inorm=jel(inorm1,inorm2)
          a(inorm)=a(inorm)+outcon(iat)
          alc(inorm)=alc(inorm)+outcon(iat)
   20   continue
   10 continue
c
      return
      end
