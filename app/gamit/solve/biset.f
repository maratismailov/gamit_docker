Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
      subroutine biset(biskut,nchi,nfix,minchi,dsave)
c
c check chi-square contrast and locate significant biases for bias-searching.
c
c biscut : contrast factor
c nchi : bias number in comparison (input)
c        significant bias number (output)

      implicit none

      include '../includes/dimpar.h'

c*     integer nfix(nchi),nchi
      integer nfix(maxbis),nchi,nbad,is,i,j
      
      real*8 minchi(10),biskut,dsave(10,20),dtemp

      logical bad
c
      bad=.false.
c first see if any bias had enough contrast.
c
      do 100 i=2,10
         is=i
         if(minchi(i).ge.biskut) go to 110
100   continue
c
c none had enough contrast
      bad=.true.
c
110   continue
c
      continue
      nbad=0
c
c modify criterian : same integer from 1 to is-1 instead of is
c           -dnd- 871026
c
      do 260 i=1,nchi
         if(bad) go to 255
         dtemp=dsave(1,i)
         do 250 j=2,is-1
            if(dtemp.ne.dsave(j,i)) go to 255
250      continue
c
c this is a good one
         go to 260
c
c bad ones here. free it and increment counter
c
255      continue
         nbad=nbad+1
         nfix(i)=11111
260   continue
c
c in bias-searching, BISET will be call again until nchi=0.
c nchi=0 means no biases were fixed
c
      nchi=nchi-nbad
c
      return
      end
