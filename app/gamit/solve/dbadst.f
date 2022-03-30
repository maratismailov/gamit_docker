c
      Subroutine DBADST(llbad,k,ngood,istat,baslen,ndex,
     .   iswtch,last_nonbias)

c     construct bad station baseline bias operator

      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'
 
      integer maxbas
      parameter (maxbas=maxsit*(maxsit-1)/2)

      integer istat(maxsit),ndex(maxbas*2),lstat(maxsit)
      real*8 baslen(maxbas),baslin

      integer llbad,ngood,iswtch,last_nonbias,kk0,ngood1
     .        ,ilen,ix1,ix2,ii1,ibs,i1,j1,j2,ii,jj,mm,k1,i,j,k

        kk0 = k
        mm = k*4
        ngood1 = ngood
        do 370 ibs = 1,llbad
            i1 = lbad(ibs)
            ngood1 = ngood1+1
            istat(ngood1) = i1
c---- calculate lengths of all baselines connected to the bad site
            jj = i1
            ix1 = 3*(jj-1)+1
            do 374 i = 1,ngood1-1
               lstat(i) = istat(i)
               if (lstat(i).lt.0) lstat(i) = -lstat(i)
               ii = lstat(i)
               ix2 = 3*(ii-1)+1
               baslen(i) = baslin(coords(ix1),coords(ix2))
 374        continue
c---- sort out all baselines in sequences of length.
      ilen = ngood1-1
      if (ilen.gt.1) call sortbl(ilen,1,baslen,lstat)
c
c---- find first live satellite for the bad site
        do 382 j = 1,nsat
           if (iuse(i1,j).gt.0) goto 384
 382    continue
 384    j1 = j
        if (j1.eq.nsat) goto 370
c---- find second live satellite for the bad site
        do 386 j = j1+1,nsat
           if (iuse(i1,j).eq.0) goto 386
           j2 = j
c---- search the nearest baseline with observation for the two
c----    satellites
           do 388 i = 1,ngood1-1
              ii1 = lstat(i)
              if (iuse(ii1,j1).eq.0.or.iuse(ii1,j2).eq.0) goto 388
c---- add a new independent baseline
              k1 = ilen*2
              ilen = ilen+1
              ix1 = min0(i1,ii1)
              ix2 = max0(i1,ii1)
              ndex(k1+1) = ix1
              ndex(k1+2) = ix2
              k = k+1
c---- fill d-operator
              call filld(ix1,ix2,j1,j2,iswtch,mm,ilen)
              goto 386
 388       continue
 386    continue
 370  continue
c---- range of bias parameters connected to bad stations
      lbfre(2) = last_nonbias+kk0+1
      lbfre(3) = last_nonbias+k
      if (l2flag.ge.2) then
          lbfre(4) = last_nonbias+k+kk0+1
          lbfre(5) = last_nonbias+k*2
      endif
c
      continue

      return
      end

