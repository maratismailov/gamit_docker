Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
c
      Subroutine DOUBLE(iphi,iobs,iuses)
c Root by Yehuda bock (1985), modified by D.N. Dong (1987)
c Reconstructured by Dong 03/16/91
c     constructs an independent set of double-difference operator
c     for a matrix of oneway phases observed simultaneously.
c
c     input:
c        iphi     flags good observations(1-good;0-bad)
c        nsite    number of total stations (solve.h)
c        nsat     number of total satellites (solve.h)
c     output:
c        D        operator matrix
c        ipntd    operator matrix row pointer
c                 it points to the column of the corresponding
c                 non-zero element in D.
c        irowd    operator matrix element per row counter
c          irowd(1)  = total number of non-zero element
c          irowd(2)  = 0 (always)
c          irowd(3+) = running count of the number of elements
c                      per row of D-matrix
c        iuses    record of oneway phases used in a double difference
c                 (single epoch record)
c        iobs     number of effective DD-observations
c
      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'
                 
      integer*4 iobs,istat,i,isum,j,ibase,isat,ibaset
      integer*4 icut,jsit,jsat,icnt,ib,ia,isat1,isat2,ll,iau
      integer*4 jj,jjj,istat1,kk,kkk,istat2,icut1,jsit1,jsit2
      integer*4 jsita,jsitb,ii,jsat2
  
      character*256 message
      logical pass
      integer ibonus(4)
      integer iphi(maxsit,maxsat),iefec(maxsit,maxsat),
     3          iuses(maxsit,maxsat),mapsat(maxsat),
     4          mapsit(maxsit)
c             
      logical debug/.false./

      call zero2i(maxsit,maxsat,nsite,nsat,iuses)
      call zero2i(maxsit,maxsat,nsite,nsat,iefec)
c
c     istat: number of effective stations
c     mapsit: mapping vector of stations
c     isat: number of effective satellites
c     mapsat: mapping vector of satellites
c
      istat = 0
      do 10 i = 1,nsite
         isum = 0
         do 5 j = 1,nsat
            if (iphi(i,j).eq.0) goto 5
            isum = isum+1
 5       continue
         if(isum.lt.1) go to 10
         istat = istat+1
         mapsit(istat) = i
         if(debug) print *,'DOUBLE i isum istat mapsit iphi'
     .                    ,i,isum,istat,mapsit(istat),(iphi(i,j),j=1,32)
 10   continue
c
      ibase = 0
      isat = 0
      do 30 i = 1,nsat
         isum = 0
         do 20 j = 1,nsite
            if(iphi(j,i).eq.1) isum = isum + 1
 20      continue
         if (isum.le.1) goto 30
         isat = isat+1
         mapsat(isat) = i
         if (isum.eq.istat.and.ibase.eq.0) ibase = isat
         if(debug) print *,'DOUBLE i isum isat mapsat ibase'
     .                    ,i,isum,isat,mapsat(isat),ibase 
 30   continue
c
c     skip when no double-difference combination
      if(istat.lt.2.and.isat.lt.2) go to 99
c
      ibaset = 0
      icut = 0
      do 40 j = 1,istat
         jsit = mapsit(j)
         isum = icut
         do 35 i = 1,isat
            jsat = mapsat(i)
            if (iphi(jsit,jsat).eq.0) goto 35
            icut = icut+1
            iefec(j,i) = icut
 35      continue
         if (icut-isum.eq.isat.and.ibaset.eq.0) ibaset = j
 40   continue
c
c==== search1 procedures :
c          find base satellite(ibase) observed by all live stations
c          if ibase does not exist, turn to search2 procedures.
c          sat1 = ibase
c          ... sat2 loop
c               ... stations loop
c                    find first site with observation to sat1 and sat2
c                    find second site with observation to sat1 and sat2
c                    fill DD-operator
c                    update istat1 with istat2
c
c     search1 always guarantes the independence of found DD-observation.
c     search1 is faster than search2, and, lucky enough, in most cases
c     search1 procedures are always available.
c
c     index:  jj (sit1,sat1)   jjj (sit1,sat2)
c     index:  kk (sit2,sat1)   kkk (sit2,sat2)
c
      icnt = 0
      iobs = 0
      if (ibase.eq.0) goto 70
      do 140 jsat2 = 1,isat
         if (jsat2.eq.ibase) goto 140
         ib = jsat2
         ia = ibase
         if (jsat2.lt.ibase) then
            ia = jsat2
            ib = ibase
         endif
         isat1 = mapsat(ia)
         isat2 = mapsat(ib)
         ll = 0
         do 150 jsit = 1,istat
            iau = iefec(jsit,jsat2)
            if (iau.eq.0) goto 150
            ll = ll+1
            if (ll.eq.1) then
               jj = iefec(jsit,ia)
               jjj = iefec(jsit,ib)
               istat1 = mapsit(jsit)
            endif
            if (ll.eq.2) then
               kk = iefec(jsit,ia)
               kkk = iefec(jsit,ib)
               istat2 = mapsit(jsit)
               call dbcnt(istat1,istat2,isat1,isat2,iuses,
     1                    icnt,jj,jjj,kk,kkk)
               iobs = iobs+1
               jj = kk
               jjj = kkk
               ll = 1
               istat1 = istat2
            endif
 150     continue
 140  continue
      goto 300
c
c==== search2 procedures :
c     if base station (ibaset) exists, same procedures as search1
c     except reverse site and sat.
c     if base station (ibaset) does not exist --
c      ... sit1 loop
c          ... sit2 loop
c               ... satellite loop
c                    find first sat with observation to sit1 and sit2
c                    find second sat with observation to sit1 and sit2
c                    check the independency of found DD observation
c                    if not independent, goto next satellite
c                    fill DD-operator
c                    update isat1 with isat2
c
 70   ib = 0
      icut1 = 0
      pass = .true.
      if (ibaset.eq.0) call chkmis(pass,istat,isat,iefec,mapsat,ibonus)
      do 80 jsit1 = 1,istat
         if (ibaset.eq.0) ib = jsit1
         if (ibaset.gt.0.and.jsit1.ne.ibaset) goto 80
         if (ibaset.eq.0.and.jsit1.eq.istat) goto 80
         do 60 jsit2 = ib+1,istat
            if (jsit2.eq.jsit1) goto 60
            jsita = jsit1
            jsitb = jsit2
            if (jsit1.gt.jsit2) jsita = jsit2
            if (jsit1.gt.jsit2) jsitb = jsit1
            istat1 = mapsit(jsita)
            istat2 = mapsit(jsitb)
            ll = 0
            do 50 jsat = 1,isat
               if (iefec(jsit1,jsat).eq.0.or.iefec(jsit2,jsat).eq.0)
     .         goto 50
               ll = ll+1
               if (ll.eq.1) then
                  isat1 = mapsat(jsat)
                  jj = iefec(jsita,jsat)
                  kk = iefec(jsitb,jsat)
               endif
               if (ll.eq.2) then
                  isat2 = mapsat(jsat)
c-----  check its independency (not regurously, depend on path)
                  if (ibaset.eq.0) then
                     if(iuses(istat1,isat1).ge.1.and.
     .                  iuses(istat1,isat2).ge.1.and.
     .                  iuses(istat2,isat1).ge.1.and.
     .                  iuses(istat2,isat2).ge.1) then
                        if (.not.pass.and.jsita.eq.ibonus(1).and.
     .                  jsitb.eq.ibonus(2).and.isat1.eq.ibonus(3)
     .                  .and.isat2.eq.ibonus(4)) then
c                           write (6,110) istat1,istat2,isat1,isat2
                           pass = .true.
                           goto 120
                        endif
                        ll = 1
                        goto 50
                     endif
                  endif
 120              jjj = iefec(jsita,jsat)
                  kkk = iefec(jsitb,jsat)
                  iobs = iobs+1
                  call dbcnt(istat1,istat2,isat1,isat2,iuses,
     1            icnt,jj,jjj,kk,kkk)
                  ll = 1
                  jj = jjj
                  kk = kkk
                  isat1 = isat2
               endif
 50         continue
 60      continue
         if (ibaset.eq.0.and.pass) call chkcut(nsite,nsat,iuses,icut1)
         if (ibaset.eq.0.and.pass.and.icut1.ge.icut) goto 300
 80   continue
c 110  format(2x,'You get a bonus DD-operator at site',2i3,' sat',2i3)
c
 300  continue
c check if there is a mismatch
      do 90 ii = 1,nsite
      do 90 jj = 1,nsat
         if(iphi(ii,jj).eq.1.and.iuses(ii,jj).eq.0) then
           write(message,88) ii,jj
 88        format (20x,'warning!'/2x,'observation between station (',
     .             i2,') and satellite (',i2,') has not been used')
           call report_stat('WARNING','SOLVE','double',' ',message,0)
         endif
         if(iphi(ii,jj).eq.0.and.iuses(ii,jj).ge.1) then
           write(message,89) ii,jj
 89        format (20x,'warning!'/2x,'observation between station (',
     .                  i2,') and satellite (',i2,') does not exist')
           call report_stat('WARNING','SOLVE','double',' ',message,0)
         endif
 90   continue

c  fill irowd-vector
      irowd(2) = 0
      icnt = 0
      do 190 i = 1,iobs
         icnt = icnt+4
         irowd(i+2) = icnt
 190  continue
      irowd(1) = icnt
c
 99   continue
c
      return
      end
c----------------------------------------------------------------------
      subroutine chkmis(pass,istat,isat,iefec,mapsat,ibonus)

      implicit none

      include '../includes/dimpar.h'

      logical pass,merge

      integer*4 iefec(maxsit,maxsat),igrp1(maxsat),igrp2(maxsat) 
      integer*4 ibonus(4),mapsat(maxsat) 
      integer*4 istat,isat,isit1,isit2,i,ll,mm,j,isum,isat2
      integer*4 isat1,iu


      isit1 = 0
      isit2 = 0
      do 10 i = 1,isat
         igrp2(i) = 0
 10   igrp1(i) = iefec(1,i)
      ll = 0
      do 20 i = 2,istat
         merge = .false.
         mm = 0
         do 30 j = 1,isat
            if (iefec(i,j).ge.1.and.igrp2(j).ge.1) mm = mm+1
 30      if (iefec(i,j).ge.1.and.igrp1(j).ge.1) merge = .true.
         if (.not.merge) then
            do 35 j = 1,isat
            if (iefec(i,j).ge.1) igrp2(j) = 1
 35         continue
         endif
         if (merge.and.mm.eq.0) then
            isum = 0
            do 45 j = 1,isat
            if (iefec(i,j).ge.1) igrp1(j) = 1
            if (igrp1(j).ge.1) isum = isum+1
 45         continue
            if (isum.ge.istat) goto 100
         endif
         if (merge.and.mm.gt.0) ll = ll+1
         if (ll.eq.1.and.merge) isit1 = i
         if (ll.eq.2.and.merge) isit2 = i
         if (ll.eq.2) then
            isat2 = 0
            isat1 = 0
            do 40 j = 1,isat
              if (igrp1(j).eq.0.and.igrp2(j).eq.0) goto 40
              if( isit1.le.0 ) call report_stat('FATAL','SOLVE','double'
     .                   ,' ','Zero subscript isit1--something wrong',0)
              if (iefec(isit1,j).ge.1.and.iefec(isit2,j).ge.1
     .            .and.igrp1(j).ge.1) isat1 = j  
              if( isit2.le.0 ) call report_stat('FATAL','SOLVE','double'
     .                   ,' ','Zero subscript isit2--something wrong',0)
              if (iefec(isit1,j).ge.1.and.iefec(isit2,j).ge.1
     .           .and.igrp2(j).ge.1) isat2 = j
              if (isat1.gt.0.and.isat2.gt.0) then
                  ibonus(1) = isit1
                  ibonus(2) = isit2
                  if (isat1.gt.isat2) then
                     iu = isat1
                     isat1 = isat2
                     isat2 = iu
                  endif
                  ibonus(3) = mapsat(isat1)
                  ibonus(4) = mapsat(isat2)
                  pass = .false.
                  goto 100
              endif
 40         continue
            ll = 1
         endif
 20   continue
 100  continue
      return
      end
c------------------------------------------------------------------
      subroutine chkcut(nstat,nsate,iuses,icut)

      implicit none

      include '../includes/dimpar.h'

      integer*4 iuses(maxsit,maxsat) 
      integer*4 nstat,nsate,i,j,icut

      icut = 0
      do 10 i = 1,nstat
         do 10 j = 1,nsate
         if (iuses(i,j).ge.1) icut = icut+1
 10   continue
      return
      end

