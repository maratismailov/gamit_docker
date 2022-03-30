      Subroutine NBIASP
c
c     fix wide-lane biases by pseudo-range data
c
      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'
      include 'parameters.h' 

      integer k1,k2,issn,ifix,ipseu,l1b,ia1,ia2,lim1,i,k

      real*8 awl,awlv,adn,xint,dev,prob,deci
c
c*      dimension junk(100),junk1(100),junk2(100)
      integer*4 junk(maxbis),junk1(maxbis),junk2(maxbis)

c      print*,'NBIASP is working ...', limitb(1),l1bias(1)
c      print *,'prdev prsig prcut ',prdev,prsig,prcut
C
      k1 = msig - 1
      issn = 0
      k2 = 0
      ifix = 0
      ipseu  =  0
         l1b = l1bias
         k1 = k1 + k2
         ia1 = issn + 1
         ia2 = issn + l1b*iband
         issn = issn + l1b
         k2 = 0
         do 2 i = ia1,ia2
            if (i.gt.issn) ipseu = ipseu + 1
            k = idxb(i)
            if(k.le.0) go to 2
            k2 = k2 + 1
            if (i.le.issn) goto 2
c---- force all continental scale biases free
            lim1 = ia1 + limitb - 1
            if(i.ge.lim1+l1b) goto 2
c---- get wl estimation
            awl = ddwl(ipseu)
            awlv = ddwv(ipseu)
            adn = adjust(k)
C---- old criteria:
c*     1. difference between two wl estimations < 0.8 cycle
c*     2. sigma of pseudo-range estimation < 0.15 cycle
c*     3. high confidence
c*
c           if (awlv.gt.0.15d0.and.noptin.eq.6) goto 2
c           if (awlv.gt.0.20d0.and.noptin.eq.7) goto 2
c*           dev = dabs(awl - adn)
c*           if (dev.gt.0.8d0.and.(noptin.eq.6.or.noptin.eq.7)) goto 2
c*           xint = dint(awl + 0.5d0*dsign(1.d0,awl))
c*           conf = dabs(xint - awl)
c*           if (noptin.eq.6.and.
c*    .          (conf.gt.3.0d0*awlv.or.conf.gt.0.3d0)) goto 2
c*           if (noptin.eq.7.and.
c*    .          (conf.gt.3.0d0*awlv.or.conf.gt.0.4d0)) goto 2
c --- new criteria:
c      1. Satisfies decision function with input values
c      2. Deviates by less than 0.5 cycles from phase estimate (if LC_HELP)
            xint = dint(awl + 0.5d0*dsign(1.d0,awl))
            dev = dabs(xint - awl)
            prob = 1.d0
            call bdeci(dev,awlv,1,prdev,prsig,prob,deci)
c            print *,'NBIASP k awl dev awlv prob deci '
c     .                     ,k,awl,dev,awlv,prob,deci
            if( deci.le.prcut) goto 2
            dev = dabs(awl - adn)
            if (dev.gt.0.5d0.and.(noptin.eq.6.or.noptin.eq.7)) goto 2
            ifix = ifix + 1
            nfix(ifix) = k1 + k2
            junk(ifix) = k
            junk1(ifix) = i
            junk2(ifix) = ipseu
            write( 6,110) k,awl,xint,adn
            write(10,110) k,awl,xint,adn
   2     continue
         if (iband.eq.2) issn = issn + l1bias

C---- Reorder array to fit BNEW
      call sort4i(ifix,junk1,junk2,junk,nfix,1)
c
         do 90 i = 1,ifix
            k = junk(i)
            k1 = junk1(i)
            k2 = junk2(i)
            awl = ddwl(k2)
c            this repeated to avoid storage:
            xint = dint(awl + .5d0*dsign(1.d0,awl))
            adn = adjust(k)
            bdev(i) = xint - adn
            adjust(k) = xint
            free(k) = 0
            idxb(k1) = -k
            nlive = nlive - 1
c*          write (6,110) k,adn,xint
c*          write (10,110) k,adn,xint
  90     continue
 110  format (2x,'(PR) Fix No.',I5,'  bias  from',f9.2,3x,'to'
     .       ,f8.1,'  (Phase est ',f9.2,' )')
c
       return
       end
