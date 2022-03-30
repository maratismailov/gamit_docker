      Subroutine BCHECK( mbias0,last_nonbias )
c
c---- count bias parameters and fill their pointors
c
      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'        
      include 'parameters.h'

      integer mbias0,last_nonbias,ntbias,icnt,nbias1,isat,istat
     .      , iparm,i1,ki,kl,i

      logical debug/.false./

         mbias=0
         msig=0
c---- count total number of bias parameters (L1 biases)
         ntbias=0
         do icnt=1,ntpart
          if(islot1(icnt).gt.2900 .and.
     .       islot1(icnt).le.7400) ntbias=ntbias+1
         enddo
c        It should = nsat * nsite
c        The l2flag option (L1-receiver) is no longer supported since the logic
c        in SOLVE is unduly complicated and unnecessary; rather cfmrg/fills1.f
c        defines L2 biases on the m-file header in all cases, leaving for SOLVE
c        the task of removing the L2-L1 biases if a single-frequency observable
c        is specified in the SOLVE batch file.
c*         if(l2flag.ne.0) ntbias=2*ntbias   
         ntbias = 2*ntbias 

c---- total number of non-bias parameters
          lpart=ntpart-ntbias     
c          print *,'BCHECK ntpart ntbias lpart ',ntpart,ntbias,lpart
c          ntpart is set in /cfmrg/fills1 and read from the m-file in /solve/read_bfl.f
c          ntbias is set in /solve/bcheck.f
c          lpart should be (3 coords)/stn + (nzen+ngrad)/stn/sesson + [optionally] 15 parm per SV
c---- count the number of live non-bias parameters
         do 10 i=1,lpart
            if (free(i).eq.0) goto 10
            msig=msig+1
  10     continue
         msig=msig+1
c
c      print *,'ntpart,lpart,ntbias',ntpart,lpart,ntbias
c number of biases this session
       nbias1=nsat*nsite
       ibias=nbias1
       if(l2flag.ne.0) ibias=2*ibias
c       print *,'ibias',ibias
c
c separate counters for bias parameters
      ibcnt1=0
      ibcnt2=0
      call zero1i(1,nbias1,ipntb1)
      call zero1i(1,nbias1,ipntb2)
c determine bias slots
      do 131 isat=1,nsat
         do 135 istat=1,nsite
            do 132 iparm=1,npartm(istat)
               ki=islot2(iparm,istat)
               kl=(islot1(ki)/100)+1
               if(kl.le.29) go to 132
c  storebias pointers
               ibcnt=ibcnt1+ibcnt2
               if(ibcnt.eq.ibias) go to 133
               if(kl.le.74) then
                  ibcnt1=ibcnt1+1
                  ipntb1(ibcnt1)=ki
               endif
c  Modified for multiple zenith delay parameters past islot1 = 21500
               if(kl.gt.74.and.kl.le.215.and.l2flag.ne.0) then
                  ibcnt2=ibcnt2+1
                  ipntb2(ibcnt2)=ki
               endif
  132       continue
c
c-----    change number of biases to fit implicit biases mode. -dnd-
         if (l2flag.ne.0.and.ibcnt2.eq.0) ibias=nbias1
  135     continue
  131  continue
  133  continue
c        print *, 'ibcnt1,ibcnt2',ibcnt1,ibcnt2   
c        print *,'IPNTB1 '
c        write(*,'(5i8)') (ipntb1(i),i=1,ibcnt1)  
c        print *,'IPNTB2'
c        write(*,'(5i8)') (ipntb2(i),i=1,ibcnt2)
c
c---- calculate total number of bias parameters
      mbias0=mbias
      do 100 i=1,ibcnt1
         i1=ipntb1(i)
         mbias=mbias+1
         idxb(mbias)=i1
         if (free(i1).eq.0) idxb(mbias)=-i1
 100  continue
      if(l2flag.ne.0) then
      do 110 i=1,ibcnt2
         i1=ipntb2(i)
         mbias=mbias+1
         idxb(mbias)=i1
         if (free(i1).eq.0) idxb(mbias)=-i1
 110  continue
      endif   
c      print *,'IDXB '
c      write(*,'(5i8)') (idxb(i),i=1,ibcnt2)

         i1=idxb(1)
         if (i1.lt.0) i1=-i1
         last_nonbias=i1-1
      iband=1
      if (mbias.eq.0) iband=0
      if (mbias.gt.nbias1) iband=2
      if (l2flag.le.0) iband=1
      if(debug) print *,'BCHECK mbias nbias1 l2flag iband last_nonbias ''
     .                 ,mbias, nbias1,l2flag,iband,last_nonbias
        
      return
      end
