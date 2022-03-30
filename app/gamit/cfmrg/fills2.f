Copyright (c) Massachusetts Institute of Technology and the University of
California at San Diego, 1995.  All rights reserved.

      Subroutine fills2(icrchk,iatchk,igcchk,iorchk,iantchk,ieochk
     $                 ,isite,itsat,nsite,norb
     $                 ,ldual,stbias,kpart,islot2,nsat
     $                 ,numzen,numgrad,is2names,is2nums)
c
c Purpose:
c     Fill the C-file-specific parameter slots in islot2
c
      implicit none
c
      include '../includes/dimpar.h'
c
      character*1  explct
      character*1  stbias(maxsit)
      character*1  upperc
      character*5  is2names(maxprm)

      integer*4 is2nums(maxprm),isite,jsite,nsite
     .      , iatchk,icrchk,igcchk,iorchk,iantchk,ieochk
     .      , kpart,lpart,lsat,norb
     .      , itsat,nsat
     .      , izen,numzen,igrad,numgrad
     .      , islot2(maxprm),k

      logical ldual(maxsit)
c
      data explct/'e'/
c
      explct = upperc(explct)

      kpart=0
      lpart=0

      if (icrchk.eq.1) then
c.....   station partials
         do  jsite=1,nsite
            if (jsite.eq.isite) then
               do  k=1,3
                  kpart=kpart+1
                  islot2(kpart)=lpart+k
                  is2names(kpart)="site"
                  is2nums(kpart) = isite
               enddo
            endif
            lpart=lpart+3
         enddo
      endif

      if (iatchk.eq.1) then
 
c.....   atmosphere  --first the average value, always used, then multiple delays
            do  jsite=1,nsite
                  if (jsite.eq.isite) then
                     kpart=kpart+1
                     islot2(kpart)=lpart+1
                     is2names(kpart)="atm"
                     is2nums(kpart) = isite
                  endif
                  lpart=lpart+1
            enddo
         if( numzen.gt.1 ) then
              do  jsite=1,nsite
                do  izen=1,numzen
                    if (jsite.eq.isite) then
                       kpart=kpart+1
                       islot2(kpart)=lpart+1
                       is2names(kpart)="atm"
                       is2nums(kpart) = isite
                    endif
                    lpart=lpart+1
                enddo
              enddo
         endif

c N/S atmospheric gradient parameter   
cd         print *,'FILLS1 numgrad ',numgrad
           do jsite=1,nsite 
             do  igrad=1,numgrad
                 if (jsite.eq.isite) then
                   kpart=kpart+1
                   islot2(kpart)=lpart+1
                   is2names(kpart)="nsatm"
                   is2nums(kpart) = isite
cd                   print *,'site igrad kpart islot2 '
cd     .                ,jsite,igrad,kpart,islot2(kpart)
               endif
               lpart=lpart+1
            enddo
           enddo
c E/W atmospheric gradient parameter
           do jsite=1,nsite
             do  igrad=1,numgrad
               if (jsite.eq.isite) then
                 kpart=kpart+1
                 islot2(kpart)=lpart+1
                 is2names(kpart)="ewatm"
                 is2nums(kpart) = isite
               endif
               lpart=lpart+1
            enddo
           enddo
      
      endif

      if (igcchk.eq.1) then
c.....   clock epoch partial (rate and acceleration no longer estimatable)
            do jsite=1,nsite
                  if (jsite.eq.isite) then
                        kpart=kpart+1
                        islot2(kpart)=lpart+1
                        is2names(kpart)="rvclk"
                        is2nums(kpart) = isite
                  endif
                  lpart=lpart+1
            enddo
      endif

      if (iorchk.eq.1) then

c.....   satellite orbital partials
         do  lsat=1,itsat
             do  k=1,norb
                islot2(kpart+k)=lpart+k
                is2names(kpart+k)="orbit"
                is2nums(kpart+k) = lsat
             enddo
             lpart=lpart+norb
             kpart=kpart+norb
         enddo
      endif      
 
     
      if (iantchk.eq.1) then
c.....   satellite antenna offsets
         do lsat=1,itsat
           do  k=1,3
              islot2(kpart+k)=lpart+k
              is2names(kpart+k)="svant"
              is2nums(kpart+k) = lsat
           enddo
           lpart=lpart+ 3
           kpart=kpart+ 3
         enddo
      endif

      if (ieochk.eq.1) then
c.....   earth orientation partials
         do  k=1,6
              islot2(kpart+k)=lpart+k
              is2names(kpart+k)="eop"
              is2nums(kpart+k) = k
         enddo
         lpart=lpart+ 6
         kpart=kpart+ 6
      endif


c.....biases
c.....   implicit or explicit l1
         do  jsite = 1, nsite
               do  lsat = 1, nsat
                  if (jsite.eq.isite) then
                     kpart = kpart+1
                     islot2(kpart) = lpart+1
                     is2names(kpart)="l1bis"
                     is2nums(kpart) = isite
                  endif
                  lpart = lpart+1
               enddo
         enddo

c.....   explicit l2
         do jsite = 1, nsite
            if ( ldual(jsite)) then
               do  lsat = 1, nsat
                  if (jsite.eq.isite) then
                     kpart = kpart+1
                     islot2(kpart) = lpart+1
                     is2names(kpart)="l2bis"
                     is2nums(kpart) = isite
                  endif
                  lpart = lpart+1
               enddo
            endif
         enddo

      return
      end
