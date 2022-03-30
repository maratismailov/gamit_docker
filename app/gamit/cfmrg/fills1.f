Copyright (c) Massachusetts Institute of Technology and the University of
California at San Diego, 1995.  All rights reserved.

      Subroutine fills1(igcchk,iatchk
     $                ,icrchk,iorchk,iantchk,ieochk
     $                ,nsite,itsat,norb,isat
     $                ,ldual,usite,fname,stbias,clabel
     $                ,mtpart,idms,islot1,alabel,totsat,nsat
     $                ,numzen,numgrad,l1_only )

c
c Purpose:  Fill islot1 (partial slots for all sessions) and label and preval values

c  Input
c    clabel (maxlab)  c*20  parameters names from c-file

c  Output
c    mtpart           i*4   number of parameters for m-file
c    islot1 (maxprm)  i*4   pointers for m-file
c    alabel (maxprm)  c*20  parameter names for m-file


C-File order   M-file order             Slot value
c----------    ----------------------   ----------
c  1            'SC GEOC. LAT. (DMS) ',     1-100
c  2            'SC LONG. (DMS)      ',  101-200
c  3            'SC RADIUS (KM)      ',  201-300
c  4            'SC ATMOS-n          ',  301-400  (average zenith delay)
c  5            'SC CLOCK-n EP SECS  ',  401-500  
c------------------------------------------------------------------------------------
c  6            'ORBIT ELEMENT 1 PNnn',  501-600
c  7            'ORBIT ELEMENT 2 PNnn',  601-700
c  8            'ORBIT ELEMENT 3 PNnn',  701-800
c  9            'ORBIT ELEMENT 4 PNnn',  801-900
c 10            'ORBIT ELEMENT 5 PNnn',  901-1000
c 11            'ORBIT ELEMENT 6 PNnn',  1001-1100
c 12            'RAD PRES DIRECT PNnn',  1101-1200
c 13            'Y AXIS BIAS     PNnn',  1201-1300
c 14            'B AXIS BIAS     PNnn',  1301-1400    (or X AXIS BIAS or Z AXIS BIAS)  
c-------------------------------------------------------------------------------------
c 15            'COS U DIRECT    PNnn',  1401-1500  |
c 16            'SIN U DIRECT    PNnn',  1501-1600  |
c 17            'COS U Y         PNnn',  1601-1700  |  optional (ECOM1/BERNE or ECOM2-extended models)
c 18            'SIN U Y         PNnn',  1701-1800  |
c 19            'COS U B         PNnn',  1801-1900  |
c 20            'SIN U B         PNnn',  1901-2000  |
c-------------------------------------------------------------------------------------   
c 21            'COS 2U DIRECT   PNnn',  2001-2100  |
c 22            'SIN 2U DIRECT   PNnn',  2101-2200  |  optional (ECOM2 model)
c 23            'COS 4U DIRECT   PNnn',  2201-2300  |                                       
c 24            'SIN 4U DIRECT   PNnn',  2301-2400  |

c-------------------------------------------------------------------------------------
c The following C-file order assumes ECOM1/BERNE or ECOM2 
c Subtract 4 for ECOM1; subtract 10 3-parameter models; subtract 23 if no orbits
c-------------------------------------------------------------------------------------
c 25            'SVANT X AXIS    PNnn',  2401-2600
c 26            'SVANT Y AXIS    PNnn',  2501-2700
c 27            'SVANT Z AXIS    PNnn',  2601-2800  
c------------------------------------------------------------------------------------
c 28            'PNnnmmmmk BIAS L1   ',  2901-7400 
c 29            'PNnnmmmmk BIAS L2-L1',  7401-11900
c------------------------------------------------------------------------------------
c                nn is the nn'th zenith delay or gradient for that site
c                s  is the session number               
c 30            'SC ATM ZEN  nn s    ', 21501-24000 Multiple zenith delays 
c 31            'SC N/S GRAD nn s    ', 24001-26500 Multiple gradient parameters
c 32            'SC E/W GRAD nn s    ', 26501-29000 Multiple gradient parameters
c-------------------------------------------------------------------------------------
c 33            'X POLE (ARCS)       ', 80001-80001
c 34            'X POLE RATE (ARCS/D)', 80002-80002
c 35            'Y POLE (ARCS)       ', 80003-80003
c 36            'Y POLE RATE (ARCS/D)', 80004-80004
c 37            'UT1-TAI (SEC)       ', 80005-80005
c 38            'UT1-TAI RATE (SEC/D)', 80006-80006
c--------------------------------------------------------------------------------------
c Note:  The C-file pointer value for a parameter, not currently used by CFMRG
c        or SOLVE, is now the first value of the range of values (islot1) 
c        written on the M-file and used by SOLVE and CVIEW/SCANDD.  rwk 980905        
c-------------------------------------------------------------------------------------

      implicit none

      include '../includes/dimpar.h'

      integer maxbin
      parameter (maxbin=38)

      logical ldual(maxsit),l1_only

      integer*4 mpos,isite,nsite
     .        , igcchk,iatchk,icrchk,iorchk,iantchk
     .        , ieochk,itsat,nsat
     .        , mtpart,norb,msat,jzen,jgrad,nn,nn1
     .        , n,i

      real*8  aprclk,aprorb,aprant,adjust,aprval,aprcrd,apreop,aprzen
     .        ,aprgrd

      integer ibin(maxbin)
      integer isat(maxsat),totsat(maxsat)
      integer islot1(maxprm)
      integer idms(maxprm)
      integer*4 numzen,izen,numgrad,igrad

      character*1   explct
      character*1   stbias(maxsit)
      character*1   upperc
      character*2   buf2
      character*4   usite(maxsit),sit4
c      character*16  blank
      character*16  fname(maxsit)
      character*20  alabel(maxprm)
      character*20  clabel(maxlab)
      character*20  rlabel(5)
      character*256 message
c
      common/mercom/aprval(maxprm),adjust(maxprm),aprcrd(maxsit,3),
     $              aprclk(maxsit,3),aprorb(maxorb,maxsat),
     $              aprant(3,maxsat),apreop(6),aprzen(maxsit),
     $              aprgrd(maxsit,2)

c.....labels for parameters not explicitly on the C-file

      data rlabel/
     $ 'PNNNMMMMK BIAS L1   ','PNNNMMMMK BIAS L1   '
     $,'PNNNMMMMK BIAS L2-L1'
     $,'SITE N/S GRAD m     ','SITE E/W GRAD m     '/

c.....bins
c              |  coords |atm|station clock | 
c               1  2   3   4   5  
      data ibin/0,100,200,300,400

c              |     satellite state  |
c                6   7   8   9   10  11
     $      ,   500,600,700,800,900,1000

c           |  satellite non-gravitational parameters    |
c             12   13   14   15   16   17   18   19   20   
     $    ,  1100,1200,1300,1400,1500,1600,1700,1800,1900
c             21   22   23   24 |
     $    ,  2000,2100,2200,2300
c
c               | SV antenna offsets |
c             25   26   27
     $    ,  2400,2500,2600
c
c    biases :   | l1 explicit | l2 explicit 
c                  28             29       
     $    ,       2900,          7400     
c
c           | multi-zenith | N/S atmos grad | E/W atmos grad |
c                  30              31               32 
     $    ,      21500,          24000,           26500
c             
c    
c    eop :   |  xp  | xp dot |  yp  | yp dot |  ut1  | ut1 dot |
c               33       34      35      36      37       38
     $    ,    80001,  80002,  80003,  80004,  80005,   80006/
c
      data explct/'e'/
c      data blank/'                '/
c
c case insensitivity
c
      explct = upperc(explct)
c
      n=0
c
c.....coordinates (1st, 2nd, and 3rd types in 'ibin')
      if (icrchk.eq.1) then
         nn=0
         do isite=1,nsite
            do i = 1, 3
               n=n+1 
               call check_limits( n,nn+i,nn+i,maxbin )
               if(i.ne.3) idms(n) = 1
               islot1(n) = ibin(nn+i) + isite
               aprval(n) = aprcrd(isite,nn+i)
               alabel(n) = clabel(nn+i)
               alabel(n) (1:4) = usite(isite)
            enddo
         enddo
      endif
c
c.....atmosphere
      if (iatchk.eq.1) then      

c average zenith delay (4th type in 'ibin')
         nn=4    
           do isite=1,nsite
              n=n+1 
              call check_limits( n,4,nn,maxbin )
              islot1(n) = ibin(nn) + isite
              aprval(n) = aprzen(isite)
              alabel(n) = clabel(4)
              alabel(n) (1:4) = usite(isite)
           enddo  

c multiple zenith delays (30th type in 'ibin') 
         if( numzen.gt.1 ) then
         nn=30
         jzen=0
            do  isite=1,nsite
              do izen=1,numzen
                jzen=jzen+1
                n=n+1  
                call check_limits( n,4,nn,maxbin )
                islot1(n) = ibin(nn) + jzen
                aprval(n) = aprzen(isite)
                alabel(n) = clabel(4)
                alabel(n) (1:4) = usite(isite)
                write (alabel(n)(16:18),'(i3)') izen
              enddo
            enddo
         endif

c N/S atmospheric gradient parameter (31st type in 'ibin')
         nn=31
         jgrad=0   
cd         print *,'fills1 nn numgrad jgrad ',nn,numgrad,jgrad
           do  isite=1,nsite      
cd             print *,'isite  ',isite
             do igrad=1,numgrad
               jgrad=jgrad+1
c.....         if no site on that session, no atmosphere
                 n=n+1 
                 call check_limits( n,4,nn,maxbin )
                 islot1(n) = ibin(nn) + jgrad
cd                 print *,'n nn jgrad islot1 ',n,nn,jgrad,islot1(n)
                 aprval(n) = aprgrd(isite,1)
                 alabel(n) = rlabel(4)
                 alabel(n) (1:4) = usite(isite)  
                 write (alabel(n)(17:18),'(i2)') igrad
c                 write (alabel(n)(20:20),'(i1)') mdy
             enddo
           enddo
c E/W atmospheric gradient parameter  (32nd type in 'ibin')
         nn=32
         jgrad=0
           do  isite=1,nsite 
             do igrad=1,numgrad   
               jgrad=jgrad+1
c.....         if no site on that session, no atmosphere
               n=n+1     
               call check_limits( n,5,nn,maxbin )
               islot1(n) = ibin(nn) + jgrad
               aprval(n) = aprgrd(isite,2)
               alabel(n) = rlabel(5)
               alabel(n) (1:4) = usite(isite)
               write (alabel(n)(17:18),'(i2)') igrad
c               write (alabel(n)(20:20),'(i1)') mdy
            enddo
           enddo
      endif


c.....station clock epoch (5th type in 'ibin')
      if (igcchk.eq.1) then
         nn=5
            do isite=1,nsite
c.....         if site on that session, no clock
                   n=n+1  
                   call check_limits( n,nn,nn,maxbin )
                   islot1(n) = ibin(nn) + isite
                   alabel(n) = clabel(5)
                   aprval(n) = aprclk(isite,1)
                   alabel(n) (1:4) = usite(isite)
c                   write (alabel(n)(12:12),'(i1)') mdy
            enddo
      endif


      if (iorchk.eq.1) then

c.....   integrated orbital parameters  (6th - 24th types in 'ibin')
         nn=5
         do msat=1,itsat
            do i=1,norb
               n=n+1    
               call check_limits( n,5+i,nn+i,maxbin )
               islot1(n) = ibin(nn+i) + msat
c              The following is to trap SVs in the input list but missing from 
c              the C-file header; this shouldn't happen unless FIXDRV is not 
c              rerun after SV changes.  A better fix would be use of 'nsat' 
c              rather than 'itsat here, or different logic in cfmrg.f, or reading 
c              of the t-file directly.
               if( i.le.6 .and. aprorb(i,msat).eq.0.d0 ) then
                 write(message,'(a,i3,a)')
     .           'Error, orbital element for satellite: ',totsat(msat),
     .           ' = 0. Values missing from C-file? '
                 call report_stat('FATAL','CFMRG','fills1',' ',
     .           message,0)
               endif         
               aprval(n) = aprorb(i,msat)
               alabel(n) = clabel(5+i)
               write(buf2,'(i2)') totsat(msat)
               if(buf2(1:1).eq.' ') buf2(1:1)= '0'
               write(alabel(n)(17:20),'(2a2)') 'PN',buf2
            enddo
         enddo
      endif   

c.....   satellite antenna offsets (25th - 27th types in 'ibin')
       if (iantchk.eq.1) then
         nn=24            
         do msat = 1,itsat
            do  i=1,3
               n=n+1  
               call check_limits( n,5+norb+i,nn+i,maxbin )
               islot1(n) = ibin(nn+i) + msat
               aprval(n) = aprant(i,msat)
               alabel(n) = clabel(5+norb+i)
               write(buf2,'(i2)') totsat(msat)
               if(buf2(1:1).eq.' ') buf2(1:1)= '0'
               write(alabel(n)(17:20),'(2a2)') 'PN',buf2
            enddo
         enddo
       endif

c
c.....   earth orientation parameters (33rd - 38th types in 'ibin')
       if (ieochk.eq.1) then
         nn=32
            do i=1,6
               n=n+1   
               call check_limits( n,5+norb+3+i,nn+i,maxbin )
               islot1(n) = ibin(nn+i)
               aprval(n) = apreop(i)
               alabel(n) = clabel(5+norb+3+i)
            enddo
       endif
c
c.....    biases (28th - 29th types in 'ibin')
c
       do isite=1,nsite
           mpos=index(fname(isite),'.')
           nn=28 
           nn1 = 1
           do msat=1,nsat
              n=n+1
* MOD TAH 190615: Added the site number (isite) back into the slot
*     number.  In 10.70, this was nfile(isite,mdy) which with single
*     session had the site number in it.  (The cfmrg.out file had the 
*     same slot number of all sites of the same bias -- strangely did
*     not seem to crash solve.  Same fix below for L2 bias.
C             islot1(n) = ibin(nn) + 100*(msat-1) + 1
              islot1(n) = ibin(nn) + 100*(msat-1) + isite
              alabel(n) = rlabel(nn1)
              write(alabel(n)(3:4),'(i2)') isat(msat)
              sit4 = fname(isite)(mpos-5:mpos-2) 
              call uppers(sit4)
              alabel(n)(5:8) = sit4
              alabel(n)(9:9) = ' '  
           enddo
c      end loop on sites for L1 biases
       enddo

       nn=29                          
       nn1 =3 
c      If single and dual frequency within the same session. assign dual-frequency
c      biases to all to ease the logic in SOLVE (at least for now)  rwk 990705
c      The same holds true for unmixed, actually, so always define L2 biases.  rwk 990906
c      See also /solve  read_bfl.f and bcheck.f
       l1_only = .false.
       do isite = 1,nsite
         if(.not.ldual(isite)) l1_only =.true.
       enddo              
c      now assign the L2 bias slots
       do isite=1,nsite
            mpos=index(fname(isite),'.')
            do msat=1,nsat
             n=n+1
* MOD TAH 190615: Added site number back to islot (see comment above associated
*     with L1 biases.
C            islot1(n) = ibin(nn) + 100*(msat-1) + 1
             islot1(n) = ibin(nn) + 100*(msat-1) + isite
             alabel(n) = rlabel(nn1)
             write(alabel(n)(3:4),'(i2)') isat(msat)
             sit4 = fname(isite)(mpos-5:mpos-2)
             call uppers(sit4)
             alabel(n)(5:8) = sit4
             alabel(n)(9:9) = ' '  
           enddo
c      end loop on sites for L2 biases
       enddo

      mtpart=n
      if (mtpart.gt.maxprm) then                  
        write (message,'(a,i5,a,a5)') 
     .   'Number of partial slots (',mtpart,') > maxprm (',maxprm,')'
        call report_stat('FATAL','CFMRG','fills1',' ',message,0)
      endif
cd     write (*,'(1x,14i5)') nsite,n

      return
      end
