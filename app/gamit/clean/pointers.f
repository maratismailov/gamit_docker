Copyright (c) Massachusetts Institute of Technology and the University of
California at San Diego, 1995. All rights reserved.

      Subroutine pointers( imsite
     .                   , coord_index,atm_index1,atm_index2,clock_index
     .                   , orb_index,svant_index,eop_index,grad_index )

C     Rewritten by R. King August 1995 to determine pointers to the first index of
c     parameter groups in the array of adjustments for a given site and satellite.

c  Input:

c   imsite                   =  index of site in solution

c   In common /postft
c     mtpart                =  number of parameters in solution
c     islot(maxprm)         =  codes for parameters adjusted in solution, from M-file

c  Output:

c    Index in adjust array for the first parameter in each group for input site and sat
c       coord_index         - station coordinates (3)
c       atm_index1          - average zenith delay parameter (1)
c       atm_index2          - time-dependent zenith delay parameters (numzen)
c       grad_index(2)       - atmospheric gradient parameters (2)
c       clock_index         - station clock parameters (3)
c       orb_index(maxsat)   - orbital parameters (norb: 9 or 15)
c       svant_index(maxsat) - satellite antenna offsets (3)
c       eop_index           - Earth rotation parameters (6)

c  Parameter codes in islot

c   GEOC. LAT. (DMS)          1-100
c   LONG. (DMS)             101-200
c   RADIUS (KM)             201-300
c   SC ATMOS-n              301-400   (Single zenith delay per station)
c   SC CLOCK-n EP SECS      401-500
c   ORBIT ELEMENT 1 PNnn    501-600
c   ORBIT ELEMENT 2 PNnn    601-700
c   ORBIT ELEMENT 3 PNnn    701-800
c   ORBIT ELEMENT 4 PNnn    801-900
c   ORBIT ELEMENT 5 PNnn    901-1000
c   ORBIT ELEMENT 6 PNnn   1001-1100
c   RAD PRES DIRECT PNnn   1101-1200
c   Y AXIS BIAS     PNnn   1201-1300
c   B AXIS BIAS     PNnn   1301-1400    (or X AXIS BIAS or Z AXIS BIAS)
c   COS DIRECT      PNnn   1401-1500  |
c   SIN DIRECT      PNnn   1501-1600  |
c   COS Y BIAS      Pnnn   1601-1700  |  optional (Berne model only)
c   SIN Y BIAS      PNnn   1701-1800  |
c   COS B BIAS      PNnn   1801-1900  |
c   SIN B BIAS      PNnn   1901-2000  |
c   SC CLOCK-n RATE        2001-2100
c   SC CLOCK-n ACC 1/D     2101-2200
c   SATL CLK-n EP SECS     2201-2300
c   SATL CLK-n RATE        2301-2400
c   SATL CLK-n ACC 1/D     2401-2500
c   PNnnmmmmk BIAS L1      2501-5500  EXPLICIT
c   PNnnmmmmk BIAS L1      5501-8500  IMPLICIT
c   PNnnmmmmk BIAS L2-L1   8501-11500 EXPLICIT
c   SC ATMOS-nn s         11501-14000 Mutiple zenith delays
c   SC   N/S ATMOS GRAD   14001-16500
c   SC   E/W ATMOS GRAD   16501-19000
c   X POLE (ARCS)         80001-80001
c   X POLE RATE (ARCS/D)  80002-80002
c   Y POLE (ARCS)         80003-80003
c   Y POLE RATE (ARCS/D)  80004-80004
c   UT1-TAI (SEC)         80005-80005
c   UT1-TAI RATE (SEC/D)  80006-80006

      implicit none

      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
c
      integer*4 imsite,coord_index,atm_index1,atm_index2,clock_index
     .        , orb_index(maxsat),svant_index(maxsat),eop_index
     .        , grad_index(2),islot_coord,islot_atm1,islot_atm2
     .        , islot_clock,islot_eop,islot_grad(2),ksat,i,k

c-----Initialize pointers

      coord_index = 0
      atm_index1  = 0
      atm_index2  = 0 
      grad_index(1)  = 0
      grad_index(2)  = 0
      clock_index = 0
      eop_index   = 0
      islot_atm1 = 0
      islot_atm2 = 0
      do i = 1, maxsat
        orb_index(i) = 0
        svant_index(i) = 0
      enddo

c-----Compute islot numbers for the first parameter of each group for this site

      islot_coord = imsite 
      if( numzen.gt.0 ) then
        islot_atm1 = 300 + imsite
        if( numzen.gt.1 ) then
          islot_atm2 = 11501 + numzen*(imsite-1)
        endif
      endif
      islot_grad(1)  = 14000 + imsite
      islot_grad(2)  = 16500 + imsite
      islot_clock = 400 + imsite
      islot_eop   = 80001


c-----Loop through all adjusted parameters

      do k = 1, mtpart

c       site coordinates

        if( islot(k).eq.islot_coord ) then
           coord_index = k


c       zenith delays

        else if( islot(k).eq.islot_atm1 ) then
           atm_index1 = k
        else if( islot(k).eq.islot_atm2 ) then
           atm_index2 = k

c       atmosphere gradient parameters

        else if( islot(k).eq.islot_grad(1) ) then
           grad_index(1) = k
        else if( islot(k).eq.islot_grad(2) ) then
           grad_index(2) = k

c       station clock

        else if( islot(k).eq.islot_clock ) then
          clock_index = k
c         check to make sure station clock rate and acceleration are next-
c**       these parameters removed to fix in SV antenna offsets rwk 980911
c**          if( islot(k+1).ne.(2000+imsite) .and.
c**     .        islot(k+2).ne.(2100+imsite) ) then
c**               write(*,'(a)') 'POINTERS: Something wrong in islot'
c**                   stop
c**          endif

c       orbital parameters

        else if( islot(k).ge.501 .and. islot(k).le.600 ) then
           ksat = mod(islot(k),100)
           if( ksat.le.maxsat) then
              orb_index(ksat) = k
           else
             write(*,'(a,i4,a)') ' POINTERS: ksat = ',maxsat,' > maxsat'
             stop
           endif

c       satellite antenna offsets 

        else if( islot(k).ge.2001 .and. islot(k).le.2100 ) then
           ksat = mod(islot(k),100)
           if( ksat.le.maxsat) then
              svant_index(ksat) = k
           else
             write(*,'(a,i4,a)') ' POINTERS: ksat = ',maxsat,' > maxsat'
             stop
           endif

c     Earth orientation parameters

        else if( islot(k).eq.islot_eop ) then
           eop_index = k

        endif

      enddo

      return
      end
