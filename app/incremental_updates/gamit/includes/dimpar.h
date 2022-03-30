C     DIMPAR.H: MAXIMUM DIMENSIONS FOR THE GAMIT PACKAGE

C     This file should be the same for all modules.
c     Last modified (except for values) by R. King 151125 

c     maximum number of parameter labels
c     3 crds +  2 clks + 6 ics + 13 rad + 3 svant + 6 eop 
      integer maxlab
      PARAMETER (MAXLAB=33)

c     maximum number of data types on x-file and used by model and solve
c     (two phase and two pseudorange): see makex.h for RINEX data types maxobt
      integer maxdat
      PARAMETER (MAXDAT=4)

* RWK 200501: add maxob11 and maxgnss for compatibility with GG11 rinex.h
*     maximum number of phase and pseurdo-range data types 
*     RINEX 3.03 defines 40 for GPS, 14 for Glonass, 18 for Beidou, 38 for Galileo, 16 for IRNSS, 28 for QZSS,
*     (two phase and two pseudorange)
      integer maxob11 
      parameter (maxob11=140)
*     maximum number of GNSS to be processed
      integer maxgnss
      parameter (maxgnss=6)
*     maximum number of channels: Used for GG 10.71 ; in GG 11 this is maxsat. (TAH 200509).
      integer maxsch   
      parameter (maxsch=200)                

c     maximum number of stations in a session
      integer maxsit
      PARAMETER (maxsit=80)

c     maximum number of satellites for t-files and GAMIT processing
      integer maxsat
      PARAMETER (maxsat=45)
                                                   
c     maximum number of satellite on an SP3 file (multi-GNSS)
      integer maxsp3sv
      PARAMETER (maxsp3sv=200)

c     maximum number of slips in AUTCLN
      integer maxslp
      PARAMETER (MAXSLP=1500)                         

c     maximum number of cycle-slip bias parameters flagged
      integer maxcsb
      PARAMETER (MAXCSB=200)

c     maximum number of orbital arcs (T-files)
      integer maxtfl
      PARAMETER (MAXTFL=7)

c     maximum number of integrated orbital parameters (ICs+rad) + SV antenna offsets
c     increased from 18 to 22 rwk 190425 to accommodate ECOM2 
      integer maxorb
      PARAMETER (MAXORB=22)
              
c     maximum number of components on T-file and equations integrated 
      integer maxytp,maxyt2
      PARAMETER (MAXYTP=3*(MAXORB+1))
      PARAMETER (MAXYT2=2*MAXYTP)

c     maximum number of eclipses for yaw rate and velocity impulse parameters  
c     max number is 2 per day x 2 planes x 4 SVs/plane = 16/day  
c     **removed from dimpar.fti and kept only in /model until we start estimating them
c      integer maxecl
c      PARAMETER (MAXECL=32)

c     maximum number of extra variables stored in extra array of C-file
      integer maxext
      PARAMETER (MAXEXT=2)

c     maximum number of extra variables stored in save array of C-file
      integer maxsav
      PARAMETER (MAXSAV=MAXSAT)

c     maximum number of extra variables stored in spare array of C-file
      integer maxspr
      PARAMETER (MAXSPR=2)


c     maximum number of text comment lines in a C-file
      integer maxtxt
      PARAMETER (MAXTXT=1000)

c     maximum number of biases
      integer maxbis
      PARAMETER (MAXBIS=2*MAXSIT*MAXSAT)
      
c     maximum number of zenith delay parameters per session
      integer maxatm
      PARAMETER (maxatm=25)
                          
c     maximum number of atmospheric gradient parameters 
c     (N/S + E/W per tabular point per session)
      integer maxgrad
      PARAMETER (maxgrad=6)
                 
c     maximum number of met rinex observations
      integer maxmetrnx
      parameter (maxmetrnx=1000)     

c     maximum number of global parameters
      integer maxeop
      PARAMETER (maxeop=6)

c     maximum number of parameters   
c     stn coords[3*maxsit], atmos[maxatm*maxsit], gradients[maxgrad*maxsit]
c     stn clock[maxsit]
c     sat orbits and SV antenna offsets [(maxorb)*maxsat], eop[maxeop]
c     biases[maxbis]
      integer maxprm
c old     .(MAXPRM =
c old     .  3*MAXNET+MAXATM*MAXSIT+3*MAXSIT
c old     . +(MAXORB+3)*MAXSAT+MAXBIS+5*MAXDAT)
      PARAMETER (MAXPRM = 3*MAXSIT + MAXATM*MAXSIT + MAXGRAD*MAXSIT +
     .           MAXSIT + MAXORB*MAXSAT + MAXBIS  +  MAXEOP )

c     maximum dimenesion of normal matrix
      integer maxnrm
      PARAMETER (MAXNRM=(MAXPRM*(MAXPRM+1)/2))

c     maximum number of observations at an epoch
      integer maxobs
      PARAMETER (MAXOBS=MAXSIT*MAXSAT)

c     maximum number of single differences at an epoch
      integer maxsng
      PARAMETER (MAXSNG=MAXSAT*(MAXSIT-1))

c     maximum number of double differences at an epoch
      integer maxdbl
      PARAMETER (MAXDBL=(MAXSAT-1)*(MAXSIT-1))

c     What's this for?
      integer maxdd
      PARAMETER (MAXDD=4*MAXDBL)

c     what's this for?
      integer maxsd
      PARAMETER (MAXSD=MAXSNG+2)

c     maximum size of a printed covariance matrix?
      integer maxcov
      PARAMETER (MAXCOV=21*MAXSAT)

      INTEGER MAXND
      PARAMETER (MAXND=MAXOBS+2)

      INTEGER MAXDM
c     maximum size number of elements in c matrix.
c     First 3 is site coords, second 3 is recv clock parameters.
c old PARAMETER (MAXDM=18*MAXOBS)
c old PARAMETER (MAXDM=(3+3+MAXATM+2+MAXORB+6)*MAXSAT*MAXNET)
      PARAMETER (MAXDM=(3+3+MAXATM+MAXGRAD+MAXORB+MAXEOP)*MAXSAT*MAXSIT)

      INTEGER MAXP2
      PARAMETER (MAXP2=MAXPRM+2)

      INTEGER MAXCRD
      PARAMETER (MAXCRD=3*MAXSIT)
              
c     maximum number of double-difference combinations
      INTEGER MAXWM1
      PARAMETER (MAXWM1=(MAXOBS*(MAXOBS+1)/2))

      INTEGER MAXWM2
      PARAMETER (MAXWM2=(MAXSNG*(MAXSNG+1)/2))
                        
c     maximum number of azimuth and elevation terms for receiver or
c     satellite antenna phase center variations
c      maxel =  46, maxaz = 181 allows 2-degree increments
c      maxel = 181, maxax = 721 allows 1 degree increments for a site-dependent model
      INTEGER MAXEL,MAXAZ              
      PARAMETER    (MAXAZ=721, MAXEL=181)


