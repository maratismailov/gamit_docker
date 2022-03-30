C     DIMPAR.H: MAXIMUM DIMENSIONS FOR THE GAMIT PACKAGE

C     This file should be the same for all modules.
c     Last modified (except for values) by R. King 151125 

c     maximum number of parameter labels
      integer maxlab
      PARAMETER (MAXLAB=29)

c     maximum number of data types on x-file and used by model and solve
c     (two phase and two pseudorange): see makex.h for RINEX data types maxobt
      integer maxdat
      PARAMETER (MAXDAT=4)

c     maximum number of stations in a session
      integer maxsit
      PARAMETER (maxsit=15)

c     maximum number of stations in a multi-session solution
c     should be larger than MAXSIT if stations not repeated each day
      integer maxnet 
      PARAMETER (MAXNET=MAXSIT)
c  ** Warning: There may be a bug in SOLVE that prohibits MAXNET>MAXSIT
c  **      PARAMETER (MAXNET=35)

c     maximum number of C files in a multi-session solution
      integer maxcfl
c     allow 4 day solutions on Apollo with dynamic memory allocation
c     PARAMETER (MAXCFL=4*MAXSIT)
c     This variable controls the size of solve.
c     At the moment, we only have 32MB of swap space, so
c     allow only one day solution
      PARAMETER (MAXCFL=1*MAXSIT)

c     maximum number of satellites for t-files and GAMIT processing
      integer maxsat
      PARAMETER (maxsat=32)
                                                   
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
      integer maxorb
      PARAMETER (MAXORB=18)
              
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
      PARAMETER (MAXBIS=2*MAXCFL*MAXSAT)
      
c     maximum number of zenith delay parameters per session
      integer maxatm
      PARAMETER (maxatm=13)
                          
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
c     stn coords[3*maxnet], atmos[maxatm*maxcfl], gradients[maxgrad*maxcfl]
c     stn clock[maxcfl]
c     sat orbits and SV antenna offsets [(maxorb)*maxsat], eop[maxeop]
c     biases[maxbis]
      integer maxprm
c old     .(MAXPRM =
c old     .  3*MAXNET+MAXATM*MAXCFL+3*MAXCFL
c old     . +(MAXORB+3)*MAXSAT+MAXBIS+5*MAXDAT)
      PARAMETER (MAXPRM = 3*MAXNET + MAXATM*MAXCFL + MAXGRAD*MAXCFL +
     .           MAXCFL + MAXORB*MAXSAT + MAXBIS  +  MAXEOP )

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
      PARAMETER (MAXDM=(3+3+MAXATM+MAXGRAD+MAXORB+MAXEOP)*MAXSAT*MAXNET)

      INTEGER MAXP2
      PARAMETER (MAXP2=MAXPRM+2)

      INTEGER MAXCRD
      PARAMETER (MAXCRD=3*MAXCFL)
              
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
      PARAMETER    (MAXAZ=181, MAXEL=46)


