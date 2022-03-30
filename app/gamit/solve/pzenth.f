c
      Subroutine PZENTH( istat,jsat,jzen,tmpart )

c     Compute partials for multiple zenith-delay parameters
c     R.W. King and T.A. Herring, partially from Y. Bock code in FILPAR  Sept-Oct 93

      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'

      integer*4 istat,jsat,jzen
      real* 8 ff,tmpart

      character*3 upperc

      logical epcontribs


c In solve.h
c-------------
c common <parts>  tpart(maxlab,maxsit,maxsat) : array of 1-way partials
c                   zenith-delay partials are in tpart(4,i,j)
c common <cutoff> istart : begining epoch of observations used
c                 iend   : ending epoch of observations used
c common <atmprms>
c  idtzen :  array (dim nzen + 1)of break-epochs for zenith parameters
c            computed in READ_BFL or READ
c            idtzen(1)         = istart
c            idtzen(nzen)    = iend
c            idtzen(nzen+1)  = iend + 1 (dummy to simplify logic)
c  zenmod :  type of zenith-delay model used
c             CON = constant (step)
c             PWL = piecewise linear

c Input
c------
c  istat  :  station index in tpart array
c  jsat   :  satellite index in tpart array
c  jzen   :  index in zenith-delay partial vector for this station and satellite
c  iepoch :  current epoch
c  ff     :  fraction of zenith-delay window of this epoch

c Output
c- -----
c tmpart  : partial derivative of this observation wrt this parameter
c           ( ipcnt in FILPAR)


c  Constant (step) model

      if( upperc(zenmod).eq.'CON' ) then

        if( ( jzen.lt.nzen .and.
     .     iepoch.ge.idtzen(jzen) .and. iepoch.lt.idtzen(jzen+1) ) .or.
     .      ( jzen.eq.nzen .and. iepoch.ge.idtzen(jzen) ) ) then
             tmpart = tpart(4,istat,jsat)
        else
             tmpart = 0.d0
        endif

c  Piecewise-linear model

      else if( upperc(zenmod).eq.'PWL' ) then

c       First determine if this epoch contributes to the partial for
c       this zenith delay tabular point.  The first and last tabular
c       point need to be treated specially becuase they have only one-
c       sided contributions.
c  NOTE: This code does not support multi-
c       session mode in the sense that the data contribution from the
c       next session is not accounted for.)
c  NOTE: The code belows assumes that the first idtzen epoch is on or
c       before the first epoch processed, and that the last idtzen is on or
c       after the last epoch processed. (Could print warning message)
c  NOTE: code assumes that PWL is not selected if there is only
c       one zenith delay parameter.
        epcontribs = .false.
        if( jzen.eq.1 ) then
            if( iepoch.le.idtzen(2) ) epcontribs = .true.
        else if ( jzen.eq.nzen ) then
            if( iepoch.gt.idtzen(nzen-1) ) epcontribs = .true.
        else
            if( iepoch.ge. idtzen(jzen-1) .and.
     .          iepoch.le. idtzen(jzen+1)   ) epcontribs = .true.
        end if

c       If this epoch contributes, compute partial based on time
c       before or after tabular point.  If it does not contribute then
c       set partial to zero.
        if( epcontribs ) then
            if( iepoch.ge. idtzen(jzen) ) then
c               epoch after the tabular point, scale partial by time
c               since tabular point
                ff = 1.d0 - dfloat(iepoch-idtzen(jzen)) /
     .                      dfloat(idtzen(jzen+1)-idtzen(jzen))
            else
c               epoch before the tabular point.  Note the epoch
c               on the tabular point is computed in the first part
c               of this if statement.)
                ff = 1.0 - dfloat(idtzen(jzen)-iepoch) /
     .                     dfloat(idtzen(jzen)-idtzen(jzen-1))
            end if
            tmpart = tpart(4,istat,jsat) * ff
        else
          tmpart = 0.d0
        end if
      endif

c     if(iepoch.eq.30.and.istat.eq.1.and.jsat.eq.2) then
c       write(*,'(a,a,26i5)') 'PZENTH: zenmod,maxatm,idtzen '
c    .                    ,   zenmod,maxatm,idtzen
c     endif
c     if(iepoch.eq.30.or.iepoch.eq.56.or.iepoch.eq.224.or.iepoch.eq.225)
c    .  then
c       if( istat.eq.1.and.(jsat.eq.2.or.jsat.eq.7)) then
c     write(*,'(a,2i5,2d9.2)') 'PZENTH: iepoch,jzen,ff,tmpart '
c    .                         ,         iepoch,jzen,ff,tmpart
c     endif
c     endif


      return
      end

