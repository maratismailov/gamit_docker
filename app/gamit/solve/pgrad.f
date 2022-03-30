c
      Subroutine PGRAD( istat,jsat,jgrad,itype,tmpart )

c     Compute partials for multiple gradient parameters
c     R.W. King from subroutine pzenth.   17 Sept 98

      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'

      integer*4 istat,jsat,jgrad,itype
      real* 8 ff,tmpart

      character*3 upperc

      logical epcontribs


c In solve.h
c-------------
c common <gparts>  gpart(maxsit,maxsat,2) : array of 1-way gradient partials
c                  at this epoch; N/S are gpart(i,j,1), E/W are gpart(i,j,2)
c common <cutoff> istart : begining epoch of observations used
c                 iend   : ending epoch of observations used
c common <atmprms>
c  idtgrad:  array (dim ngrad + 1)of break-epochs for gradient parameters
c            computed in READ_BFL or READ
c            idtgrad(1)          = istart
c            idtgrad(ngrad)    = iend
c            idtgrad(ngrad+1)  = iend + 1 (dummy to simplify logic)
c  gradmod :  type of gradient-delay model used
c             CON = constant (step)
c             PWL = piecewise linear

c Input
c------
c  istat  :  station index in tpart array
c  jsat   :  satellite index in tpart array
c  jgrad  :  index in gradient-delay partial vector for this station and satellite
c  iepoch :  current epoch
c  itype  :  N/S = 1,   E/W = 2
c  ff     :  fraction of gradient-delay window of this epoch

c Output
c- -----
c tmpart  : partial derivative of this observation wrt this parameter
c           ( ipcnt in FILPAR)


c  Constant (step) model

      if( upperc(gradmod).eq.'CON' ) then

        if( ( jgrad.lt.ngrad .and.
     .     iepoch.ge.idtgrad(jgrad) .and. iepoch.lt.idtgrad(jgrad+1) ) 
     .   .or. ( jgrad.eq.ngrad .and. iepoch.ge.idtgrad(jgrad) ) ) then
             tmpart = gpart(istat,jsat,itype)
        else
             tmpart = 0.d0
        endif

c  Piecewise-linear model

      else if( upperc(gradmod).eq.'PWL' ) then

c       First determine if this epoch contributes to the partial for
c       this gradient delay tabular point.  The first and last tabular
c       point need to be treated specially becuase they have only one-
c       sided contributions.
c  NOTE: This code does not support multi-
c       session mode in the sense that the data contribution from the
c       next session is not accounted for.)
c  NOTE: The code belows assumes that the first idtgrad epoch is on or
c       before the first epoch processed, and that the last idtgrad is on or
c       after the last epoch processed. (Could print warning message)
c  NOTE: code assumes that PWL is not selected if there is only
c       one gradient delay parameter.
        epcontribs = .false.
        if( jgrad.eq.1 ) then
            if( iepoch.le.idtgrad(2) ) epcontribs = .true.
        else if ( jgrad.eq.ngrad ) then
            if( iepoch.gt.idtgrad(ngrad-1) ) epcontribs = .true.
        else
            if( iepoch.ge. idtgrad(jgrad-1) .and.
     .          iepoch.le. idtgrad(jgrad+1)   ) epcontribs = .true.
        end if

c       If this epoch contributes, compute partial based on time
c       before or after tabular point.  If it does not contribute then
c       set partial to zero.
        if( epcontribs ) then
            if( iepoch.ge. idtgrad(jgrad) ) then
c               epoch after the tabular point, scale partial by time
c               (fraction of window) since last tabular point
                ff = 1.d0 - dfloat(iepoch-idtgrad(jgrad)) /
     .                      dfloat(idtgrad(jgrad+1)-idtgrad(jgrad))
            else
c               epoch before the tabular point.  Note the epoch
c               on the tabular point is computed in the first part
c               of this if statement.)
                ff = 1.0 - dfloat(idtgrad(jgrad)-iepoch) /
     .                     dfloat(idtgrad(jgrad)-idtgrad(jgrad-1))
            end if
            tmpart = gpart(istat,jsat,itype) * ff
        else
          tmpart = 0.d0
        end if
      endif

c     if(iepoch.eq.30.and.istat.eq.1.and.jsat.eq.2) then
c       write(*,'(a,a,26i5)') 'PGRAD: gradmod,maxatm,idtgrad '
c    .                    ,   gradmod,maxatm,idtgrad
c     endif
c     if(iepoch.eq.30.or.iepoch.eq.56.or.iepoch.eq.224.or.iepoch.eq.225)
c    .  then
c       if( istat.eq.1.and.(jsat.eq.2.or.jsat.eq.7)) then
c     write(*,'(a,2i5,2d9.2)') 'PGRAD: iepoch,jgrad,ff,tmpart '
c    .                         ,         iepoch,jgrad,ff,tmpart
c     endif
c     endif


      return
      end

