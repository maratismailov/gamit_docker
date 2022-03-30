      Subroutine rate_est ( iepoch, iprn, isite, pl1, ierrf, pp
     .                    , c , npnts, debug )

c     Estimate the SV clock polynomial coefficients using Eqn. 6 in
c     Feigl et al 1991.  Hardwired to fit 3rd-order polynomial from 3 pts
c     Feigl/King/Dong   1991-June 1992

c       Input:  iepoch = epoch number
c               iprn   = satellite prn index
c               isite  = input C-file index
c               pl1    = L1 phase residuals
c               ierrf  = error flags for pl1
c               pp     = partial derivatives (3 epochs, 3 coefficients)
c               nn     = number of data points to use
c               debug  = logical for debug (T/F)
c       Output: c      = polynomial coefficients (Feigl et al., eqn.  )
c               npnts  = number of good points found (-1 if bad inversion)

      implicit none

      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/errflg.h'

      integer*4 ierrf(maxepc,maxsat,maxsit),info,iprn,iepoch,isite
     .        , mm,nn,npnts,i,j,k,l,m

c     must be changed in all routines:  MAKEX, RATE_EST, ALLANV, AVGCLK
      parameter(mm=3,nn=3)

      real*8 pl1(maxepc,maxsat,maxsit),x1(mm),b1(mm),a1(mm,mm)
     .     , c(mm,maxsit),pp(mm,nn),scale(mm),determ(2),rcond

      logical lgood,lbias,debug

         if( debug ) print *,' RATE_EST:'

c        Initialize the output estimates

        do k=1,mm
           c(k,isite) = 0.d0
        enddo


c        Zero the normal equations

         do k=1,mm
             x1(k) = 0.0d0
             b1(k) = 0.0d0
             do l = 1,mm
                a1 (l,k) = 0.0d0
             enddo
          enddo


c        Fill the normal equations

          npnts = 0

          j = 0
          do 10 k = iepoch-nn/2, iepoch+nn/2
            j = j + 1
            if( debug ) print *,'  k, iprn, isite, ierrf: '
     .                ,k,iprn,isite,ierrf(k,iprn,isite)
c           use the point if its not unweighted and not flagged for a bias
            if ( lgood(ierrf(k,iprn,isite)) .and.
     .           .not. lbias(ierrf(k,iprn,isite)) ) then
               npnts = npnts + 1
c              fill normal equations
c              upper triangle of left hand side
               do m = 1,mm
                  do l = 1,m
                     a1(l,m) = a1(l,m) + pp(l,j)*pp(m,j)
                  enddo
               enddo

c             right hand side
              do m = 1,mm
                  b1(m) = b1(m) + pp(m,j)*pl1(k,iprn,isite)
              enddo

            endif

  10     continue


c       Solve the normal equations if there are three good points

      if( npnts.eq.nn ) then

            info = 0
c           scale the normal matrix
            do k = 1,mm
               scale(k) = 1.0d0/dsqrt(a1(k,k))
            enddo

            do l = 1,mm
               do k = 1,mm
                  a1(k,l) = scale(k)*scale(l)*a1(k,l)
               enddo
            enddo


c            invert normal matrix A1 using linpack (blas 1) routines
c            assuming that the matrix is positive definite (and thus
c            symmetric). DPODI returns inverse in upper triangle, so
c            the full matrix must be filled in

            call dpoco(a1,mm,mm,rcond,x1,info)

c           We need a good inversion and a delay rate to proceed:
            if (info                 .eq. 0
     .         .and. rcond+1.0d0     .ne. 1.0d0  ) then

               call dpodi(a1,mm,mm,determ,11)

c              fill in lower half and rescale
               do k = 1,mm
                  do l = 1, k
                     a1(l,k) = scale(k)*scale(l)*a1(l,k)
                     a1(k,l) = a1(l,k)
                  enddo
               enddo

c              multiply inv(A1) * B1 = X1
               call dgemv ('N', mm, mm
     .                    , 1.d0, a1, mm
     .                    , b1, 1
     .                    , 0.d0, x1, 1)
               do i=1,3
                 c(i,isite)=x1(i)
               enddo

            else
c              bad inversion
               if( debug ) then
                  print *,' Bad inversion of polynomial fit in RATE_EST'
                  print *,'  iepoch, iprn, isite = ',iepoch, iprn, isite
                  print *,'  a1 = ',a1
                endif
                npnts = -1
                return
             endif

      else
c       less than 3 good points-- do not use
        if( debug) then
           print *,' Less than 3 good points in RATE_EST'
           print *,'   iepoch, iprn, isite, npnts = '
     .            , iepoch, iprn, isite, npnts
        endif
      endif

      return
      end
