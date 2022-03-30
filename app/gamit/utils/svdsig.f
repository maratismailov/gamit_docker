      program svdsig

*=+=+=+=+=+=+=+=+=+=+ MODIFICATIONS by M.Burc Oral =+=+=+=+=+=+=+=+
** the foll line is added (M.Burc Oral)
** input : X(meters)    Y(meters)
**output : a(mm)   b(ppm)

**** sh_globk_repetability_gmt printouts: added fake entries, first and last,
**** so that we have a nice line plotted

*           trying  with different damping parameters
***  rather creates almost an infinite loop
*** it is better not to use this feature for such a simple estimation...
***  I commented it out and put a bound for damping

*** keep the warning off the print out! we want to create plots!
*** commented out!
c         print*,'B negative, using mean'

*scaled the uncertaineties with nrms === sqrt (chi**2 * cov / dof )


*=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+


***   Estimate the precision of the measurements as a function
*     of baseline length according to the general model:
*           sigma^2 = a^2 + b^2*L^2
*     Estimate a and b given sigma(L), where L=length
*     This is a non-linear estimation problem.
*     Units: a (mm), b (mm/km), L (km)
*
*     Note that iteration decisions are somewhat empirically derived.
*     Here's the philosophy.  The cutoff criterion, defined by 'eps',
*     is the maximum allowable difference between successive model
*     adjustments.  It is currently set at 1.d-06, which corresponds
*     to meters and ppm.
*
*     The estimation problem becomes ill-posed
*     when b is near zero (no length dependence), a somewhat common
*     situation, especially for the vertical component.  To get around
*     this, the a priori model is set so that 'a' is the mean value of
*     of the RMS's and 'b' is assumed to be 1 ppb.  The program will
*     assign a stochastic damping parameter if the iterations appear
*     to not be converging, and will increase the damping until it
*     does stabilize.  In most cases, the horizontal components will
*     converge without need of damping in less than 10 iterations.
*     Ill-posed cases will not converge, so maximum iterations is set
*     to 20.  It the number of iterations is less than maxiter and
*     the adjustments are greater than eps, it then checks to make
*     sure the adjustment is not much greater than the previous
*     adjustments, which would indicate divergence; then maximum
*     allowable adjustment is 1.8 times the previous.  If the
*     number of iterations is greater than 10 and all the above
*     tests are negative, the damping parameter is increased by 0.01 and
*     the iterations are started over.  Some particularily degenerate
*     cases may require many iterations, and a large damping parameter.
*     This is a strong indication of non-length dependence.
*
      integer i,numdat,numpar,maxdat,maxpar,iter,maxiter,info
      parameter (maxdat=1000,maxpar=2)
      real*8  amat(maxdat,maxpar),d(maxdat)
      real*8  m(maxpar),chisq,covm(maxpar,maxpar)
      real*8  a,b,a0,b0,eps,sigmak,avesig,damp,prechi,ddot
      real*8  l(maxdat),sigma(maxdat)
      real*8  u(maxdat,maxpar),v(maxpar,maxpar),s(maxpar),e(maxpar)
      real*8  work(maxdat),scale(maxpar),w1,w2,trace,doff
      logical debug

*     minimum iteration steps
      data maxiter,eps,damp/20,1.d-06,0.d0/

      debug = .false.

***   read in data

      do i = 1, maxdat
         read(5,*,end=10) l(i),sigma(i)
         sigma(i) = sigma(i)
         l(i) = l(i)*1.d-06
      enddo

10    numdat = i-1

***   apriori model: a = mean(sigma), b = 0

      avesig = 0.d0
      do i = 1, numdat
         avesig = avesig + sigma(i)
      enddo

      numpar = 2
      a0 = avesig/dble(numdat)
      b0 = 1.d-03

15    continue
      w1 = 999.d0
      w2 = 999.d0
      a  = a0
      b  = b0
      m(1) = 0.d0
      m(2) = 0.d0

***   iteration loop : mk' = mo + Gk~[d - g(mk) + Gk(mk-mo)]

      iter = 0
20    iter = iter+1

*     update partials and data
      do i = 1, numdat
         sigmak = dsqrt(a**2 + b**2 * l(i)**2)
         amat(i,1) = a/sigmak
         amat(i,2) = b/sigmak * l(i)**2
         d(i) = sigma(i) - sigmak
      enddo

      call dgemv('N', numdat, numpar
     .          , 1.d0, amat, maxdat
     .          , m, 1
     .          , 1.d0, d, 1)

*     calculate prefit chi-square
      prechi = ddot(numdat,d,1,d,1)

*     singular value decomposition of partial derivitive matrix
      call dsvdc(amat, maxdat, numdat, numpar, s, e, u, maxdat
     .          , v, maxpar, work, 21, info)

      if (info.ne.0) then
         write(6,'(a,i4)')
     .       'SVDSIG: all singular values are not correct: ',info
         stop
      endif

***   estimate model

*     u'd
      call dgemv('T', numdat, numpar
     .          , 1.d0, u, maxdat
     .          , d, 1
     .          , 0.d0, e, 1)

*     scale by significant singular values
      do i = 1, numpar
         scale(i) = s(i)*s(i) + damp

*        scale the denominator if non-zero
         if (1.d0 + scale(i) .gt. 1.d0) scale(i) = s(i)/scale(i)

*        minimum singular value cutoff
         if (scale(i) .le. eps) scale(i) = 0.d0

         e(i) = scale(i)*e(i)
      enddo

*     vsu'd
      call dgemv('N', numpar, numpar
     .          , 1.d0, v, maxpar
     .          , e, 1
     .          , 0.d0, m, 1)

*     output individual interation information
      if (debug) then
         write(6,'(i4,1p9e11.3)') iter,damp,s(1),s(2)
     .            ,m(1),m(2),a-a0-m(1),b-b0-m(2),a,b
      endif

*     check iteration
      if (iter.le.maxiter
     .    .and. (abs(a-a0-m(1)).gt.eps
     .     .or.  abs(b-b0-m(2)).gt.eps)) then

         if (   (abs(a-a0-m(1)).lt.w1*1.8d0
     .     .and. abs(b-b0-m(2)).lt.w2*1.8d0
     .     .and. iter.lt.maxiter)
     .     .or. iter.lt.20) then

*           keep iterating
            w1 = abs(a-a0-m(1))
            w2 = abs(b-b0-m(2))
            a = a0 + m(1)
            b = b0 + m(2)
            goto 20
         else

*           try with different damping parameter
*** this rather creates almost an infinite loop
*** it is better not to use this feature for such a simple estimation...
*** well, if this feature is not used we have some problems, too.
            damp = damp + 1.d-02
c            goto 15
** lets bound this infinite loop
            if ( damp .lt. 1.d0 ) goto 15
         endif
      endif

      if (b.lt.0.d0) then
*** keep the warning off the print out! we want to create plots!
c         print*,'B negative, using mean'
         b = 0.d0
         a = a0
      endif

***   covariance

*     scale the columns of right singular vectors by the singular values
      do i = 1, numpar
         call dscal(numpar,scale(i),v(1,i),1)
      enddo

*     covm = vv'
      call dgemm('N', 'T', numpar, numpar, numpar
     .          , 1.d0, v, maxpar
     .          , v, maxpar
     .          , 0.d0, covm, maxpar)

***   chisq = d'd - (u'd)'(u'd)

*     scaling for the left singular vectors
      trace = 0.d0
      do i = 1, numpar
         if (1.d0 + scale(i) .gt. 1.d0) then
            scale(i) = dsqrt(s(i)**4 + 2.d0*s(i)*s(i)*damp) /
     .                 (s(i)*s(i) + damp)
            trace = trace + scale(i)*scale(i)
         endif

         call dscal(numdat,scale(i),u(1,i),1)
      enddo

*     u'd
      call dgemv('T', numdat, numpar
     .          , 1.d0, u, maxdat
     .          , d, 1
     .          , 0.d0, e, 1)

      chisq = ddot(numdat,d,1,d,1) - ddot(numpar,e,1,e,1)

***   output results
      doff = numdat - trace
***   output results
*     single line summary
      write (6,'(i3,a,f8.2,a,f8.5,a,f8.5,a,2f15.5,a,2f15.5)')
c     .   iter,' dof: ',numdat-trace
     .   iter,' dof: ',doff
     .,       ' SDUW: ',dsqrt(chisq),' Damp: ',damp
c     .,       ' a (mm) : ', a*1000.d0, dsqrt(chisq*covm(1,1))*1000.d0
c     .,       ' b (ppm): ', b, dsqrt(chisq*covm(2,2))
     .,       ' a (mm) : ', a*1000.d0, dsqrt(chisq*covm(1,1)/doff)*1000.d0
     .,       ' b (ppm): ', b, dsqrt(chisq*covm(2,2)/doff)



*     output data and estimated data
*** this is for repeatability plots
*** a point way out of the scale of x-axis
      write(6,'(3(f20.5,1x))') 0.d0, a*1.d03,  a*1.d03
      write(6,'(i4)') numdat

      do i = 1, numdat
         write(6,'(3(f20.5,1x))') l(i)*1.d03,sigma(i)*1.d03
     .         ,dsqrt(a**2 + (b*l(i))**2)*1.d03
      enddo

*** this is for repeatability plots
*** a point way out of the scale of x-axis
         write(6,'(3(f20.5,1x))') 2. * l(numdat)*1.d03
     .         ,dsqrt(a**2 + (b*2.*l(numdat))**2)*1.d03
     .         ,dsqrt(a**2 + (b*2.*l(numdat))**2)*1.d03



      end


