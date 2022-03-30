      subroutine lsq_soln(npar)
c
c     perform normal matrix inverse and get least_square solution
c    
      implicit real*8(a-h,o-z)
      include 'solvem.fti'
c
      integer i1,i2,ii,k,ier,npar,i
c
c     solve the normal equation
      print*,' Begin to solve the normal equations ...'
      if (pmode.eq.3) then
         do ii = 1,nlive
         i1 = ii*(ii-1)/2
         print*,'bnorm: ',ii,bnorm(ii)
         write(*,'(6(1pe14.6,1x))') (anorm(k),k=i1+1,i1+ii)
         enddo
      endif

c     copy bnorm for chi square calculation
      call copy1d(1,nlive,0,bnorm,aprm)
      print*,' Scaling the normal matrix ...'
      call nrmscl(anorm,bnorm,scale,nlive,1,1)
      print*,' Cholesky for ',nlive,' parameters ...'
      call cholsk(anorm,bnorm,3,nlive,ier)

c     report problem if ier .ne. 0
      if (ier.ne.0) print*, 'Error: ier=',ier,'  IN LSQ_SOLN'
      print*,' Re-scaling the normal matrix ...'
      call nrmscl(anorm,bnorm,scale,nlive,2,1)
c
c     calculate chi squares
      do ii = 1,nlive
         chi2 = chi2-aprm(ii)*bnorm(ii)
      enddo
      call zero1d(1,nlive,aprm)
      write(*,'(a,f8.2)') ' Calculated chi squares = ',chi2
c     print*,' Calculated chi squares = ',chi2
c
      if (pmode.eq.3) then
      do ii=1,nlive
      i1 = (ii+1)/2
      i2 = ii-i1*2
      f = 1.0d0/finv
      e2 = 2.0d0*f-f*f
      as = radius/(1.0d0-e2*dsin(slat(i1)))
      bn = bnorm(ii)/3.6d3*dtor*as*1.0d3
      if (i2.eq.-1) bn = bn*dcos(slat(i1))
         if (ii.gt.6.and.ii.le.10) bn = bnorm(ii)
         print*,'solution: ',ii,bn
         enddo
         if (smode.gt.0) stop
      endif
c
c     check the matrix quality
c      print*,' check the matrix quality'
c      call chknrm(npar,anorm,scale)
      
      return
      end
