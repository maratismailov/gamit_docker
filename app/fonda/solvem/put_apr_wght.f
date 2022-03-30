      subroutine put_apr_wght(apr_val,id_sit)
c
c     put apriori weighting matrix
c
c     unit:  
c        corodinate: meter
c        velocity:   meter/year (temporary)
c     frame of apriori values
c        1 = xyz
c        2 = enu
c
      implicit real*8(a-h,o-z)
      include 'solvem.fti'
c
      integer i1,is,jd,id,isit1,isit2,id_sit,ier
      integer nrow,ncol
      dimension apr_val(6),cjac1(9),tmp1(9),tmp2(9)
      dimension tmp3(9)
c
c     assign apriori covariance
      if (id_sit.eq.0) then
         isit1 = 1
         isit2 = nsit
      else
         isit1 = id_sit
         isit2 = id_sit
      endif
      nrow = nlive
      ncol = 3
c 
      do 20 is = isit1,isit2
         i1 = jtoi(is)
         if (i1.le.0) goto 20
         if (id_frame.ne.1) then
            a1 = slat(is)
            a2 = slon(is)
            a3 = srad(is)
            call getjac(a2,a1,a3,cjac1,6)
         endif
c     
         id = (i1-1)*6+1
         jd = id
         call zero1d(1,6,tmp1)
         tmp1(1) = apr_val(1)**2
         tmp1(3) = apr_val(2)**2
         tmp1(6) = apr_val(3)**2
         if (tmp1(1).lt.1.0d-6) tmp1(1) = 1.0d-6
         if (tmp1(3).lt.1.0d-6) tmp1(3) = 1.0d-6
         if (tmp1(6).lt.1.0d-6) tmp1(6) = 1.0d-6
         if (id_frame.eq.1) then
            call cholsk(tmp1,tmp2,1,3,ier)
            call fill_apr_mtx(nrow,ncol,id,jd,3,3,aprm,tmp1,2)
         else
            call atwa(3,3,cjac1,tmp1,tmp2,tmp3,1)
            call cholsk(tmp2,tmp3,1,3,ier)
            call fill_apr_mtx(nrow,ncol,id,jd,3,3,aprm,tmp2,2)
         endif
c        
         id = (i1-1)*6+4
         jd = id
         call zero1d(1,6,tmp1)
         tmp1(1) = apr_val(4)**2
         tmp1(3) = apr_val(5)**2
         tmp1(6) = apr_val(6)**2
         if (tmp1(1).lt.1.0d-12) tmp1(1) = 1.0d-12
         if (tmp1(3).lt.1.0d-12) tmp1(3) = 1.0d-12
         if (tmp1(6).lt.1.0d-12) tmp1(6) = 1.0d-12
         if (id_frame.eq.1) then
            call cholsk(tmp1,tmp2,1,3,ier)
            call fill_apr_mtx(nrow,ncol,id,jd,3,3,aprm,tmp1,2)
         else
            call atwa(3,3,cjac1,tmp1,tmp2,tmp3,1)
            call cholsk(tmp2,tmp3,1,3,ier)
            call fill_apr_mtx(nrow,ncol,id,jd,3,3,aprm,tmp2,2)
         endif

c         print 100,sname(is),(apr_val(k),k=1,6)
 100     format (1x,'PUT_APR_WGHT: ',a,6(1x,1pe10.4))
 20   continue

      return
      end

