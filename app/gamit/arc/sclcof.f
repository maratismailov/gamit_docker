Copyright (c) Massachusetts Institute of Technology, 1986. All rights reserved.
      Subroutine SCLCOF

c     Unscale the zonal and tesseral harmonic coefficients to match the
c     formulas used in SBFN and SBFN1, from Ash's LL TN 1972-5.
c     We calculate the scale factors (sclcf) for Legendre functions for
c     normalization to four pi, n is the degree and ih is the order
c     Rick Abbot - November 1984

c     Original version had changed the sign of both zonals and tesserals.
c     Modification by R. King, April, 1995, to retain input signs.
c     Zonals will thus be C's (as input), and there is a sign-change
c     wrt Ash's Eq. 103 for variable czone(i) since Ash assumes J's
c     but we have C's.

      implicit none

      include '../includes/dimpar.h'   
      include '../includes/arc.h'

             
      real*8 sclcf,f 
                                                       
      real*8 temp1, temp2,temp3,temp4
      integer*4 n,m,jh,ih,ff,i
            
      do 100 n=2,nczone
        if (czhar(n-1).eq.0.d0) go to 100
        sclcf=dsqrt(dble(2*n+1))
        czhar(n-1) = czhar(n-1)*sclcf
  100 continue

      jh=0
      do 400 n=2,nctess
         do 300 ih=1,n
            jh=jh+1
            if (cchar(jh).eq.0.d0) go to 300
            ff = n - ih
            f = ff + 1
            m = 2 * ih
            do i=2,m
               f = f * dble(ff+i)
            enddo
c**            sclcf = - dsqrt(dble(4*n+2)/f)   Don't change sign any more
            sclcf = dsqrt(dble(4*n+2)/f)
            cchar(jh) = cchar(jh)*sclcf
            cshar(jh) = cshar(jh)*sclcf
  300    continue
  400 continue

      return
      end

