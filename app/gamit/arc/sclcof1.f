c-------------------------------------------------------------------------------------
Copyright (c) Massachusetts Institute of Technology, 1986/2017. All rights reserved.
      Subroutine SCLCOF1( dir,ztflag,n,m,harcofs ) 
                                          
c     Scale or unscale zonal or tesseral hamronic coefficients for
c     Legendre functions for nomalization to 4pi and assuming the
c     zonal sign convention of C20 (not J2).   This is a modification
c     of subroutine SCLCOF but now accepts the hamronic coefficients
c     as input arguments so that it can be applied to tidal perturbations
c     as well as the static geopotential.  

c     Note: Dimensions currently hard-wired to 12th degreeand order

c     R. King March 10, 2017

c     Input
c       dir  : -1  unscale 
c            :  1  scale 
c       ztflag: 0  harmonics are for zonals
c               1  harmonics are for tesserals
c       n    : degree 
c       m    : order
c      Input and output
c         harcofs(77) : coefficients to degree and order 12 

      implicit none

c**      include '../includes/dimpar.h'   
c**      include '../includes/arc.h'

                              
      integer*4 dir,ztflag,n,m,ff,ih,jh,i,j,mm
      real*8 sclcf,harcofs(77),f
             
cd      print *,'dir ztflag n  m ',dir,ztflag,n,m
                                          
      if( ztflag.eq.0 ) then
c       scale the zonals
        do i=2,n            
          sclcf = dsqrt(dble(2*i+1))           
cd          print *,'i-1 sclcf ',i-1,sclcf
          if( dir.eq.-1 ) harcofs(i-1) = harcofs(i-1) * sclcf
          if( dir.eq. 1)  harcofs(i-1) = harcofs(i-1) / sclcf 
        enddo   
                    
      elseif ( ztflag.eq.1 ) then
c       scale the tesserals 
        jh = 0 
        do i=2,n
          do ih=1,i
            jh=jh+1
            ff = i - ih
            f = ff + 1
            mm = 2 * ih
cd            print *,'i jh ff mm ',i,jh,ff,mm
            do j=2,mm         
              f = f * dble(ff+j)                        
cd              print *,'f ',f
            enddo
            sclcf = dsqrt(dble(4*i+2)/f)
cd            print *,'i sclcf ',i,sclcf
            if( dir.eq.-1 ) harcofs(jh) = harcofs(jh)*sclcf
            if( dir.eq. 1 ) harcofs(jh) = harcofs(jh)/sclcf
          enddo
        enddo
      endif

      return
      end



