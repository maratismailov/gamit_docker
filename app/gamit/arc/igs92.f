Copyright (c) Massachusetts Institute of Technology and University of
California at San Diega, 1986-1994. All rights reserved.

      Subroutine IGS92

c     Initialize constants according to IERS/IGS standards
c     R. King and Y. Bock  May 1994

      implicit none       

      include '../includes/dimpar.h'   
      include '../includes/arc.h'

      integer*4 i

      do i=1,10
        czhar(i)=0.d0
      enddo
      do i=1,54
        cchar(i)=0.d0
        cshar(i)=0.d0
      enddo

c        IERS Standards for Constants - coded from NOAA iers.f and
c        then correct GM from .4418 to .4415 for consistency with other
c        IGS centers. Ray/King 950920

      ltvel= 299792.458D0
      aunit= 1.4959787061D8
      ertrad= 6378.1363D0  
      monrad = 1738.d0
      sunrad= 696000.D+00
C               GM OF THE EARTH
      gm(1)= 398600.4415D+0
cold  gm(1)= 398600.4418D+0
C               GM OF THE MOON
      gm(2)= 4902.7989D+0
C               GM OF THE SUN
      gm(3)= 1.32712440D+11
C
      crad= ertrad

c-------------------------------------------------------------------------
c         ** former MIT code - effectively equivalent
c         LTVEL=299792.458D0
c         AULTSC=499.00478370D0
c         ERTRAD=6378.137D0
c         SUNRAD=6.96D+05
c         CONST=(0.01720209895D0)**2*(LTVEL)**3*(AULTSC)**3/(86400.D0)**2
c               GM OF THE EARTH
c         Value used until 9103012: SEM = 328900.553  => GM = 398600.4445
c              (Proper MERIT value: SEM = 328900.5501 << GM = 398600.448)
c         New value from CSR8901  : SEM = 328900.5564 << GM = 398600.4401
c         IERS/IGS 1992 stds value: SEM = 328900.5552 << GM = 398600.4418
c         GM(1)= CONST*81.3005883D0/(1.D0+81.3005883D0)*1.D0/328900.55523496d0
c               GM OF THE MOON
c         GM(2) = GM(1)/81.3005883D0
c               GM OF THE SUN
c         GM(3) = (GM(1)+GM(2))*328900.55523496d0
c----------------------------------------------------------------------------

C         Gravity field is GEM-T3
c         scaling between 2-body GM and GM for GEM-T3
c**     no, this appears not to be used
c**     gmfct = 398600.436d0/gm(1)
      crad = ertrad

C     Zonals  J2,J3,J4,J5,J6,J7,J8
C
      czhar(1)=-484.1650994D-06
      czhar(2)=  +0.9572011D-06
      czhar(3)=  +0.5395212D-06
      czhar(4)=  +0.0683433D-06
      czhar(5)=  -0.1495135D-06
      czhar(6)=  +0.0913009D-06
      czhar(7)=  +0.0488832D-06
      czhar(8)=  +0.0268624D-06

C     Sectorials and Tesserals

C     (2,1), (2,2)

      cchar(1)=  -0.17D-09
      cshar(1)=  +1.19D-09
      cchar(2)=  +2.4390658D-06
      cshar(2)=  -1.4000946D-06

C        (3,1), (3,2), (3,3)

      cchar(3)=  +2.0277142D-06
      cshar(3)=  +0.2492171D-06
      cchar(4)=  +0.9044707D-06
      cshar(4)=  -0.6194477D-06
      cchar(5)=  +0.7203425D-06
      cshar(5)=  +1.4138845D-06
C
C        (4,1), (4,2), (4,3), (4,4)
C
      cchar(6)=  -0.5361511D-06
      cshar(6)=  -0.4734360D-06
      cchar(7)=  +0.3502181D-06
      cshar(7)=  +0.6630152D-06
      cchar(8)=  +0.9909337D-06
      cshar(8)=  -0.2009274D-06
      cchar(9)=  -0.1887706D-06
      cshar(9)=  +0.3094237D-06
C
C        (5,1), (5,2), (5,3), (5,4), (5,5)
C
      cchar(10)= -0.0582802D-06
      cshar(10)= -0.0960839D-06
      cchar(11)= +0.6527110D-06
      cshar(11)= -0.3238637D-06
      cchar(12)= -0.4523301D-06
      cshar(12)= -0.2152958D-06
      cchar(13)= -0.2955841D-06
      cshar(13)= +0.0496903D-06
      cchar(14)= +0.1737635D-06
      cshar(14)= -0.6689070D-06
C
C        (6,1), (6,2), (6,3), (6,4), (6,5), (6,6)
C
      cchar(15)= -0.0768942D-06
      cshar(15)= +0.0269984D-06
      cchar(16)= +0.0487345D-06
      cshar(16)= -0.3740131D-06
      cchar(17)= +0.0572032D-06
      cshar(17)= +0.0093728D-06
      cchar(18)= -0.0868265D-06
      cshar(18)= -0.4713064D-06
      cchar(19)= -0.2673304D-06
      cshar(19)= -0.5367802D-06
      cchar(20)= +0.0096846D-06
      cshar(20)= -0.2371348D-06
C
C        (7,1), (7,2), (7,3), (7,4), (7,5), (7,6), (7,7)
C
      cchar(21)= +0.2748687D-06
      cshar(21)= +0.0974659D-06
      cchar(22)= +0.3277950D-06
      cshar(22)= +0.0932467D-06
      cchar(23)= +0.2512201D-06
      cshar(23)= -0.2152927D-06
      cchar(24)= -0.2755610D-06
      cshar(24)= -0.1237672D-06
      cchar(25)= +0.0013262D-06
      cshar(25)= +0.0186200D-06
      cchar(26)= -0.3588314D-06
      cshar(26)= +0.1517387D-06
      cchar(27)= +0.0009703D-06
      cshar(27)= +0.0240836D-06
C
C        (8,1), (8,2), (8,3), (8,4), (8,5), (8,6), (8,7), (8,8)
C
      cchar(28)= +0.0236282D-06
      cshar(28)= +0.0588472D-06
      cchar(29)= +0.0775985D-06
      cshar(29)= +0.0660087D-06
      cchar(30)= -0.0177852D-06
      cshar(30)= -0.0863470D-06
      cchar(31)= -0.2463398D-06
      cshar(31)= +0.0701796D-06
      cchar(32)= -0.0250411D-06
      cshar(32)= +0.0894628D-06
      cchar(33)= -0.0649237D-06
      cshar(33)= +0.3091226D-06
      cchar(34)= +0.0674622D-06
      cshar(34)= +0.0750948D-06
      cchar(35)= -0.1241984D-06
      cshar(35)= +0.1201722D-06
C
      return
      end













































