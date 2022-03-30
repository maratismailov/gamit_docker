      subroutine addaux(iterm,coef,dphi1,dphi2,dlam1,dlam2,dfac,mode)
c
c     form coefficients related to CD parameters 
c
c     iterm : size of total coefficients
c     coeff : calculated coefficients
c     mode  : 0 -- 
c             1 --
c
      include 'solvem.fti'

      real*8 coef,dphi1,dphi2,dlam1,dlam2,dfac
      integer iterm,mode
      dimension coef(iterm)
c
      if (iterm.le.12.or.jaux.eq.0) goto 100
c
      coef(13) = (coef(2)*dlam1+coef(8)*dlam2)*dfac
      coef(14) = (coef(2)*dphi1+coef(8)*dphi2)*dfac
      coef(15) = (coef(1)*dlam1+coef(7)*dlam2)*dfac
      coef(16) = (coef(1)*dphi1+coef(7)*dphi2)*dfac
c
 100  continue
      return
      end
c
