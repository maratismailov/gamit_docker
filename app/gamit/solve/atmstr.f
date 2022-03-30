      real*8 function atmstr(rho)
c
c     J.L. Davis  870305
c
c     Routine to calculate the value of the atmospheric structure
c     function for the zenith delay at a separation rho.  The
c     value is returned in units of cm^2.
c
c     Input variables
c     ---------------
c     rho             The separation distance, in km.
c
c     Formula
c     -------
c     D(rho) = (C^2 h^2) rho^(2/3) / [1 + (rho/L)^(2/3)]
c
c     where C is a climate and weather dependent constant (see below
c     for value used in this routine), h is a troposphere "length"
c     (also see below), and rho is the separation, and L is the
c     decorrelation length.
c
c     Constants used
c     --------------
c     The following values are for southern California (see reference).
c     C is given in m^(-1/3) and h in m.  For rho in m, D will be
c     in m^2.  Alpha is the exponent for 2-D Kolmogorov turbulence.
c
      real*8 C, h, alpha, rho, L
c
      data C       / 2.4D-07 /
      data h       / 1.0D+03 /
      data alpha   / 0.666666667D+00 /
      data L       / 3000.0D+00 /
c
c     Reference
c     ---------
c     Treuhaft, R.N., G.E. Lanyi, The effect of the dynamic wet tropo-
c         sphere on radio interferometric measurements, to appear in
c         Radio Science, spring 1987.
c
c.... Calculate D in m^2 for rho in km
      atmstr = (C * h) ** 2 * (1.0D+03 * rho) ** alpha
      atmstr = atmstr / (1.0D+00 + (rho / L) ** alpha)
c
c.... Convert to cm^2
      atmstr = atmstr * 1.0D+04
c
      end
