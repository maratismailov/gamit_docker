C Subroutine to calculate the third-order ionospheric correction - editing to add the 2nd order correction.
C EJP May 2007
C EJP Jul 2009 Edited to interpolate Nmax using VTEC instead of STEC
    
      Subroutine iono_corr(VTEC,STEC,vlight,freq,mag_vec,uni_vec
     .                  ,phdel1,phdel2,phtwodel1,phtwodel2,debug)

C In: STEC - slant TEC in TECU
C     vlight - velocity of light in m/s
c     freq(2) -- L1 and L2 frequencies
C Out: IRC - return code

      implicit none

C Shape factor for correction (iono_corr.f)
      Real*8 fac
      Parameter (fac=0.66d0)

C Slant TEC/TECU,Nmax, phase delay, frequency/Hz (L1, L2),slant TEC in m-2
      double precision STEC, Nmax,phdel1,phdel2,freq(2),tec,VTEC,nmaxtec
C Second order phase delay, lightspeed in km/s, magnetic vector components 
C in Tesla, satellite to site unit vector,dotproduct of magnetic vector & unit vector
C Magnitude of magnetic vector
      double precision phtwodel1, phtwodel2, vlight, mag_vec(3),
     . uni_vec(3), vecdot, dot, magsize
      Integer*4 IRC,i,j,debug

C-------------------------------------------------------------------
C Third order phase correction - NB: sign to fit convention of a positive delay
C-------------------------------------------------------------------
C Debug - print shape factor, signal frequency
        if(debug.ge.3) then
      Print*,'MODEL\iono_corr Shape factor',fac
C      Print*,'L1 freq:',freq(1),'L2 freq:',freq(2)
        endif



C Convert slant TEC from TECU to m-2 (1TECU =1x16electrons/m2)
      tec=STEC*1.0D16
C Convert VTEC from TECU to m-2 (1TECU =1x16electrons/m2)
      nmaxtec=VTEC*1.0D16
        if(debug.ge.3) then
      Print*,'MODEL\iono_corr  nmaxtec',nmaxtec
        endif
C Calculate Nmax (peak electron density along signal path/in m-3) with linear
C interpolation, after Fritsche et al. GRL, 2005 
C Oct 2007 Adapt to slightly different version following pers. comm. from Fritsche
C to interpolate properly.

C      Nmax=((20.0-6.0)*1.0E12)/((4.55-1.38)*1.0E18)*tec
      Nmax=(14.0D12/3.17D18)*(nmaxtec-4.55D18)+20.0D12
C For VTEC < ~ 3TECU, Nmax is negative from this interpolation so set to 1 in this case.
C      Nmax = 0.8d0
      if (Nmax.le.1.d0) Nmax = 1.d0

C Calculate phase delay (812.045= (1/8)(Ap)sq), equivalent to 2437/3 from B&H
C Checked with most up-to-date values for constants: use 812.37519. EJP Oct 07

      phdel1=-(812.37519D0*Nmax*fac*tec)/((freq(1))**4)
      phdel2=-(812.37519D0*Nmax*fac*tec)/((freq(2))**4)

        if(debug.ge.2) then
C      Print*,'MODEL\iono_corr  phdel1',phdel1
C      Print*,'MODEL\iono_corr  phdel2',phdel2
      Print*,'MODEL\iono_corr  Nmax: ',Nmax,' TEC/m-2: ',tec
        endif
C-------------------------------------------------------------------
C Second order phase correction - NB: sign to fit convention of a positive delay
C-------------------------------------------------------------------

C Debug - Magnitude of magnetic vector 
      magsize = sqrt(((mag_vec(1))**2)+((mag_vec(2))**2)
     . +((mag_vec(3))**2))

C      find dot product of sat to site unit vector and magnetic field vector
       vecdot=dot(mag_vec,uni_vec)

C Debug
        if(debug.ge.3) then
         Print*, 'MODEL\iono_corr vecdot:',vecdot
         Print*, 'MODEL\iono_corr magsize/teslas: ',magsize
c*  this debug in Liz's dipole version but costheta not used
c        Print*, 'MODEL\iono_corr costheta : ', vecdotother/magsize
         Print *,'MODEL\iono_corr costheta : ',vecdot/magsize
        endif
C 7527 is the value of the constants. From Bassiri & Hajj (1993) but checked.
      phtwodel1=-(7527*vlight*vecdot*tec)/(2*((freq(1))**3))
      phtwodel2=-(7527*vlight*vecdot*tec)/(2*((freq(2))**3))

C Debug
C         Print*, 'MODEL\iono_corr phdeltwo1: ',phtwodel1
C         Print*, 'MODEL\iono_corr phdeltwo2: ',phtwodel2
      End
