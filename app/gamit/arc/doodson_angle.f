      Subroutine Doodson_angle(tjd,Doodson_number,fund_arg,gmst,theta)

c     Calculate the tide-argument angle from Brown's fundamental 
c     arguments (returned from Subroutine 'tide_angles') and the
c     Doodson number.  R. King 170223
                                       
c     In:                              
c       TJD r*8  True Julian date in TDT  
c       Doodson_numeber c*7 ( e.g. '166.554')
c       fund_arg : Brown's fundamental arguments from subroutine
c                  tide_angles (l l' F D Om GMST+pi)  (radians)
c       gmst     : More precise evaluation of GMST 
c     Out:
c       theta    : tide-angle argument (radians)

      implicit none
                                         
      integer*4 Darg(6),i
      real*8 tjd,fund_arg(6),gmst,beta(6),theta,pi,twopi,norm_angle
      character*7 Doodson_number  
      data pi/3.14159265358979d0/

c  Convert the Doodson number to integer coefficients of the Doodson angles
      read(Doodson_number,'(3i1,1x,3i1)') (Darg(i),i=1,6)
      do i=2,6
        Darg(i) = Darg(i) - 5
      enddo            
   
c  Compute the Doodson angles from the Brown angles
 
      call tide_angles( tjd,fund_arg )  
    
c   Brown (returned by tide_angles )
c F(1) el   Moon lon - per   s - w  = mean anomaly l 
c F(2) elp  Sun  lon - per   s' -w' = mean anmalyy l'
c F(3) Moon lon - Asc   s - Omega  = w + l - W  = arg latitude
c F(4) D    Moon elon from Sun  =
c F(5) Om   Moon Asc = Omega 
c F(6_ gmst+pi 
c** Note: We use the input argument GMST rather than fund_arg(6) because the latter
c         has been evaluated with tdj terrestrial dynamical time (TDT), which
c         I think is correct for Brown's theory.  However, to get GMST correctly
c         from this routine, we would have needed to use Universal time (UTC+[UT1-UTC),
c         so for precise use, GMST value (sidtm) from routine rotsnp should be used.
c         rwk 170413

c Doodson from Brown (from calc_acc.cc but with subscripts 1-6 vs 0-5)
c    beta[2] = F[2] + F[4]; // Moon's mean longitude
c    beta[3] = beta[1] - F[3]; // Sun's mean longitude
c    beta[4] = beta[1] - F[0]; // Longitude of Moon's mean perigee
c    beta[5] = -F[4]; // Negative longitude of Moon's mean node
c    beta[6] = beta[1] - F[3] - F[1]; // Longitude of Sun's mean perigee
c    beta[1] = GMST + DPI - beta[2]; // Time angle in lunar days reckoned from lower transit

c Brown     GMST+pi - 2F - 2Om = GMST+pi -2(s-Om) - 2Om  
c                              = GMST+pi -2s 
c Doodson   GMST    - s        = GMST+pi-s  -s = GMST+pi -2s

c Mapping ORB 0-5 to ARC 1-6:
c ORB                                                              ARC 
c beta                                                             beta
c 0 = GMST + pi -beta(1)             GMST+pi-f2-f4 = GMST+pi-F-Om   1    GMST+pi-f3-f5
c 1   f2 + f4        moon long       F+Om                           2    f3+f5
c 2   beta1 - f3     sun mean long   f2+f4-f3 = F+Om-D              3    f3+f5-f4
c 3   beta1 -f0      moon perigee    f2+f4-f0 = F+Om-l              4    f3+f5-f1
c 4   -f4            neg moon node   -Om                            5    -f5
c 5   beta1 - f3 -f1 sun perigee  f2+f4-f3-f1 = F+Om-D-lp           6    f3+f5-f4-f2 

c Example: K1 
c Brown             0 0 0 0  0 1   GMST+pi    
c Doodson 165.555   1 1 0 0  0 0   GMST+pi -s  + s
                  
c Example: O1
c  Brown             0  0 -2 0 -2 1 
c  Doodson 145.555   1 -1  0 0  0 0 
            
      beta(1) = gmst + pi - fund_arg(3) - fund_arg(5)
      beta(2) = fund_arg(3) + fund_arg(5)
      beta(3) = fund_arg(3) + fund_arg(5) - fund_arg(4)
      beta(4) = fund_arg(3) + fund_arg(5) - fund_arg(1)   
      beta(5) = -fund_arg(5)
      beta(6) = fund_arg(3) + fund_arg(5) -fund_arg(4) - fund_arg(2)
                    
cd      print *,'DEBUG fund_arg  beta ' ,fund_arg,beta 
c     Normalize the beta angles
      do i=1,6
        beta(i) = norm_angle(beta(i))
      enddo               

c  Compute the tide argument
      theta = 0.d0
      do i = 1,6
        theta = theta + beta(i)*Darg(i)
      enddo 

cc      if(doodson_number.eq.'255.555'.or.
cc    .   doodson_number.eq.'245.655') then
cc         print *,'arg beta theta ',doodson_number
cc    .      ,(darg(i),beta(i),i=1,6),theta
cc      endif 

      return
      end


c-------------------------------------------------------------------------------------

      Function norm_angle(angle)

      real*8 norm_angle,angle,pi,twopi

      data pi/3.14159265358979d0/
                                                        
      twopi = 2.d0*pi 
      norm_angle = dmod(angle,twopi)  
      if( norm_angle.lt.0.d0 ) norm_angle = norm_angle + twopi

      return
      end 





