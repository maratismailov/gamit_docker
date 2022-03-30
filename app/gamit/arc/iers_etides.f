      Subroutine iers_etides  

c     Set the Doodson numbers and coefficents for the frequency-dependent
c     corrections to the second-order Stokes coefficients.  Values from 
c     IERS Conventions (2010), Table 6.5, retaining only the constituents
c     for which the in-phase (ip) or out-of-phase (op) values are greater 
c     than 1.e-12.  Includes the zonal (C20) and diurnal (C21/S21) tides.
c     The single semi-diurnal term (M2) is coded directly in sbfn.f.
c     R. King 11 April 2017 

      implicit none           
                      
      include '../includes/dimpar.h'
      include '../includes/arc.h' 

c      in common/ef2tides/ in arc.h:
c Stokes harmonic coefficients for frequency-dependent 2nd-degree zonal 
c and (diurnal) solid-Earth tides
c      common/ef2tides/ztidip(6),ztidop(6),dtidip(11),dtidop(11)
c     .               , ztid_doodson(6),dtid_doodson(11)
c      real*8 ztidip,ztidop,dtidip,dtidop
c      character*7 ztid_doodson,dtid_doodson

c     Corrections to C20 
c       unnamed
      ztid_doodson(1) = ' 55.565'
      ztidip(1) = 16.6d-12
      ztidop(1) = -6.7d-12
c       Ssa
      ztid_doodson(2) = ' 57.565'  
      ztidip(2) = -5.5d-12
      ztidop(2) =  4.3d-12
c       Mm 
      ztid_doodson(3) = ' 65.455'      
      ztidip(3) = -1.2d-12
      ztidop(3) =  3.7d-12
c       Mf
      ztid_doodson(4) = ' 75.555'      
      ztidip(4) = 0.6d-12
      ztidop(4) = 6.3d-12
c       unnamed
      ztid_doodson(5) = ' 75.565'   
      ztidip(5) = 0.2d-12
      ztidop(5) = 2.6d-12
c       Mtm 
      ztid_doodson(6) = ' 85.455'  
      ztidip(6) = 0.4d-12
      ztidop(6) = 1.1d-12 


c     Corrections to C21 and S21 (diurnal)
c       unnamed
      dtid_doodson(1) = '145.545'
      dtidip(1) = -1.3d-12
      dtidop(1) =  0.1d-12 
c       O1
      dtid_doodson(2) = '145.555' 
      dtidip(2) = -6.8d-12
      dtidop(2) =  0.6d-12
c       No1
      dtid_doodson(3) = '155.655'   
      dtidip(3) = 1.3d-12
      dtidop(3) = 0.1d-12
c       P1
      dtid_doodson(4) = '163.555' 
      dtidip(4) = -43.4d-12
      dtidop(4) =   2.9d-12
c       unnamed
      dtid_doodson(5) = '165.545'
      dtidip(5) =  -8.8d-12
      dtidop(5) =   0.5d-12
c       K1
      dtid_doodson(6) = '165.555'  
      dtidip(6) =  470.9d-12
      dtidop(6) =  -30.2d-12
c       unnamed 
      dtid_doodson(7) = '165.565'
      dtidip(7) =  68.1d-12
      dtidop(7) =  -4.6d-12
c       unnamed 
      dtid_doodson(8) = '163.575'    
      dtidip(8) =  -1.6d-12
      dtidop(8) =   0.1d-12
c       psi1
      dtid_doodson(9) = '166.554' 
      dtidip(9) =  -20.6d-12
      dtidop(9) =   -0.3d-12
c       phi1
      dtid_doodson(10)= '167.555' 
      dtidip(10) = -5.0d-12
      dtidop(10) =  0.3d-12
c       J1
      dtid_doodson(11) = '175.455' 
      dtidip(11) =   -2.1d-12
      dtidop(11) =    0.1d-12

                                               
      return
      end                           


                                
      

 


