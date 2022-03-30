C E Petrie 2007
C Program that will calculate the coordinates of the 'pierce point'
C through a 'single layer' ionosphere, along a line between a GPS satellite
C and a receiver. Coordinates should be cartesian.
C Adapted to calculate the unit vector from sat to site as well. 

      Subroutine pierce_pt (satcoord, recoord, ionht, eradius, 
     1 ppcoord,unvec,IRC,debug)

C    In:  satcoord (satellite coordinates, geocentric cartesian/km)
C         reccoord (reciever coordinates, geocentric cartesian/km)
C         ionht    (modelled height of single layer ionosphere/km)
C         eradius  (value for average radius of earth /km)
C    Out: ppcoord  (calculated pierce point coordinates, cartesian/km
C         IRC      (returned value, 1 = ok, 0 = error)
C	  unvec	   (unit vector satellite to site geocentric cartesian)

      implicit none

C Declare passed variables
      Real*8 satcoord(3), recoord(3), ionht, eradius, ppcoord(3),
     . unvec(3)
      Integer*4 IRC,debug
      
C Declare used variables
C lambda: fraction of total satellite-receiver vector
C coordiff: satellite coordinates-receiver coordinates
C recsat: |satellite coordinates - receiver coordinates|
C dtty:  dot product of recoord with coordiff
C modrec: |receiver coordinates|
C dot:   variable for function dot which calculates dot product


       Real*8 lambda, coordiff(3), recsat, dtty, modrec, dot
       Integer*4 k,j,i
       Real*8 lambda1, lambda2, lambda3

C Variables for function "mod"

      Real*8 vector(3), modulus


C Calculate satellite-receiver vector
 10   Do i=1,3
      coordiff(i)=satcoord(i)-recoord(i)
      End do
C Debug
        if(debug.ge.3) then
      Print*, 'MODEL\pierce_pt  coordiff', coordiff
      Print*, 'MODEL\pierce_pt  eradius ',eradius, ' ionht: ',ionht 
        endif
C Calculate modulus (size) of coordiff

      recsat=modulus(coordiff)
C Debug
        if(debug.ge.3) then
      Print*, 'MODEL\pierce_pt  mod coordiff ie sat-site dist ', recsat
        endif
C Calculate lambda for height

      dtty=dot(recoord,coordiff)
C Debug
C      Print*,'MODEL\pierce_pt dtty',dtty

      modrec=modulus(recoord)
C Debug
C      Print*, 'MODEL\pierce_pt modulus recoord',modrec

      lambda=(-dtty+sqrt((dtty*dtty)-((recsat*recsat)
     1 *((modrec*modrec)-
     2 (eradius+ionht)*(eradius+ionht)))))/(recsat*recsat)

C Now calculate coordinates at pierce point using lambda
      Do j=1,3
      ppcoord(j)=recoord(j)+lambda*coordiff(j)
      Enddo

C Calculate unit vector sat-site

      Do k=1,3
       unvec(k)=(-coordiff(k))/recsat
      Enddo
C Debug    
      if( debug.ge.2 ) then
       Print*, 'MODEL\pierce_pt unit vector', unvec
       Print*, (unvec(1)**2)+(unvec(2)**2)+(unvec(3)**2)
C (tests that magnitude of unit vector =1)
       endif
      End

      Function modulus(vector)
      Real*8 modulus, vector(3)
      modulus=sqrt((vector(1)*vector(1))+(vector(2)*vector(2))+
     1(vector(3)*vector(3)))




      End



