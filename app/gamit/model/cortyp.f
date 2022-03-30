      subroutine cortyp ( pos,offl1,offl2,semi,finv,shft
     .                  , evec0,gdlat,gdhgt,sitepart0 )
C
C     Calculate the geodetic latitude and height from Cartesian coordinates 
c     and the Jacobian  D(GEOCENTRIC)/D(INPUT) 

c     Input

c       pos(3)             R*8 : coordinates of monument (x y z) (m)
c       offl1(3)           R*8 : antenna offsets for L1 (U N E) (m)
c       offl2(3)           R*8 : antenna offsets for L2 (U N E) (m)
c       semi,finv, shft(3) R*8 : parameters of gdetic datum   (m)
c        
c     Output

c       evec0(3,2)         R*8 : Cartesian coordinates of L1 and L2 phase center (km)
c       gdlat              R*8 : geodetic latitude in radians
c       gdhgt              R*8 : geodetic height (km)
c       sitepart0(3,3)     R*8 : partial derivatives (dX,Y,Z wrt spherical lat long rad)

      implicit none

      integer*4 i,j

      real*8 pos(3),fin(3,2),out(3,2),shft(3),shftkm(3)
     .     , cormat(3,3),part(3,3),sitepart0(3,3),evec0(3,2)
     .     , carg,sarg,semi,finv,cosl,sinl,tmprad,tmplat,tmplon
     .     , gdlat,gdhgt,offl1(3),offl2(3)

 
c    Get the spherical coordinates to calculate the geodetic coordinates and partials
        
      tmprad = dsqrt(pos(1)**2 + pos(2)**2 + pos(3)**2 ) / 1000.d0 
      tmplon = datan2(pos(2),pos(1))    
      tmplat = dasin(pos(3)/(tmprad*1.d3))  
                                      
c      write(*,'(a,3f16.7)') 'CORTYP pos ',pos
c      print *,' tmprad tmplon tmplat',tmprad,tmplon,tmplat
c        Need to convert to geodetic because we need the height for
c        the atmospheric model and because we want to add antenna
c        offsets in geodetic coordinates
        
      do i=1,3
        shftkm(i) = shft(i)/1.d3
      enddo
 
c    Loop over L1 and L2

      do i=1,2

        fin(1,i)=tmplat
        fin(2,i)=tmplon
        fin(3,i)=tmprad

C       HEIGHT calls GDETIC to convert from geocentric to geodetic
        call height( fin(1,i),out(1,i),semi/1.d3,finv,shftkm ) 
c        print *,'Aft HEIGHT fin out ',(fin(j,i),j=1,3),(out(j,i),j=1,3) 
c        print *,'offl1 offl2 ',offl1,offl2

C       Apply the antenna offsets in geodetic coordinates
C       Convert horizontal offsets to km

        if (i.eq.1) then
          out(1,I)=out(1,I) + offl1(2)*1.D-3/tmprad
          out(2,I)=out(2,I) + offl1(3)*1.D-3/tmprad/dcos(tmplat)
          out(3,I)=out(3,I) + offl1(1)*1.D-3
        else
          out(1,I)=out(1,I) + offl2(2)*1.D-3/tmprad
          out(2,I)=out(2,I) + offl2(3)*1.D-3/tmprad/dcos(tmplat)
          out(3,I)=out(3,I) + offl2(1)*1.D-3
        endif 

c       Convert back to geocentric
                   
cd        print *,'Bef GDETIC out fin ',(out(j,i),j=1,3),(fin(j,i),j=1,3) 
cd        print *,'semi finv shftkm ',semi,finv,shftkm
        call gdetic( out(1,i), fin(1,i), cormat(1,1)
     .             , semi/1.d3, finv, shftkm )    
cd        print *,'AFt GDETIC out fin ',(out(j,i),j=1,3),(fin(j,i),j=1,3)

      enddo

      do i=1,3
        do j=1,3
          cormat(i,j) = 0.d0
        enddo
        cormat(i,i) = 1.d0
      enddo
      gdlat = out(1,1)
      gdhgt = out(3,1)

c ----------------------------------------------------------------------

c        Compute the Cartesian Coordinates
        
      do i = 2, 1, -1
         sinl=dsin(FIN(1,I))
         cosl=dcos(FIN(1,I))
         sarg=dsin(FIN(2,I))
         carg=dcos(FIN(2,I))
         evec0(1,I)  =  fin(3,I)*cosl*carg
         evec0(2,I)  =  fin(3,I)*cosl*sarg
         evec0(3,I)  =  fin(3,I)*sinl   
      enddo    

c        Compute the Jacobian from the L1 phase center

C     PART(I,J)    I=X,Y,Z   J=LAT,LON,R

      part(1,1)=-sinl*carg*fin(3,1)
      part(2,1)=-sinl*sarg*fin(3,1)
      part(3,1)=cosl*fin(3,1)

c     New coordinate convention (East is positive)
      part(1,2)=-cosl*sarg*fin(3,1)
      part(2,2)=cosl*carg*fin(3,1)
      part(3,2)=0.0D0
      part(1,3)= cosl*carg
      part(2,3)= cosl*sarg
      part(3,3)= sinl

      call matmpy(part,cormat,sitepart0(1,1),3,3,3)

      return
      end

