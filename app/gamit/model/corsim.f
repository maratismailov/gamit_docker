      subroutine corsim ( pos,offl1,offl2
     .                  , semi,finv,shft,dispneu,simvec0 )
C    
c     R. King 981210, from sb CORTYP; modified for apr-type files and new model.h
c          R. King 030509

c     Calculate the Earth-fixed Cartesian vector for simulated observations
c     from L-file coordinates and a displacement in E,N,U
      
c     Input

c       pos(3)             R*8 : coordinates of monument (x y z) (m)
c       offl1(3)           R*8 : antenna offsets for L1 (U N E)  (m)
c       offl2(3)           R*8 : antenna offsets for L2 (U N E)  (m)
c       semi,finv, shft(3) R*8 : parameters of gdetic datum (m) 
c       dispneu(3)         R*8 : displacement in E, N, U (m)

c     Output

c       simvec0(3,2)       R*8 : coordinates of L1 and L2 phase center (km)


      implicit none
                    
      integer*4 i

      REAL*8 pos(3),FIN(3),OUT(3),SHFT(3),shftkm(3)           
     .     , offl1(3),offl2(3)
     .     , CORMAT(3,3),dispneu(3),simvec0(3,2)
     .     , carg,sarg,semi,finv,cosl,sinl,tmprad,tmplat,tmplon

C-----------------------------------------------------------------------------

c    Get the spherical coordinates to calculate the geodetic coordinates and partials
        
      tmprad = dsqrt(pos(1)**2 + pos(2)**2 + pos(3)**2 ) / 1000.d0 
      tmplon = datan2(pos(2),pos(1))    
      tmplat = dasin(pos(3)/(tmprad*1.d3))  

C     Need to convert to geodetic because we want to add the 
c     displacements in geodetic coordinates
        
      do i=1,3
        shftkm(i) = shft(i)/1000.d0
      enddo
 
c     Loop over L1, L2 antenna offsets

      do i=1,2

        FIN(1)=TMPLAT
        FIN(2)=TMPLON
        FIN(3)=TMPRAD

C       HEIGHT calls GDETIC to convert from geocentric to geodetic
        CALL HEIGHT( FIN,OUT,SEMI/1.d3,FINV,shftkm )

C       Apply the displacements and antenna offsets in geodetic coordinates
C       (need to convert offsets from m to km
C                   
        if( i.eq.1 ) then
          OUT(1)=OUT(1) + (dispneu(1)+offl1(2))*1.D-3/TMPRAD
          OUT(2)=OUT(2) + (dispneu(2)+offl1(3))*1.D-3/TMPRAD
          OUT(3)=OUT(3) + (dispneu(3)+offl1(1))*1.D-3
        else
          OUT(1)=OUT(1) + (dispneu(1)+offl2(2))*1.D-3/TMPRAD
          OUT(2)=OUT(2) + (dispneu(2)+offl2(3))*1.D-3/TMPRAD
          OUT(3)=OUT(3) + (dispneu(3)+offl2(1))*1.D-3
           endif

C      Convert back to geocentric

        CALL GDETIC( OUT,FIN,CORMAT,SEMI/1.d3,FINV,shftkm )

C      Compute the Cartesian Coordinates

        SINL=DSIN(FIN(1))
        COSL=DCOS(FIN(1))
c       coordinate convention East is positive)
        SARG=DSIN(FIN(2))
        CARG=DCOS(FIN(2))
        simvec0(1,i)  =  FIN(3)*COSL*CARG
        simvec0(2,i)  =  FIN(3)*COSL*SARG
        simvec0(3,i)  =  FIN(3)*SINL

      enddo
c     end loop on L1, L2

      RETURN
      END
