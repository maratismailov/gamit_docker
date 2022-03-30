Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
      SUBROUTINE ROTCRD( ROT,ROTD,VIN,VOUT,NVEL,IPAR )
C
C Written by Yehuda Bock and Robert King
C
C Rotate state vector and partial derivatives
c
C         NVEL= 0   Position only
C             = 1   Position and velocity
c
C         IPAR= 1   Coordinates only
C             > 1   Coordinates and IPAR-1 partial derivatives
C
      implicit none

      integer*4 nvel,ipar,indx,i,j,ii,jj
      real*8 rot,rotd,vin,vout,temp1,temp2
      dimension rot(3,3),rotd(3,3),vin(*),vout(*)
     1          ,temp1(3),temp2(3)    

      logical debug/.false./

*  Rotate position and partial derivative coordinates

      if(debug) then 
        print *,'ROTCRD nvel,ipar ',nvel,ipar
        WRITE(*,200) ((ROT(II,JJ),JJ=1,3),II=1,3)
  200   FORMAT(' Rotation Matrix :',/,3(3d15.8,/))
      endif 

       DO J=1,IPAR
         INDX=3*(J-1) + 1
         CALL MATMPY( ROT,VIN(INDX),VOUT(INDX),3,3,1 )
         if(debug) then 
           WRITE(6,300) (vin(ii),ii=indx,indx+2),
     .                (vout(ii),ii=indx,indx+2)
  300      FORMAT(' Input position : ',3d15.8,/,
     .          ' Rotated position : ',3d15.8/)
         endif 
       enddo

*  Rotate velocity
 
      IF( NVEL.GT.0 ) THEN
        CALL MATMPY( ROT,VIN(4),TEMP1,3,3,1)
        CALL MATMPY( ROTD,VIN(1),TEMP2,3,3,1)
CD       write(6,100) (temp1(i),i=1,3)
CD       write(6,101) (temp2(i),i=1,3)
CD 100   format(' TEMP1 : ',3d22.14)
CD 101   format(' TEMP2 : ',3d22.14)
        DO 20 I=1,3
        VOUT(I+3)=TEMP1(I)+TEMP2(I)
   20   CONTINUE
      ENDIF
c
c
      RETURN
      END
