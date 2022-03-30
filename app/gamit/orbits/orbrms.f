      subroutine orbrms(lurms,iref,i2nd,iepoch,nprn,iprn,orb,delt,
     .                  use_runstring, refepoch)
c
c purpose:	compute delta-x/y/z and delta-radial/along/cross at
c		each epoch for each matching satellite, and then obtain
c		the rms' of delta-x/y/z and of delta-along/cross/radial
c		optionally prepare a set of delta serieses for ploting
c
c by Simon McClusky,   JAN, 94
c
c variable: iepoch - number of epoches of the orbits (ref. and 2nd)+1
c		     nprn - number of satellites of the orbits
c		     iprn - prn#'s of the orbits
c	         orb(orbit_A/orbit_B interpolated to orbit_A times,
c		     number of satellites,
c		     number of epoches,
c		     x/y/z and xdot,ydot,zdot components)
c		     delt - time interval (sec)
c		     dab - delta_x/y/z between orbit_A and orbit_B
c		     dac - velocity components of orbit A
c		     dbc - velocity components of orbit A
c		     sigxyz - rms of delta_x/y/z between ref_orbit and 2nd_orbit
c		     sigarc - rms of delta_along/cross/radial between ref & 2nd
c		     sum - sum of delta_x/y/z and delta_along/cross/radial
c		     sum2 - sum of delta_x/y/z squared and delta_a/c/r squared
c		     drab/dalon/dcros - delta-radial/along/cross between ref & 2nd
c		     prms - combined positional errors
* MOD TAH 180106: Added 3-D RMS (See new variables below).
* MOD TAH 200331: Addded passing refepoch (PEPJD) so that MJD can be output on
*       residual lines.
c
        implicit none
c
      include '../includes/dimpar.h'
c
      character*28 pltf
      character*20 ext
      character*3 in_sys, out_sys
      integer lurms,luplt,nprn(2),iprn(2,85),iepoch
* MOD TAH 190625: Made max Pos XYZ and Vel XYZ real*8
      INTEGER l, m,iref,i2nd
      integer iop,isav,i,ii,k,kk,j,maxepc,indx,indx_end
      parameter (maxepc=673)
      real*8 orb(4,maxsat,maxepc,6),delt(2)
      real*8 MAXDAB(MAXSAT,3), MAXDAC(MAXSAT,3)
      real*8 maxdrab(maxsat),maxdalon(maxsat),maxdcros(maxsat)
      real*8 dab(maxsat,maxepc,3),dac(maxsat,maxepc,3)
      real*8 dbc(maxsat,maxepc,3)
      real*8 sigxyz(maxsat,3),sigacr(maxsat,3)
      real*8 sum(2,maxsat,3),sum2(2,maxsat,3)
      real*8 drab(maxsat,maxepc),dalon(maxsat,maxepc)
      real*8 dcros(maxsat,maxepc),prms(2),big,dt
      real*8 sprodr,sproda,sprodc,radlen,vellen,croslen,c1,c2,c3
      real*8 xyzsat_pos(3), llhsat_pos(3), dxyz(3), out_comp(3)
      real*8 orbllh(1,maxsat,maxepc,3), dneu(maxsat,maxepc,3),pi
      real*8 rotmat(3,3),az_diff(maxsat,maxepc),len_diff(maxsat,maxepc)
      real*8 maxdneu(maxsat,3)

* MOD TAH 180106: Added 3-D RMS
      real*8 sum23D(maxsat)   ! Sum of 3D difference squured
      real*8 rms3D(maxsat)    ! Computed 3D RMS
      real*8 sum2mean(8)      ! Sum squared for total, XYZ, ACR, and 3D RMS
      integer*2 num_mean       ! Number of satellites in mean values

* MOD TAH 200331: Added refepoch and mjd
      real*8 refepoch   ! PEPJD of time of epoch 1
      real*8 mjd        ! Computed MJD for time of residual

* MOD TAH 180106: Added batch mode
      logical use_runstring   ! True if runstring is used: PASSED
      character*10 runstring


      data big/999.99/
      parameter ( pi            = 3.1415926535897932D0 )
c
c see if plot file making is required
c
      if( use_runstring ) then
* MOD TAH 201020: Initialize the ext string 
          ext = 'NoPlots'
          call rcpar(6,runstring)
          if( runstring(1:1).ne.' ' ) then
              read(runstring,* ) iop
              if( iop.gt.0 ) then
                 call rcpar(7,ext)
              endif
          else
              iop = 0
          endif
          write(*,'(a,1x,i2,1x," Ext ",a)') 
     .             ' Plot Options',iop, trim(ext)
      else  
         print *, ('Enter the following options')
         print *, ('0  to return to skip plot file making')
         print *, ('1  to make a set of plot files of dx,dy,dz')
         print *, ('2  to make a set of plot files of dr,da,dc')
         print *, ('3  to make a set of plot files of vx,vy,vz')
         print *, ('4  dN,dE,dU, and satellite ground track positions')

        read *, iop

c       ask for a plot file extension
        if(iop.ne.0)then
           write(*,'(a)')
     .             'Enter a plot file extension (ret for no extn):'
           read(*,'(a)')ext
         endif
      end if

      isav=i2nd
c
C CLEAN OUT MAX DIFFERENCE ARRAYS
c
      DO l=1,NPRN(iref)
        MAXDRAB(L)=0
        MAXDALON(L)=0
        MAXDCROS(L)=0
          DO m=1,3
            MAXDAB(L,M)=0
            MAXDAC(L,M)=0
            MAXDNEU(L,M)=0
          ENDDO
      ENDDO

* MOD TAH 180107: Clear summation variables for mean vallues
      num_mean = 0    ! Number of statellites reported
      do k = 1,8      ! Total, XYZ, ACR, 3-D (8 values)
         sum2mean(k) = 0
      end do
c
c loop over each satellite, ii is the index for 2nd orbit
c
      do i=1,nprn(iref)
      ii=1
c
c advance to next satellite until both ref. and 2nd satellites match
c
      do while (iprn(isav,ii).ne.iprn(iref,i).and.ii.le.nprn(isav))
         ii=ii+1
      enddo
c
c clean sum array
c
      do k=1,3
        do kk=1,2
        sum(kk,i,k)=0.
        sum2(kk,i,k)=0.
        enddo
      enddo
      sum23D(i) = 0.0d0  ! Initialize sum**2
c
c see if the satellites of ref. and 2nd match, if not, then skip
c
      if (iprn(iref,i).eq.iprn(isav,ii)) then
c
c loop over each epoch
c
      do j=1,(iepoch-1)
c
c get x/y/z differences between ref. and 2nd orbits, then get the sums
c
      do k=1,3
         dab(i,j,k)=orb(i2nd,ii,j,k)-orb(iref,i,j,k)
c
C GET MAXIMUM X/Y/Z ORBIT DIFFERENCES
c
         IF (abs(DAB(I,J,K)).GT.MAXDAB(I,K)) THEN
            MAXDAB(I,K)=abs(DAB(I,J,K))
         ENDIF

         sum(1,i,k)  = sum(1,i,k)  +dab(i,j,k)
         sum2(1,i,k) = sum2(1,i,k) +dab(i,j,k)*dab(i,j,k)

         sum23D(i)   = sum23D(i)  + dab(i,j,k)**2  ! i is sat number, j epoch
c
c fill satellite velocity arrays
c
         dac(i,j,k)=orb(iref,i,j,k+3)
c
C GET MAXIMUM velocity DIFFERENCES
c
         IF (abs(DAC(I,J,K)).GT.MAXDAC(I,K)) THEN
            MAXDAC(I,K)=abs(DAC(I,J,K))
         ENDIF

         dbc(i,j,k)=orb(i2nd,ii,j,k+3)
      enddo
c
c ***************S. McClusky calculation of radial difference.
c
c calculation of the length of the ref sat position vector.
c
        radlen=dsqrt(orb(iref,i,j,1)**2+orb(iref,i,j,2)**2+
     #              orb(iref,i,j,3)**2)
c
c calculation of the scalar product of the (ref and 2nd orbit difference
c vectors) with the position vector of the ref sat.
c
        sprodr=dab(i,j,1)*orb(iref,i,j,1)+dab(i,j,2)*orb(iref,i,j,2)+
     #         dab(i,j,3)*orb(iref,i,j,3)
c
c radial component differences given by the above scalar product divided by
c the length of the ref sat position vector.
c
        drab(i,j)=sprodr/radlen
c        print*,'radial drab(i,j) = ',drab(i,j)
c
c ****************S. McClusky calculation of the along track difference.
c
c calculation of the length of the ref sat velocity vector.
c
        vellen=dsqrt(dac(i,j,1)**2+dac(i,j,2)**2+dac(i,j,3)**2)
c       if (mod(j,20).eq.0) print*, 'epoch vel-leng ',j,vellen
c
c calculation of the scalar product of the (ref and 2nd orbit difference
c vectors) with the velocity vector of the ref sat.
c
        sproda=dab(i,j,1)*dac(i,j,1)+dab(i,j,2)*dac(i,j,2)
     #           +dab(i,j,3)*dac(i,j,3)
c
c along track differences given by the above scalar product divided by
c the length of the ref sat velocity vector.
c
       dalon(i,j)=sproda/vellen
c       print*,'along dalon(i,j) = ',dalon(i,j)
c
c ****************S. McClusky calculation of cross track difference.
c
c calculation of cross product of ref sat position and
c interpolated ref sat velocity.
c
         c1=orb(iref,i,j,2)*dac(i,j,3)-orb(iref,i,j,3)*dac(i,j,2)
         c2=orb(iref,i,j,3)*dac(i,j,1)-orb(iref,i,j,1)*dac(i,j,3)
         c3=orb(iref,i,j,1)*dac(i,j,2)-orb(iref,i,j,2)*dac(i,j,1)
c
c calculation of the length of the pos/vel cross product vector
c
        croslen=dsqrt(c1**2+c2**2+c3**2)
c       if (mod(j,20).eq.0) print*, 'epoch cross-leng ',j,croslen
c
c calculation of the scalar product of the (ref and 2nd orbit difference
c vectors) with the pos/vel cross product vector.
c
        sprodc=dab(i,j,1)*c1+dab(i,j,2)*c2+dab(i,j,3)*c3
c
c cross track differences given by the above scalar product divided by the
c length of the the pol/vel cross product vector.
c
        dcros(i,j)=sprodc/croslen
c        print*,'cross dcros(i,j) = ',dcros(i,j)

c sum up the deltas of radial, along track, and cross track AND
C DETERMINE THE MAXIMUM RADIAL/ALONG/CROSS ORBIT DIFFERENCES.
c
c PT 950425: use real value rather than integer in order to compute correct plot range

c Original way
c        IF (IABS(INT(DRAB(I,J))).GT.MAXDRAB(I)) THEN
        IF (dabs(DRAB(I,J)).GT.MAXDRAB(I)) THEN
c          MAXDRAB(I)=IABS(INT(DRAB(I,J)))
          maxdrab(i) = dabs(drab(i,j))
        ENDIF

        sum(2,i,1) =sum(2,i,1) +drab(i,j)
        sum2(2,i,1)=sum2(2,i,1)+drab(i,j)*drab(i,j)

c        IF (IABS(INT(DALON(I,J))).GT.MAXDALON(I)) THEN
c          MAXDALON(I)=IABS(INT(DALON(I,J)))
        IF (dabs(DALON(I,J)).GT.MAXDALON(I)) THEN
          MAXDALON(I)=dABS(DALON(I,J))
        ENDIF
        sum(2,i,2) =sum(2,i,2) +dalon(i,j)
        sum2(2,i,2)=sum2(2,i,2)+dalon(i,j)*dalon(i,j)

        IF (dABS(DCROS(I,J)).GT.MAXDCROS(I)) THEN
          MAXDCROS(I)=dABS(DCROS(I,J))
        ENDIF

        sum(2,i,3) =sum(2,i,3) +dcros(i,j)
        sum2(2,i,3)=sum2(2,i,3)+dcros(i,j)*dcros(i,j)
c
c compute satellite position in lat, long, height, convert delta
c XYZ's between orbits into delta NEU's, and compute azimuth and length
c of orbit difference vector. S McClusky..........
c
       if (iop.eq.4) then
c
c make temporary storage for reference satellite position, and delta X,Y,Z between orbits
c
         do l = 1,3
           xyzsat_pos(l) = orb(iref,i,j,l)
           dxyz(l) = dab(i,j,l)
         enddo
         in_sys = 'XYZ'
         out_sys = 'NEU'
c
c convert satellite position X,Y,Z to Lat, Long, Height, and delta X,Y,Z, to delta N,E,U...
c
         call rotate_geod(dxyz,out_comp,in_sys,out_sys,
     .                   xyzsat_pos,llhsat_pos,rotmat)
c
c convert Lat, Long, Height of satellite from radians/meters to degrees/kilometres.
c
         do l = 1,3
           if (l.eq.1) then
           orbllh(iref,i,j,l) = 90.d0 - (llhsat_pos(l)*180.d0/pi)
           elseif(l.eq.2) then
           orbllh(iref,i,j,l) = llhsat_pos(l)*180.d0/pi
           elseif(l.eq.3) then
           orbllh(iref,i,j,l) = llhsat_pos(l)/1000.d0
           endif
           dneu(i,j,l) = out_comp(l)
c
c get max NEU differences
c
           IF (IABS(INT(DNEU(I,J,L))).GT.MAXDNEU(I,L)) THEN
             MAXDNEU(I,L)=IABS(INT(DNEU(I,J,L)))
           ENDIF
         enddo
c
c compute Azimuth and Length of horizintal orbit difference components......
c
c Can't use this line with g77 compiler as the datan2d intrinsic is not supported 
c We have to output the computation in radians and convert to degrees.....
c        az_diff(i,j) = datan2d(out_comp(2),out_comp(1))  
         az_diff(i,j) = datan2(out_comp(2),out_comp(1)) 
         az_diff(i,j) = az_diff(i,j)*180.d0/pi
         if ( az_diff(i,j).lt.0.d0) az_diff(i,j) = 360.d0 + az_diff(i,j)
         len_diff(i,j) = dsqrt( out_comp(2)**2 + out_comp(1)**2 )
       endif
c
c end of epoch loop
c
        enddo
c
c get rms' of x, y, z and of radial, along, cross
c
       do k=1,3
          sigxyz(i,k)=min(big,dsqrt(sum2(1,i,k)/(iepoch-1)))
          sigacr(i,k)=min(big,dsqrt(sum2(2,i,k)/(iepoch-1)))
       enddo
       rms3D(i) = min(big,sqrt(sum23D(i)/(iepoch-1)))
c
c get combined positional rms
c
       prms(1)=min(big,dsqrt((sigxyz(i,1)*sigxyz(i,1)+
     *		sigxyz(i,2)*sigxyz(i,2)+
     *		sigxyz(i,3)*sigxyz(i,3))/3))

       prms(2)=min(big,dsqrt((sigacr(i,1)*sigacr(i,1)+
     *		sigacr(i,2)*sigacr(i,2)+
     *		sigacr(i,3)*sigacr(i,3))/3))

c
c this test may be failed if the ref and 2nd orbits are significantly
c different thus having significantly different velocities.
c The along track difference this case may appear to be in error
c by several centimetres.
c
       if (abs(prms(1)-prms(2)).gt.1.e-2) then
       print *, ("X/Y/Z => Radial/along/cross may be WRONG !")
       endif
c
c write out rms' for current satellite
c MOD TAH 190625: Moved 1x in format make consistent with write_summary.f
c     use in orbfit.
      write(lurms,910) iprn(iref,i),prms(1),
     .     (sigxyz(i,k),k=1,3),(sigacr(i,k),k=1,3), rms3D(i)
910   format(i3,1x, 8(f9.5,1x))

* MOD TAH 180107: Accumulate statitics for MEAN values
       num_mean = num_mean + 1
       sum2mean(1) = sum2mean(1) + prms(1)**2
       do k = 1,3 
          sum2mean(k+1) = sum2mean(k+1) + sigxyz(i,k)**2
          sum2mean(k+4) = sum2mean(k+4) + sigacr(i,k)**2
       end do
       sum2mean(8) = sum2mean(8) + rms3D(i)**2

c
c see if plot file making is needed
c
       if (iop.eq.1.or.iop.eq.2.or.iop.eq.3.or.iop.eq.4) then
c
c convert time interval from day into min
c
       dt=delt(iref)/60
c
c prepare plot output file
c
       luplt=20+i
       if(ext.ne.'                    ')then
         indx = 1
c.... Find next nonblank character
         do while (ext(indx:indx) .eq. ' ' .and. indx .lt. 20)
            indx = indx+1
         enddo
c....   Find the next blank character
         indx_end = indx
         do while (ext(indx_end:indx_end) .ne. ' '
     .            .and. indx_end .lt. 20)
c
c....     Increment end indx
           indx_end = indx_end + 1
c
         end do

c we now have the start and end indexes for the name extension
         write(pltf,'("plt_",a,".",i2.2)')ext(indx:(indx_end-1))
     .         ,iprn(iref,i)
       else
         write(pltf,'("plt.",i2.2)') iprn(iref,i)
       endif
c       if (pltf(5:5).eq." ") write(pltf(5:5),'("0")')
       open(unit=luplt,file=pltf,status='unknown')
c
c prepare the label text
c
       if (iop.eq.1) then
         write(luplt,'("Delta X (m)")')
         write(luplt,'("Delta Y (m)")')
         write(luplt,'("Delta Z (m)")')
       endif
       if (iop.eq.2) then
          write(luplt,'("D-radial (m)")')
          write(luplt,'("D-along (m)")')
          write(luplt,'("D-cross (m)")')
       endif
       if (iop.eq.3) then
         write(luplt,'("V-X (km/s)")')
         write(luplt,'("V-Y (km/s)")')
         write(luplt,'("V-Z (km/s)")')
       endif
       if (iop.eq.4) then
         write(luplt,'("Latitude  (deg)")')
         write(luplt,'("Longitude (deg)")')
         write(luplt,'("Height (km)")')
         write(luplt,'("Delta N (m)")')
         write(luplt,'("Delta E (m)")')
         write(luplt,'("Delta U (m)")')
         write(luplt,'("VECT AZ (deg)")')
         write(luplt,'("VECT LEN (m)")')
       endif
       write(luplt,'("Epoch Number (Time Interval=",f6.2," min)")') dt
       write(luplt,'(i4)') (iepoch-1)
c
c write out delta series
c MOD TAH 190526: Updated formats
       do j=1,(iepoch-1)

* MOD TAH 200331: Compute the MJD of this epoch
          mjd = (refepoch-2400001) + (j-1)*delt(1)/86400.d0

          if (iop.eq.1) then
              IF (J.LE.1) THEN
                  WRITE( LUPLT, '(F9.5," Max X")') MAXDAB(I,1)
                  WRITE( LUPLT, '(F9.5," Max Y")') MAXDAB(I,2)
                  WRITE( LUPLT, '(F9.5," Max Z")') MAXDAB(I,3)
              ENDIF
     	      write(luplt,'(i4,3f11.5,1x,F13.5)') 
     .                     j,(dab(i,j,k),k=1,3), mjd
           endif
           if (iop.eq.2) then
              IF (J.LE.1) THEN

c PT 950425: change I4 format to f7.4
                  WRITE( LUPLT, '(f9.5," MAx R")') MAXDRAB(I)
                  WRITE( LUPLT, '(f9.5," Max A")') MAXDALON(I)
                  WRITE( LUPLT, '(f9.5," MAx C")') MAXDCROS(I)
              ENDIF
     	      write(luplt,'(i4,3f11.5,1x,F13.5)') 
     .                      j,drab(i,j),dalon(i,j),dcros(i,j), mjd
           ENDIF
           if (iop.eq.3) then
               IF (J.LE.1) THEN
                  WRITE( LUPLT, '(F9.5," Max VX")') MAXDAC(I,1)
                  WRITE( LUPLT, '(F9.5," Max VY")') MAXDAC(I,2)
                  WRITE( LUPLT, '(F9.5," Max VZ")') MAXDAC(I,3)
              ENDIF
     	      write(luplt,'(i4,3f11.5,1x,F13.5)') 
     .                      j,(dac(i,j,k),k=1,3), mjd
           endif
           if (iop.eq.4) then
              IF (J.LE.1) THEN
                 WRITE( LUPLT, '(f7.4)') MAXDNEU(I,1)
                 WRITE( LUPLT, '(f7.4)') MAXDNEU(I,2)
                 WRITE( LUPLT, '(f7.4)') MAXDNEU(I,3)
              ENDIF
     	      write(luplt,'(i4,8f14.7,1x,F13.5)') 
     .                      j,(orbllh(iref,i,j,k),k=1,3),
     .                     (dneu(i,j,k),k=1,3),az_diff(i,j),
     .                      len_diff(i,j), mjd
           endif
       enddo

       close (luplt)
c
c end of plot file making if
c
       endif

       else	
       write(*,'("Warning: PRN",i2," not in 2nd Orbit")') iprn(iref,i)
c
c end of matching satellite if
c
       endif
c
c end of satellite loop
c
       enddo

* MOD TAH 180107: Write out the mean values
c MOD TAH 190625: Moved 1x in format make consistent with write_summary.f
c     use in orbfit.
       write(lurms,'(84("-"))')
       write(lurms,920) (sqrt(sum2mean(k)/num_mean),k=1,8)
920    format('MEAN',8(f9.5,1x))
       write(lurms,'(84("-"))')

       return
       end
