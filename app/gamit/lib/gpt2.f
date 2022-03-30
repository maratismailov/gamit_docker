      subroutine gpt2 (dmjd,dlat,dlon,hell,nstat,it,p,T,dT,e,ah,aw,undu)
*  RWK MOD 201216: File name no longer needed (see below)
*     .                ,gpt_filnam)
          
c  TU-WIEN subroutine from J Boehm 2012-12-11; mods by RWK:
c     - single station only (vectors become scalars)
c     - change name and unit number for grid file 
c     - define dmjd1 (wrt 2000.0) to avoid changing value of dmjd 
*     - more mods 201216 to accommodate moving the file opening
*       and initial reading to the calling routine 


c% This subroutine determines pressure, temperature, temperature lapse rate, 
c% water vapour pressure, hydrostatic and wet mapping function coefficients 
c% ah and aw, and geoid undulation for specific sites near the Earth 
c% surface. It is based on a 5 x 5 degree external grid file ('gpt2_5.grd') 
c% with mean values as well as sine and cosine amplitudes for the annual and 
c% semiannual variation of the coefficients.
c
c% input parameters:
c%
c% dmjd:  modified Julian date (scalar, only one epoch per call is possible)
c% dlat:  ellipsoidal latitude in radians [-pi/2:+pi/2] 
c% dlon:  longitude in radians [-pi:pi] or [0:2pi] 
c% hell:  ellipsoidal height in m 
c% nstat: number of stations in dlat, dlon, and hell
c%        maximum possible: not relevant for Matlab version
C% it:    case 1: no time variation but static quantities
C%        case 0: with time variation (annual and semiannual terms)
c% 
c% output parameters:
c%
c% p:    pressure in hPa  
c% T:    temperature in degrees Celsius 
c% dT:   temperature lapse rate in degrees per km 
c% e:    water vapour pressure in hPa 
c% ah:   hydrostatic mapping function coefficient at zero height (VMF1) 
c% aw:   wet mapping function coefficient (VMF1) 
c% undu: geoid undulation in m 
c%
c% The hydrostatic mapping function coefficients have to be used with the
c% height dependent Vienna Mapping Function 1 (vmf_ht.f) because the
c% coefficients refer to zero height.
c%
c% Example 1 (Vienna, 2 August 2012, with time variation):
c%
c% dmjd = 56141.d0
c% dlat(1) = 48.20d0*pi/180.d0
c% dlon(1) = 16.37d0*pi/180.d0
c% hell(1) = 156.d0
c% nstat = 1
c% it = 0
c%
c% output:
c% p = 1002.56 hPa
c% T = 22.12 deg Celsius
c% dT = -6.53 deg / km
c% e = 15.63 hPa
c% ah = 0.0012647
c% aw = 0.0005726
c% undu = 44.06 m
c%
c% Example 2 (Vienna, 2 August 2012, without time variation, i.e. constant values):
c%
c% dmjd = 56141.d0
c% dlat(1) = 48.20d0*pi/180.d0
c% dlon(1) = 16.37d0*pi/180.d0
c% hell(1) = 156.d0
c% nstat = 1
c% it = 1
c%
c% output:
c% p = 1003.49 hPa
c% T = 11.95 deg Celsius
c% dT = -5.47 deg / km
c% e = 9.58 hPa
c% ah = 0.0012395
c% aw = 0.0005560
c% undu = 44.06 m
c%
c% Klemens Lagler, 2 August 2012
c% Johannes Boehm, 6 August 2012, revision
c% Klemens Lagler, 21 August 2012, epoch change to January 1 2000
c% Johannes Boehm, 23 August 2012, adding possibility to determine constant field
c% ---
        
      implicit none

                 
      integer*4 nstat,it,indx(4),ipod,ipod1,ilon,ilon1,ibilinear,ix,l,n
     .        , ioerr,rcpar,len
                                         
      real*8 dmjd,dlat,dlon,hell,p,T,dT,e,ah,aw,undu
     .     , gm,dMtr,Rg,pi,dmjd1,cosfy,coshy
     .     , sinfy,sinhy,ppod,plon,diffpod,difflon,hgt,heli,redh,Tv,c
     .     , T0,p0,Q,Hs1,dnpod1,dnpod2,dnlon1,dnlon2,R1,R2
     .     , vec(34),pgrid(2592,5),Tgrid(2592,5),Qgrid(2592,5)
     .     , dTgrid(2592,5),u(2592),Hs(2592),ahgrid(2592,5)
     .     , awgrid(2592,5)
     .     , undul(4),Ql(4),dTl(4),Tl(4),pl(4),ahl(4),awl(4)

      character*80 line, prog_name      
* RWK MOD 201216: No longer needed
*       character*16 gpt_filnam
                            
c     get calling program name report_stat
      len = rcpar(0,prog_name)

c%  mean gravity in m/s**2
      gm = 9.80665d0
c% molar mass of dry air in kg/mol
      dMtr = 28.965d-3;
c% universal gas constant in J/K/mol
      Rg = 8.3143d0

      pi = 3.1415926535d0

c% change the reference epoch to January 1 2000
      dmjd1 = dmjd-51544.5;

c% factors for amplitudes
      if (it.eq.1) then  ! constant parameters
        cosfy = 0.d0
        coshy = 0.d0
        sinfy = 0.d0
        sinhy = 0.d0
      else 
        cosfy = dcos(dmjd1/365.25*2*pi)
        coshy = dcos(dmjd1/365.25*4*pi)
        sinfy = dsin(dmjd1/365.25*2*pi)
        sinhy = dsin(dmjd1/365.25*4*pi)
      end if
      
* RWK MOD 201216: Move the file opening and reading of the MIT-added
*                 first line to the calling program
c% read gridfile
c      open(11,file='gpt2_5.grd')  
*      open(41,file='gpt.grid',iostat=ioerr)
*      if( ioerr.ne.0 )  call report_stat('FATAL','MODEL','gpt2'
*     .    ,'gpt.grid','Error opening GPT2 grid file: ',ioerr)

c  read comment lines: so far only two: 
c    line 1: % [filename]  (added by MIT, so far 'gpt2_5.grd' 
c    line 2: % column headers in original file
*      gpt_filnam = ' ' 
*      read (41,'(a80)') line                 
*      if( line(3:4).ne.'gp' ) then
*        call report_stat('WARNING',prog_name,'lib/gpt2',' '
*     .       ,'GPT file name missing from first line of gpt.grid',0)  
*      else
*        gpt_filnam = line(3:18)
*      endif
* The second (column-header) line still read here
      read (41,'(a80)') line       

c% loop over grid points
      do n = 1,2592
        
        !% read data line
        read (41,*) vec
        pgrid(n,1:5)  = vec(3:7)           !% pressure in Pascal
        Tgrid(n,1:5)  = vec (8:12)         !% temperature in Kelvin
        Qgrid(n,1:5)  = vec(13:17)/1000.d0 !% specific humidity in kg/kg
        dTgrid(n,1:5) = vec(18:22)/1000.d0 !% temperature lapse rate in Kelvin/m
        u(n)          = vec(23)            !% geoid undulation in m
        Hs(n)         = vec(24)            !% orthometric grid height in m
        ahgrid(n,1:5) = vec(25:29)/1000.d0 !% hydrostatic mapping function coefficient, dimensionless
        awgrid(n,1:5) = vec(30:34)/1000.d0 !% wet mapping function coefficient, dimensionless

      end do
* RWK MOD 201216: Move the close to the calling program
*      close (41)

c     !% only positive longitude in degrees
      if (dlon.lt.0.d0) then
        plon = (dlon + 2.d0*pi)*180.d0/pi
      else
        plon = dlon*180.d0/pi
      end if
c     transform to polar distance in degrees
      ppod = (-dlat + pi/2.d0)*180.d0/pi 

c     find the index (line in the grid file) of the nearest point
      ipod = floor((ppod+5.d0)/5.d0) 
      ilon = floor((plon+5.d0)/5.d0)
    
c     normalized (to one) differences, can be positive or negative
      diffpod = (ppod - (ipod*5.d0 - 2.5d0))/5.d0
      difflon = (plon - (ilon*5.d0 - 2.5d0))/5.d0
    
      !% get the number of the corresponding line
      indx(1) = (ipod - 1)*72 + ilon
    
      !% near the poles: nearest neighbour interpolation, otherwise: bilinear
      ibilinear = 0
      if ((ppod.gt.2.5d0).and.(ppod.lt.177.5d0)) then 
        ibilinear = 1          
      end if         
    
      !% case of nearest neighbourhood
      if (ibilinear.eq.0) then

        ix = indx(1)
        
        !% transforming ellipsoidial height to orthometric height
        undu = u(ix)
        hgt = hell-undu
            
        !% pressure, temperature at the heigtht of the grid
        T0 = Tgrid(ix,1) + 
     .       Tgrid(ix,2)*cosfy + Tgrid(ix,3)*sinfy + 
     .       Tgrid(ix,4)*coshy + Tgrid(ix,5)*sinhy
        p0 = pgrid(ix,1) + 
     .       pgrid(ix,2)*cosfy + pgrid(ix,3)*sinfy + 
     .       pgrid(ix,4)*coshy + pgrid(ix,5)*sinhy
         
        !% specific humidity
        Q = Qgrid(ix,1) + 
     .      Qgrid(ix,2)*cosfy + Qgrid(ix,3)*sinfy + 
     .      Qgrid(ix,4)*coshy + Qgrid(ix,5)*sinhy
            
        !% lapse rate of the temperature
        dT = dTgrid(ix,1) + 
     .          dTgrid(ix,2)*cosfy + dTgrid(ix,3)*sinfy + 
     .          dTgrid(ix,4)*coshy + dTgrid(ix,5)*sinhy 

        !% station height - grid height
        redh = hgt - Hs(ix)

        !% temperature at station height in Celsius
        T = T0 + dT*redh - 273.15d0
        
        !% temperature lapse rate in degrees / km
        dT = dT*1000.d0

        !% virtual temperature in Kelvin
        Tv = T0*(1.d0 + 0.6077d0*Q)
        
        c = gm*dMtr/(Rg*Tv)
        
        !% pressure in hPa
        p = (p0*exp(-c*redh))/100.d0
        
        !% water vapour pressure in hPa
        e = (Q*p)/(0.622d0 + 0.378d0*Q);
            
        !% hydrostatic coefficient ah 
        ah = ahgrid(ix,1) + 
     .       ahgrid(ix,2)*cosfy + ahgrid(ix,3)*sinfy + 
     .       ahgrid(ix,4)*coshy + ahgrid(ix,5)*sinhy
            
        !% wet coefficient aw
        aw = awgrid(ix,1) + 
     .       awgrid(ix,2)*cosfy + awgrid(ix,3)*sinfy + 
     .       awgrid(ix,4)*coshy + awgrid(ix,5)*sinhy           
                    
      else !% bilinear interpolation
        
        ipod1 = ipod + int(sign(1.d0,diffpod))
        ilon1 = ilon + int(sign(1.d0,difflon))
        if (ilon1.eq.73) then
          ilon1 = 1
        end if
        if (ilon1.eq.0) then
          ilon1 = 72
        end if
        
        !% get the number of the line
        indx(2) = (ipod1 - 1)*72 + ilon   !% along same longitude
        indx(3) = (ipod  - 1)*72 + ilon1  !% along same polar distance
        indx(4) = (ipod1 - 1)*72 + ilon1  !% diagonal
        
        do l = 1,4
                
          !% transforming ellipsoidial height to orthometric height:
          !% Hortho = -N + Hell
          undul(l) = u(indx(l));
          hgt = hell-undul(l);
        
          !% pressure, temperature at the heigtht of the grid
          T0 = Tgrid(indx(l),1) + 
     .         Tgrid(indx(l),2)*cosfy + Tgrid(indx(l),3)*sinfy + 
     .         Tgrid(indx(l),4)*coshy + Tgrid(indx(l),5)*sinhy;
          p0 = pgrid(indx(l),1) + 
     .         pgrid(indx(l),2)*cosfy + pgrid(indx(l),3)*sinfy + 
     .         pgrid(indx(l),4)*coshy + pgrid(indx(l),5)*sinhy

          !% humidity 
          Ql(l) = Qgrid(indx(l),1) + 
     .            Qgrid(indx(l),2)*cosfy + Qgrid(indx(l),3)*sinfy + 
     .            Qgrid(indx(l),4)*coshy + Qgrid(indx(l),5)*sinhy
 
          !% reduction = stationheight - gridheight
          Hs1 = Hs(indx(l))
          redh = hgt - Hs1

          !% lapse rate of the temperature in degree / m
          dTl(l) = dTgrid(indx(l),1) + 
     .             dTgrid(indx(l),2)*cosfy + dTgrid(indx(l),3)*sinfy + 
     .             dTgrid(indx(l),4)*coshy + dTgrid(indx(l),5)*sinhy 

          !% temperature reduction to station height
          Tl(l) = T0 + dTl(l)*redh - 273.15d0

          !% virtual temperature
          Tv = T0*(1.d0+0.6077d0*Ql(l))  
          c = gm*dMtr/(Rg*Tv)
            
          !% pressure in hPa
          pl(l) = (p0*exp(-c*redh))/100.d0
            
          !% hydrostatic coefficient ah
          ahl(l) = ahgrid(indx(l),1) + 
     .             ahgrid(indx(l),2)*cosfy + ahgrid(indx(l),3)*sinfy + 
     .             ahgrid(indx(l),4)*coshy + ahgrid(indx(l),5)*sinhy
            
          !% wet coefficient aw
          awl(l) = awgrid(indx(l),1) + 
     .             awgrid(indx(l),2)*cosfy + awgrid(indx(l),3)*sinfy + 
     .             awgrid(indx(l),4)*coshy + awgrid(indx(l),5)*sinhy
            
      end do
            
        dnpod1 = abs(diffpod)  !% distance nearer point
        dnpod2 = 1.d0 - dnpod1 !% distance to distant point
        dnlon1 = abs(difflon) 
        dnlon2 = 1.d0 - dnlon1
        
        !% pressure
        R1 = dnpod2*pl(1)+dnpod1*pl(2)
        R2 = dnpod2*pl(3)+dnpod1*pl(4)
        p = dnlon2*R1+dnlon1*R2
            
        !% temperature
        R1 = dnpod2*Tl(1)+dnpod1*Tl(2)
        R2 = dnpod2*Tl(3)+dnpod1*Tl(4)
        T = dnlon2*R1+dnlon1*R2
        
        !% temperature in degree per km
        R1 = dnpod2*dTl(1)+dnpod1*dTl(2)
        R2 = dnpod2*dTl(3)+dnpod1*dTl(4)
        dT = (dnlon2*R1+dnlon1*R2)*1000.d0
            
        !% humidity
        R1 = dnpod2*Ql(1)+dnpod1*Ql(2)
        R2 = dnpod2*Ql(3)+dnpod1*Ql(4)
        Q = dnlon2*R1+dnlon1*R2
        e = (Q*p)/(0.622d0+0.378d0*Q)
            
        !% hydrostatic
        R1 = dnpod2*ahl(1)+dnpod1*ahl(2)
        R2 = dnpod2*ahl(3)+dnpod1*ahl(4)
        ah = dnlon2*R1+dnlon1*R2
           
        !% wet
        R1 = dnpod2*awl(1)+dnpod1*awl(2)
        R2 = dnpod2*awl(3)+dnpod1*awl(4)
        aw = dnlon2*R1+dnlon1*R2
        
        !% undulation
        R1 = dnpod2*undul(1)+dnpod1*undul(2)
        R2 = dnpod2*undul(3)+dnpod1*undul(4)
        undu = dnlon2*R1+dnlon1*R2
                    
      end if

      end subroutine
  
