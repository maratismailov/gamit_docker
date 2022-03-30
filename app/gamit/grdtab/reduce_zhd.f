      subroutine reduce_zhd(slat,slon,dmjd,nod_ht,gps_ht,dzhd)

c  subroutine to correct a ZHD value for a height
c  difference. That is, given a ZHD at a particular height, calculate
c  the ZHD for a different height knowing the height difference and
c  the station latitude. Heights in this subroutine are in metres.
c
c  The subroutine will invoke the GPT to get pressure and temperature
c  at the height of the original ZHD value and also at the 
c  required station height. The returned ZHD value will be corrected
c  for the difference, as computed using formulae of Hopfield and 
c  Saastamoinen. The bilinear interpolation of the four ZHD heights
c  will then be done back in get_map_grid.f at the height of the
c  GPS station.
c
c  P. Tregoning
c  13 December 2006

c  INPUT
c  slat        station latitude  (radians)
c  slon        station longitude (radians)
c  dmjd        modified Julian date
c  gps_ht      station ellipsoidal height
c  dzhd        height correction to zhd value
c  nod_ht      height to which ZHD value refers
c
c  OUTPUT
c  dzhd        ZHD correction required because of height difference

      implicit none

      real*4 slat,slon,nod_ht,gps_ht,dzhd

c  local variables
      real*8 Tk(2),Psur(2),nod_zhd,gps_zhd,gm0,gm,R

c  variables for calling GPT
      real*8 dlat,dlon,hgt(2),dmjd,undu,pi

      pi = 3.141592653d0

      dlat = slat*1.d0
      dlon = slon*1.d0
      hgt(1) = nod_ht*1.d0
      hgt(2) = gps_ht*1.d0

c  get the surface pressure and temperature from the GPT for
c  the height of the grid node and the GPS site
      call gpt(dmjd,dlat,dlon,hgt(1),Psur(1),Tk(1),undu )
c      print*,'pressure and temp from GPT at node are:',Psur(1),Tk(1)    
      call gpt(dmjd,dlat,dlon,hgt(2),Psur(2),Tk(2),undu )
c      print*,'pressure and temp from GPT at stn are:',Psur(2),Tk(2)    

c  now, at each height, calculate the ZHD using the pressure mapping
c  of  Saastamoinen [1972].  These formulas should match those in
c  model/atmdel.f. 
c        gm is at the centroid of the atmsphere; gm0 is adjusted so that 
c        height in the equation is station height  (km/sec**2)
      gm0 = 9.784d-3
c     R  is the specific gas constant for dry air (km**2/sec**2/degC = 28.9644d9 kg/kmol)
      R = 2.89644d-4
c PT 0601229 using old model routines
c      nod_zhd = saaszd(Psur(1),Tk(1),0.d0,"R",dlat,hgt(1)/1.d3)
c      gps_zhd = saaszd(Psur(2),Tk(2),0.d0,"R",dlat,hgt(2)/1.d3) 
c      dzhd = (nod_zhd - gps_zhd)*1.d3   ! output in mm
c      print*,'PT zhd at node and gps are:',nod_zhd,gps_zhd,dzhd
c RWK 070124 code from model/atmdel:** with this from atmdel   
c   (exact match in test)
      gm = gm0*(1.d0 - 0.266D-02 * cos(2.0D0 * dlat) 
     .                    - 0.28D-03 * hgt(1)/1.d3 ) 
      nod_zhd = 0.2277d-2 * Psur(1) * (gm0/gm)   
c      print *,'RWK gm nod_zhd ',gm,nod_zhd
      gm = gm0*(1.d0 - 0.266D-02 * cos(2.0D0 * dlat) 
     .                    - 0.28D-03 * hgt(2)/1.d3 )     
      gps_zhd = 0.2277d-2 * Psur(2) * (gm0/gm)   
c      print *,'RWK gm gps_zhd ',gm,gps_zhd
      dzhd = (nod_zhd - gps_zhd)*1.d3   ! output in mm
c      print *,'dzhd ',dzhd
      return
      end

