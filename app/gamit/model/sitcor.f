Copyright (c) Massachusetts Institute of Technology and the University of
California, San Diego. 1995.   All rights reserved.

      subroutine sitcor( jd,t, sun, moon, evec, sitepart,pnsmat
     .                 , eop,polepart,ut1part,simvec,ipass )
C
C     Compute the site coordinates.  If not earth-fixed, rotate to 1950.0
C     R. W. King    2 June 1987
C     modified 1-Oct-93 by peterm to accept sidereal time, sidtm,
c     from sidmat by way of roturt as it is needed in the
c     computation of the earth tides done in etide.
C
c     P. Tregoning  27th March, 1995
c     sun and moon positions as interpolated by solred, lunred are passed in
c     for use in etide to replace the approximate estimates originally
c     coded by PJM

c     P. Tregoning 9th may 1995
c     inertial frame (eg J2000) and precession model (eg IAU76) used in creating
c     the tfile for this model run are passed to routines to compute the
c     terrestrial/inertial rotations    

c     R. King 24th Feb 2000
c     added the ocean tide components to the argument list so that they
c     can be passed from setup (via model) to etide.

c     R. King 16th Feb 2001
c     Added 'ipass' to the argument list to skip tide calculations during the
c     first pass for clock corrections.

c     P. Tregoning 7 January 2004
c     added atmospheric loading displacements to the argument list so that
c     they can be passed from setup (via model) to etide. If no atm loading
c     has been requested then the input loading displacements will have
c     been passed in here as zero.

      implicit none                        
  
      include '../includes/dimpar.h'   
      include '../includes/units.h' 
      include '../includes/global.h'
      include '../includes/model.h'

      logical eop,first_call

      integer*4 jd,jds,ipass,year,doy,nyr
     .        ,ndoy,nt,start_year,nydays,i,j, idir
                           
      character*4 sitcd
      character*256 message

      real*4 tmpatmlod,exttim

      real*8 evec(6,2)
     .     , sitepart(3,3),t,ts,fjd,fjdeop,diff,sidtm
     .     , polepart(3,4),ut1part(3,2),tdtgpst,tdtutc,gpstutc,taiutc
     .     , prec(3,3),rnut(3,3),pnsmat(3,3),pnsdmt(3,3),srot(3,3)
     .     , sdrot(3,3),sidten(3,3),prot(3,3),sun(6),moon(6)
     .     , eqe,xpole,ypole,simvec(3,2)
     .     , tidevec(3),atmvec(3),doy_decimal
* MOD TAH/SCM 190918: Add era variables
      real*8 era,erarot(3,3) 

      logical debug / .false. / 
                                                                 

      data first_call/.true./

 
c     Calculate the precession, nutation, sidereal rotation, and polar motion
c     for rotation from earth-fixed to inertial, and for partials
c
c PT 950418: pass value of iau76 to pnrot,srotat to use appropriate precession
c PT 950509: pass inertial frame and precession model to pnrot,srotat
      tdtgpst = 32.184d0 + 19.d0

* MOD TAH/SCM: Added for SOFIA routines and IAU0A, IAU0C, IAU06 models
      if( precmod .eq. 'IAU76' ) then

         call pnrot( inut,jd,t,tdtgpst,eqe,prec,rnut,frame,precmod)
c        srotat (unique in GAMIT) expects UTC
         jds = jd
         ts = t
         gpstutc = taiutc(jd) - 19.d0
         call timinc(jds,ts,-gpstutc)

c        catch possible leap second
         if( jds.ne.jd ) then
           gpstutc = taiutc(jds) - 19.d0
           jds = jd
           ts = t
           call timinc(jds,ts,-gpstutc)
         endif
         tdtutc = taiutc(jd) + 32.184d0
         call srotat( jds,ts,tdtutc,eqe,isptide,srot,sdrot,sidtm
     .              , xpole,ypole,sidten,prot,iut1,ipole,precmod )
         call pns( prec,rnut,srot,pnsmat )
         call pns( prec,rnut,sdrot,pnsdmt )
      else
*        Use the SOFA routines

         idir = 1
         call rotsnp_sofa( idir,jd,t,tdtgpst,isptide
     .                    , frame,precmod,iut1,ipole,inut
     .                    , pnsmat,pnsdmt,sidtm,xpole,ypole
     .                    , prec,rnut,prot,sidten,era,erarot )   

         if ( precmod .ne. 'IAU76' .and. precmod .ne. 'IAU0A' ) then
c Is this what is required?
           sidtm = era
           sidten = erarot
cd          print*,'CIO based CRS: ',precmod,' : Copying ERA > SIDTM'             
         endif

      endif 

c Debug output
      if (debug) then
        write(*,100) 'SITCOR - prec: ',precmod,
     .             ((prec(i,j),j=1,3),i=1,3) 
        write(*,100) 'SITCOR - rnut: ',precmod,
     .             ((rnut(i,j),j=1,3),i=1,3) 
        write(*,100) 'SITCOR - prot: ',precmod,
     .             ((prot(i,j),j=1,3),i=1,3) 
        write(*,100) 'SITCOR - sidten: ',precmod,
     .             ((sidten(i,j),j=1,3),i=1,3)
        write(*,*) 'SITCOR - sidtm,xp,yp: ',precmod,sidtm,xpole,ypole
        write(*,100) 'SITCOR - pnsmat: ',precmod,
     .             ((pnsmat(i,j),j=1,3),i=1,3) 
        write(*,100) 'SITCOR - pnsdmt: ',precmod,
     .             ((pnsdmt(i,j),j=1,3),i=1,3) 
100     format(a,1x,a,1x,/,3(1x,3D22.14,/))
      endif

c     Calculate Earth-orientation partials

      if (eop) then
c     calculate the time of observation relative to the midpoint of the
c     data span, the reference time for EOP values
        fjd = jd + t/86400.d0 - 0.5d0 + tdtgpst/86400.d0
        fjdeop = jdmidpt + tmidpt/86400.d0 - 0.5d0 + tdtgpst/86400.d0
        diff = fjd - fjdeop
        call eopart( xpole, ypole, sidtm, prec, rnut, prot, sidten
     .             , diff, evec0, polepart, ut1part )
      endif


c     Rotate the coordinates and velocity from earth-fixed to inertial

      do i = 1, 2
         call matmpy( pnsmat,evec0(1,i),evec(1,i),3,3,1 )
         call matmpy( pnsdmt,evec0(1,i),evec(4,i),3,3,1 )  
         if( simulation ) 
     .     call matmpy( pnsmat,simvec0(1,i),simvec(1,i),3,3,1)
      enddo

      if (debug) then
        write(*,110) 'SITCOR - evec: ',precmod,jd,t,evec*1000.d0 
110     format(a,1x,a5,1x,i10,1x,f16.8,1x,2(1x,3f16.5,1x,3f11.5/))
      endif

c     Rotate the partial matrix (Jacobian) to mean of 1950.0

c     x,y,z wrt latitude
      call matmpy(pnsmat,sitepart0(1,1),sitepart(1,1),3,3,1)
c     x,y,z wrt longitude
      call matmpy(pnsmat,sitepart0(1,2),sitepart(1,2),3,3,1)
c     x,y,z wrt radius or height
      call matmpy(pnsmat,sitepart0(1,3),sitepart(1,3),3,3,1)


c         Correct the site coordinates for loading deformations

cd      print*,'SITCOR: ipass and atmlod',ipass,atmlod
      if( ipass.eq. 2 ) then
                         
c      First compute the non-tidal loading from the u-file values (common model.h)
               
         call dayjul(jd,year,doy)   
cd         print *,'SITCOR jd year doy ',jd,year,doy  
c        prevent errors at boundaries due to round-off
         call even_minute(year,doy,t,nyr,ndoy,nt)    
         doy_decimal = float(ndoy) + nt/86400.d0
cd         print *,'doy t nyr ndoy nt doy_decimal '
cd     .       ,doy,t,nyr,ndoy,nt,doy_decimal   
c        check if interpolating beyond the year boundary  
         if( first_call ) then    
           start_year = nyr     
           first_call = .false.
         endif             
cd         print *,'year start_year doy_decimal '
cd     .      ,year,start_year,doy_decimal
         if( nyr.gt.start_year.and.doy_decimal.lt.2. ) then
           doy_decimal = doy_decimal + dfloat(nydays(start_year))
         endif
         do i=1,3
           atmlod(i) = 0.
         end do
         if( ntatml.gt.0) then   
cd           print *,'SITCOR ntatml time values ',ntatml 
cd           do i=1,ntatml  
cd             print *,atml_time(i),(atml_val(i,j),j=1,3)
cd           enddo
           call lininterp( maxatml,3,atml_time,atml_val,doy_decimal
     .                   , ntatml,3,atmlod,exttim )     
cd          print *,'doy_decimal ntatml atml_time exttim '
cd    .        ,doy_decimal,(atml_time(i),i=1,ntatml),exttim
           if( exttim.ne.0. ) then      
              write(message,'(a,f7.3,a,2f8.3,a)') 'Obs time ('
     .        ,doy_decimal,') outside range of ATML values on u-file ('
     .         ,atml_time(1),atml_time(ntatml),' )'
              call report_stat('FATAL','MODEL','sitcor',' ' ,message,0)
           endif
c          change from UNE to NEU for sb etide and c-file   
           tmpatmlod = atmlod(1)
           atmlod(1) = atmlod(2)
           atmlod(2) = atmlod(3)
           atmlod(3) = tmpatmlod
cd           print *,' atmlod NEU ',(atmlod(i),i=1,3)     
c          transformation to inertial XYZ (atmvec) is done in ETIDE
         endif

c      Hydrological loading not yet supported but values on c-file
         hydrlod(1) = 0.
         hydrlod(2) = 0.
         hydrlod(3) = 0.              

c      The variable ietide is a binary coded variable that was read in routine
c      SETUP. It contains the value for the type of tide that needs to be applied.
c      the appropriate values are:
c              1 = Solid earth tides
c              2 = K1 frequency dependant earth tide (not needed for IERS2003 model)
c              4 = Pole tide
c              8 = Ocean Tide  
c             16 = Remove mean from pole tide
c             32 = Atmospheric tide
c
c      There is a  LOGICAL function KBIT(ivariable,ibit) in the library which
c      return the value .true. if ibit is set.  
       call etide( jd,t,sidtm,tdtgpst
     .           , evec,sun,moon,tidevec,atmvec ) 
cd       write(*,'(f10.5,3f9.4)')23.0+t/86400.d0,(atmlod(i)*1.d3,i=1,3)   
cd       print *,'Aft ETIDE evec    ',evec
cd       print *,'          tidevec ',tidevec
cd       print *,'          atmvec  ',atmvec
       do j=1,2
         do i=1,3
           evec(i,j)  = evec(i,j)  + tidevec(i) + atmvec(i)
           if( simulation ) then
             simvec(i,j)  = simvec(i,j)  + tidevec(i) + atmvec(i)
           endif
         enddo
       enddo  
      endif

      continue

      return
      END





