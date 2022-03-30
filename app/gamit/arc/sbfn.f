Copyright (c) Massachusetts Institute of Technology and the University of
California at San Diego, 1994.  All rights reserved.
      Double Precision Function sbfn(k,j,s)
c
c     Evaluation of right side of satellite-probe differenial equations
c     Ash/Amuchastegui/friedman - June 1969
c     Rick Abbot  October 1984, modification for gps satellites
c     Documentation and improved structure by Mark Schenewerk incoporated by R. King - April 199
c     Solid Earth tides added by Jim Ray and R. King - April 1995
c     Paul Tregoning May 1995, allow choice of Inertial frame and precession model
c     Add 3rd and 4th degree terms to the solid-Earth tides, add the soild-Earth
c       and ocean pole tides, and add the ocean tides from an external file. R. King April 2017
c
      implicit none
      include '../includes/dimpar.h'        
      include '../includes/units.h'   
      include '../includes/global.h'
      include '../includes/arc.h'

      dimension fn(2)

c     Input Parameters:
c     -----------------
c           k  i   equation number
c           j  i   iteration number
c           s  i   time (seconds past epoch of ics)
c     In common:
c         frame    inertial frame
c         precmod  precession model

c     Output Parameters:  Function value (rhs of eqn) only
c     ------------------

c     Global variables and constants:
c     -------------------------------
c     bcor()      coordinates sun -> earth satellite
c     ccor()      coordinates sun -> earth
c     kount       number of partial derivatives = n/6-1
c     pbcor()     coordinates earth satellite -> moon
c     pccor()     coordinates earth -> moon  
c     vccor()     coordinates earth -> venus
c     vbcor()     coordinates earth satellite -> venus
c     jccor()     coordinates earth -> jupiter
c     jbcor()     coordinates earth satellite -> jupiter
c     rb          magnitude sun -> earth satellite
c     rb2         magnitude squared sun -> earth satellite
c     rb3         magnitude cubed sun -> earth satellite
c     rc          magnitude sun -> earth
c     rc2         magnitude squared sun -> earth
c     rc3         magnitude cubed sun -> earth
c     rpb         magnitude earth satellite -> moon
c     rpb2        magnitude squared earth satellite -> moon
c     rpb3        magnitude cubed earth satellite -> moon   
c     rpc         magnitude earth -> moon
c     rpc2        magnitude squared earth -> moon
c     rpc3        magnitude cubed earth -> moon
c     rsb         magnitude earth -> earth satellite
c     rsb2        magnitude squared earth -> earth satellite
c     rsb3        magnitude cubed earth -> earth satellite
c     rvb         magnitude earth satellite -> venus
c     rvb3        magnitude cubed earth satellite -> venus
c     rjb         magnitude earth satellite -> jupiter
c     rjb3        magnitude cubed earth satellite -> jupiter
c     rsbh()      earth harmonic quantities
c     sbcor()     coordinates earth -> earth satellite

c RWK Note 180302: With the addition of Venus and Jupiter, it would be more
c logical to replace the null for sun and 'p' (planet) for moon with 's',
c and 'm', or to make the 'cor' arrays doubly dimensioned with 1 for sun, 2
c for moon, 3 for venus, and 4 for jupiter.  However, to minimize coding
c errors, I'm keeping the original names for now.  Additionally,  we could
c sacrifice some computational efficiency in favor of fewer variables by
c retaining in common only the four primary vectors in SBFN and recomputing 
c the body-differences and scalars in SBFN1, reusing the secondary variable
c names when computing the accelerations from each body.  

c     Local Variables and Constants:
c     ------------------------------
c     dum1        dummy variable
c     dum2        dummy variable
c     dum3        dummy variable
c     dum4        dummy variable
c     i           loop counter
c     l           linear array index
c     idx         linear array index 
c     idum1       dummy variable
c     idum2       dummy variable
c     vec1()      dummy vector
c     vec2()      dummy vector
c     vec3()      dummy vector
c     vec4()      dummy vector
c     vec5()      dummy vector
c     vec6()      dummy vector

c
c    This module called by:  start_int, adam
c
c    This module calls:      ertorb, legnd2, legndr, lunred, roturt,
c                            solred
c
c    Include files used:
c
c    common blocks used:     const, coraux, ertaux, hrmaux, harcof,
c                            incon, output, paraux, polmot, timxtr
c
c    References:  Ash, M.E., Deterination of Earth Satellite Orbits
c                   MIT Lincoln Lab Tech. Note 1972-5, 19 April, 1972.

      character*6 lowerc
                                                 
c       These reserved for time, equation and interation numbers:
      integer*4 j,k,kk,kkk       
      real*8 s 
c       These local variables:
      integer*4 idum1,idum2,nctp1,nczp1,ntopc,imode,ispeed,ll1
     .        , jd,jds,idir,iut1pol,iarg,ih,i,i1,l,idx
      real*8 dummy,rpc2,fn,rc2,cclatr,rsbh1,t,taiutc,fjd,fjdtmp
     .     , rpb,tcur,ts,sidtm,gpstutc,tdtgpst
     .     , slatm,clatm,slats,clats
     .     , xeqm,xeqs,gmor3m,gmor3s,gmst
     .     , clons,slons,slonm,clonm,coef2,coef20,coef21,coef22
     .     , C_mlon_m(4),S_mlon_m(4),C_mlon_s(4),S_mlon_s(4)
     .     , Pn_mon(3),Pnm_mon(5),Pn_sun(3),Pnm_sun(5) 
      real*8 fund_arg(6),tdj
      real*8 delc20,delc21, dels21,delc22,dels22 
      real*8 cmoon2,cmoon3,csun2,csun3,erad3,erad4,rc4,rpc4,rvec(3)
      real*8 m1,m2,convds                                  

c ertorb,sbfn,sbfn1,shadow 
c Sun, Moon, Venus, Jupiter and satellite coordinates - ertorb,sbfn,sbfn1,shadow 
c         (rsbh dimensioned for a 12-degree gravity field)
c     common/coraux/sbcor(6),bcor(3),ccor(6),ccor3(3),pccor(6)
c    .    , pccor3(3),pbcor(3),vccor(3),vbcor(3),jccor(3),jbcor(3)
c    .    , rsb,rsb2,rsb3,rb,rb2,rb3,rc,rc3,rpc,rpc3,rpb2,rpb3
c    .    , rvc,rvb,rjc,rjb,rsbh(12)
c     real*8 sbcor,bcor,ccor,ccor3,pccor,pccor3,pbcor,vccor,vbcor
c    .    , jccor,jbcor,rsb,rsb2,rsb3,rb,rb2,rb3,rc,rc3,rpc,rpc3
c    .    , rpb2,rpb3,rvc,rvb.rjc,rjb, rsbh

c** For debugging 
      logical kbit,debug/.false./,norad/.false./
     .      , noetide/.false./,nootide/.false./,nodiurnal/.false./
     .      , nopermtide/.true./,nopoletide/.false./
*     .      , nopermtide/.true./,nopoletide/.true./
     .      , fcheck 
      real*8 brown,theta,tjd,tempdb 
      character*7 Doodson_number   
      real*8 klove20, klove21, klove22
      real*8 delc21old(6),dels21old(6) 
      logical firstcall/.true./                                             
      logical diag /.false./   ! New diagnostic logical to output some debug
      integer*2 ndiag / 0 /    ! Only allow 1 count of diagonotic output 
      save ndiag                             
c**  
c     local for zero-ing out tesserals: values are the index of the 1st degree higher than requested
      integer*4 tessindx(11)/1,3,6,10,15,21,28,35,45,55,66/

      data convds/4.8481368110953599d-6/
             
cd      print *,'Entering SBFN k j s ',k,j,s 
cd      if( s.eq.-75.00 ) then
cd        debug = .true.
cd        print *,'  set debug T  s j k ',s,j,k 
cd      else
cd        debug = .false.  
cd      endif

* DEBUG TAH 171226: set diag based on epoch                          
* RWK 190501: Don't print this unless 'diag' is initialized 'true'
      if(diag) then 
        if( s.eq.0.d0 .and. ndiag.lt.4 ) then
           diag = .true.
           ndiag = ndiag + 1
        else
           diag = .false.
        endif
      endif 
                
c         Explicit initialization

      fn(1) =0.0d0
      fn(2) =0.0d0

c         Determine time quantities

      jd = jde
      tcur = te
      call timinc( jd,tcur,s )
      fjd = jd + tcur/86400.d0
c*old fjd    =  jde + (te+s)/86400.d0
      if( debug ) print *,'SBFN 1 fjd ',fjd


c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c         Compute the quantities needed for accelerations
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      if ( k.eq.4 ) then
c         k = 4  => first velocity component-- rhs is acceleration
c         Evaluate quantities for all rhs components each iteration

        if( j.ne.2 ) then
c         j = 1  => first iteration of integrator
c         Need to evaluate the perturbing body and Earth-rotation quantities only once per step

c         Get PEP JD for planetary ephemerides interpolation
          fjd    =  jde + (te+s)/86400.d0 
          if( debug ) print *,'SBFN 2 fjd s/m tdtoff  ',fjd,tdtoff

c
c         Determine Earth relative to Sun
          ispeed=1
          if( fcheck('nbody') ) then 
            call ephtrp(fjd+tdtoff/86400.d0,3,ispeed,ccor)
            if( debug ) print *,'called EPHTRP 3 t ccor '
     ,          ,fjd+tdtoff/86400.d0,ccor
          else
            call solred(ispeed,fjd+tdtoff/86400.d0,ccor)
            if( debug ) print *,'called SOLRED t ccor '
     .           ,fjd+tdtoff/86400.d0,ccor
          endif 
          rc2=ccor(1)**2+ccor(2)**2+ccor(3)**2
          rc =dsqrt(rc2)
          rc3=rc2*rc
          do i=1,3
            ccor3(i)=ccor(i)/rc3
          enddo

c         Determine Moon relative to Earth
          ispeed=0
          if( fcheck('nbody') ) then 
            call ephtrp(fjd+tdtoff/86400.d0,10,ispeed,pccor)
            if( debug ) print *,'called EPHTRP 10 t pccor '
     ,          ,fjd+tdtoff/86400.d0,pccor
          else
            call lunred(ispeed,fjd+tdtoff/86400.d0,pccor)  
            if( debug ) print *,'called LUNRED t pccor '
     .           ,fjd+tdtoff/86400.d0,pccor
          endif 
          rpc2  =   pccor(1)**2+pccor(2)**2+pccor(3)**2
          rpc   =   dsqrt(rpc2)
          rpc3  =   rpc2*rpc
          do i=1,3
            pccor3(i)=pccor(i)/rpc3
          enddo     

c         Determine rotation matrix for earth
c         --- Don't invoke sub-daily EOP terms since the effect is very small
C         iut1pol = 0 
c         iut1pol = 11
* MOD TAH 200505: Read the sestbl. to get the correct value   
C         iut1pol = 0   ! Original code (reverse note above about not using). 
          call get_iut1pol( iut1pol )

          idir = -1
c         ROTSNP expects GPST so convert if ARC running in UTC
c         No need to account for leap seconds since ARC doesn't honor them
          jds = jd
          ts = tcur
          if( time_type.eq.'UTC ') then
            gpstutc = taiutc(jd) - 19.d0
            call timinc(jds,ts,gpstutc)
            tdtgpst = tdtoff - gpstutc
            if( debug)   print *,'SBFN jds ts gpstutc tdtgpst '
     .                            ,jds,ts,gpstutc
          else
            tdtgpst = tdtoff
            if( debug )  print *,'time_type ', time_type
          endif
          if(debug)  write(*,'(a,1x,a5,1x,a5)')' SBFN fr model '
     .                       ,frame,precmod
          call rotsnp( idir,jds,ts,tdtgpst,iut1pol
     .               , frame,precmod,iut1,ipole,inut
     .               , cntrot,cntrtd,sidtm,xpole,ypole ) 
          fjdtmp = fjd+tdtoff/86400.d0
          if( debug ) then
            print *,'SBFN 3 fjd ',fjd 
            print *,'frame,precmod ',frame,precmod
            print *,'jd,tcur,tdtgpst,sidtm,cntrot: '
     .       , jd,tcur,tdtgpst,sidtm,cntrot
          endif
        endif
c       end quantities calculated only on the first iteration
        if(debug) print *,'Start of SBFN s j y(1,1-4) '
     .     ,s,j,(y(1,i),i=1,4) 
 

c       Calculate satellite position quantities each interation

c       Determine earth satellite relative to earth
        do i=1,6
          sbcor(i)=y(i,j)
        enddo
        rsb2=sbcor(1)**2+sbcor(2)**2+sbcor(3)**2
        rsb =dsqrt(rsb2)
        rsb3=rsb2*rsb

c       Determine earth satellite relative to sun
        do i=1,3
c          bcor = (earth wrt sun) + (satellite wrt earth)
           bcor(i)=ccor(i)+sbcor(i)
        enddo
        rb2=bcor(1)**2+bcor(2)**2+bcor(3)**2
        rb=dsqrt(rb2)
        rb3=rb2*rb   
        if(debug) print *,'SV-Sun ',bcor           


c       Determine moon relative to earth satellite
        do i=1,3
c         pbcor = (moon wrt earth) - (satellite wrt earth)
          pbcor(i)=pccor(i)-sbcor(i)
        enddo
        rpb2  =   pbcor(1)**2+pbcor(2)**2+pbcor(3)**2
        rpb   =   dsqrt(rpb2)
        rpb3  =   rpb2*rpb               
        if(debug) print *,'SV-Moon ',pbcor    


c        --latitude of satellite in body-fixed coordinates
        cslat= (cntrot(3,1)*sbcor(1) + cntrot(3,2)*sbcor(2) +
     1        cntrot(3,3)*sbcor(3) )/ rsb
        cclat = dsqrt(1.0d0-cslat**2)
        do  i=1,3
          cslat1(i) = cntrot(3,i)-sbcor(i)*cslat/rsb
c** trial for s-zonal terms
          cclat1(i) = cntrot(3,i)-sbcor(i)*cclat/rsb 
        enddo
        cclatr = cclat*rsb
        rsbh1 = rsb/crad
        rsbh(1) = rsbh1**2
        ntopc = max0(nczon1,nctes1) 
        do  i=2,ntopc
          rsbh(i) = rsbh(i-1)*rsbh1
        enddo      

c       --longitude
        cslng(1) = (cntrot(2,1)*sbcor(1) + cntrot(2,2)*sbcor(2) +
     1            cntrot(2,3)*sbcor(3) )/ cclatr
        cclng(1) = (cntrot(1,1)*sbcor(1) + cntrot(1,2)*sbcor(2) +
     1            cntrot(1,3)*sbcor(3) )/ cclatr       
        do i=2,nctess
          cslng(i)= cslng(i-1)*cclng(1) + cclng(i-1)*cslng(1)
          cclng(i)= cclng(i-1)*cclng(1) - cslng(i-1)*cslng(1)
        enddo    
        do i=1,3
          clng1(i)=( cntrot(2,i)*cclng(1)-cntrot(1,i)*cslng(1) )/cclat
        enddo     
               

         if(debug) print *,'SV coords ',sbcor 
c         Determine Venus and Jupiter relative to Earth and the satellite
          if( fcheck('nbody').and.lbody.gt.0 ) then   
            ispeed = 0 
            call ephtrp(fjd+tdtoff/86400.d0,2,ispeed,rvec)
            if(debug) print *,'SBFN 4 fjd ',fjd 
            if( debug ) print *,'called EPHTRP 2 t rvec ',rvec
            do i=1,3
              vccor(i) = rvec(i) - ccor(i) 
              vbcor(i) = rvec(i) - bcor(i)
            enddo
            rvc = dsqrt(vccor(1)**2+vccor(2)**2+vccor(3)**2)        
            rvb = dsqrt(vbcor(1)**2+vbcor(2)**2+vbcor(3)**2) 
            if( debug ) then
               print *,'Venus wrt Earth ',vccor,rvc
               print *,'Venus wrt SV    ',vbcor,rvb 
            endif 
            call ephtrp(fjd+tdtoff/86400.d0,5,ispeed,rvec)
            if(debug) print *,'SBFN 5 fjd ',fjd
            if( debug ) print *,'called EPHTRP 5 t rvec ',rvec 
            do i=1,3
              jccor(i) = rvec(i) - ccor(i) 
              jbcor(i) = rvec(i) - bcor(i) 
            enddo                
            rjc = dsqrt(jccor(1)**2+jccor(2)**2+jccor(3)**2)  
            rjb = dsqrt(jbcor(1)**2+jbcor(2)**2+jbcor(3)**2)  
           if( debug ) then
               print *,'Jupiter wrt Earth ',jccor,rjc
               print *,'Jupiter wrt SV    ',jbcor,rjb 
            endif 
          endif   


c       Determine the unnormalized Legendre functions and the partials for the 
c       gravity field terms in sbfn1:  p(n), p'(n), p(n,h), p'(n,h)
        call legndr(cslat,cclat,nczone,nctess,cleg,cleg1,cgleg,cgleg1)
        if( debug ) then
          print *,'cslat cleg',cslat,(cleg(i),i=1,3)
          print *,'cleg1',(cleg1(i),i=1,4)
        endif    

c       Determine quantities for solid-Earth tides
                                                
c          sine and cosine of geocentric latitudes of Moon, Sun in body-fixed coordinates
        slatm = (cntrot(3,1)*pccor(1) + cntrot(3,2)*pccor(2) +
     .           cntrot(3,3)*pccor(3) )/ rpc
        clatm = dsqrt(1.d0-slatm**2)
        slats = -(cntrot(3,1)*ccor(1) + cntrot(3,2)*ccor(2) +
     .            cntrot(3,3)*ccor(3) )/ rc
        clats = dsqrt(1.d0-slats**2)

c          sine, cosine of geocentric east longitudes of Moon, Sun in body-fixed coordinates
        xeqm = rpc * dsqrt( 1.0d0 - slatm**2 )
        xeqs = rc * dsqrt( 1.0d0 - slats**2 )
        slonm = (cntrot(2,1)*pccor(1) + cntrot(2,2)*pccor(2) +
     .           cntrot(2,3)*pccor(3) )/ xeqm
        slons = - (cntrot(2,1)*ccor(1) + cntrot(2,2)*ccor(2) +
     .              cntrot(2,3)*ccor(3) )/ xeqs
        clonm = (cntrot(1,1)*pccor(1) + cntrot(1,2)*pccor(2) +
     .            cntrot(1,3)*pccor(3) )/ xeqm
        clons = - (cntrot(1,1)*ccor(1) + cntrot(1,2)*ccor(2) +
     .             cntrot(1,3)*ccor(3) )/ xeqs
        if( debug ) then
           write(*,'(a,4d12.4)') 'clonm slonm clons slons '
     .           ,clonm,slonm,clons,slons
        endif                


c       Determine the normalized Legendre functions used in the tide code 
  
        call legndr_p( slatm,clatm,3,3,Pn_mon,Pnm_mon )
        call sclcof1( -1,0,3,3,Pn_mon )
        call sclcof1( -1,1,3,3,Pnm_mon )
        call legndr_p( slats,clats,3,3,Pn_sun,Pnm_sun )
        call sclcof1( -1,0,3,3,Pn_sun )
        call sclcof1( -1,1,3,3,Pnm_sun )
        if( debug ) then
           write(*,'(a,4d12.4)') 'Pn_mon',(Pn_mon(i),i=1,3)
           write(*,'(a,7d12.4)') 'Pnm_mon',(Pnm_mon(i),i=1,7)   
           write(*,'(a,4d12.4)') 'Pn_sun',(Pn_sun(i),i=1,3)
           write(*,'(a,7d12.4)') 'Pnm_sun',(Pnm_sun(i),i=1,7)
        endif 
                   
c          Intermediate products
        erad3  = ertrad * ertrad * ertrad
        erad4 = erad3*ertrad 
        rpc4 = rpc3*rpc
        rc4 = rc3*rc 
        gmor3m = gm(2) / rpc3
        gmor3s = gm(3) / rc3
        C_mlon_m(1) = 1.d0
        S_mlon_m(1) = 0.d0
        do i=2,4
          C_mlon_m(i) = clonm * C_mlon_m(i-1) - slonm * S_mlon_m(i-1)
          S_mlon_m(i) = slonm * C_mlon_m(i-1) + clonm * S_mlon_m(i-1)
        enddo
        C_mlon_s(1) = 1.d0
        S_mlon_s(1) = 0.d0
        do i=2,4
          C_mlon_s(i) = clons * C_mlon_s(i-1) - slons * S_mlon_s(i-1)
          S_mlon_s(i) = slons * C_mlon_s(i-1) + clons * S_mlon_s(i-1)
        enddo   
        if( debug ) then
          write(*,'(a,3d12.4)') 'C_mlon_m ',(C_mlon_m(i),i=2,4)
          write(*,'(a,3d12.4)') 'C_mlon_s ',(C_mlon_s(i),i=2,4)
        endif 
 
c           Step 1:  Stokes coefficients with frequency-independent Love numbers

c       see IERS Standard 2010 10 August 2012 revision, Section 6.2                     
        cmoon2 = (erad3/rpc3) * (gm(2)/gm(1))
        csun2  = (erad3/rc3)  * (gm(3)/gm(1))
        cmoon3 = (erad4/rpc4) * (gm(2)/gm(1))
        csun3  = (erad4/rc4)  * (gm(3)/gm(1))
c          C20 
        cztid(1) = (cmoon2/5.d0) * Pn_mon(1) * k2mr(1) 
     .         +   (csun2/5.d0)  * Pn_sun(1) * k2mr(1)   
c          C21, C22, S21, S22
        do i=1,2
          cctid(i) = (cmoon2/5.d0) * Pnm_mon(i) * 
     .              (k2mr(i+1)*C_mlon_m(i+1) + k2mi(i+1)*S_mlon_m(i+1))
     .            +   (csun2/5.d0)  * Pnm_sun(i) * 
     .              (k2mr(i+1)*C_mlon_s(i+1) + k2mi(i+1)*S_mlon_s(i+1))  
           cstid(i) = (cmoon2/5.d0) * Pnm_mon(i) * 
     .               (k2mr(i+1)*S_mlon_m(i+1) - k2mi(i+1)*C_mlon_m(i+1))
     .           +    ( csun2/5.d0) * Pnm_sun(i) * 
     .               (k2mr(i+1)*S_mlon_s(i+1) - k2mi(i+1)*C_mlon_s(i+1)) 
        enddo        
c          C20
        cztid(2) = (cmoon3/7.d0) * Pn_mon(2) * k3m(1) 
     .           + (csun3 /7.d0) * Pn_sun(2) * k3m(1)
c          C3m, S3m  m=1,3
        do i=1,3    
          idx = i + 2                        
          i1 = i + 1 
          cctid(idx)=(cmoon3/7.d0)*k3m(i1) * Pnm_mon(idx) * C_mlon_m(i1)
     .            +  (csun3/7.d0) *k3m(i1) * Pnm_sun(idx) * C_mlon_s(i1)
          cstid(idx)=(cmoon3/7.d0)*k3m(i1) * Pnm_mon(idx) * S_mlon_m(i1)
     .            +  (csun3/7.d0) *k3m(i1) * Pnm_sun(idx) * S_mlon_s(i1)
        enddo        
c          C40                         
        cztid(3) = (cmoon2/5.d0)*k2mp(1) * Pn_mon(1) 
     .           + (csun2 /5.d0)*k2mp(1) * Pn_sun(1) 
c         C4m, S4m   m=1,4
* COM TAH 1712228: The degree 4 terms due to degree 2 term.  Only m=0,1,2
*       terms effected and k2m(+) Love number used.
        do i=1,2                          
          i1= i + 1
          idx = i + 5
          cctid(idx)= (cmoon2/5.d0)*k2mp(i1) * Pnm_mon(i) * C_mlon_m(i1)
     .           +    (csun2 /5.d0)*k2mp(i1) * Pnm_sun(i) * C_mlon_s(i1)
          cstid(idx)= (cmoon2/5.d0)*k2mp(i1) * Pnm_mon(i) * S_mlon_m(i1)
     .           +    ( csun2/5.d0)*k2mp(i1) * Pnm_sun(i) * S_mlon_s(i1)
        enddo
        if( debug .or. diag ) then   
           call print_csnm('SETIDE','(normalized', 4, cztid, 
     .                      cctid, cstid)
        endif  

c           Step 2: Correct for deviations of the k20 and k21 Love number 
c                   for several constituent tides in the diurnal and 
c                   semi-diurnal bands
             
c       Get Brown's astronomical angular arguments plus GMST+pi
        call tide_angles( fjd-0.5d0+tdtoff/86400.d0, fund_arg )
        if(debug) print *,'SBFN 7 fjd ',fjd
c       rwk 170414: tide_angles has been called with dynamical time, but
c       GMST should use Universal time (UTC + [UT1-UTC]), so subsitute
c       the value obtained from routine rotsnp above. 
c          gmst = fund_arg(6) - (twopi/2.d0)           
        gmst = sidtm 
        fund_arg(6) = gmst + (twopi/2.d0)
        tjd = fjd-0.5d0+tdtoff/86400.d0
        if(debug) write(*,'(a,f15.7,7f12.7)') 
     .    'E-tides epoch fund_arg gmst ',tjd, (fund_arg(i),i=1,6),gmst 
                                      
c         Zonal tides (k20, C20)  - 6 constituents > 1.e-12 (see iers_etides.f)
        delc20 = 0.d0   
        do i=1,6
          call doodson_angle(tjd,ztid_doodson(i),fund_arg,gmst,theta)
          delc20= delc20 + ztidip(i)*dcos(theta) - ztidop(i)*dsin(theta)
        enddo   
        if( nodiurnal ) then 
           delc20 = 0.d0
           if( firstcall ) print *
     .        ,'SBFN: Dirunal and semi-diurnal tides omitted '
        endif
        if( debug ) print *,'freq dep zonal delc20 ',delc20
        cztid(1) = cztid(1) + delc20 
c         Tesseral tides (k21, C21, S21) 11 constituents > 1.e-12 (see iers_etides.f)
        delc21 = 0.d0
        dels21 = 0.d0     
        if( debug ) print *,'C21/S21 freq-dependent normalized '
        do i=1,11 
          call doodson_angle(tjd,dtid_doodson(i),fund_arg,gmst,theta)
          delc21 =delc21 + dtidip(i)*dsin(theta) + dtidop(i)*dcos(theta)
          dels21 =dels21 + dtidip(i)*dcos(theta) - dtidop(i)*dsin(theta)
          if( debug ) then                                      
            print *,'delc21 for freq-dep ',dtid_doodson(i)
     .                ,dtidip(i)*dsin(theta) + dtidop(i)*dcos(theta)
            print *,'dels21 for freq-dep ',dtid_doodson(i)
     .               ,dtidip(i)*dcos(theta) - dtidop(i)*dsin(theta)
          endif
        enddo 
        if (nodiurnal ) then
          delc21 = 0.d0
          dels21 = 0.d0
        endif           
        if(debug) print *,'diurnal delc21 dels21 ',delc21,dels21 
        cctid(1) = cctid(1) + delc21
        cstid(1) = cstid(1) + dels21                     
c         Semi-diurnal tide (k22, C22, S22) only 2 constituents (M2 and N2) > 2.e-13
        call doodson_angle(tjd,'255.555',fund_arg,gmst,theta)
        delc22 = -1.2d-12*dcos(theta)
        dels22 =  1.2d-12*dsin(theta)                        
        if( debug ) print *,'255.555 ',theta,' 1.2'
        call doodson_angle(tjd,'245,655',fund_arg,gmst,theta)
        delc22 = delc22 - 0.3d-12*dcos(theta)
        dels22 = dels22 + 0.3d-12*dsin(theta)      
        if( debug ) print *,'245.655 ',theta,' 0.3 '
        if( debug ) print *,'semi-diurnal delc22 dels22 ',delc22,dels22
        cctid(2) = cctid(2) + delc22
        cstid(2) = cstid(2) + dels22 
        if( debug ) then 
          print *,'C20 w/ freq-dependence normalized ',cztid(1)             
          print *,'C21 w/ freq-dependence normalized ',cctid(1)
          print *,'S21 w/ freq-dependence normalized ',cstid(1)
          print *,'C22 w/ freq-dependence normalized ',cctid(2)
          print *,'S22 w/ freq-dependence normalized ',cstid(2)
        endif
                      

c          Permanent tide   (C20)
 
        if( zero_tide .or. .not.nopermtide ) then
c         From NGS orb with m1,m2 in arc-seconds:
          cztid(1) = cztid(1) +  4.4228d-8*(-0.31460)*k2mr(1) 
          if( debug ) then 
             print *,'Permanant tide C20 correction '
     .           , 4.4228d-8*(-0.31460)*k2mr(1)      
             print *,'C20 w/ permanent normalized ',cztid(1) 
          endif    
        else  
          if( debug.and.firstcall ) print *
     .        ,'SBFN: Permanent tide not applied'
        endif   
             
                
c       Solid-Earth pole tide (IERS Section 6.4)
                    
        if( .not.nopoletide ) then 
c           dC21 =  -1.333d-9*(m1 + 0.0115d0*m2) 
c           dS21 =  -1.333d-9*(m2 - 0.0115d0*m1    
c           where m1 = xp-xpm, m2 = -(yp-ypm) in arc-seconds
          m1 =  (xpole-xpm)/casr 
          m2 = -(ypole-ypm)/casr     
          cctid(1) = cctid(1) -1.333d-9*(m1 + 0.0115d0*m2) 
          cstid(1) = cstid(1) - 1.333d-9*(m2 - 0.0115d0*m1)
          if(debug) then 
            print *,'xpole xpm m1 ypole ypm m2'
     .            , xpole,xpm,m1,ypole,ypm,m2
            print *,'SE   pole-tide dC21 dS21 ',
     .       -1.333d-9*(m1 + 0.0115d0*m2), - 1.333d-9*(m2 - 0.0115d0*m1)
            print *,'C21, S21 w/ pole tide normalized '
     .             ,cctid(1),cstid(1)
          endif 
        endif

*       Debug:
        if( diag ) then
            print *,'NoPole tide: ',nopoletide
            print *,'xpole xpm m1 ypole ypm m2'
     .            , xpole,xpm,m1,ypole,ypm,m2
            print *,'SE   pole-tide dC21 dS21 ',
     .       -1.333d-9*(m1 + 0.0115d0*m2), - 1.333d-9*(m2 - 0.0115d0*m1)
            print *,'C21, S21 w/ pole tide normalized '
     .             ,cctid(1),cstid(1)
        endif    

c       convert to unnormalized coefficients 
        call sclcof1(-1,0,2,2,cztid)               
        call sclcof1(-1,1,2,2,cctid )
        call sclcof1(-1,1,2,2,cstid )
       
c       Turn off etides for debugging
        if( noetide ) then
          do i=1,nczon1   
           cztid(i)=  0.d0
          enddo               
          do i=1,nctes2 
            cstid(i) = 0.d0
            cctid(i) = 0.d0
          enddo
          if(firstcall ) print *,'SBFN: solid-Earth tides off '
        endif                
c**     replace the above code with the following: zeroes out coefficients
c       higher than the degree requested (all if etidedeg=0). This 'backward'
c       logic necessary because the various terms are lumped above 
        if( etidedeg.lt.4 ) then
          if( etidedeg.lt.2 ) then
            i1 = 1
          else
            i1 = etidedeg 
          endif                                                     
          do i=i1,nczon1
            cztid(i) = 0.d0
          enddo    

          if( etidedeg.lt.2 ) then
            i1 = 1
          else
            i1 = tessindx(etidedeg)
          endif       
          do i=i1,nctes2
            cstid(i) = 0.d0
            cctid(i) = 0.d0
          enddo 
        endif 

c       Determine quantities for the ocean tides           
 
        do i=1,nczon1                
          czcotid(i) = 0.d0
          czsotid(i) = 0.d0
        enddo
        do i=1,nctes2
          ctcotid(i) = 0.d0
          ctsotid(i) = 0.d0
          ctcotid(i) = 0.d0
          ctsotid(i) = 0.d0
        enddo         
c            Using the mapping in egm08.f, we have
c            degree 2-12 --> cz[cs]otid(1-11)
c            2, 1-2  ct[cs]otid(1-2)
c            3, 1-3  ct[cs]otid(3-5)
c            4, 1-3  ct[cs]otid(6-9)
c            5, 1-5  ct[cs]otid(10-14)
c            6, 1-6  ct[cs]otid(15-20)
c            7, 1-7  ct[cs]otid(21-27)
c            8, 1-8  ct[cs]otid(28-35)
c            9  1-9  ct[cs]otid(36-44)
c           10 1-10  ct[cs]otid(45-54)
c           11 1-11  ct[cs]otid(55-65)
c           12 1-12  ct[cs]otid(66-77)

        if( debug ) then  
          call print_csnm('NEWTIDE','normalized',4,cztid,cctid,cstid)
       endif

c       Ocean pole tide 
c         dC21 =  -2,1778d-10*(m1 - 0.01724d0*m2) 
c         dS21 =  -1.7232d-10*(m2 - 0.03365d0*m1)
        ctcotid(1) = ctcotid(1) -2.1778d-10*(m1 - 0.01724d0*m2)
        ctsotid(1) = ctsotid(1) -1.7232d-10*(m2 - 0.03365d0*m1)
        if(debug) then
          print *,'Ocean pole-tide dC21 dS21 '
     .           , -2.1778d-10*(m1 - 0.01724d0*m2) 
     .           , -1.7232d-10*(m2 - 0.03365d0*m1)
          print *,'C21, S21 w/ ocean pole tide normalized '
     .           ,ctcotid(1),ctsotid(1)
        endif
c         Ocean tides for 18 waves
        call tide_angles( fjd-0.5d0+tdtoff/86400.d0, fund_arg )
        if( debug ) print *,'SBFN 8 fjd ',fjd 
        do iarg=1,18
          call Doodson_angle(tjd,otid_doodson(iarg),fund_arg,gmst,theta) 
* MOD TAH 171226: Replaced code based on analysis of EOT11a 
*         do i=1,nczon1     
*           czcotid(i) = czcotid(i) + ozcp(iarg,i)*cos(theta) 
*    .                              + ozsp(iarg,i)*sin(theta)       
c    Tom's formula has these all zero:
c            czsotid(i) = czsotid(i) + ozsp(iarg,i)*cos(theta)
c     .                              - ozcp(iarg,i)*sin(theta) 
*           if( debug .and. i.eq.1 ) then
*             print *,'Ocean C20 iarg doodson theta ozcp ozsp dH ',iarg
*    .         ,otid_doodson(iarg),theta,ozcp(iarg,i),ozsp(iarg,i)
*    .         ,ozcp(iarg,i)*cos(theta)+ozsp(iarg,i)*sin(theta)        
* 
*             print *,'      S20 dH   '
*    .          , ozsp(iarg,i)*cos(theta)-ozcp(iarg,i)*sin(theta)        
*           endif 
*         enddo                      
*         do i=1,nctes2
c            ctcotid(i) = ctcotid(i) +
c     .                 (otcp(iarg,i)+otcm(iarg,i))*cos(theta)
c     .               + (otsp(iarg,i)+otsm(iarg,i))*sin(theta) 
c            ctsotid(i) = ctsotid(i) + 
c     .                 (otsp(iarg,i)-otsm(iarg,i))*cos(theta)
c     .               - (otcp(iarg,i)-otcm(iarg,i))*sin(theta) 
c    Altenative formula for fes2004_Cnm-Snm.dat suggested by Tom:
*           ctcotid(i) = ctcotid(i) +
*    .                   otcp(iarg,i)*cos(theta)+otsp(iarg,i)*sin(theta)
*           ctsotid(i) = ctsotid(i) +
*    .                   otcm(iarg,i)*cos(theta)+otsm(iarg,i)*sin(theta)
*          enddo   
*          Zonals first.  The IERS 2010 fes2004_Cnm-Snm.dat file, the 
*          DelC+ and DelS+ are the zonal Cos and Sin for Cnm; Snm values
*          are zero.
           do i=1,nczon1     
              czcotid(i) = czcotid(i) + ozcp(iarg,i)*cos(theta) 
     .                                + ozsp(iarg,i)*sin(theta)
cd              if( debug .and. i.le.8 ) 
cd     .            print *,'Ocean Zonal iarg doodson theta ozcp ozsp ',
cd     .                     i+1, iarg ,otid_doodson(iarg),theta,
cd     .                     ozcp(iarg,i), ozsp(iarg,i)
              if( debug .and. iarg.eq.10 ) then   ! O1 
                  print *,'Zonal Doodsin C ', i+1, iarg,
     .                 ozcp(iarg,i), ozsp(iarg,i)     
              endif 
           end do

*          Now do the terrereral terms.  Here we sum and difference
*          the + and - terms.
*          Colmns are DelC+ DelS+ DelC- DelS-
*                     otcp  otsp  otcm  otsm   ! GAMIT naes
*                       1     2     3     4
           do i=1,nctes2
              ctcotid(i) = ctcotid(i) +
     .                 (otcp(iarg,i)+otcm(iarg,i))*cos(theta)
     .               + (otsp(iarg,i)+otsm(iarg,i))*sin(theta) 
              ctsotid(i) = ctsotid(i) + 
     .                 (otsp(iarg,i)-otsm(iarg,i))*cos(theta)
     .               - (otcp(iarg,i)-otcm(iarg,i))*sin(theta)  
    
              if( debug .and. (iarg.eq.10 .or. i.le.5 )) then   ! O1 coefficients/2nd degree
*                   i <= 5 will be 2,1 2,2 3,1 3,2 3,3
                  write(*,899) i, iarg,  
     .                 otcp(iarg,i), otsp(iarg,i),     
     .                 otcm(iarg,i), otsm(iarg,i),
     .                 (otcp(iarg,i)+otcm(iarg,i))*cos(theta)
     .               + (otsp(iarg,i)+otsm(iarg,i))*sin(theta), 
     .                 (otsp(iarg,i)-otsm(iarg,i))*cos(theta)
     .               - (otcp(iarg,i)-otcm(iarg,i))*sin(theta) 
 899              format('Non-Zonal ',2I4,' dCS+ ',2E15.6,' dCS- ',
     .                    2E15.6,' dCnm dSnm ', 2E15.6) 


              endif 
           end do                                                    
        enddo                         
c       Turn off ocean tides for debugging
        if( nootide ) then
          do i=1,nczon1                
            czcotid(i) = 0.d0
            czsotid(i) = 0.d0
          enddo
          do i=1,nctes2
            ctcotid(i) = 0.d0
            ctsotid(i) = 0.d0
            ctcotid(i) = 0.d0
            ctsotid(i) = 0.d0
          enddo       
          if( firstcall ) print *,'SBFN: ocean tides off '
        endif       
c**     replace the above code with the following, zeroes out coefficients
c       higher than the degree requested (all if otidedeg=0)
        if(debug) print *,'nootide otidedeg nczone '
     .                    ,nootide,otidedeg,nczone
        if( otidedeg.lt.nczone ) then
          if( otidedeg.lt.2 ) then
            i1 = 1
          else
            i1 = otidedeg 
          endif       
          do i=i1,nczon1                
            czcotid(i) = 0.d0
            czsotid(i) = 0.d0
          enddo           
          if( otidedeg.lt.2 ) then
            i1 = 1
          else
            i1 = tessindx(otidedeg)
          endif       
          do i=i1,nctes2
            ctsotid(i) = 0.d0
            ctcotid(i) = 0.d0
          enddo 
        endif

        if( debug ) then  
           call print_csnm('OTIDE','normalized',4,czcotid, 
     .                                   ctcotid, ctsotid)
        end if

c       Convert to unnormalized coefficients
        call sclcof1(-1,0,nczone,nctess,czcotid)    
        call sclcof1(-1,0,nczone,nctess,czsotid)
        call sclcof1(-1,1,nczone,nctess,ctcotid)
        call sclcof1(-1,1,nczone,nctess,ctsotid)
        if( debug ) then  
           print *,'Sample unormalized ocean tide coefficients '
           do i = 2, 5   ! Don't do all 
              do i1 = 0,i
                 if( i1.eq.0 ) then  ! Zonal
                     write(*,'("UNNOR C/S ",I2,1x,I2,1x,2e15.6)') 
     .                         i,i1, czcotid(i-1)
                 else    ! Non-Zonal
                     write(*,'("UNNOR C/S ",I2,1x,I2,1x,2e15.6)') 
     .                         i,i1, ctcotid((i-1)*i/2+i1-1),
     .                               ctsotid((i-1)*i/2+i1-1)
                 endif
              end do
           enddo

        endif
 
c       Determine non-gravitational accelerations
         if( .not.norad ) then
           imode = -1       
           if(debug) print *,'calling ertorb imode fjd ',imode,fjd
           call ertorb(imode,idum1,idum2,fjd,dummy)
         else
           if(firstcall) print *,'SBFN: Radiation pressure off'
         endif                 

       elseif ( k.eq.7 ) then
c          k = 7  => first partial component

c       Determine partial derivatives of motion relative to central body
c          kount  = n/6-1 = number of partial derivatives
        do ll1=1,kount
          l=ll1*6
          do i=1,3
            l=l+1
            dsbcor(i,ll1)=y(l,j)  
cd            print *,'SBFN ll1 l i j y dsbcor '
cd     .        ,ll1,l,i,j,y(l,j),dsbcor(i,ll1)
          enddo
        enddo

c       Determine central body harmonic quantities for partials
c          this should be in sbset    (next 2 lines)
cc?      nczp1=min0(nczone,nczonp+1)
cc?      nctp1=min0(nctess,nctesp+1)
        nczp1=2
        nctp1=1
        call legnd2(cslat,cclat,nczp1 ,nctp1 ,cleg,cleg1,cleg2,
     1            cgleg,cgleg1,cgleg2)


      endif
c     end quantities for accelerations calculated once per step

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c         Evaluate right-hand-side for each equation component
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c     Determine the class of equation
c          k   = eqn number ( 1 -> neq )
c          kk  = component number for each parameter
c                  ( -2 -> 0 for position ; 1 -> 3 for velocity )
c          kkk = eqn group  ( 0 for motion, 1 -> 9 for partials )

      kkk   = (k-1)/6
      kk    =  k-kkk*6-3
      if(debug) print *,'SNFN bottom j k kkk kk y(k+3,j) '
     .   ,j,k,kkk,kk,y(k+3,j)
      if ( kk.le.0 ) then
        fn(1) =y(k+3,j)
      else   
      if( debug ) print *,'SBFN bottom s j k kk kkk  sbcor '
     .     ,s,j,k,kk,kkk,sbcor(kk)     
c       flag for printing tide turn-off messages on the first call 
        firstcall = .false.
        call sbfn1(k,kk,kkk,s,fn)  
        if(debug) then
           print *,'SBFN bottom s j k fn ',s,j,k,fn(1)
           if(k.gt.3 ) stop 1
        endif 
      endif                                               
      sbfn = fn(1)
      return
      end                                                                     






