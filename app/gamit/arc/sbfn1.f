Copyright (c) Massachusetts Institute of Technology and the University of
California at San Diego, 1994.   All rights reserved.

      Subroutine SBFN1 ( k,kk,kkk,s,fn )  

c     Evaluate right side of motion and partials
c     ash/amuchastegui - June 1969
c     Revised 1977 April - J.F.Chandler
c     Rick Abbot November 1984, modification for GPS satellites
c     King / Ho Sept 93 - clean up sign convensions
c     King May 94 - Add option for separate GM for harmonics
c     King Sep 94 - Add solid-earth tides using Rotchacher vector formulation
c     King Apr 95 - Add solid-earth tides using harmonic formulation from
c                   IERS standards as coded by Jim Ray.

c     Input arguments:
c     ----------------
c           k   int*4  equation number (1 --> neq)
c           kk  int*4  component number (-2 -> 0 position; 1 -> 3 velocity
c           kkk int*4  eqn group (0 for motion, 1 -> 9 for partials)     
c     Output argument:
c     ---------------  
c           fn  r*8  acceleration  (km/s**2)
c           s   r*8  time in seconds past epoch of ICs   (currently used only for debug)


      implicit none
                    
      include '../includes/dimpar.h'     
      include '../includes/global.h'
      include '../includes/arc.h'

      integer*4 nczonp,ih,jh,iq,jq,kk,kkk,i1,i,j,k

      real*8 sumtc,fn,czone,szone,termg,sumhc,tg,zd,fnfb,tessc1
     .     , tessc2,zgaa,di2,di3,di4,zgax,sump,tm1,sums,tm2,sbforc,zgxx
     .     , csnh1,csnh2,djh,dot,dadxc,dadxh,sum,dadxp,dadxs,termd
     .     , rotdot_e2i,angvelten,angvel,Me,Ie,grelacc,s,sumv,sumj 
     .     , fjd 

      dimension fn(2)
* Set for 12x12 field
      dimension czone(12),szone(12),tessc1(77),tessc2(77)
     .        , csnh1(77),csnh2(77)
      dimension dadxc(3,3),dadxp(3,3),dadxs(3,3),dadxh(3,3)
      dimension rotdot_e2i(3,3),angvelten(3,3),angvel(3),grelacc(3)
      equivalence (dadxc,dadx)
                                                       
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

c for i=1,2,3 we have in this routine
c sbcor(i)  coordinates of Earth satellite relative to Earth (now contains velocities)
c bcor(i)   coordinates of Earth satellite relative to Sun
c ccor(i)   coordinates of Earth relative to Sun (should be i=1,6)
c pccor(i)  coordinates of Moon relative to Earth
c pbcor(i)  coordinates of Moon relative to Earth satellite
c vccor(i)  coordinates of Earth relative to Venus
c vbcor(i)  coordinates of Earth satellite relative to Venus      
c jccor(i)  coordinates of Earth relative to Jupiter
c jbcor(i)  coordinates of Earth satellite relative to Jupiter

c  G*mass for Venus and Jupiter (km**3/s**2)
      real*8 gmv/3.24858598826000000d5/,gmj/1.26712767857796000d8/
                 
c  For debug
      logical debug/.false./,debug2/.false./,norad/.false./
     .    ,accelprt/.false./
      real*8 tempdb 
c  Temporary for testing new and old code togethr
      logical fcheck 

                                
c Class of Equation
c          k   = eqn number ( 1 -> neq )
c          kk  = component number for each parameter
c                  ( -2 -> 0 for position ; 1 -> 3 for velocity )
c          kkk = eqn group  ( 0 for motion, 1 -> 9 for partials
    
 
      if (kkk.eq.0) then

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c           Determine right side of equations of motion
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c          Effect of the Moon on the equations of motion
c          sump = (moon wrt sat) - (moon wrt earth)
        sump = pbcor(kk)/rpb3-pccor3(kk)
        if( debug ) then
            print *,'kk start sbcor fn ',kk,sbcor(kk),fn(1)
            print *,'Moon gm pbcor/rpb3 pccor3 sump accel'
     .           ,gm(2),pbcor(kk)/rpb3,pccor3(kk),sump,gm(2)*sump
        endif    
        fn(1) = fn(1)+gm(2)*sump
        if(debug) print *,'Moon fn ',gm(2)*sump,fn(1)

c          Effect of the Sun on the equations of motion
c          sums = (earth wrt sun)-(satellite wrt sun)
        sums = ccor3(kk)-bcor(kk)/rb3     
        if( debug ) print *,'Sun gm ccor3 bcor/rb3 sums accel '
     .           ,gm(3),ccor3(kk),bcor(kk)/rb3,sums,gm(3)*sums
        fn(1) =fn(1) + gm(3)*sums
        if(debug) print *,'Sun fn ',gm(3)*sums,fn(1) 
                  
c          Effect of Venus and Jupiter on the equations of motion
        if( fcheck('nbody').and.lbody.gt.0 ) then
c         sumv = (Earth wrt Venus)-(satellite wrt Venus) 
          sumv = vccor(kk)/rvc**3 - vbcor(kk)/rvb**3    
c         reverse since acceleration should be negative
          if(debug) then          
            print *,'Venus vccor vbcor ',vccor,vbcor
            print *,'Venus gm vccor/rvc3 vbcor/rvb**3 sumv accel '
     .            ,gmv,vccor(kk)/rvc**3,vbcor(kk)/rvb**3,sumv,gmv*sumv
          endif
          fn(1) =fn(1)+gmv*sumv                           
          if(debug) print *,'Venus fn ',gmv*sumv,fn(1)
c         sumj = (Earth wrt Jupiter) - (satellite wrt Jupiter)   
c         reverse since accelration should be negative
         sumj = jccor(kk)/rjc**3 - jbcor(kk)/rjb**3  
          if(debug) then
            print *,'Jupiter jccor jbcor ',jccor,jbcor 
            print *,'Jupiter gmj jccor/rjc3,jbcor/rjb**3 sumj accel'
     .           ,gmj,jccor(kk)/rjc**3,jbcor(kk)/rjb**3,sumj,gmj*sumj 
          endif
          fn(1) =fn(1) + gmj*sumj         
          if(debug) print *,'Jupiter fn ',gmj*sumj,fn(1)
        endif
                                                          
c          Effect of Earth zonal harmonics and tides on equations of motion

cd         if( rdebug ) then 
cd           write(*,*) 'cleg    ',(cleg(j),j=1,nczon1)
cd           write(*,*) 'cleg1   ',(cleg1(j),j=1,nczon1)
cd           write(*,*) 'czhar   ',(czhar(j),j=1,nczon1)
cd           write(*,*) 'cztid   ',(cztid(j),j=1,nczon1)
cd           write(*,*) 'czcotid ',(czcotid(j),j=1,nczon1)
cd         endif
cd        if( debug ) then
cd          print *, ' SBFN1 '   
cd          print *,'C20 EGM08 alone ',czhar(1)
cd          print *,'C21 EGM08 alone ',cchar(1)
cd          print *,'S21 EGM08 alone ',cshar(1)
cd          print *,'C22 EGM08 alone ',cchar(2)
cd          print *,'S22 EGM08 alone ',cshar(2)
cd          print *,'C20 solid tide ',cztid(1)
cd          print *,'C21 solid tide ',cctid(1)
cd          print *,'S21 solid tide ',cstid(1)
cd          print *,'C22 solid tide ',cctid(2)
cd          print *,'S22 solid tide ',cstid(2)        
cd          print *,'C20 ocean tide ',czcotid(1)
cd          print *,'S20 ocean tide ',czsotid(1)
cd          print *,'C21 ocean tide ',ctcotid(1)
cd          print *,'S21 ocean tide ',ctsotid(1)
cd          print *,'C22 ocean tide ',ctcotid(2)
cd          print *,'S22 ocean tide ',ctsotid(2)
cd        endif
        sumhc = 0.0d0                                              
c  RWK 170420:  Need to incorporate the out-of-phase zonal terms for the
c               ocean tides (czsotid). 
        do i=1,nczon1
          di2 = dble(i+2)
c         zonal term is negative if coefficients are C's rather than J's
c         see Ash p. 71, Eq. 103.
          czone(i) = - ( di2 * sbcor(kk)*cleg(i)/rsb-cleg1(i)*cslat1(kk)
     .              ) / rsbh(i) / rsb2                        
c** trial S-zonal terms - I'm guessing on the form of this term 
          szone(i) =  ( di2 * sbcor(kk)*cleg(i)/rsb-cleg1(i)*cclat1(kk)
     .              ) / rsbh(i) / rsb2 
          sumhc = sumhc + (czhar(i) + cztid(i) + czcotid(i) )*czone(i)           
c** add S-zonals for ocean tides
          sumhc = sumhc + czsotid(i)*szone(i)    

         if(debug) print *,'di2 sbcor cleg rsb cleg1 cslat1 rsbh rsb2'
     .     ,di2,sbcor(kk),cleg(i),rsb,cleg1(i),cslat1(kk),rsbh(i),rsb2
          if(debug) 
     .       print *,'i czhar cztid czcotid czsotid,czone rrszone sumhc'
     .                     ,i,czhar(i),cztid(i),czcotid(i),czsotid(i)
     .                     ,czone(i),szone(i),sumhc
        enddo      
        if( debug ) then
          print *,'Coords: sbcor rsb rsb2 rsbh ',sbcor,rsb,rsb2,rsbh
           print *,'Zonals sumhc fn ',sumhc,fn(1)
          print *,'  cleg cleg1 cslat,cslat1  '
     .           ,   cleg,cleg1,cslat,cslat1
          print *,'  nczon1,czhar ',nczon1,czhar
          print *,'  czone ',czone
          print *,'  nctes1 ',nctes1
        endif

c          Effect of Earth tesseral harmonics and tides on equations of motion
         
cd       if( debug ) then 
cd          write(*,*) 'cgleg   ',(cgleg(j) ,j=1,nctess)
cd          write(*,*) 'cgleg1  ',(cgleg1(j),j=1,nctess)
cd          write(*,*) 'cchar   ', (cchar(j),j=1,nctess)
cd          write(*,*) 'cctid   ', (cctid(j),j=1,nctess)
cd          write(*,*) 'cstid   ', (cstid(j),j=1,nctess)
cd          write(*,*) 'ctcotid  ', (ctcotid(j),j=1,nctess)
cd          write(*,*) 'ctsotid  ', (ctsotid(j),j=1,nctess)
cd        endif

        ih=0
        do i=1,nctes1
           di2=dble(i+2)
           sumtc = 0.0d0
           i1=i+1
           do jh=1,i1
              djh = dble(jh)
              ih=ih+1   
              tessc1(ih) =
     .            (-di2*sbcor(kk)*cgleg(ih)/rsb + cgleg1(ih)*cslat1(kk))
     .            / rsbh(i)/ rsb2
              tessc2(ih) =
     .             cgleg(ih)*clng1(kk)/rsbh(i) / rsb2
              if(kk.eq.1) then
c                only once per iteration
                 csnh1(ih) =   ( cchar(ih) 
     .                 + cctid(ih) + ctcotid(ih) ) * cclng(jh)
     .                       + ( cshar(ih) 
     .                 + cstid(ih) + ctsotid(ih) ) * cslng(jh)
                 csnh2(ih) = - ( cchar(ih) 
     .                 + cctid(ih) + ctcotid(ih) ) * cslng(jh)
     .                       + ( cshar(ih) 
     .                 + cstid(ih) + ctsotid(ih) ) * cclng(jh)
              endif
              sumtc = sumtc + csnh1(ih) * tessc1(ih)
     .                      + djh * csnh2(ih) * tessc2(ih)
              if(debug) print * 
     .           ,'jh ih  cchar cshar cctid cstid ctcotid ctsotid sumtc'
     .           ,jh,ih,cchar(ih),cshar(ih),cctid(ih),cstid(ih)
     .           , ctcotid(ih),ctsotid(ih),sumtc
cd             if(debug) then
cd                 print *,'jh ih cchar cshar cctid cstid ctcotid ctsotid'
cd     .                  ,jh,ih,cchar(ih),cshar(ih),cctid(ih),cstid(ih)
cd     .                  ,ctcotid(ih),ctsotid(ih)
cd                 print *,'  cclng cslng tessc1 tessc2 csnh1 csnh2 sumtc'
cd     .              ,cclng(jh),cslng(jh),tessc1(ih),tessc2(ih),csnh1(ih)
cd     .              ,csnh2(ih),sumtc
cd              endif  
           enddo
c          end inner loop (order)        
          sumhc = sumhc + sumtc
        enddo
c       end outer loop (degree)

cd        if( debug ) then
cd          print *,'  cgleg cgleg1 ',cgleg,cgleg1
cd          print *,'  i1 nctes1 cchar cshar ',i1,nctes1,cchar,cshar
cd          print *,'  cclng cslng ',cclng,cslng
cd          print *,'  clng1 ',clng1
cd          print *,'  tessc1 tessc1 ',tessc1,tessc2
cd          print *,'  csnh1 csnh2 ',csnh1,csnh2
cd        endif

c       gmhfct accounts for different gm for harmonics and 2-body accelerations
        if( debug ) then 
          print *,'fn before harmonics ',fn(1)
          print *,'Harmonics sumhc gm gmhfct fn '
     .        ,sumhc,gm(1),gmhfct,gm(1)*sumhc*gmhfct
        endif
        fn(1)= fn(1) + gm(1)*sumhc*gmhfct
        if( debug ) print *,'Gravity fn ',gm(1)*sumhc*gmhfct
        if( debug) print *,'fn ',fn(1)

c          Additional accelerations for earth satellite
        if( .not.norad ) then
          if(debug) print *,' SBFN1 calling ERTORB k ',k
          if( debug ) tempdb = fn(1)        
          call ertorb(k,kk,kkk,fjd,fn(1)) 
          if( debug2 ) print *,'SBFN1 ERTORB s acc ',s,fn(1) - tempdb 
        endif           
        if(debug) print *,'After ertorb fn ',fn(1)                
        

c          Effect of general reativity due to motion of the Earth Satellite
                      
c         Temporarily set EGM08 = EGR08 if general relativity to be turned on
           if( gravmod.eq.'EGR08' ) then

c          Convert cntrtd rotation matrix which is i2e to e2i using the transpose operator.
        call transp(cntrtd, rotdot_e2i, 3, 3)
cd        call printmat(cntrot,3,3,'cntrot ')
cd        call printmat(cntrtd,3,3,'cntrtd ')
cd        call printmat(rotdot_e2i,3,3,'rotdot_e2i ')

c          Compute the angular velocity tensor of the terrestrial frame
        call matmpy(rotdot_e2i, cntrot, angvelten, 3, 3, 3)
cd        call printmat(angvelten,3,3,'angvelten ')

c          Angvelten is the angular velocity tensor, with the rotation angles about x (3,2), y(1,3) and z (2,1).
c          It also contains the -ives of them, and has zeroes down the diagonal.

        angvel(1) = angvelten(3,2)
        angvel(2) = angvelten(1,3)
        angvel(3) = angvelten(2,1)
cd        print *,'angvel(1:3) ',angvel

c          Some constants needed that should be moved to arc.h?
c          Mass of the Earth, Me (kg)
        Me = 5.9736e24

c          Moment of inertia of Earth around polar axis, (kgkm^2)
        Ie = 8.0365d31

        call generalrel(czhar(1), sbcor(1:3), sbcor(4:6), ccor, angvel
     .  , gm, ertrad, Me, Ie, ltvel, grelacc)
cd         print *,'kk, ccor sbcor fn(1) grelacc '
cd     .       ,kk, ccor, sbcor, fn(1), grelacc

c          Add general relativity to other perturbing accelerations
        fn(1) = fn(1) + grelacc(kk)        
        if(debug) print *,'general rel acc kk ',grelacc(kk),kk  
        if(debug) print *,'genrel fn ',grelacc(kk),fn(1)
           
c end of test to turn on general relativity
        endif

c          Effect of central force on motion of the Earth satellite
        sbforc= sbcor(kk)/rsb3
        sum=-gm(1)*sbforc   
        if(debug) print *,'kk,sbcor(kk),gm sbforc central_acc '
     .                    ,kk,sbcor(kk),gm(1),sbforc,sum

c          Add central force to other perturbing forces
        fn(1) =   fn(1) + sum
        if( accelprt ) then
          write(*,'(f7.0,i3,d22.14)') s,kk,fn(1)
        endif 
cd        if( debug.and.kk.eq.3  ) then 
cd            print *,'Final FN stop ',fn(1)
cd           stop
cd         endif

      elseif (kkk.gt.0 ) then

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c           Determine right side of equations for partial derivatives
c           with respect to parameter alpha(kkk)
c           (kkk goes from 1 to number of partials).

c           Need to evaluate most quantities only for the first
c           equation element (k=10)
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        if( k.eq.10 ) then

c          Effect of the Moon on partial derivatives
          do iq=1,3
            do jq=1,3
              dadxp(iq,jq)=0.0d0
            enddo
          enddo
          termd= gm(2)/rpb3
          tg  = -termd*3.d0/rpb2
          do iq=1,3
            do jq=1,iq
              termg=tg*pbcor(iq)*pbcor(jq)
               dadxp(iq,jq)=dadxp(iq,jq)+termg
               if(iq.ne.jq) dadxp(jq,iq)=dadxp(jq,iq)+termg
            enddo
            dadxp(iq,iq)=dadxp(iq,iq)+termd 
          enddo

c           Effect of Sun on partial derivatives
          do iq=1,3
            do jq=1,3
              dadxs(iq,jq)=0.0d0
            enddo
          enddo
          termd= gm(3)/rb3
          tg=-gm(3)*3.d0/rb3
          do iq=1,3
            do jq=1,iq
              termg=tg*bcor(iq)*bcor(jq)/rb2
              dadxs(iq,jq)=termg
              if(iq.ne.jq) dadxs(jq,iq)=termg
            enddo
            dadxs(iq,iq)=dadxs(iq,iq)+termd  
          enddo

c          Effect of Earth zonal harmonics on partial derivatives
          do iq=1,3
            do jq=1,3
              dadxh(iq,jq)=0.0d0
            enddo
          enddo
c       gmhfct accounts for different gm for harmonics and 2-body accelerations
          tg= -gm(1)*gmhfct/rsb3
          zgaa=0.0d0
          zgax=0.0d0
          zgxx=0.0d0
          zd  =0.0d0
          di3=3.d0
          di4=4.d0
          do 3040 i=1,nczon1
            nczonp=1
            if( i.gt.nczonp) go to 3040
            di2=di3
            di3=di4
            di4=dble(i+4)
            tm2= czhar(i)/rsbh(i)
            zgaa=zgaa+cleg2(i)*tm2
            tm1=tm2*(di2*cleg(i)+cslat*cleg1(i))
            tm2=tm2*(di3*cleg1(i)+cslat*cleg2(i))
            zgax=zgax+tm2
            zgxx=zgxx+di4*tm1+cslat*tm2
            zd=zd+tm1
 3040     continue

          zgaa=zgaa*tg
          zgax=zgax*tg/rsb
          zgxx=zgxx*tg/rsb2
          zd=-zd*tg

          do iq=1,3
            do jq=1,3
              termg = zgxx*sbcor(iq)*sbcor(jq)
     .           - zgax*(sbcor(iq)*cntrot(3,jq)+sbcor(jq)*cntrot(3,iq) )
     .           + zgaa*cntrot(3,iq)*cntrot(3,jq)
              dadxh(iq,jq) = termg
              if(iq.eq.jq) go to 3060
              dadxh(jq,iq) = termg
            enddo
 3060       dadxh(iq,iq) =dadxh(iq,iq) + zd
          enddo

c           Effect of central force on partial derivatives
          termd= -gm(1)/rsb3
          tg= -termd*3.d0/rsb2
          do iq=1,3
            do jq=1,iq
              termg=tg*sbcor(iq)*sbcor(jq)
              dadxc(iq,jq)=termg
              if(iq.ne.jq) dadxc(jq,iq)=termg
            enddo
            dadxc(iq,iq)=dadxc(iq,iq)+termd
          enddo

c          Form dadx each iteration
          do iq=1,3
            do jq=1,3
              dadx(iq,jq) = dadxc(iq,jq) - dadxh(iq,jq)
     .                    - dadxp(iq,jq) - dadxs(iq,jq)
            enddo
          enddo
c       note - dadx is symmetric except for dadxl (drag contribution)
                           
cd        print *,'j dadxc ',j,((dadxc(iq,jq),iq=1,3),jq=1,3)
cd        print *,'  dadx  ',((dadx(iq,jq),iq=1,3),jq=1,3)


        endif
c       end branch for quantities evaluated only once (k=10)

c           Effect of additional accelerations for earth satellite on partials
        if(debug) print *,'SBFN1 calling ERTORB kkk k ',kkk,k
        if (kkk.gt.6) call ertorb(k,kk,kkk,fjd,fn(1))
c           
cd        print *,'j dadxc ',j,((dadxc(iq,jq),iq=1,3),jq=1,3)
cd        print *,'  dadx  ',((dadx(iq,jq),iq=1,3),jq=1,3)
cd        print *,'kk kkk dsbcor ',kk,kkk,((dsbcor(iq,jq),iq=1,3),jq=1,15)
        fnfb= dot( dadx(1,kk),dsbcor(1,kkk) )  
cd        print *,'SBFN1 bottom k kk kkk dsbcor dadx fn fnfb '
cd     .        ,  k,kk,kkk,dsbcor(1,kkk),dadx(1,kk),fn(1),fnfb
        fn(1)=fn(1)+fnfb

      endif
c     end of branch on kkk ( =0 motion ; >0 partials )
            
cd      if( debug ) stop 
      return
      end
