Copyright (c) Massachusetts Institute of Technology, 1986. All rights reserved.
      Subroutine adam
c
c     Adams-Moulton integration subroutine
c     W.B. Smith   Oct 1968
c     Rick Abbot  Nov  1984, modification to do GPS satellites only

      implicit none

      include '../includes/dimpar.h'    
      include '../includes/arc.h'


      integer*2 bakfor
      dimension dydtnu(maxyt2),yyy(maxyt2)
      dimension ytemp(2)

      integer*4 kstop,mm,np,idim,iter,nstmid,nsigns,l2
     .       , nst,j,k,l

      real*8 yyy,xsign,ytemp,ttemp1,ttemp2,sbfn,dydtnu,tstopt,s

c     common/adams/
c            npredt  i  = number of predictor terms
c            ncoret  i  = number of corrector terms
c            neq     i  = number of equations
c            nh      i  = positive, step size is nh seconds
c            nh      i  <= zero, step size is 2**nh seconds
c            inthmx  i  = positive, tabular interval is inthmx seconds
c            inthmx  i  <= zero, tabular interval is 2**inthmx seconds  
c            neq     i  = number of equations to integrate
c            meq     i  = number of equations for accuracy test (=6)

      logical debug/.false./,debug2/.false./,debug3/.false./

      open(unit=3,status='scratch',form='unformatted')

      np=npredt+1
      bakfor=0
      ttemp1=trun0
      ttemp2=trunf
      tint0=0.d0    
      tstop=trunf
      if (trun0.lt.0.d0.and.trunf.lt.0.d0) go to 7
      go to 8
    7 tstop=trun0
      trun0=trunf
      trunf=tstop
    8 if (trunf.eq.0.d0) nsign=-1
      if (trunf.eq.0.d0) go to 9
      if (trun0/trunf.lt.0.d0) bakfor=1
      if (bakfor.eq.1) nsign=-1
c
      if (nsign.lt.0.and.bakfor.eq.1) go to 9
      go to 10
    9 tstop=trun0
      trun0=0.d0
   10 continue 
      if(debug) print *,'ADAM 10 nsign hc trun0 trunf tstop '
     .       ,nsign,hc,trun0,trunf,tstop
   99 continue
      if(debug2) then
         print *,'ADAM 99 s kstop nsign np hc hmx tstop'
     .     , s,kstop,nsign,np,hc,hmx,tstop
      endif
      if(tstop-tint0)97,999,98
   97 kstop=1
      go to 96
   98 kstop=2
c          
      if(debug)  print *,'ADAM start trun0 trunf bakfor nsign kstop '
     .                  , trun0,trunf,bakfor,nsign,kstop


c        set up control constants for starting procedure using start_int
c        to calculate np derivatives prior to tint0

   96 k=1
      tstopt=tstop
      xsign=nsign
      hmx=dnh*nsign
      tstop=tint0-np*hmx
      hmx=-hmx
      nsigns=nsign
      nsign=-nsign 
      if(debug2) print *
     . ,'ADAM calling start_int tstop tstopt nsign hmx l2 s '
     .                , tstop,tstopt,nsign,hmx,l2,s 
      if(debug) print *,'  y(14,1-4) ',(y(14,mm),mm=1,4)

c        calculation of the np derivatives

  182 continue 
      if(debug)  print *,'ADAM calling START_INT s l2 ',s,l2
      call start_int( s,l2 )   
      if(debug) print *,'after start_int tstop tstopt nsign hmx hc l2 s'
     .        ,tstop,tstopt,nsign,hmx,hc,l2,s   
      if(debug) print *,'  y(14,1-4) ',(y(14,mm),mm=1,4)

      if(l2.gt.0)goto 999
      k=k+1
      idim = 3
      do 183 l=1,neq
      dydt(l,k)=sbfn(l,idim,s)
      if(debug) print *,'ADAM s k l dydt ',s,k,l,dydt(l,k)
      if(debug) print *,'  y(14,1-4) ',(y(14,mm),mm=1,4)

  183 continue     
      if(debug) print *,'ADAM 183 s np k dydt(1,k) ',s,np,k,dydt(1,k)
      if(k-np)182,184,184

c        set up of control constants for continuing proceedure using
c        adams-moulton method

  184 do 185 l=1,neq
  185 y(l,4)=satics(l)
      hc=-hmx
      nsign=nsigns
cc      if(inthmx)190,190,191
cc  190 hmx=2.d0**inthmx*xsign*86400.d0
cc  190 hmx = stdint/(-nh)*xsign
cc      go to 192
cc  191 hmx=inthmx*nsign
      hmx=dinthmx*nsign
cc  192 tstop=tstopt   
      tstop=tstopt
      nstmid=hmx/hc+1.d-15
      s=tint0
      if(debug) print *,'ADAM bef 124 bakfor s trun0 tint0 '
     .     ,bakfor,s,tint0,trun0
      if (bakfor.eq.0.and.trun0.eq.tint0) call sbout
  124 continue          
cd      if(s.lt.0.d0 ) debug2 = .true. 
      if(debug) print *,'ADAM aft 124 tstop nsign hmx hc s '
     .                 ,tstop,nsign,hmx,hc,s 
c*    debug stop when two tabular points into the backward integation
      if(debug) then 
        if(s.lt.-1800.d0 ) then
           print *,'DEBUG 124 stop at s=-1800.'
           stop
        endif
      endif 
      do 900 nst=1,nstmid   
         if(debug) print *,'nstmid nst ',nstmid,nst 
c        prediction of y(n+1)
cd      if(s.eq.0.0d0) then
cd        debug3 = .true.
cd      else
cd        debug3 = .false. 
cd      endif
      do l=1,neq
        ytemp(1)=0.d0
        do k=1,np
          ytemp(1)=ytemp(1)+cadams(k)*dydt(l,k)  
          if(debug) print *,'ADAM prediction i k cadams dydt ytemp '
     .        ,l,k,cadams(k),dydt(l,k),ytemp(1)
          if(debug) print *,'  y(14,1-4) ',(y(14,mm),mm=1,4)
        enddo
        y(l,1)=y(l,4)+hc*ytemp(1)
        if(debug3) print *,'l y1 y4 hc ytemp1 '
     .     ,l,y(l,1),y(l,4),hc,ytemp(1)
      enddo
      j=1      
      if(debug) print *,'ADAM prediction s neq np hc y(14,1-4) '
     .                   ,s,neq,np,hc,(y(14,l),l=1,4)
      do 165 iter=1,2
c calculation of dy/dt at y(n+1)
        do l=1,neq
          dydtnu(l)= sbfn(l,j,s+hc) 
          if(debug) print *,'s iter j l hc dydtnu y(14,1-4) '
     .        ,s,iter,j,l,hc,dydtnu(l),(y(14,mm),mm=1,4)
        enddo
        j=j+1     
c        corrector equation being applied twice
        do 160 l=1,neq
          ytemp(1) =dadams(1)*dydtnu(l)
          go to (150,155),iter
 150      yyy(l)=0.d0
          do  k=1,ncoret
            yyy(l)=yyy(l)+dadams(k+1)*dydt(l,k)
          enddo
 155      y(l,j)=y(l,4)+hc*(ytemp(1) +yyy(l))     
          if(debug3) print *
     .     ,'l j hc temp1 ncoret dadams12 dydt12 yyy y(l,j)', l,j,hc
     .    ,ytemp(1),ncoret,dadams(12),dydt(l,12),yyy(l),(y(l,mm),mm=1,4)
 160    continue
       if(debug) print *,'correction s iter ncoret dydt(1.nc) y(14,1-4)'
     .       ,s,iter,ncoret,dydt(1,ncoret),(y(14,l),l=-1,4)
 165  continue

c        update things, preparatory to the next step.
      do  l=1,neq
        y(l,4)=y(l,3)
        do  k=2,np
          mm=np-k+2
          dydt(l,mm)=dydt(l,mm-1)
        enddo
        dydt(l,1)=dydtnu(l)
      enddo       
      if(debug) print *,'Update dydt(1,1) ',dydt(1,1)
      continue
 900  s = s + hc
c
c        call output routine     
      go to (1110,1109) kstop
 1109 if (s.ge.tint0) go to 1101
      go to 1012
 1110 if (s.le.tint0) go to 1101
      go to 1012
 1101 npstep=npstep+1          
      if(debug) print *,'ADAM 1101 s npstep nstout tint0 calling sbout'
     .     ,s,npstep,nstout,tint0
      if (npstep.eq.nstout.and.(dabs(s).ge.(dabs(trun0)-1.d-08))) then
        call sbout  
        if(debug2) then
          print *,'Tab pt s y ',s
          do j=1,24
            print *,j,y(j,4) 
          enddo 
        endif
      else
        if(debug2) print *,'       s y1 ',s,y(1,4)
      endif 
      if (npstep.eq.nstout) npstep=0
 1012 continue
c
c        test for end
      go to (1013,1014) kstop
 1013 if(s-tstop)1015,1015,124
 1014 if(s-tstop)124,1015,1015
 1015 if(l1.lt.0) go to 99
c
      if (nsign.gt.0) go to 999
c        backwards integration, reverse output
      call revrec
c
c        now check whether to continue with a forward integration
      if (nsign.lt.0.and.bakfor.eq.1) go to 899
      go to 999
  899 tint0=0.d0
      tstop=trunf
      bakfor=0
      nsign=+1
      npstep=0
      l1=0
      if(debug2) print *,'Finished back s bakfor nsign npstep l1 tstop '
     .       , s,bakfor,nsign,npstep,l1,tstop 
      go to 99
c
  999 trun0=ttemp1
      trunf=ttemp2
c
      close (unit=3)
c
      return
      end
