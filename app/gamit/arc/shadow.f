Copyright (c) Massachusetts Institute of Technology and the University of
California at San Diego, 1994.  All rights reserved.

      Subroutine SHADOW ( lambda,ecltyp )
c
c    Computation of fraction (lambda) of solar disk seen by spacecraft
c           Beebe, King, Reasonberg, Preston   June 19971
c
      implicit none  
   
      include '../includes/dimpar.h'   
      include '../includes/arc.h'
  
      integer i
      real*8 lambda  
      character*1 ecltyp 

c     function
      logical kbit 
  
c        quantities internal to the routine 
      
      real*8 rsbx,rs,rp,sep,sepp(3),ubcor(3)
     .     , pscor(3),rps,sbpcor(3),dot,area1,area2,area3,ari,hgt
     .     , phi,pi,r1,r2,thet

c          shadow computation - geometric model
c          lambda = 1 - no shadow
c          lambda = 0 - no sunlight
c          0 .lt. lambda .lt. 1 - partial shadow

c        no consideration is given to the change of direction associated
c        with partial shadow.

* Initialization
      lambda=1.0d0
      ecltyp = 'N' 

c ** Check for both eclipses of the Sun by both the Earth and the Moon

c    First the Earth
       

cc      write (6,55555) rb,rc
      if(rb.le.rc) go to 200
c
c        get the unit vector of the satellite relative to the sun
      do  i=1,3
       ubcor(i)=bcor(i)/rb
      enddo     
               

cRWK 040512: old code to adjust the positions for the light-time delay
c            has been removed; effect is < 1.d-5 m for the GPS case (?)

         call cross(sbcor,ubcor,sepp) 
c        rsbx is the projection of sbcor along bcor
         rsbx = dot(sbcor,ubcor)

c        rs, rp are apparent (from satellite) radii of sun and earth
c        sep is apparent separation of their centers
      rs=sunrad/rb
      rp=ertrad/rsbx
      sep=dsqrt(sepp(1)**2+sepp(2)**2+sepp(3)**2)/rsbx        
      call get_lambda( rs,rp,sep,lambda,idbprt,idebug )  
cd       print *,'rs rp sep lambda ',rs,rp,sep,lambda  
      if( kbit(idbprt,3) ) then
        if ( kbit(idbprt,1).or.(lambda.lt.1.d0.and.kbit(idbprt,2))) then
         write(idebug,*) 'SHADOW sbcor rsb bcor rb ccor rc '
     .                    , sbcor,rsb,bcor,rb,ccor,rc 
         write(idebug,*) 'SHADOW lambda rs rp sep ',lambda,rs,rp,sep
        endif
      endif
                            
c     If no Earth eclipse, check the Moon  
             

  200 if( lambda.lt.1.d0 ) then
         ecltyp = 'E'
         return   
      
       else    

c*** Temporary:  allow override of lunar shadowing to make model consistent with old h-files
c         if( ecl_flag.eq.'nolunecl' ) return
         if( kbit(idbprt,12)) return
              
c                                Vector    Distance
c        Moon wrt Earth          pccor      rpc
c        Earth wrt Sun           ccor       rc
c        Moon wrt Sun            pscor      rps   
c        Satellite wrt Earth     sbcor      rsb  
c        Satellite wrt Sun       bcor       rb 
c        Satellite wrt Moon      sbpcor     rsbp

      do i=1,3
       pscor(i) = pccor(i) + ccor(i) 
       sbpcor(i) = sbcor(i) - pccor(i)
      enddo  
      rps = dsqrt(pscor(1)**2+pscor(2)**2+pscor(3)**2)   
c      print *,'M  wrt E pccor rpc ',pccor,rpc
c      print *,'E  wrt S ccor  rc  ',ccor,rc
c      print *,'M  wrt S pscor prs ',pscor,rps
c      print *,'SV wrt E sbcor rsb ',sbcor,rsb  
c      print *,'SV wrt S bcor  rb  ',bcor,rb   
c      print *,'SV wrt M sbpcor rsbp ',sbpcor,rsbp

      if(rb.le.rps) return

c       unit vector of SV wrt Sun
      do i=1,3
        ubcor(i)=bcor(i)/rb
      enddo

c        rsbx is the projection of sbcor along bcor
      rsbx=dot(sbpcor,ubcor)     
c      print *,'rb sunrad rsbx monrad ',rb,sunrad,rsbx,monrad

c        rs, rp are apparent (from satellite) radii of sun and moon
c        sep is apparent angular separation of their centers
      rs=sunrad/rb
      rp=monrad/rsbx    
      call cross(sbpcor,ubcor,sepp)  
      sep=dsqrt(sepp(1)**2+sepp(2)**2+sepp(3)**2)/rsbx  
c      print *,'SHADOW moon rs rp sep ',rs,rp,sep  
 
      call get_lambda( rs,rp,sep,lambda, idbprt, idebug ) 
      if( lambda.lt.1.d0 ) ecltyp = 'M'

      endif

      return
      end


c-----------------------------------------------------------------------------------

      Subroutine get_lambda( rs, rp, sep, lambda, idbprt, idebug )

c     Calculate lambda

c      rs  : apparent radius of sun as viewed from satellite (radians)
c      rp  : apparent radius of eclipsing body as viewed from satellite (radians)
c      sep : apparent separation of the center of the Sun and eclipsing body (radians)

c      lambda : fraction of Sun's disk visible (1.0 = no eclipse; 0, = total eclipse)

      implicit none  
                         
      integer*4 idbprt,idebug

      real*8 rs,rp,sep,lambda
     .     , r1,r2,phi,hgt,thet,area1,area2,area3,ari,pi

c     function
      logical kbit
       
      pi =  3.14159265359d0

      if (rs+rp.le.sep) then
c       no eclipse 
        return
      elseif( rp-rs.ge.sep ) then 
c       full eclipse
        lambda = 0.d0
        return
      else
c       partial eclipse, do the calculations

        if(sep.le.rs-rp)go to 361
c        set r1 = smaller disc, r2 = larger
        if(rs.gt.rp)go to 371
        r1=rs
        r2=rp
        go to 370
  371   continue
        r1=rp
        r2=rs
  370   continue
c         phi = 1/2 angle subtended in disc 1 by arc of intersection
        phi = dacos((r1*r1+sep*sep-r2*r2)/(2.0d0*r1*sep))
        if (phi.lt.0.d0) phi = pi + phi
        if(r2/r1.gt.5.0d0)go to 365
c        thet = 1/2 angle subtended in disc 2 by arc of intersection
c        hgt  = 1/2 linear distance between ends of arc of intersection
        hgt=r1*dsin(phi)
        thet=dasin(hgt/r2)
        area2=sep*hgt
        area3=thet*r2**2   
        go to 366
  365   continue
c         one disc much bigger - treat boundary as a straight line
        hgt=dsqrt(r1**2-(sep-r2)**2)
        area2=hgt*(sep-r2)
        area3=0.0d0
  366 continue
        area1=(pi-phi)*r1**2
c         ari = area of non-overlapped portion of small disc
        ari=area1+area2-area3  
        if ( kbit(idbprt,1).or.(lambda.lt.1.d0.and.kbit(idbprt,2)))  
     .   write(idebug,*) 
     .      'SHADOW r1 r2 phi thet hgt area1 area2 area3 ari '
     .            , r1,r2,phi,thet,hgt,area1,area2,area2,ari
        area1=pi*rs**2
        if(rs.gt.rp)go to 362
c         sun is small disc
        lambda=ari/area1
        go to 999
  362   continue
c         eclipsing body is small disc
        area2=pi*rp**2
        lambda=(area1+ari-area2)/area1
        go to 999
  361  continue
c         eclipsing body lies within sun's disc - what fraction of sun's disk is blocked
        lambda=(rs**2-rp**2)/rs**2   

      endif
                 
  999 continue

      return
      end

