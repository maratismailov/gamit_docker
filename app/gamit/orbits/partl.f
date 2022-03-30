      Subroutine partl( jd,t,iepoch,ksat,rvec1 )

c     Calculate the partial derivatives of orbital coordinates with 
c     respect to global (translation/rotation/scale) and orbital parameters

      implicit none
                
c     Input:  jd,t           epoch of satellite position   
c             iepoch         epoch number in reference T-file
c             ksat           satellite index in array for reference orbit
c             rvec1(maxyt2)   vector of position and velocity plus partials for satellite  
c             Stored in common /nrmcom/:
c             nparam         number of parameters to be estimated
c             islot(mxoprm)  pointers to parameters to be estimated (see below)

c     Output: part(3,mxoprm,maxsat,maxepc)   partials of coordinates wrt parameters,
c                                     stored in common /omcpart/ 

c     Pointers to parameters (islot values)
c       1-3      translation
c       4-6      inertial rotation
c       7        scale
c       8-10     terrestrial rotation
c       11-13    terrestrial rotation rate   
c       n01-n15  satellite parameters (partials interpolated from ephemeris)
c                  where n is the satellite index on the reference T-file
c                  (n=1-nsat)
            
      include '../includes/dimpar.h'
      include '../includes/orbits.h'
      include 'orbfit.h'  
           
      integer*4 jd,iepoch,ksat
      real*8      rvec1(maxyt2)
                              
      logical eop,debug
      
      character*256 message

      integer*4 iut1pol,jds,ioerr,i,j,jj,kk,idir

      real*8 rvec_in(6),rvec_out(6),xrot(3,3),xdtrot(3,3)
     .     , rot(3,3),rotdot(3,3)
     .     , t,ts,diff,sidtm,eqe,radsec,casr
     .     , polepart(3,4),ut1part(3,2),tdtgpst,tdtutc,gpstutc,taiutc
     .     , prec(3,3),rnut(3,3),snpmat(3,3),snpdmt(3,3),srot(3,3)
     .     , sdrot(3,3),sidten(3,3),prot(3,3),xpole,ypole 
     .     , tran_part_t(3,7),tran_part_i(3,7)
     .     , pi,twopi
     .     , era, erarot(3,3)
     
c  Initialize UT1/pole flag and constants

      debug = .false.
* MOD TAH 200505: Read the sestbl. to get the correct value   
C     iut1pol = 3 
      call get_iut1pol( iut1pol )

      ioerr   = 0
      pi = 4.d0*datan(1.d0)
      twopi = pi * 2.d0
      radsec = twopi/86400.d0
      casr = pi/180.d0/3600.d0
         
c   Store the position vector

      do i=1,6
        rvec_in(i) = rvec1(i)
      enddo
 

c   Calculate the global partials
                                         
c     see if terrestrial rotation partials are needed
      eop = .false.
      do i=1,nparam 
         if( islot(i).ge.8.and.islot(i).le.13) eop=.true.
      enddo
c     if so, compute the rotation matrices
      if( eop ) then
        tdtgpst = 32.184d0 + 19.d0
c Old Inertial reference frame calculations 
        if ( precmod .eq. 'IAU76' .or. precmod .eq. 'IAU68' ) then                
          call pnrot( inut,jd,t,tdtgpst,eqe,prec,rnut,frame,precmod)
c       srotat (unique in GAMIT) expects UTC
          jds = jd
          ts = t
          gpstutc = taiutc(jd) - 19.d0
          call timinc(jds,ts,-gpstutc)
c       catch possible leap second
          if( jds.ne.jd ) then
            gpstutc = taiutc(jds) - 19.d0
            jds = jd
            ts = t
            call timinc(jds,ts,-gpstutc)
          endif
          tdtutc = taiutc(jd) + 32.184d0
c       print*, 'jds,ts ',jds,ts
          call srotat( jds,ts,tdtutc,eqe,iut1pol,srot,sdrot,sidtm
     .               , xpole,ypole,sidten,prot,iut1,ipole,precmod )
          call snp( prec,rnut,srot,snpmat )
          call snp( prec,rnut,sdrot,snpdmt )
          do j=1,3
            do i=1,3
              rot(i,j)    = snpmat(i,j)
              rotdot(i,j) = snpdmt(i,j)
            enddo
          enddo
          
        else 
        
          idir = -1
          call rotsnp_sofa( idir,jd,t,tdtgpst,iut1pol
     .                    , frame,precmod,iut1,ipole,inut
     .                    , rot,rotdot,sidtm,xpole,ypole
     .                    , prec,rnut,prot,sidten,era,erarot )   

          if ( precmod .ne. 'IAU76' .and. precmod .ne. 'IAU0A' ) then
c Is this what is required?
            sidtm = era
            sidten = erarot
cd            print*,'CIO based CRS: ',precmod,' : Copying ERA > SIDTM'             
          endif
            
        endif 
c Debug output
         if (debug) then
           write(*,100) 'PARTL - prec: ',precmod,
     .             ((prec(i,j),j=1,3),i=1,3) 
           write(*,100) 'PARTL - rnut: ',precmod,
     .             ((rnut(i,j),j=1,3),i=1,3) 
           write(*,100) 'PARTL - prot: ',precmod,
     .             ((prot(i,j),j=1,3),i=1,3) 
           write(*,100) 'PARTL - sidten: ',precmod,
     .             ((sidten(i,j),j=1,3),i=1,3)
           write(*,*) 'PARTL - sidtm,xp,yp: ',precmod,sidtm,xpole,ypole
           write(*,100) 'PARTL - rot: ',precmod,
     .             ((rot(i,j),j=1,3),i=1,3) 
           write(*,100) 'PARTL - rotdot: ',precmod,
     .             ((rotdot(i,j),j=1,3),i=1,3) 
100        format(a,1x,a,1x,/,3(1x,3D22.14,/))
         endif
                  
c       rotate the coordinates and velocity from inertial to earth fixed
        call matmpy( rot,rvec_in(1),rvec_out(1),3,3,1 )
        call matmpy( rotdot,rvec_in(4),rvec_out(4),3,3,1 )
        if (debug) then
          print*, 'PARTL rvec_in,rvec_out ',rvec_in,rvec_out
        endif
c       calculate the time of observation relative to the midpoint of the
c       data span, the reference time for EOP values
        diff = (jd + t/86400.d0) - (jde + te/96400.d0)   
        call eopart( xpole, ypole, sidtm, prec, rnut, prot, sidten
     .              , diff, rvec_out,  polepart, ut1part )
      endif
      

c     Calculate 7 parameter terrestrial transformation partials
      call tran_part( rvec_out,tran_part_t)

c     Calculate 7 parameter inertial transformation partials
      call tran_part( rvec_in,tran_part_i)

c  Assign the partials
                              
      do j=1,nparam
         if( islot(j).ge.1.and.islot(j).le.7 )  then  
c          inertial translation , rotation, or scale
             do i = 1,3
               part(i,j,ksat,iepoch) = tran_part_i(i,islot(j))
             enddo  
         elseif (islot(j).ge.8.and.islot(j).le.13) then  
c        terrestial rotation and rate
           do i = 1,3
             if( islot(j).eq.8 ) part(i,j,ksat,iepoch) = polepart(i,1)
             if( islot(j).eq.9 ) part(i,j,ksat,iepoch) = polepart(i,3)
             if( islot(j).eq.10) part(i,j,ksat,iepoch) = ut1part(i,1)
             if( islot(j).eq.11) part(i,j,ksat,iepoch) = polepart(i,2)
             if( islot(j).eq.12) part(i,j,ksat,iepoch) = polepart(i,4)
             if( islot(j).eq.13) part(i,j,ksat,iepoch) = ut1part(i,2) 
           enddo
         elseif( islot(j).gt.100 ) then
c          satellite parameters
           kk=islot(j)/100
           if( kk.eq.ksat) then
             jj=islot(j) - kk*100
             do i=1,3
               part(i,j,ksat,iepoch) = rvec1(6+(jj-1)*3+i)
             enddo
           endif   
         else
           write(message,'(a,i3,a,i5)') 'Invalid islot(',j,')=',islot(j)
          call report_stat('FATAL','ORBFIT','orbits/orbfit',' '
     .                 ,message,0)
         endif
      enddo   
c      print *,'tran_part_t ',tran_part_t
c      print *,'tran_part_i ',tran_part_i
c      print *,'polepart ',polepart
c      print *,'ut1part  ',ut1part 

      return
      end


