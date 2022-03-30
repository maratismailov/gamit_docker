Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.

      Subroutine LSQDO1( covkep,free_fix )
c
c     Write results to O-file, Q-file and screen
C
C     Written by Yehuda Bock and Sergei Gourevitch
C     Modified to add local system (north,east,up) - yb 8/2/87
c     Modified and cleaned up -- K. Feigl May 8, 1990
c     Moved inside of write_soln --  R. King 93/9/21
c
c flow chart :
c
c     LSQDO1
c        |
c        +---> (get postfit values)
c        |
c        +---> (page loop)
c           |
c           +---> QHEAD3                        } output parameter adjustments & sigmas
c        |
c        +---> (baseline loop)
c           |
c           +---> COVBLN                        } get covariance of baselines
c           |
c           +---> (output)                      } baseline vector and uncertainties
c           |
c           +---> OLINE                         } one-line output to O-file
c        |
c        +---> UPDAT3                           } questions of updating M-, L-, I-, G-files
c

      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'
      include 'parameters.h'

      integer jparm,ilive,ictyp,istat,jstat,ijk,i1,i2,j,k

      real*8 coslat,sigtmp,dummy,baseln,dbase

      character*1 symbl
c      character*1 ff
      character*4 free_fix
c**   character*8 stru
cc    character*4 buf4a,buf4b,upperc
cc    character*256 message

c     1 = lat, 2 = lon, 3 = rad
      integer kcoord

      logical iflg

      real*8  sigs(3,3),siglcl(3,3),covkep(maxcov)
      integer id(maxsit)

cc rwk 190509       equivalence (dr,postvl)
CC      equivalence (id,work1)
          
      logical debug/.false./

      if(ntpart.eq.0) call report_stat('FATAL','SOLVE','lsqdo1',' '
     .                   , 'No parameters estimated',0)

*  Apply adjustments to the prefits
      do  j=1,ntpart
         postvl(j)=preval(j) + adjust(j)
      enddo

*  Initialize the live parameter counter
      ilive=0

*  Loop over all parameters
                    
      do jparm = 1,ntpart         
  
*        Print the title line, but only once for o-files
         if(ioflag.eq.0) write(15,'()')
         if( jparm.eq.1 .and. logprt )   write(6,3001)
         if (iqflag.eq.1 .and. jparm.eq.1 ) write(10,3001)
         if (ioflag.eq.1 .and. jparm.eq.1)  write(15,3001)

         if(debug) print *,'LSQDO1 ntpart jparm ',ntpart,jparm
*        Is the paramter a lat, lon or a rad?
         j = islot1(jparm)
          kcoord = (j-1)/100 + 1
         if(debug)  print *,'LSQDO1 jparm islot1(jparm) j kcoord '
     .                        , jparm,islot1(jparm),j,kcoord

c        compute cos(lat) for the station to use in converting the
c        longitude adjust and sigma to meters (assumes the next parameter
c        is longitude for this same station--should always work, and in any
c        case applies only to q- and o-file printout)
         if( kcoord.eq.1 ) coslat = dcos(preval(jparm))
         if (kcoord.gt.3.or.kcoord.lt.1) kcoord = 0

         k=free(jparm)
         if(k.eq.0) then 

*          Parameter is not adjusted
           symbl=' ' 
           call qhead3( jparm,symbl,dummy,kcoord,coslat,1)
         
         else

*          Parameter is adjusted
           symbl='*'
           ilive = ilive+1 
           sigtmp = sigma(ilive)    
           call qhead3(jparm,symbl,sigtmp,kcoord,coslat,2)

         
        endif 

      enddo

c--Formats for adjustment printout
 3001  format (//,7X,'Label (units)',12X,'a priori',5X,'Adjust (m) ',
     $   2X,' Formal  Fract',5x,'Postfit')
  
*  Write to the q-file a warning if the adjustments are larger 
*  than the tolerance or larger than twice their aprior constraints

      if( iqflag.eq.1 ) then     
        call check_adjust( rlabel,postvl,free,sigma )  
      endif      

*  Now do the baseline estimates and covariances

*     check if station partials are present; if not skip baseline covariance 
*     propagation. We assume that station coordinates are first in the list
      if(islot1(1).le.300) then 

*       ictyp=1 : spherical coordinates (geodetic no longer supported)
        ictyp=1
        do istat = 1,nsite
          id(istat) = istat
        enddo
              
*       Propagate to cartesian covariances
        do  istat=1,nsite-1 
          i1=(istat-1)*3+1
          do  jstat=istat+1,nsite   
            i2=(jstat-1)*3+1  
            call covbln(istat,jstat,sigs,baseln,dbase,iflg,siglcl)
            if( logprt ) write(6,1895) rlabel(i1)(1:4),istat
     .          ,rlabel(i2)(1:4),jstat,(sigs(ijk,1),ijk=1,3),baseln
            if(iqflag.eq.1)  write (10,1895) rlabel(i1)(1:4),istat
     .          , rlabel(i2)(1:4),jstat,(sigs(ijk,1),ijk=1,3),baseln
*           write a single - line record to o-file (x,y,z)
            if(ioflag.eq.1)
     .      call oline(id,istat,jstat,free_fix,sigs,baseln,dbase,1)
            if(iflg) then
*             cartesian (x,y,z) covariances
              if( logprt ) write(6,1896)
     .          (sigs(ijk,2),ijk=1,3),dbase,(sigs(ijk,3),ijk=1,3)
              if(iqflag.eq.1) write(10,1896)
     .          (sigs(ijk,2),ijk=1,3),dbase,(sigs(ijk,3),ijk=1,3)   
              if( logprt ) write(6,2896)(siglcl(ijk,2),ijk=1,3),dbase
     .                               , (siglcl(ijk,3),ijk=1,3) 

              if(iqflag.eq.1) then 
                write (10,2895) (siglcl(ijk,1),ijk=1,3),baseln
                write(10,2896)
     .            (siglcl(ijk,2),ijk=1,3),dbase,(siglcl(ijk,3),ijk=1,3)
              endif 
              if(ioflag.eq.1)
     .         call oline(id,istat,jstat,free_fix,siglcl,baseln,dbase,2)
            else
*             local covariances
*             note that format is different for o and q files
              if( logprt ) write(6,2895) (siglcl(ijk,1),ijk=1,3),baseln
              if(iqflag.eq.1)  write (10,2895) 
     .          (siglcl(ijk,1),ijk=1,3),baseln
            endif
          enddo
        enddo
      endif                      


*  Formats for the baseline printout
 1895 format (/, ' Baseline vector (m ): ', a4,8x, '(Site', i2, ') to ', 
     +  a4,8x, '(Site', i2, ')'/' X', f15.5, ' Y(E)', f15.5, ' Z', 
     +  f15.5, '  L', f15.5)
 1896 format (6x, '+-', f9.5, 3x, 2(6x, '+-', f9.5), 7x, '+-', f9.5, 
     +  '  (meters)', /, ' Correlations (X-Y,X-Z,Y-Z) = ', 3f10.5)
 2895 format (' N', f15.5, ' E   ', f15.5, ' U', f15.5, '  L', f15.5)
 2896 format (6x, '+-', f9.5, 3x, 2(6x, '+-', f9.5), 7x, '+-', f9.5, 
     +  '  (meters)'/' Correlations (N-E,N-U,E-U) = ', 3f10.5)

*  Write the zenith-delay and gradient summaries to the o-file
             
      if(ioflag.eq.1) then
        call zenout( free_fix,id )
        call gradout(free_fix,id )
      endif
      

*  Update the M-, L-, I-, and G-files
    
      call updat3( covkep )          

      return
      end
