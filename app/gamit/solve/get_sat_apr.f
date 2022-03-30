Copyright (c) Massachusetts Institute of Technology and the University of
California at San Diego, 1995. All rights reserved.

C     Get a priori constraints for satellite parametrs

      subroutine get_sat_apr( indorb )

      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'       
      include 'parameters.h'

      real*4 temp(maxorb)
      real*8 factor
      character*256 message
      character*256 wcmd
      character*16 code
      character*4 snam,upperc
      integer j,i1,i2,i0,ib,ic,type,indorb,mchkey,i
      integer lcmd,count_arg,lift_arg,ioerr
c
c     Input
c          rlabel(maxprm) (parameters.h)   character array of parameter labels
c          indorb                          index in rlabel array of first orbit parameter
c
c      Output (in common /constr/ in solve.h)
c          sat_apr(maxsat,maxorb)   real*8      array of tight constaints
c          sat_apr2(maxsat,maxorb)  real*8      array of loose constraints
c          satwt_type               character*3 constraints are Keplerian (KEP) or Cartesian (CAR)

c      Determine type of orbital constraints (Keplerian or Cartesian)

      satwt_type = 'KEP'
      call getcmd(5,'aprior',wcmd,lcmd,2)
      call getcmd(5,'apr_orbit_type',wcmd,lcmd,3)
      if (lcmd.gt.0 ) then
         satwt_type = wcmd(1:3)
         call uppers(satwt_type)
       endif
     
c     Set the defaults but if no input for any satellite, keep satwgt = .false. 

      do i=1,nsat
        if( satwt_type.eq.'KEP' ) then  
c       --for Keplerian elements, convert ppm to km for semi-major axis, 
c         dimensionless for angular elements;  
c         convert percent to dimensionless for radiation pressure parameters   
             factor = 1.d-6
             sat_apr(i,1)  = 10.d0*factor*26000.d0 
             sat_apr2(i,1) = 10.d0*factor*26000.d0
           do j=2,6
             sat_apr(i,j)  = 10.d0*factor
             sat_apr2(i,j) = 10.d0*factor
           enddo 
           do j=7,norb
c            defaults for non-gravitational force parametes are 1000%
             sat_apr(i,j) = 1000.d0*.01d0
             sat_apr2(i,j) = 1000.d0*.01d0 
           enddo                                                     
           do j=norb+1,norb 
             sat_apr(i,j) = 10.d0*.001d0
             sat_apr2(i,j) = 10.d0*.001d0
           enddo
        elseif ( satwt_type.eq.'CAR' ) then
           do j=1,3
             sat_apr(i,j) = 0.2d0
             sat_apr2(i,j) = 0.2d0 
           enddo
           do j=4,6
             sat_apr(i,j) = 3.d-5
             sat_apr2(i,j) = 3.d-5 
           enddo
           do j=7,norb
             sat_apr(i,j) = 1000.d0*.01d0
             sat_apr2(i,j) = 1000.d0*.01d0 
           enddo    
           do j=norb+1,norb 
             sat_apr(i,j) = 10.d0*.001d0
             sat_apr2(i,j) = 10.d0*.001d0
           enddo
        else
            call report_stat('FATAL','SOLVE','get_sat_apr',' '
     .                      ,'Illegal constraint type',0)
        endif
       enddo   
c      SV antenna offsets (units m)
       do i=1,nsat
         do j=norb+1,norb+3 
          if( j.gt.maxorb ) then
           write(message,'(a,i2,a,i2)')
     .        '# SV parameters =',j,' exceeds maxorb=',maxorb
           call report_stat('FATAL','SOLVE','get_sat_apr',' ',message,0)
          endif    
          sat_apr(i,j) = 10.d0
          sat_apr2(i,j) = 10.d0
         enddo
       enddo
            

c      Get a priori for tightly constrained solution

c     Orbit parameters from the T-file (1,norb)
      i2 = 0
      call getcmd(5,'aprior',wcmd,lcmd,2)
      do 120 i0 = 1,1000
         call getcmd(5,'tight_apr_or',wcmd,lcmd,3)
         if (lcmd.le.0) goto 130
         i2 = i2+1
c        decompose command line
         ic = count_arg(wcmd)
c        not enough arguments
         if (ic.le.1) goto 120
c        pointer to
         ib = lift_arg(wcmd,code,1)
         if (ib.le.0) goto 120
         type = 1
         if (upperc(code(1:4)).eq.'ALL_') type = 2
         if (type.eq.1) snam = upperc(code(1:4)) 
         read (wcmd(ib+1:lcmd),*,iostat=ioerr) (temp(i),i=1,norb) 
c*       read (wcmd(ib+1:lcmd),*) (temp(i),i=1,3) 
         if( ioerr.ne.0 ) then
          write(message,'(2a,i2,a)') 'Error reading orbital constraints'
     .     ,' (too few for norb=',norb,' ?)'
          call report_stat('FATAL','SOLVE','get_sat_apr',' ',message,0)
         endif
         do i=1,norb
           if( temp(i).eq.0.d0 ) then
               call report_stat('FATAL','SOLVE','get_sat_apr',' ',
     .           'Orbital tight constraints are zero in batch file',0)
               goto 100
           endif
         enddo            
 100  do 110 i=1,nsat     
            i1=norb*(i-1)+indorb
            if (type.eq.1) then     
               j = mchkey(rlabel(i1),snam,20,4)
               if (j.le.0) goto 110
            endif
            do j=1,norb
              sat_apr(i,j) = temp(j)
            enddo    
 110  continue  
 120  continue
 130  if (i2.gt.0) then
         satwgt = .true.
         if( satwt_type.eq.'KEP' ) then
           if( logprt ) write(6,135)
           write(10,135)
 135       format(/          
     .   ,' Keplerian a priori orbital errors (dimensionless except'
     .   ,' semi-major axis (km))',/
     .   ,' Sat#             Semiaxis   Eccen.  Inclin.   Asc.node  '
     .   ,'Perigee   M.anom.     rad1      rad2      rad3      rad4'
     .   ,'      rad5      rad6      rad7      rad8      rad9'
     .   ,'      rad10     rad11     rad12     rad13',/)
           factor=1.d-06
           do i=1,nsat
             do j=1,norb
c              convert semi-major axis error units from ppm to km
               if(j.eq.1) sat_apr(i,j)=sat_apr(i,j)*factor*26000.d0
c              convert angular elements from ppm to radians
               if(j.ge.2.and.j.le.6) sat_apr(i,j)=sat_apr(i,j)*factor
c              convert non-gravitational parameter units from % to fractions
               if(j.gt.6.and.j.le.norb) sat_apr(i,j)=sat_apr(i,j)*0.01d0
             enddo
             if(logprt) write( 6,140) i,isprn(i),(sat_apr(i,j),j=1,norb)
             write(10,140) i,isprn(i),(sat_apr(i,j),j=1,norb)
  140        format(i3,2x,'PRN ',i2,2x,f13.4,1p30e10.1) 
           enddo
        elseif (satwt_type.eq.'CAR' ) then
           if( logprt ) write(6,145)
           write(10,145)
  145      format(/,
     .' Cartesian a priori orbital errors (km km/s) + force paramters'
     .,' (fraction))',/,
     .' Sat#      X          Y          Z          Xdot       Ydot     '
     .,'   Zdot     rad1      rad2      rad3      rad4      rad5'
     .,'      rad6      rad7      rad8      rad9 '
     .,'    rad10     rad11    rad12     rad13',/)
           do  i=1,nsat
             do j=7,norb
c              convert rad parameter units from percent to fraction
               if(j.gt.6.and.j.le.norb) sat_apr(i,j)=sat_apr(i,j)*0.01d0
             enddo
             if( logprt ) write( 6,150) i,(sat_apr(i,j),j=1,norb)
             write(10,150) i,(sat_apr(i,j),j=1,norb)
  150        format(i3,6f6.3,1p9e10.1)    
           enddo
        endif
      endif
                    
c     SV antenna offsets (norb+1,norb+3)

      if( svantest ) then

      i2 = 0
      call getcmd(5,'aprior',wcmd,lcmd,2)
      do 170 i0 = 1,1000
         call getcmd(5,'tight_apr_sv',wcmd,lcmd,3)
         if (lcmd.le.0) goto 180
         i2 = i2+1
c        decompose command line
         ic = count_arg(wcmd)
c        not enough arguments
         if (ic.le.1) goto 170
c        pointer to
         ib = lift_arg(wcmd,code,1)
         if (ib.le.0) goto 170
         type = 1
         if (upperc(code(1:4)).eq.'ALL_') type = 2
         if (type.eq.1) snam = upperc(code(1:4))
         read (wcmd(ib+1:lcmd),*,iostat=ioerr) (temp(i),i=1,3) 
         if( ioerr.ne.0 ) call report_stat('FATAL','SOLVE','get_sat_apr'
     .                  ,' ','Error reading SV ant tight constraints',0)
         do i=1,3
           if( temp(i).eq.0.d0 ) then
               call report_stat('FATAL','SOLVE','get_sat_apr',' ',
     .          'SV antenna tight constraints are zero in batch file',0)
               goto 160
           endif
         enddo
 160  do 165 i=1,nsat   
c           i1 is index of first SV ant offset in parameters array 
            i1=3*(i-1)+ (indorb+norb*nsat)
            if (type.eq.1) then
               j = mchkey(rlabel(i1),snam,20,4)
               if (j.le.0) goto 165
            endif
            do j=1,3                 
              sat_apr(i,norb+j) = temp(j)
            enddo
 165  continue
 170  continue
 180  if (i2.gt.0) then
         if( logprt ) write( 6,185)
         write(10,185)
 185     format(/,' Apriori SV antenna offset errors (m)',/,
     .          ' Sat#    dX       dY       dZ ',/)
         do i=1,nsat
           if( logprt ) 
     .         write( 6,'(i4,3x,3f8.3)') i,(sat_apr(i,norb+j),j=1,3)
           write(10,'(i4,3x,3f8.3)') i,(sat_apr(i,norb+j),j=1,3)
         enddo
      endif

      endif

  
c      Get a priori for loosely constrained solution

      i2 = 0
      call getcmd(5,'aprior',wcmd,lcmd,2)
      do 240 i0 = 1,1000
         call getcmd(5,'loose_apr_or',wcmd,lcmd,3)
         if (lcmd.le.0) goto 250
         i2 = i2+1
c        decompose command line
         ic = count_arg(wcmd)
c        not enough arguments
         if (ic.le.1) goto 240
c        pointer to
         ib = lift_arg(wcmd,code,1)
         if (ib.le.0) goto 240
         type = 1
         if (upperc(code(1:4)).eq.'ALL_') type = 2   
         if (type.eq.1) snam = upperc(code(1:4))
         read (wcmd(ib+1:lcmd),*) (temp(i),i=1,norb)
         do i=1,norb
           if( temp(i).eq.0.d0 ) then
               call report_stat('FATAL','SOLVE','get_sat_apr',' ',
     .           'Orbital loose constraints are zero in batch file',0)
               goto 210
           endif
         enddo
  210    do 230 i=1,nsat
           i1=norb*(i-1)+indorb
           if (type.eq.1) then
              j = mchkey(rlabel(i1),snam,20,4)
              if (j.le.0) goto 230
           endif
           do j=1,norb
             sat_apr2(i,j) = temp(j)
           enddo
           do  j=1,norb
             if( satwt_type.eq.'KEP') then
               factor = 1.d-06
c              convert semi-major axis error units to km
               if(j.eq.1) sat_apr2(i,j)=sat_apr2(i,j)*factor*26000.d0
c              convert other element error units to radians
               if(j.ge.2.and.j.le.6) sat_apr2(i,j)=sat_apr2(i,j)*factor
               if(j.gt.6.and.j.le.norb) 
     .           sat_apr2(i,j)=sat_apr2(i,j)*0.01d0
             elseif( satwt_type.eq.'CAR') then
c              convert rad parameter units from percent to fraction
c              convert non-gravitational parameter error units to fractions
               if(j.gt.6.and.j.le.norb) 
     .           sat_apr2(i,j)=sat_apr2(i,j)*0.01d0
             endif
           enddo
 230    continue
 240  continue
 250  continue

      
c     SV antenna offsets (norb+1,norb+3)   
              
      if( svantest ) then

      i2 = 0
      call getcmd(5,'aprior',wcmd,lcmd,2)
      do 270 i0 = 1,1000
         call getcmd(5,'loose_apr_sv',wcmd,lcmd,3)
         if (lcmd.le.0) goto 280
         i2 = i2+1
c        decompose command line
         ic = count_arg(wcmd)
c        not enough arguments    
         if (ic.le.1) goto 270 
c        pointer to
         ib = lift_arg(wcmd,code,1)
         if (ib.le.0) goto 270
         type = 1
         if (upperc(code(1:4)).eq.'ALL_') type = 2   
         if (type.eq.1) snam = upperc(code(1:4))
         read (wcmd(ib+1:lcmd),*,iostat=ioerr) (temp(i),i=1,3) 
         if( ioerr.ne.0 ) call report_stat('FATAL','SOLVE','get_sat_apr'
     .         ,' ','Error reading SV ant loose constraints',0)
         do i=1,3
           if( temp(i).eq.0.d0 ) then
               call report_stat('FATAL','SOLVE','get_sat_apr',' ',
     .          'SV antenna loose constraints are zero in batch file',0)
               goto 260
           endif
         enddo
 260  do 265 i=1,nsat   
c           i1 is index of first SV ant offset in parameters array 
            i1=3*(i-1)+ (indorb+norb*nsat)
            if (type.eq.1) then
               j = mchkey(rlabel(i1),snam,20,4)
               if (j.le.0) goto 265
            endif
            do j=1,3                 
              sat_apr2(i,norb+j) = temp(j)
            enddo
 265  continue
 270  continue
 280  continue 

      endif
              
      return
      end
