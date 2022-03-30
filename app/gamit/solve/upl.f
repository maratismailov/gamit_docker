C
      Subroutine UPL

C     Updates an L-file

C     Modified from UPSIT -- Yehuda Bock 11/25/91
c     Modified by T Herring 990624 to clean up writing so that blanks are not written
c         also added comment to line so history can be traced. 
c     Modified by R King 020807 to allow for updating an apr (rather than l-) file 
c     
 
      implicit none

      include '../includes/dimpar.h'  
      include 'solve.h' 
      include 'parameters.h'
                            
      character*1 cvalue
      CHARACTER*4 SNAME(MAXSIT),STRBK,lowerc  
      character*8 site8
      CHARACTER*16 TMPNAM
      CHARACTER*16 LSTR(2*MAXSIT)
      character*256 line

      LOGICAL FCHECK,lupdate(maxsit)

      integer*4 istat,lstat,iparm,ilive,indx,indx2
     .        , jd0,ioerr,kfflg,i,j,k

* nblen  -- Gamit function to return non-blank length of string
      integer nblen

      real*8 erad,coslat,lcoords(6),newpos(3),alat,alon,radius,t0
     .     , epoch,epoch0,rvalue,prelatr,prelonr,prerad,prexyz(3)
                
c     Function
      integer*4 julday 

      data erad/6378137.d0/


      if( logprt ) write(6,'(/,2a,/,2a,/,a,f8.3,a)') 
     .                  'Updating coordinate file : ',linf
     .                , 'New coordinate file      : ',loutf
     .                , 'Coordinate tolerance     : ',coord_upd_tol,' m'
      if( iqflag.eq.1 )
     .  write(10,'(/,2a,/,2a,/,a,f8.3,a)' ) 
     .                  'Updating coordinate file : ',linf
     .                , 'New coordinate file      : ',loutf
     .                , 'Coordinate tolerance     : ',coord_upd_tol,' m'
 


C     LOOP THRU ALL STATIONS AND UPDATE COORDINATES
      ilive = 0
      do istat = 1,nsite  
c        Do not update a station's coordinates if the adjustment is less
c        than the input threshold (coord_upd_tol, nominally 30 cm) or
c        the sigma exceeds 10 m.
         lupdate(istat) = .false. 
c          latitude (radians) 
         iparm = 3*(istat-1) + 1     
         sname(istat) = rlabel(iparm)(1:4)  
cd         print *,'istat sname ',istat,sname(istat)
         if( free(iparm).gt.0 ) then
           ilive = ilive + 1
cd           print *,'istat iparm ilive adj sig  ',istat,iparm,ilive
c     .            , adjust(iparm),sigma(ilive)   
           if( dabs(adjust(iparm))*erad.gt.coord_upd_tol .and.
     .              sigma(ilive)*erad.lt.10.d0 ) then
             lupdate(istat) = .true.
cd             print *,'lupdate ',lupdate(istat)   
           endif
         endif

c          longitude (radians)
         iparm = 3*(istat-1) + 2   
         if( free(iparm).gt.0 ) then
           ilive = ilive + 1
           coslat = dcos(postvl(iparm-1))
cd           print *,'istat iparm ilive adj sig ',istat,iparm,ilive
c     .            , adjust(iparm), sigma(ilive) 
           if( dabs(adjust(iparm))*erad*coslat.gt.coord_upd_tol .and.
     .            sigma(ilive)*erad*coslat.lt.10.d0 ) then
              lupdate(istat) = .true.  
cd              print *,'lupdate ',lupdate(istat)           
           endif
         endif
c
c          radius  
         iparm = 3*(istat-1) + 3  
         if( free(iparm).gt.0 ) then  
          ilive = ilive + 1
cd           print *,'istat iparm ilive adj sig ',istat,iparm,ilive
cd     .            ,adjust(iparm), sigma(ilive)   
           if( dabs(adjust(iparm))*1.d3.gt.coord_upd_tol .and.
     .            sigma(ilive)*1.d3.lt.10.d0 ) then
             lupdate(istat) = .true.      
cd             print *,'lupdate ',lupdate(istat)    
           endif
         endif
c        Convert from N positive radians to string
         indx = 2*(istat-1)
         iparm = 3*(istat-1) + 1
         call wdms(1,postvl(iparm),lstr(indx+1))
c        Convert from E positive radians to string
         iparm = 3*(istat-1) + 2
         call wdms(2,postvl(iparm),lstr(indx+2))
c        take out the colons for safety.
         do  k=1,16
            if (lstr(indx+1)(k:k) .eq. ':') lstr(indx+1)(k:k) = ' '
            if (lstr(indx+2)(k:k) .eq. ':') lstr(indx+2)(k:k) = ' '
         enddo  
      enddo
      

C     LOOP THRU L or APR FILE   

      TMPNAM=linf
      call lowers(tmpnam)
      if (fcheck(tmpnam)) then
         OPEN (UNIT=16,FILE=tmpnam,STATUS='OLD',IOSTAT=IOERR)
         IF (IOERR.NE.0) THEN
           call report_stat('WARNING','SOLVE','upl',tmpnam
     .          , 'Cannot open coordinate file--no update',ioerr)
            RETURN
         ENDIF
      else
         call report_stat('WARNING','SOLVE','update',tmpnam
     .          , 'Cannot find old coordinate file--no update',ioerr)
         RETURN
      endif              
c     determine if l- or apr file 
      kfflg = 0
      call crd_file_type( linf,kfflg)   
cd      print *,'CRD_FILE kfflg ',kfflg                   
c     if an apr file, get midpoint epoch 
      if( kfflg.eq.1 ) then
        jd0 = julday(it0(1),it0(2),it0(3))    
        t0 =  3600.d0*t00(1) + 60.d0*t00(2) + t00(3) +
     .       dfloat((nepoch-1)*inter)/2.d0     
c       jd_to_decyrs wants true JD, not PEP JD  
        call jd_to_decyrs( dfloat(jd0) -0.5d0  + t0/86400.d0, epoch )
      endif   
cd      print *,'it0 t00 jd0 t0 epoch ',it0,t00,jd0,t0,epoch

      OPEN (UNIT=17,STATUS='SCRATCH',ERR=601)

C     SEARCH FOR ALL OCCURRENCES OF STATION
      do 220 j=1,32000
         READ(16,'(a)',END=209) line
         if( kfflg.eq.0 ) then
            STRBK = line(1:4) 
         else            
c          if apr-file, skip a comment line
           if( line(1:1).ne.' ' ) then 
             strbk = '   '
           else  
c            site name can start in any column 
             indx = 2 
             call read_line( line,indx,'CH',ioerr,rvalue,site8 )
             strbk = site8(1:4)
           endif 
         endif  
cd         print *,'strbk ',strbk

C        LOOK FOR MATCHING 4-CHARACTER STATION NAME  
         DO 500 ISTAT=1,nsite
cd            print *,'strbk,sname',strbk,sname(istat)
            IF(lowerc(STRBK).NE.lowerc(SNAME(ISTAT))) GO TO 500  
C           WE FOUND IT 
            lstat = istat
            INDX =2*(ISTAT-1)
            INDX2=3*(ISTAT-1)                    
            GO TO 501
 500     CONTINUE

C        NO MATCH
         WRITE(17,'(A)') line(1:max(1,nblen(line)))
         go to 220

 501     CONTINUE

C        A MATCH - Update the coordinates if this station was set to be updated 
c                  by the adjustments exceeding the threshold test

         if( lupdate(lstat) ) then 
cd           print *,'Updating station ',lstat
               
c          Convert radius to meters
           radius=postvl(indx2+3)*1000.d0  
           if( kfflg.eq.0 ) then
c            find and update the fields of the L-file; since the format is fixed
c            we can overwrite by column number the part of the string with coordinates
             WRITE(line(17:62),'(A16,1X,A16,1X,F12.4)')
     .               LSTR(INDX+1),LSTR(INDX+2),RADIUS
c            add a comment string at the end of the line  
             write(line(63:),'(2x,2a)') 'Updated from ', linf 
           elseif( kfflg.eq.1 ) then 
c            find and update the fields of the apr file; since the format is free,
c            we need to find the positions of the values 
             iparm = 3*(istat-1) + 1 
             alat = postvl(iparm) 
             prelatr = preval(iparm)
             iparm = 3*(istat-1) + 2    
             alon = postvl(iparm)   
             prelonr = preval(iparm)   
             iparm = 3*(istat-1) + 3
             prerad = preval(iparm)*1.d3
             newpos(1) = radius*dcos(alat)*dcos(alon)
             newpos(2) = radius*dcos(alat)*dsin(alon)
             newpos(3) = radius*dsin(alat)  
cd             print *,'newpos ',newpos     
c            Compute the apriori cartesian coordinates for check multiple apr-file entries below   
             call sph2xyz( prelatr,prelonr,prerad,prexyz )      
cd             print *,'prelatr prelonr prerad prexyz '
cd     .              , prelatr,prelonr,prerad,prexyz
             indx= 2 
             call read_line( line,indx,'CH',ioerr,rvalue,site8 )
             do i = 1,6                             
               lcoords(i) = 0.d0
               call read_line( line,indx,'R8',ioerr,lcoords(i),cvalue ) 
             enddo  
             call read_line( line,indx,'R8',ioerr,epoch0,cvalue) 
cd             print *,'coords from l-file midpt ',lcoords,epoch0,epoch
             do i=1,3
                lcoords(i) = lcoords(i) + lcoords(i+3)*(epoch-epoch0)
             enddo 
cd             print *,'corrected coords ',(lcoords(i),i=1,3)

c ** rwk 080610: Update only if coordinates match w/in 1 mm (surrogate for 8-character names)
cd          print *,'strbk site8 coords prexyz ',strbk,site8,lcoords,prexyz
             if( dabs(lcoords(1)-prexyz(1)).lt.1.d-3 .and.
     .           dabs(lcoords(2)-prexyz(2)).lt.1.d-3 .and.
     .           dabs(lcoords(3)-prexyz(3)).lt.1.d-3 ) then
cd                print *,'updating coords epoch ',lcoords,epoch
c            it's too complicated to try to match the existing format--and may lose precision
c             rwk 060524: change vel output from f10.5 to f11.5 to allow for motions on ice
               write( line,'(1x,a8,3f15.5,3f11.5,f10.4,2x,2a)') 
     .             site8,(newpos(i),i=1,3),(lcoords(i),i=4,6),epoch
     .            ,'Updated from ', linf       
             endif 
           endif
         endif     
         WRITE(17,'(a)') line(1:max(1,nblen(line)))
      
220   continue
C     End global search on L-file

209   CLOSE(UNIT=16)
      REWIND 17
      TMPNAM = loutf
      CALL LOWERS(TMPNAM)
      OPEN (UNIT=18,FILE=TMPNAM,STATUS='UNKNOWN')
C     COPY SCRATCH FILE TO NEW S-FILE
      DO  I=1,32000
         READ(17,'(a)',END=601) line
         WRITE(18,'(a)') line(1:max(1,nblen(line)))
      enddo

  601 CLOSE(UNIT=18)
      CLOSE(UNIT=17)  
      RETURN
      END
