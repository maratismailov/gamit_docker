      SUBROUTINE PUTSPH( IDMS,NUMSIT,X,SITNAM,OFILE)
C
C       Write out a set or file of spherical coordinates.
c23456789012345678901234567890123456789012345678901234567890123456789012
c       re coded by pjm April 98

      implicit none

      include '../includes/tform.h'

      integer*4 ofile,idms,numsit,ltdeg,ltmin,lgdeg,lgmin,i
c*      integer izero 

      real*8 x(3,nsdim),alat,along,rad,altsec,algsec

      character *1 nors,eorw   
      character *8 siteid  
      character*12 sname
      character*16 sitnam(nsdim) 
      character*60,lfmt1
      character*63 lfmt2

      data lfmt1/     
     .'(a4,1x,a12,a1,i2,1x,i2,1x,f8.5,1x,a1,i3,1x,i2,1x,f8.5,f13.4)'/  
      data lfmt2/
     .'(1x,a4,1x,a12,a1,i2,1x,i2,1x,f8.5,1x,a1,i3,1x,i2,1x,f8.5,f13.4)'/
 

c     'sitnam' stores a 4-character id followed by a 12-character descriptor.
c     For apr-type files, these are output as an 8-character name SITE_GPS;
c     for l-files, the output is 'SITE description '
                 
                           
c     Allowed formats
c     ---------------
c     if screen output (one site): free-format decimal degrees, no site name

c     if file output with decimal degrees:  8-character site-name beginning in column 2, 
c                                           free-format coordinates
c     if file output with deg/min/sec:  L-file format, 4-char id in column 1, name in columns 5-16
                
c
C       Write the output file

         IF( IPRNT.GT.0 ) WRITE(IPRNT,'(1X)')
      do i=1,numsit
         CALL CARSPH(X(1,I),ALAT,ALONG,RAD)
         IF( IDMS.EQ.1 ) THEN     
            call getsname(siteid,sname,sitnam(i),2)
            IF( OFILE.GT.0) WRITE(OFILE,'(1x,a8,f14.9,1x,f16.9,f13.4)')
     .                          siteid,ALAT,ALONG,RAD
            IF( IPRNT.GT.0) WRITE(IPRNT,'(1x,a8,f14.9,1x,f16.9,f13.4)')
     .                          siteid,ALAT,ALONG,RAD
         ELSEIF( IDMS.EQ.2 ) THEN   
            call degdms(alat,nors,ltdeg,ltmin,altsec) 
            if( nors.eq.'-' ) then
               nors = 'S'
            else
               nors = 'N'
            endif       
            call degdms(along,eorw,lgdeg,lgmin,algsec)   
            if( eorw.eq.'-' ) then
                eorw = 'W'
            else
                eorw = 'E'
            endif
            call getsname(siteid,sname,sitnam(i),1)
            if( ofile.gt.0 )
     .         write(ofile,lfmt1) siteid,sname,nors,ltdeg,ltmin
     .                          , altsec,eorw,lgdeg,lgmin,algsec,rad    
            if( iprnt.gt.0 )
     .         write(iprnt,lfmt2) siteid,sname,nors,ltdeg,ltmin
     .                          , altsec,eorw,lgdeg,lgmin,algsec,rad 
         endif
      enddo 

      return
      end


