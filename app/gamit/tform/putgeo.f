      Subroutine PUTGEO( idms,numsit,x,sitnam,datum,epoch,ofile)

C       Write out a set or file of geodetic coordinates.
c       recoded by pjm April 98 to fix equatorial lat problem  
c       Recoded by rwk October 2002 to handle new format and ITRF->NAD83 conversions

c       If called by TFORM, 'datum' will be blank, to be entered interactively
c       If called by PUTGEO, 'datum' will be used, with no questions asked

      implicit none

      include '../includes/tform.h'

      integer ofile
      integer*4 ioerr,idms,numsit,ltdeg,ltmin,lgdeg,lgmin,nblen,i

      real*8 x(3,nsdim),semi,finv,shft(3),shftdot(3),rot(3),rotdot(3)
     .     , scale,scaledot,refepoch,epoch,alat,along,ht,altsec,algsec
     .     , geodrad,pi,convmas
     
      character *1 nors,eorw   
      character *5 datum  
      character *8 siteid
      character*12 sname   
      character*16 sitnam(nsdim)
      character*60,lfmt1
      character*63 lfmt2   
      character*256 message,line 

      data lfmt1/     
     .'(a4,1x,a12,a1,i2,1x,i2,1x,f8.5,1x,a1,i3,1x,i2,1x,f8.5,f13.4)'/  
      data lfmt2/
     .'(1x,a4,1x,a12,a1,i2,1x,i2,1x,f8.5,1x,a1,i3,1x,i2,1x,f8.5,f13.4)'/
 

      logical fcheck,eof
                         
c     'sitnam' stores a 4-character id followed by a 12-character descriptor.
c     For apr-type files, these are output as an 8-character name SITE_GPS;
c     for l-files, the output is 'SITE description '
                 
                           
c     Allowed formats
c     ---------------
c     if screen output (one site): free-format decimal degrees, no site name

c     if file output with decimal degrees:  8-character site-name beginning in column 2, 
c                                           free-format coordinates
c     if file output with deg/min/sec:  L-file format, 4-char id in column 1, name in columns 5-16
     

           
c       Compute constants

      pi = 4.d0*datan(1.d0) 
      convmas = 4.84813481d-9

              
c       Open the datum file to display its contents, then close and call 
c         library routine read_gdatum to get the values.  

      if (fcheck('gdetic.dat')) then
        open(unit=idatum,file='gdetic.dat',iostat=ioerr,status= 'OLD')
        if (ioerr .ne. 0)  then
          call report_stat('FATAL','TFORM','getgeo',' '
     .                             ,'Error opening gdetic.dat',ioerr)    
        else  
           eof = .false.
           if( iscrn.gt.0) write(iscrn,'(/,a)') 
     .            'Datums available from gdetic.dat:'
           do while (.not.eof )           
             read(idatum,'(a)',iostat=ioerr) line
             if( ioerr.eq.0 ) then 
               if( iscrn.gt.0 ) 
     .           write(iscrn,'(1x,a)') line(1:nblen(line))
             elseif( ioerr.eq.-1 ) then
               eof = .true.
             else
               call report_stat('FATAL','TFORM','getgeo',' '
     .                          ,'Error reading gdetic.dat',ioerr)
             endif  
          enddo
        endif    
        close( idatum ) 
        if( iscrn.gt.0 .and.iterm.gt.0 ) then   
   5      write(iscrn,'(a)') 
     .      'Enter 5-character code for datum (case insensitive):' 
          read(iterm,'(a5)',iostat=ioerr) datum   
          call uppers(datum) 
          if( ioerr.ne.0 ) go to 5    
        endif
        call read_gdatum( idatum,datum,semi,finv,shft,shftdot
     .                  , rot,rotdot,scale,scaledot,refepoch ) 
        if( iscrn.gt.0 ) write(iscrn,'(a,a5,a,f11.3,a,f13.9)') 
     .       'Selected ',datum,' a=',semi,' 1/f=',finv  
        if( refepoch.eq.0.d0 ) then
          if( iscrn.gt.0)  write(iscrn,'(/,14x,a,3f9.3)') 'dX=',shft
        else     
          if( refepoch.lt.1900.d0 .or. refepoch.gt.2100.d0 ) then
             write(message,'(a,f20.1)') 'Invalid epoch from gdetic.dat:'
     .               ,refepoch
             call report_stat('FATAL','TFORM','putgeo',' ',message,0) 
          endif
          if( iterm.gt.0 .and. epoch.eq.0.d0 ) then 
    6       write(iscrn,'(a)') 
     ,       'Enter epoch for transformation (decimal yr, e.g. 1997.0):'
            read(iterm,*,iostat=ioerr) epoch
            if( ioerr.ne.0 ) then
              write(iscrn,'(a,f20.0)') 'Error reading epoch: ',epoch
              goto 6                                   
            elseif( refepoch.lt.1900.d0 .or. refepoch.gt.2100.d0 ) then
              write(iscrn,'(a,f20.0)') 'Invalid epoch: ',epoch
              goto 6
            else
              write(iscrn,'(/,2(a,3f9.3),a)') 
     .          '  dX=',shft,'  dXdot=',shftdot,'(m m/yr)'
              write(iscrn,'(2(a,3f9.3),a)') 
     .              '  R =',(rot(i)/convmas,i=1,3)
     .             ,'  Rdot =',(rotdot(i)/convmas,i=1,3),'(mas mas/yr)'
              write(iscrn,'(2(a,d9.2),a)') '  Scale=',scale*1.d9
     .        ,'  Sdot=',scaledot*1.d9,' (ppb ppb/yr)'
              write(iscrn,'(a,f7.2)')     '  Ref Epoch=',refepoch    
            endif
          endif
          shft(1) = shft(1) + shftdot(1)*(epoch-refepoch)
          shft(2) = shft(2) + shftdot(2)*(epoch-refepoch)
          shft(3) = shft(3) + shftdot(3)*(epoch-refepoch) 
          rot(1)  = rot(1)  + rotdot(1) *(epoch-refepoch) 
          rot(2)  = rot(2)  + rotdot(2) *(epoch-refepoch) 
          rot(3)  = rot(3)  + rotdot(3) *(epoch-refepoch) 
          scale = scale + scaledot *(epoch-refepoch)   
          if( iscrn.gt.0 ) then
            write(iscrn,'(/,a,f7.2)') 'Transformation at epoch ',epoch
            write(iscrn,'(/,a,3f12.3,a)') '  dX=',shft,' m'
            write(iscrn,'(a,1p3d12.4,a)') '  dR=',rot,' radians'
            write(iscrn,'(a,d12.4)') '  dS=',scale 
          endif
        endif 

      else     
        call report_stat('WARNING','TFORM','putgeo',' '
     .  ,'Missing gdetic.dat, use NAD83 with time-dependent corrections'
     .       ,0)
        semi=6378137.d0
        finv=298.257222101d0   
        scale = 0.d0 
        scaledot = 0.d0
        refepoch = 0.d0
        do i=1,3      
          shft(i) = 0.d0
          shftdot(i) = 0.d0
          rot(i) = 0.d0 
          rotdot(i) = 0.d0 
        enddo
      endif
        

c      Convert and output the coordinates
             
      if( idms.eq.1 ) then   
         continue
      elseif( idms.eq.2 ) then     
c       Deg min sec in L-file format: no header or format statement
        continue             
      else
        call report_stat('FATAL','TFORM','putgeo',' ' 
     .        ,'Illegal type of coordinate (idms)',0)   
      endif
c     skip a line
      if( iprnt.gt.0 ) write(iprnt,'(1x)')


c     Convert to geodetic coordinates
                  
      if( refepoch.ne.0.d0.and.iprnt.gt.0 ) 
     .    write(iprnt,'(a)') 
     .      'Cartesian coordinates after applying datum shift: '
      do i=1,numsit     

         x(1,i) = shft(1) + (1.d0+scale)*x(1,i) 
     .                 + rot(3)*x(2,i) - rot(2)*x(3,i) 
         x(2,i) = shft(2) + (1.d0+scale)*x(2,i)
     .                 - rot(3)*x(1,i) + rot(1)*x(3,i)
         x(3,i) = shft(3) + (1.d0+scale)*x(3,i)
     .                 + rot(2)*x(1,i) - rot(1)*x(2,i)  
         if( refepoch.ne.0.d0.and.iprnt.gt.0 ) then
             write(iprnt,'(3f14.4)') x(1,i),x(2,i),x(3,i)  
            if( numsit.eq.1 ) write(iprnt,'(1x)')
         endif
         call geoxyz( semi,finv,alat,along,ht,geodrad
     .              , x(1,i),x(2,i),x(3,i),2 ) 
         alat = alat*180.d0/pi
         along = along*180.d0/pi   
         if( idms.eq.1 ) then   
            call getsname(siteid,sname,sitnam(i),2)
            if( iprnt.gt.0 ) write(iprnt,'(1x,a8,2f14.8,f14.4)')
     .               sitnam(i)(1:8),alat,along,ht
            if( ofile.gt.0 ) write(ofile,'(1x,a8,2f14.8,f14.4)')
     .               sitnam(i)(1:8),alat,along,ht   
         elseif( idms.eq.2 ) then  
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
     .        write(ofile,lfmt1) siteid,sname,nors,ltdeg,ltmin
     .                         , altsec,eorw,lgdeg,lgmin,algsec,ht    
           if( iprnt.gt.0 )
     .        write(iprnt,lfmt2) siteid,sname,nors,ltdeg,ltmin
     .                         , altsec,eorw,lgdeg,lgmin,algsec,ht 
         endif

      enddo   
c     end of loop on sites

      return
      end
