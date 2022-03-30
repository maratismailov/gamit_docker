      Subroutine gmake (prog,tname,gname,itsat )

c   Program gmake to create a G-file by from a T-file, either by
c   reading the header or interpolating

c   R. King  19 Jul 94
c     Roots:
c         catot - yb 20 Jun 88
c         gmaket- rwk 3 Apr 87
          
c  P Tregoning May 1995: Write the inertial frame and the precession 
c                        model onto the gfile header   
c  RWK November 2015: Remove interactive options

      implicit none
c
      include '../includes/dimpar.h'

      integer*4 julday,jde,jdb,jdstp,iye,iyear2,ime,ide,ihe,imine,isece
     .        , nics,ksat,lsat,nsat,nepchs,itsat(maxsat)
     .        , ispeed,ji0,jil,jnow,jlast,nintrs,iendf,idoye
     .        , iut,iscrn,iprnt,iug,iflag,iy1,iy2,i,j,k
     .        , ioerr,len,rcpar,idoy

      real*8 te,tb,tstp,sece,rvec,doye,trun
     .     , ytrp,sdelt,yy,satics,tend

      character*1  lowerc,gnss
      character*4  end,icsnam(maxorb),g_time_flag
      character*5  precmod,nutmod,gravmod,frame,srpmod,eradmod,antradmod
      character*6  prog
      character*16 tname,gname,satnam(maxsat),buf16
      character*19 dattim
      character*80 head,prog_name

      dimension rvec(maxyt2),ytrp(5,2,maxyt2,maxsat),
     1          yy(10,maxytp,maxsat),satics(maxorb,maxsat)
c
c     common/orb/nsat,nintrs,sdelt,iendf,ji0,jil,iy1,iy2,jlast,jnow

      data end/'END '/
      data iut/16/,iug/17/,iscrn/6/,iprnt/0/ 

c
c     Get the program name calling gmake
      len = rcpar(0,prog_name)
c
c     Open the T-file

      call lowers(tname)
      open(unit=iut,file=tname,status='old', form='unformatted'
     .,iostat=ioerr)
      if(ioerr .ne. 0 ) then
        call report_stat('FATAL',prog_name,'orbits/gmake',tname,
     .  'Error opening T-file: ',0)
      endif      

c     Open the G-file

      gname= tname
      gname(1:1)='g'
      call lowers (gname)
      open( unit=iug,file=gname,status='unknown',form='formatted'
     .              , iostat=ioerr)
      if( ioerr.ne.0 ) then
        call report_stat('WARNING',prog_name,'orbits/gmake',gname,
     .     'Error, G-file name exists, overwriting. ',ioerr)    
        open(unit=iug,file=gname,status='old',form='formatted'
     .               ,iostat=ioerr)    
        if( ioerr.ne.0 ) 
     .    call report_stat('FATAL',prog_name,'orbits/gmake',gname,
     .        'Error opening G-file',ioerr)
       endif

c     Read the T-file header

   30  frame = 'INERT'
c      G-files must be inertial
       call thdred ( iut,iscrn,iprnt,nsat,gnss,itsat,satnam
     .             , jdb,tb,jdstp,tstp,sdelt,nepchs,jde,te
     .             , nics,satics,nintrs,icsnam
     .             , precmod,nutmod,gravmod,frame,srpmod
     .             , eradmod,antradmod )

c          If ICs not to be taken from header must interpolate
c**rwk:    To avoid problem with dynamically inconsistent positions and velocities
c          (or possible bug in rotating velocites in /orbits/rotcrd.f), always
c          interpolate even if the IC epoch is not to change

c            If interactive (TTOG or optionally BCTOT) display the epoch
c            and change if desired

c       compute the time after start for the requested epoch
        trun= (jde-jdb)*86400.d0 + (te-tb)
        if( (trun-5.d0*sdelt) .le.0. ) then
          call report_stat('FATAL',prog_name,'orbits/gmake',' ',
     .    'Tabular ephemeris doesn''t start early enough ',0)
        endif
        tend= (jdstp-jdb)*86400.d0 + (tstp-tb)
        if( (trun+5.d0*sdelt) .ge. tend ) then
          call report_stat('FATAL',prog_name,'orbits/gmake',' ',
     .    'Tabular ephemeris ends too early ',0)
        endif
c       setup interpolation indices
        ji0= 0
        jil= 0
        iy1= 0
        iy2= 0
        jlast= 0
        iendf= 0
        ispeed=0
        do ksat=1,nsat
          lsat=ksat
          call gsatel( 2,trun,iut,lsat,rvec,ytrp,yy,nsat,sdelt
     .               , nintrs,ji0,jil,iy1,iy2,jlast,jnow,iendf,nepchs)
          do i=1,6
            satics(i,ksat)= rvec(i)
          enddo
        enddo
c     end interpolation code

c**     if not interactive, accept IC epoch and values from T-file header
c**     no, accept epoch but not values
c**      endif

c        Write the output G-file

c*      call mdyjul( ime,ide,iye,jde )
      call dayjul( jde,iye,idoy )
      call monday( idoy,ime,ide,iye )
      call daynum(ime,ide,iye,doye)
      idoye=idint(doye)
      call ds2hms( iye,ide,te,ihe,imine,sece )
      isece = sece
      g_time_flag = 'GPST' 
      iyear2 = mod(iye,100)  

      write(iug,'(i4,1x,i3,1x,i2,1x,i2,1x,i2,18x,a4,7(1x,a5))')
     .      iye,idoye,ihe,imine,isece,g_time_flag
     .     ,frame,precmod,srpmod,nutmod,gravmod,eradmod,antradmod
      write(iug,'(i2,15(1x,a4))') nics,(icsnam(i),i=1,nics)
      head = ' '
      head(1:48) = 'G-file generated from '//tname//' by '//prog
      call runtim( dattim )
      head(49:69) = '  '//dattim
      write(iug,'(a80)') head
      write(iug,'(a4)') end
      do 70 j=1,nsat
        write(iug,'(a16)') satnam(j)
        if( nics.eq.6 ) then
c       this is the default number of ic's for a g-file
           nics= 9
           do k=1,nsat
             satics(7,k)= 1.d0
             satics(8,k)= 0.d0
             satics(9,k)= 0.d0
           enddo
        endif
        do i=1,nics
          write(iug,'(d20.14)') satics(i,j)
        enddo
   70 continue
      write(iug,'(a4)') end
      close (iug)
      call report_stat('STATUS',prog_name,'orbits/gmake',gname
     .,'Successfully wrote G-file: ',0)

      return
      end
