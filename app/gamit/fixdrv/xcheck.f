      subroutine xcheck( xfile, lxfil, lstnfo, exprmt, nsat, ixsat
     .                 , iy, idoy, gnss )

c     Determine the X-file format by reading the header or guessing
c     from the type of experiment ('Old' if static, 'Expanded' if kinematic)
c     R. King   9 April 1992
c MOD TAH 200203: Added AntDAZ token for antenna  Alignment from True N

      implicit none

      include '../includes/dimpar.h'

      logical fcheck, old_stinf

      integer*4 lxfil,lstnfo,lscrn,nsat, ixsat(maxsat)
     .        , lprnt,nchan,ischan(maxsat),nepoch,inter,mtime
     .        , iy,im,id,ihr,min,ndat
     .        , lambda(maxsat,maxdat),dattyp(maxdat)
     .        , latd,latm,lond,lonm
     .        , ntext,ircint,span
     .        , icall,idoy,isessn,istarts(5),istops(5),isod,i
                                   
      real*4 swver
      real*8 sec,offarp(3),seclat,seclon
     .      , height, anthgt,offstn,offste
      real*8 antdaz  ! Antenna aligment from True N (deg).

       real*4 swverx

       parameter( lscrn=6, lprnt=8 )

       character* 1  latflag,lonflag,gnss 
       character* 3  rcvrsw,rxobtyp(maxdat)
       character* 4  sitcod
       character* 5  hgtcod,exprmt,radome
       character* 6  rcvcod,antcod
       character*16  sitnam,xfile,satnam(maxsat)
       character*20  rctype,rcvnum,anttyp,antnum
       character*80  text(maxtxt)
       character*256 message

                              
           if( fcheck(xfile) ) then

c            There is an X-file to check
             open( lprnt, status='scratch' )
             open( lxfil, file=xfile, status='old' )
             call  xhdred( lxfil, lprnt, lscrn, nepoch, inter, ircint,
     .                 isessn, mtime, iy, im,id, ihr, min, sec,
     .                 nchan, ischan, satnam,
     .                 ndat, dattyp, rxobtyp, lambda,
     .                 offarp, sitnam, rcvrsw, swverx, antcod,
     .                 rctype,rcvnum,anttyp,antnum,
     .                 latflag, latd, latm, seclat,
     .                 lonflag, lond, lonm, seclon, height,
     .                 ntext, text,  gnss )  
c**  rwk 070104: What was the purpose of this  rwk :
c             if( iy.gt.1900 ) iy = iy - 1900     
             close( lxfil )
             close( lprnt )
            do i=1,nsat
              if( ischan(i).ne.ixsat(i) ) then
                 write(message,'(a,i2,a,a16,a)')
     .                  'PRN ',ischan(i),' on file ',xfile
     .                 ,' does not match SV list'
                 call report_stat('FATAL','FIXDRV','xcheck',' '
     .                      ,message,0)
               endif
             enddo

           else if ( fcheck('station.info') ) then

              sitcod = xfile(2:5)  
              span = 0
* MOD TAH 200203: Added AntDAZ to list of values from station.info
              call rstnfo( lstnfo,sitcod,iy,idoy,isod
     .                   , span, sitnam, anthgt, offstn, offste, antdaz
     .                   , rcvcod, antcod
     .                   , hgtcod,radome,swver,rctype,rcvnum,antnum
     .                   , istarts, istops )  
             call report_stat('WARNING','FIXDRV','xcheck',' '
     .        , 'No X-file to check format for MODEL batch file',0)
           endif
         return
         end
