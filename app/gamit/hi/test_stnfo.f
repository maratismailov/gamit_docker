c     Program to test new rstnfo for new-style station.info
c     R. King 29 April 2001
c     Mod Feb 2010 for changed argument list for rstnfo
c     Mod Sep 2010 for changed argument list for rstnfo

      implicit none


c       station.info variables
                                
      character*4 sitcod
      character*5 htcod,radome  
      character*6 rcvcod,antcod     
      character*16 sname  
      character*20 rcvers,rcvrsn,antsn  
      real*4 swver
      real*8 anth,antn,ante 
      real*8 antdaz  ! Antenna aligment from True N (deg).
      integer*4 isessn,start(5),stop(5),span

c** not used
c      character*15 anttyp
c      character*20 rctype 

                      
c        other variables

      integer*4 luold,lunew,iyr,idoy,isod
     .        , icall,istart(3),istop(3),ioerr  

      data luold/10/,lunew/11/

c Open the files

      open(unit=luold,file='station.info.old',status='old',iostat=ioerr) 
      open(unit=lunew,file='station.info.new',status='old',iostat=ioerr)
       
c Set the test entries
                        
      sitcod = '2353'
      iyr = 1990   
      idoy = 337
      isod = 86400  
      isessn = 1
        
c Read the old file and print the results
c** no longer works: rstnfo obsolete  rwk 080509

c      prject = '    '
c      orbit =  '    '   
c      icall = 3
c      call rstnf1( luold,prject,orbit ) 
c      call rstnfo( luold, icall, prject, orbit
c     .            , sitcod, iyr, idoy, isessn, isod
c     .            , sname, anth, antn, ante, rcvcod, antcod
c     .            , htcod, swver, istart, istop )
c      print *,'Old:'
c      print *,' prject orbit trkcod sitcod ',prject,orbit,trkcod,sitcod
c      print *,' sname ',sname    
c      print *,' anth antn ante htcod ',anth,antn,ante,htcod
c      print *,' rcvcod swver antcod ',rcvcod,swver,antcod
c      print *,' istart istop ',istart,istop

c Read the new file and print the results
* MOD TAH 200203: Added AntDAZ to list of values from station.info
                         
      call rstnfo(lunew, sitcod, iyr, idoy, isod
     .                 , span, sname, anth, antn, ante, antdaz
     .                 , rcvcod, antcod
     .                 , htcod, radome, swver, rcvers, rcvrsn, antsn
     .                 , start, stop )   
      print *,'New:'
      print *,' sitcod ',sitcod
      print *,' span sname ',span,sname    
      print *,' anth antn ante htcod ',anth,antn,ante,htcod
      print *,' AntDAZ ',antdaz
      print *,' rcvcod swver antcod ',rcvcod,swver,antcod
      print *,' rcvers rcvrsn antsn ',rcvers,rcvrsn,antsn
      print *,' istart istop ',istart,istop

      stop
      end
