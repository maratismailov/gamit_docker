      Subroutine GET_SVCOORDS

*     Read the t-file to get inertial coordinates at the epochs of the y-file
*     and transform these to Earth-fixed coordinates

      implicit none

      include '../includes/dimpar.h'
      include 'orbex.h'
            
* T-file variables and interpolation indices
      integer*4 ntsat,ntinter,ntepoch,itsat(maxsat)
     .        , jdb,jdf,nepcht,jde,nics,nintrs
      real*8 sdelt,tb,tf,te,tsatic(maxorb,maxsat),satarg,svec(6)
     .     , ytrp(5,2,maxyt2,maxsat),yy(10,maxytp,maxsat)
      character*1 tgnss   
      character*4 icsnam(maxorb)  
      character*5 precmod,nutmod,gravmod,frame,srpmod,eradmod,antradmod
      character*16 satnam(maxsat)
      integer*4 ji0,jil,iy1,iy2,jlast,jnow,iendf
                      
* Input to t-file read for printout -- not used in call or thdred
      integer*4 iscrn,iprnt                                      
     
* Quantities for rotating from inertial to Earth-fixed
      integer*4 iut1pol
      real*8 tdtgpst,rot(3,3),rotdot(3,3),sidtm,xpole,ypole
     .     , svec_out(6)
                  
* Local
      integer*4 i,isat,iepc 
      character*80 message 

* Function
      real*8 timdif

* Debug flag
      logical debug/.true./

* TDT-GPST for inertial-to-Efix rotations
      tdtgpst = 32.184d0 + 19.d0
             
* Read the t-file header to determine start/stop times of tfile
            
      iscrn = 6               
      iprnt = 0
      frame = 'UNKWN'
      call thdred ( lut,iscrn,iprnt,ntsat,tgnss,itsat,satnam
     .            , jdb,tb,jdf,tf,sdelt,nepcht,jde,te
     .            , nics,tsatic,nintrs, icsnam
     .            , precmod,nutmod,gravmod,frame,srpmod,eradmod
     .            , antradmod ) 

* Check for t-/y- files compatibility
           
      if(tgnss.ne.gnss) then 
        write(message,'(a,a1,a,a1)') 'T-file GNSS ',tgnss
     .    ,' differs from y-file ',gnss
        call report_stat('FATAL','Y2ORBEX','get_svcoords',tfiln
     .     ,message,0)
      endif
      if(nysat.ne.ntsat) then
        write(message,'(a,a1,a,a1)') '# SVs on t-file ',ntsat
     .    ,' differs from y-file ',nysat
        call report_stat('FATAL','Y2ORBEX','get_svcoords',tfiln
     .     ,message,0)                                     
      endif          
      do isat=1,ntsat
        if(itsat(isat).ne.iprn(i)) then
           write(message,'(a,i2,a,i2,a,i2)') 'For isat ',isat
     .    ,' PRN on t-file ',itsat(isat),' differs from PRN on y-file '
     .    , iprn(isat)
        endif
      enddo           

* Set the default ORBEX day from the t-file IC epoch 
*  (override in y2orbex if the start and stop times are input)

      obx_jdstart = jde

* MOD TAH 200505: Read the sestbl. to get the correct value
*     Old code had no value set at all.  
      call get_iut1pol( iut1pol ) 
                             
* Loop over all epochs of the y-file and interpolate the t-file coordinates

      do iepc = 1,nyepoch
        satarg= timdif(jdy(iepc),ty(iepc),jdb,tb)   
        do isat = 1,ntsat
cd          if(debug) print*,' satarg isat sdelt yy11 '
cd     .       ,satarg,isat,sdelt,yy(1,1,1)
          call gsatel( 1,satarg,lut,itsat,svec,ytrp,yy
     .               , ntsat,sdelt,nintrs,ji0,jil,iy1,iy2,jlast
     .               , jnow,iendf,nepcht)
cd          if(debug) print *,'From gsatel svec: ',svec
*         Rotate to Earth-fixed 
          call rotsnp( -1,jdy(iepc),ty(iepc),tdtgpst,iut1pol
     .               , frame,precmod,luut1,lupole,lunut
     .               , rot,rotdot,sidtm,xpole,ypole )
          call rotcrd( rot,rotdot,svec,svec_out,0,1 ) 
          do i=1,3
            satcrd(i,isat,iepc) = svec_out(i)
          enddo             
          if(debug.and.iepc.eq.271.and.isat.eq.1) then 
            write(*,*) 'GET_SVCOORDS epoch isat jdy ty  '
     .                         , iepc,isat,jdy(iepc),ty(iepc)
            print *,'  rot ',rot
            print *,'  svec svec_out ',svec,svec_out
          endif 
        enddo
      enddo                                  
      return
      end


 
      
