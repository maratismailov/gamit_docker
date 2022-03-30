      Subroutine check_gmodels( gfframe,gfprecmod,gfnutmod,gfgravmod
     .                        , gfsrpmod,gferadmod,gfantradmod
     .                        , gftime_type,srpmod_mismatch ) 

c     Compare the input models with those on the g-file.  If the frame is
c     different, rotate the g-file.  If the models are different, warn the
c     user and trust that the fitting (sh_sp3fit) or an iteration (gamit/globk)
c     will result in a consistent orbit.   Code moved from arc.f R. King 21 May 2012

      implicit none

      include '../includes/dimpar.h' 
      include '../includes/units.h'
      include '../includes/global.h'
      include '../includes/arc.h'
                                     
      integer*4 ioerr,i
        
      character*4 srpnam(13),gftime_type
      character*5  gfframe,gfprecmod,gfnutmod,gfgravmod
     .            ,gfsrpmod,gferadmod,gfantradmod
      character*80 nut_header 
      character*256 message 

      real*8 dt

      logical srpmod_mismatch,moveold,fcheck
              
c functions
      integer mchkey
      real*8 taiutc 

      if( frame.ne.gfframe.or.precmod.ne.gfprecmod ) then
         moveold = .true.
         call rot_gfile(gfname,frame,gfframe,precmod,gfprecmod
     .                 ,moveold,0,iug,'ARC ',te,jde,tdtoff
     .                 ,inut,iut1,ipole )
      endif 
c       check whether the nutation model and gravity model are compatible
      if( .not.fcheck('nbody') ) then 
        read(inut,'(a)',iostat=ioerr) nut_header
        if(ioerr.ne.0) call report_stat('FATAL','ARC','check_gmodels'
     .     ,' ','Error reading first line of nutation file to get model'
     .     , ioerr)
        rewind(inut)
        if(mchkey(nut_header,'IAU20',80,5).gt.0) then
          nutmod='IAU00'
        else
          nutmod='IAU80'
        endif    
      else            
c       nutations to be computed by MHB_2000
        if (precmod . eq. 'IAU76') then
          nutmod='IAU00'                      
        else 
          nutmod = precmod
        endif
      endif 
      if( nutmod.ne.gfnutmod ) then   
        write(message,'(a,a5,a,a5,a)') 'Available nutation model ('
     .     ,nutmod,') differs from g-file (',gfnutmod,')'  
        call report_stat('FATAL','ARC','check_gmodels',gfname,message,0)
      endif
      if( gravmod.ne.gfgravmod ) then 
        write(message,'(a,a5,a,a5,a)') 'Requested gravity model ('
     .    ,gravmod,') differs from g-file (',gfgravmod,')'
        call report_stat('WARNING','ARC','check_gmodels',gfname,message
     .       ,0)
      endif 
      if( srpmod.ne.gfsrpmod ) then
        write(message,'(a,a5,a,a5,a)') 
     .    'Requested radiation-pressure model (',srpmod
     .     ,') differs from g-file (',gfsrpmod,')'
        call report_stat('WARNING','ARC','check_gmodels',gfname,message
     .     ,0)
        call report_stat('WARNING','ARC','check_gmodels',' '
     .       ,'Use requested model, set RAD(1)=1.0, others = 0.0',0)
        call assign_srpnames(srpmod,nics,srpnam)
        do i=7,nics
          icsnam(i)=srpnam(i-6)
        enddo 
        srpmod_mismatch = .true.
      else
        srpmod_mismatch = .false.  
      endif
      if( eradmod.ne.gferadmod ) then
        write(message,'(a,a5,a,a5,a)') 
     .    'Requested Earth radiation-pressure model (',eradmod
     .     ,') differs from g-file (',gferadmod,')'
        call report_stat('WARNING','ARC','check_gmodels',gfname,message
     .     ,0)
      endif
      if( antradmod.ne.gfantradmod ) then
        write(message,'(a,a5,a,a5,a)')   
     .    'Requested antenna thurst model (',antradmod
     .     ,') differs from g-file (',gfantradmod,')'
        call report_stat('WARNING','ARC','check_gmodels',gfname,message
     .    ,0)
      endif
      if( time_type.eq.'GPST'.and. gftime_type.eq.'UTC ' ) then
         dt = taiutc(jde) - 19.d0
         call timinc( jde,te,dt )
      elseif( time_type.eq.'UTC '.and. gftime_type.eq.'GPST' ) then
         dt = 19.d0 - taiutc(jde)
         call timinc( jde,te,dt )
      endif
   
      return
      end    

