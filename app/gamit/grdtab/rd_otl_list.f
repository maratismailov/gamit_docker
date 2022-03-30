      Subroutine rd_otl_list

c     Read the header of a Scherneck-style station file for ocean tidal loading
c     Written by R. King from old utils/octtab  
c     24 July 2004

      implicit none  

      include '../includes/grdtab.h'
       
c Input from grdtab.h
c   luprnt     i*4 unit number of print file ('grdtab.out')
c   luotll     i*4 unit number for station list file ('olt.list')   
c   maxlsit    i*4 maximum number of sites allowed on station list file
 

c  Output to grdtab.h      
c   otllmod   c*8  model name    
c   otllsrc   c*3  institutional source of tides (OSO or NAO),sets order and #
c   notll     i*4  number of tidal components read from list file 
c                  (11 or OSO, 54 for NAO)
c Note: Unlike the Vienna mapping function list files, the Scherneck OTL files
c       do not have a list of sites at the top, nor do they have on the station
c       entries themselves unambiguous 4-character ids correspondingto IGS 
c       sites.  Hence, the logic to match site names is not used here, although
c       the slots for the stations and coordinates are putinto grdtab.h.

c  Local

      integer*4 ioerr

c   File format : Lines begin with $$ or two blanks   

c        For Scherneck (OSO) file:
c       $$ Ocean loading displacements
c       $$  Comments: variable number of lines beginning with $$
c       $$  
c        $$ COLUMN ORDER:  M2  S2  N2  K2  K1  O1  P1  Q1  MF  MM SSA
c        For Matsumoto (NAO) file:  
c       $$ COLUMN ORDER:  M2 S2 K1 O1 N2 P1 K2 Q1 M1 J1 OO1 2N2 Mu2 Nu2 L2 T2
c       $$                2Q1 Sigma1 Q1' Rho1 O1' Tau1 M1' Kai1 Pi1 P1' K1' K1'
c       $$                Psi1 Phi1 Theta1 J1' OO1' Eps.2 Mu2' N2' Nu2' Gamma2 
c       $$                Alpha2 M2' Beta2 Delta2 LAmbda2 S2' R2 K2' Zeta2 Eta2
c       $$                Eta2' Mtm Mf Mm Ssa Sa 
                                     

          
c   Open the station list file
                         
      open(luotll,file='otl.list',iostat=ioerr,status='old')  
      if (ioerr.eq. 0 ) then    
         call report_stat('STATUS','GRDTAB','rd_otl_list','otl.list'
     .                   ,' Opened station ocean tide table',0)  
      else
         call report_stat('FATAL','GRDTAB','rd_otl_list'
     .            ,'otl.list','Error opening ocean tide table : ',ioerr)   
      endif  
                  


c   Read the first line to get the GAMIT 8-character ID of the model

      read(luotll,'(a8)',iostat=ioerr) otllmod   
      if(ioerr.ne.0) call report_stat('FATAL','GRDTAB','get_otl_list'
     .  ,' ','Error reading model from 1st line of station table',ioerr)
 
c   Determine the source of the file (OSO or NAO)
                               
c     Should be evident from the model name
      if( otllmod(1:3).eq.'NAO' ) then 
         otllsrc = 'NAO'   
         notll = 54
      else
         otllsrc = 'OSO'  
         notll = 11
      endif

      return

c** This code no longer needed    
c      otllsrc = ' ' 
c      otllmod = ' '
c      found_source = .false.
c      found_model = .false.
c      do while( .not.found_source .or. .not.found_model )
c        read(luotll,'(a60)',iostat=ioerr) line(1:60)  
c        if( ioerr.eq.0 ) then 
c          if( index(line,'COLUMN ORDER').ne.0 ) then
c             if( index(line,'Mu2').gt. 0 ) then
cc               call report_stat('STATUS','GRDTAB','rd_otl_list'
cc     .               ,'otl.list'
cc     .              ,' Using Matsumoto NAO model for stations in ',0)  
c               otllsrc = 'NAO' 
c               notll = 54  
c               otllmod = 'NAO     '  
c               found_source = .true.
c               found_model = .true.
c             else
cc               call report_stat('STATUS','GRDTAB','rd_otl_list'
cc     .             ,'otl.list'
cc     .          ,' Using Scherneck OSO model for stations in ',0)
c               otllsrc = 'OSO'  
c               notll = 11      
c               found_source = .true.
c             endif
c          endif    
c          if( index(line,'Ocean tide model:').ne.0 ) then  
c           newer versions of the OSO files have the model listed
c            read(line(22:28),'(a7)',iostat=ioerr) otllmod(1:7)
c            if( ioerr.ne.0 )  call report_stat('FATAL','GRDTAB'
c     .         ,'rd_otl_list','otl.list'
c     .         ,'Error reading model from OSO header',ioerr)
c            found_model=.true.
c           endif
c        elseif( ioerr.eq.-1 ) then 
c          call report_stat('FATAL','GRDTAB','rd_otl_list'
c     .    ,'otl.list','EOF before finding COLUMN ORDER on file',ioerr)
c        else
c          call report_stat('FATAL','GRDTAB','rd_otl_list'
c     .     ,'otl.list'
c     .    ,'Error reading station file header for COLUMN ORDER',ioerr)  
c        endif 
c      enddo   
c      return
 
      end
                                   
