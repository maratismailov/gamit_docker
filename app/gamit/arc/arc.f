Copyright (c) Massachusetts Institute of Technology and the University of
California at San Diego, 1995. All rights reserved.

      Program ARC

      implicit none

      include '../includes/dimpar.h' 
      include '../includes/units.h'  
      include '../includes/global.h'
      include '../includes/arc.h'           


      logical found,moveold,fcheck,srpmod_mismatch

      integer*4 ngsats,ngics,isat,norder,ioerr,iun,ivel,i,j
                
      real*8 gfsatics(maxorb,maxsat),fjd
                         
      character*4  gftime_type
      character*5  gfframe,gfprecmod,gfnutmod,gfgravmod
     .            ,gfsrpmod,gferadmod,gfantradmod
      character*7 stat
      character*16 gfsatnam(maxsat)
      character*22 buff22
      character*76 buf76
      character*40 versn
      character*80 junk,buff80
      character*256 message                  

      logical debug/.false./

c   File definitions
c
c     unit number    file type
c
c        6           screen
c        8           iarh    (archive print file)
c        9           iyawtmp (temporary output ascii yaw file created by arc)
c       10           iyaw    (ascii yaw file output by arc, read by yawtab)  
c       11           idebug  (debug output, controlled by idbprt )
c       14           iyawtab (binary yaw table output from yawtab and read by arc
c                             only for UCL radiation pressure models; if not availabe, set iyawtab=0)
c       17           r or x  (for getting obs epochs - not implemented)
c       18           t       (output ephemeris)
c       29           g       (input cartesian initial conditions)
c       30           ut1     (for rotating earth for gravity field)
c       31           nut     (nutation, for rotating earth)
c       32           pole    (pole position, for rotating earth)
c       33           ocean tides  
c       34           lunar ephemeris (old code)
c       35           solar ephemeris (old code)
c       36           planetary ephemeris (new code)
c       40+          scratch files, one per sat  
c       69           opened in lib/svnav_read, change unit number eventually to avoid conflicts
c       99           UCL SRP grid files - names set in filopn, opened in rdsrpgrd 
c       5678         UCL SRP fourier series files - names assigned and opened in srpfmod 
c         
                                     
c   Binary-code options for debug printout and models (hidden in cols 14-16 of afname line  
        
c      Bit #  Value
c        1     1   print at all epochs (at least time and total accel in sbfn)
c        2     2   print only during eclipses
c        3     4   print shadowing values 
c        4     8   print sun, moon, satellite positions
c        5    16   print radiation-pressure quantities        
c        6    32   ** temporary: switch on earth radiation quantities (use model name)
c        7    64   print SRP quantities
c        8   128   ** temporary: switch on antenna thrust (use model name)
c        9   256   print earth radiation and antenna thrust quantities 
c       10   512   ** temporary: switch on arcyaw  (use rad pressure model name)
c       11  1024   turn off solid-Earth tides
c       12  2048   turn off lunar eclipese
c       13  4096   turn off ocean tides 

c       Skip ARC if a previous step has failed

      if( fcheck('GAMIT.fatal') ) 
     .  call report_stat('FATAL','ARC','arc',' '
     .                  ,'GAMIT.fatal exists: ARC not executed',0)

c       Get the ARC version 

      call aversn(versn)
      found = .false.
                         


c  Read the batch file to get the SVs, file names, models, and integration controls

      call read_arc_batch
C     print *,'ARC lbody ',lbody            

            
c  Open the print file (arcout.DDD) 
   
      iarh  = 8  
      open( unit=iarh,file=afname,status='unknown'
     .    , form='formatted',iostat=ioerr)
      if ( ioerr .ne. 0 ) then
        call report_stat('FATAL','ARC','filopn','afname',
     .  'Error opening ARC print file',ioerr)
      endif
      write (iarh,'(/,1x,a,a40,//)')  'ARC Version ',versn


c  Set up the model parameters
                
      call set_model_controls


c  Open most of the files
  
      call filopn( versn )        


c  Read the ocean tide coefficients 
                             
      if( iotide.ne.0 ) call read_otides
                                                                        
c  Open and read the g-file

      write(iarh,'(a,a16)') ' Input ICs file: ',gfname
      call read_gfile( gfname,iug,jde,te,gfframe,gfprecmod, gfnutmod
     .               , gfgravmod,gfsrpmod,gferadmod,gfantradmod 
     .               , gftime_type 
     .               , ngsats,nics,icsnam,gfsatnam,gfsatics ) 
          
cd       print *,'ARC aft read_gfile ',gfname,iug,jde,te 
cd       print *,   gfframe,gfprecmod,gfnutmod,gfgravmod,gfsrpmod
cd     .   ,gferadmod,gfantradmod
cd       print *,ngsats,nics,icsnam
cd       write(*,'(35a16)') (gfsatnam(i),i=1,ngsats)
cd       do i=1,ngsats
cd         write(*,'(15d8.2)') (gfsatics(j,i),j=1,15)
cd       enddo
                              
                                 

c  Check the compatibility of the ARC frame and models with the g-file frame and models 
                  
      call check_gmodels( gfframe,gfprecmod,gfnutmod,gfgravmod
     .                  , gfsrpmod,gferadmod,gfantradmod,gftime_type
     .                  , srpmod_mismatch )

                                        
c  Open and read then planetary ephemeris file 
c     GPST or UTC ok for determining what records to read into storage
      fjd =  dfloat(jde) + te/86400.d0 
c     temorary: if nbody exsits use it, otherwise revert to the old code
      if( fcheck('nbody')) then 
c       file unit number          
        ibody = 36      
c       read Earth, Moon, Venus, and Jupiter
c** rwk 181205: lbody now read from batch file (default is 1) 
c        lbody = 1    
cd       Earth and Moon only for testing
cd        lbody = 0 
c       don't need velocities
        ivel = 0 
        call ephred(ibody,fjd,lbody,ivel) 
      endif 

c  Write the control and models to the print file    

      call write_arcouthd
                      
c-----------------------------------------------------------------------------
c  Loop through all the selected satellites, incrementing the scratch file #
                            
      yaw_entries = 0
      isunit=39
      do  isat=1,nsats   
        isunit=isunit+1        

c          Initialize eclipse print variables
        first = .true.
        pjdlast = 0.d0
        lamprt = .false.
        neclipse = 0   

c         Open the scratch files 
        call oscrat(isunit,jdb)      
                       

c         Get satellite characteristics and ICs   
cd        print *,'ngsats gfsatnam ',ngsats,(gfsatnam(i),i=1,ngsats)
cd        print *,'nsats  satnam   ',nsats, (satnam(i),i=1,nsats) 
        do i=1,ngsats
          if( gfsatnam(i).eq.satnam(isat) ) then  
            do j=1,6
              satics(j) = gfsatics(j,i)
            enddo
            if( srpmod_mismatch ) then
               satics(7) = 1.d0
               do j=8,nics 
                 satics(j) = 0.d0
               enddo
             else
               do j=1,nics
                 satics(j) = gfsatics(j,i)
               enddo
            endif  
            call get_sat_info( isat ) 
            sname = satnam(isat) 
          endif
         enddo
                    
c         Initialize integration control parameters. Mod TAH 030304: Passed
c         IC epoch (PEP MJD), jde, into init for time variable gravity.
cd        print *,'DEBUG before init jde te jdb tb ',jde,te,jdb,tb
        call init( isat )  
cd        print *,'DEBUG past init jdb tb ',jdb,tb
                    

c         Write T-file header                     
        call wrthed  
                    
c          Open and read headers of lunar and solar ephemeris files and load the
c          Everett interpolation coefficients and factors (old code only)

        if( .not.fcheck('nbody') ) then
cd          print *,'ARC ',frame,precmod,gravmod,nutmod,speopmod 
          call ephdrd(fjd)
          call evrtcf   
        endif          

c         Compute the binomial coefficients required by calcof
        norder=12
        call binoml(norder)

c         Compute the coefficients required for the Adams-Moulton integrator
        call calcof        
                         
c         Perform the numerical integration
        write (message,300) isat,sname,iprn
  300   format (1x,' Integrating satellite ',i2,2x,a16,2x,'PRN ',i2)
         call report_stat('STATUS','ARC','arc',' ',message,0)
        call adam  
cx        print *,'stop after adam '
cx        stop 

                        
c         Write the eclipse summary for this satellite
        call eclout
c         Rewind and close the files

c       rewind the scratch unit in preparation for writing the T-file
        rewind isunit
c       close the lunar and solar ephemeris files
        close (unit=ilun)
        close (unit=isun)
c     --end of loop on satellites
      enddo   
c-----------------------------------------------------------------------------------

c  Rewind the yaw file and write the number of yaw entries on the 2nd line of header
      rewind(iyawtmp)

c       Read 2 lines of header
      read(iyawtmp,'(a80)') buff80
      write(iyaw,'(a80)') buff80
      read(iyawtmp,'(a22)') buff22
c       Add the number of yaw entries to the start of the 2nd line
      write(buff22(1:4),'(i4)') yaw_entries
      write(iyaw,'(a22,a41)') buff22
     .      ,' - explanation of file in read_yaw.f'
c       Now read in and write out the rest of the file
      do while (ioerr.eq.0)
          read(iyawtmp,'(a80)',iostat=ioerr) junk
          write(10,'(a80)')junk
      enddo
c        Close the yaw files
      close(iyawtmp)
      close(iyaw)       
      write(iarh,'(a,a16)')' Output yaw file: ',yfname
      

c Write the T-file from the scratch files
      call arcmrg  
     
c Finished!

      call report_stat('STATUS','ARC','arc',tfname,
     .'Normal stop in ARC',ioerr)
      stop
      end
