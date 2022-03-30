c     Program to convert receiver or satellite phase-center-variation tables between
c     IGS ANTEX, NGS ant_info, and GAMIT antmod.dat formats.  Executation will convert
c     either an entire table or a specified antenna entry.  The program will not yet 
c     handle azimuthal variations.  In the main program, units of all corrections are mm.

c     R. King 11 April 2003 
c     R. King 31 May 2017: Modified to add atxfrq (and 'found' variables) to calling 
c         argument for read_antex.  HOWEVER, atxfrq is not defined here, so this routine 
c         probably doesn't work anymore. That's ok since it probably has no use.
c     T. Herring 12 May 2020: Modified atxfrq to be 4 characters and 3 elements to be
c         conistent with read_antex call.  

      
      implicit none
            
      include '../includes/dimpar.h'

      integer*4 nblen,luin,luout,ioerr,iclarg,iarg,ncomments
     .        , jd,i,j
             
      real*4 antex_vers

      real*8 offsetl1(3),offsetl2(3),elvtabl1(maxel),elvtabl2(maxel)
     .     , tablel1(maxel,maxaz),tablel2(maxel,maxaz)
     .     , dazi,zen1,zen2,dzen,offtmp
                                           
      character*1 ans,pcvtype,buf1,gnss
* MOD TAH 200512: Changed atxfrq declaration. 3 entries, 4 characters
      character*4 atxfrq(3)
      character*5 intype,outtype,upperc,radome                            
      character*6 gamit_code,buf6
      character*10 sinex_code
      character*16 amodel
      character*20 infile,outfile,antenna,anttype,antsn
     .           , refant,buf20   
      character*60 comments(maxtxt)   
      character*80 line
      character*256 message
         
      logical fcheck,finished,one_ant,found_ant,found_f1,found_f2
     .       , antidline,ant_unknown
   

c** Initialize the arrays and constants
      
      infile = ' '
      outfile = ' ' 
      intype = ' ' 
      outtype = ' ' 
      pcvtype = ' ' 
      refant = ' ' 
      radome = ' '
      luin = 1
      luout = 2   
      antex_vers = 1.0  
      dazi = 0.d0 
      do i=1,3
       offsetl1(i) = 0.d0
       offsetl2(i) = 0.d0
      enddo
      do i=1,maxel
        elvtabl1(i) = 0.d0   
        elvtabl2(i) = 0.d0
        do j=1,maxaz  
          tablel1(i,j) = 0.d0
          tablel2(i,j) = 0.d0
        enddo
      enddo

c** Read the command line

      iarg = iclarg(1,infile)
      if ( iarg.le.0 ) then 
         write(*,100)
  100    format(/,'CONVERT_ANTPCV: Program to convert antenna tables'
     . ,/,'Runstring: '
     . ,',  convert_antpcv  <in file>  <in type>  <out file> '
     .       ,' <out type>  <antenna> <pcvtype>'
     . ,/,'  where <in file>  input file to be converted'
     . ,/,'        <in type>  format of input file (antext ngs gamit)'
     . ,/,'        <out file> output file'
     . ,/,'        <out type> format of output file (antex nsg gamit)'
     . ,/,'        <antenna>  specific antenna to be converted '
     . ,/,'                   name must match colums 1-20, 1-15, or 1-6'
     . ,/,'                   if omitted, do entire table'
     . ,/,'        <pcvtype>  A or R  absolute (default) or relative',/)
        stop
      endif
      iarg = iclarg(2,intype)
      if( iarg.le.0 .or. nblen(intype).gt.6 ) then   
         call report_stat('FATAL','CONVERT_ANTPCV','convert_antpcv'
     .      ,' ','Missing or incorrect input file type',0)
      else
         call lowers(intype)
      endif
      iarg = iclarg(3,outfile)
      if( iarg.le.0 ) call report_stat('FATAL','CONVERT_ANTPCV'
     .    ,'convert_antpcv',' ','Missing output file name',0)
      iarg = iclarg(4,outtype)
      if( iarg.le.0 .or. nblen(outtype).gt.6 )  then
        call report_stat('FATAL','CONVERT_ANTPCV','convert_antpcv'
     .     ,' ','Missing or incorrect output file type',0)
      else
        call lowers(outtype)
      endif               
      one_ant = .false.
      iarg = iclarg(5,antenna)    
c*      print *,'antenna iarg intype(1:3) ',antenna,iarg,intype(1:3)
      if( iarg.gt.0) then 
         one_ant = .true.
         if ( intype.eq.'antex' ) then
           anttype  = antenna                  
         elseif( intype(1:3).eq.'ngs' ) then
           anttype = antenna 
c*           print *,'anttype ',anttype
         elseif( intype.eq.'gamit' ) then
           gamit_code = antenna(1:6)
         else
           call report_stat('FATAL','CONVERT_ANTPCV','convert_antpcv'
     .        ,' ','Invalid input file type',0)
         endif
       endif 
       iarg = iclarg(6,pcvtype)
       if( iarg.gt.0.and.pcvtype.eq.'R') then
         refant = 'AOAD/M_T'
       elseif( iarg.eq.0 ) then
         pcvtype = 'A'
       endif


c** Open the files

      open( unit=luin,file=infile,status='old',iostat=ioerr)
      if( ioerr.ne.0 ) then
        write(message,'(2a)') 'Error opening input file ',infile
        call report_stat('FATAL','CONVERT_ANTPCV','convert_antpcv'
     .     ,' ',message,ioerr)
      else
         write(*,'(4a)') 'Opened input file ',infile,' Type='
     .       ,upperc(intype)     
      endif  
      if( fcheck(outfile) ) then
        write(*,'(a)') 'Output file exists, overwrite? (y/n)'
        read(*,'(a)') ans
        if( ans.eq.'y' ) then  
          open( unit=luout,file=outfile,status='old',iostat=ioerr) 
        else
          write(*,'(a)') 'Stop. File not overwritten'
        endif 
      else
        open( unit=luout,file=outfile,status='new',iostat=ioerr) 
      endif 
      if( ioerr.ne.0 ) then 
         write(message,'(2a)') 'Error opening output file ',outfile
         call report_stat('FATAL','CONVERT_ANTPCV','convert_antpcv'
     .       ,' ',message,ioerr)
      endif
                
c** Read the header lines of the input file 

      if( intype.eq.'antex' ) then      
        call read_antex_head( luin,antex_vers,pcvtype,refant
     .                      , sinex_code,ncomments,comments)

      elseif( intype.eq.'ngs  ' ) then  
c      not sure whether file id and blank line is present, so read
c      to first prescribed line and then read the  format rigid, simply skip the header lines 
        antidline =.false.
        do while (.not.antidline )
          read(luin,'(a)',iostat=ioerr)  line
          if( ioerr.ne.0 ) then 
           write(message,'(a)') 'Error finding start of NGS file header'     

           call report_stat('FATAL','CONVERT_ANTPCV','convert_antpcv'
     .       ,' ',message,ioerr)
          elseif( line(1:10).eq.'ANTENNA ID' ) then
            antidline =.true.
          endif
        enddo
        do i=1,8      
          read(luin,'(a)',iostat=ioerr)
          if( ioerr.ne.0 ) then 
            write(message,'(a)') 'Error reading header of NGS file'   
            call report_stat('FATAL','CONVERT_ANTPCV','convert_antpcv'
     .       ,' ',message,ioerr)
          endif
        enddo
        zen1 = 0. 
        zen2 = 90.
        dzen = 5.      

      elseif( intype.eq.'gamit' ) then
c       format rigid, header lines are all comments
        zen1 = 0. 
        zen2 = 90.
        dzen = 5.      
c       assume for now that all GAMIT files are relative
        pcvtype = 'R'
        refant = 'AOAD/M_T'
      endif
          

c** Write the header lines of the output file

      if( outtype.eq.'antex' ) then
        write(luout,'(f8.1,12x,a1,39x,a20)')
     .    antex_vers,'G','ANTEX VERSION / SYST' 
        write(luout,'(a1,19x,a20,20x,a20)')
     .    pcvtype,refant,'PCV / TYPE / REFANT '
        write(luout,'(a27,a20,13x,a)') 
     .    ' Converted from input file ',infile,'COMMENT'
        write(luout,'(60x,a)') 'END OF HEADER' 

      elseif( outtype.eq.'ngs  ') then   
        call report_stat('FATAL','CONVERT_ANTPCV','convert_antpcv'
     .       ,' ','Output of NGS-format file not yet coded',0)
      elseif( outtype.eq.'gamit' ) then
        write(luout,'(a)')'# Antenna phase center models for GAMIT ' 
        write(luout,'(a,a)') '# File converted from ',infile 
        write(luout,'(a)') '#'  
      
      endif

c** Begin the loop over antennas
                          

      finished = .false.
      ant_unknown = .false.
      do while (.not.finished)    

        found_ant = .false.  
        ant_unknown = .false.  
        if( .not.one_ant ) then
c         for complete file copy, reset the codes each time through e
          anttype = ' ' 
          antsn = ' ' 
          gamit_code = ' ' 
        endif                 

        if( intype.eq.'antex' ) then 
c        this subroutine to be used in /lib later
c        if antenna type blank (here only) read one record each call
c        set gnss arbitrarily to 'G'since we're unlikely to use this routine
         gnss = 'G'
         call read_antex( luin,antex_vers,anttype,antsn,gnss,atxfrq,jd
     .                  , found_ant,found_f1,found_f2
     .                  , sinex_code,dazi,zen1,zen2,dzen
     .                  , offsetl1,offsetl2,elvtabl1,elvtabl2
     .                  , tablel1,tablel2 )   
c        keep ANTEX order (N E U or X Y Z) of offsets
         call read_rcvant(2,1,gamit_code,anttype,radome,buf6,buf20,buf1)
c        temporary assignment until a coversion table is constructed
         if( sinex_code(1:6).eq.'IGS_05') amodel(1:8) = 'I1_IGS05'  
         amodel(9:16) = '       '
        elseif( intype.eq.'ngs  ' ) then 
c*        print *,'Before READ_NGS anttype found_ant ',anttype,found_ant  
         call read_ngs( luin,anttype
     .                , found_ant,zen1,zen2,dzen
     .                , offsetl1,offsetl2,elvtabl1,elvtabl2 )  
c*         write(*,'(a,l1,1x,a20,3f8.1)') 
c*     .     'Aft READ_NGS found_ant anttype offsetl1 '
c*             ,found_ant,anttype,offsetl1
         dazi = 0.                   
         radome(1:4) = anttype(17:20)
         if( anttype(1:5).ne.'NONE ' )     
     .   call read_rcvant(2,1,gamit_code,anttype,radome,buf6,buf20,buf1)
c*          write(*,'(a,a5,1x,a20 )') 
c*     .      'After read_rcvant radome anttype ',radome,anttype
c        assume the NGS version is '4'
         amodel(1:8) = 'N4_NGS04' 
         amodel(9:16)= '      '     
c        GAMIT code returned blank if not found (warning issued)
         if( gamit_code.eq.'      ' ) ant_unknown = .true.   
        elseif( intype.eq.'gamit' ) then            
         call read_antmod( luin,gamit_code
     .                   , found_ant,amodel,dazi,zen1,zen2,dzen
     .                   , offsetl1,offsetl2,elvtabl1,elvtabl2
     .                   , tablel1,tablel2 )     
c        swap offsets to get ANTEX/NGS order (N E U) 
         offtmp = offsetl1(2)
         offsetl1(2) = offsetl1(3)
         offsetl1(3) = offsetl1(1)
         offsetl1(1) = offtmp
         offtmp = offsetl2(2)
         offsetl2(2) = offsetl2(3)
         offsetl2(3) = offsetl2(1)
         offsetl2(1) = offtmp             
         call read_rcvant(1,1,gamit_code,anttype,radome,buf6,buf20,buf1)
        endif   
        if ( found_ant .and. .not.ant_unknown ) then  
           if( anttype(1:5).eq.'NONE ') then
             continue
           elseif( outtype.eq.'antex' ) then  
c            tempoarary conversion for PCV model codes
             if( amodel(1:8).eq.'I1_IGS01') then
                sinex_code = 'IGS_01    '
             elseif( amodel(1:8).eq.'N4_NGS04') then
                sinex_code = 'NGS_04    '
             else
                sinex_code = 'UNKNOWN   '
             endif 
             call write_antex( luout,anttype,sinex_code
     .                       , dazi,zen1,zen2,dzen
     .                       , offsetl1,offsetl2,elvtabl1,elvtabl2
     .                       , tablel1,tablel2 )
           elseif( outtype.eq.'ngs  ' ) then   
              call report_stat('FATAL','CONVERT_ANTPCV','convert_antpcv'
     .           ,' ','Output of NGS-format file not yet coded',0)
           elseif( outtype.eq.'gamit' ) then 
c            ANTEX/NGS are N E U, but GAMIT expects U N E
             offtmp = offsetl1(1)
             offsetl1(1) = offsetl1(3)
             offsetl1(3) = offsetl1(2)
             offsetl1(2) = offtmp
             offtmp = offsetl2(1)
             offsetl2(1) = offsetl2(3)
             offsetl2(3) = offsetl2(2)
             offsetl2(2) = offtmp                              
c debug             print *,'amodel ',amodel
             if( amodel.eq.' ' ) amodel = 'CONVERTED      '  
c debug             print *,'amodel ',amodel
             call write_antmod(luout,gamit_code,amodel,dazi,zen1,zen2
     .                        ,dzen, offsetl1,offsetl2,elvtabl1,elvtabl2
     .                        ,tablel1,tablel2 ) 
           endif
        endif                                    
               
        if( one_ant ) then
          if( .not.found_ant)  then                               
            write(message,'(a,a20,a)') 'Requested antenna (',anttype
     .           , ') not found on input PCV file'
            call report_stat('FATAL','CONVERT_ANTPCV','convert_antpcv'
     .        ,' ',message,0) 
          else
            call report_stat('STATUS','CONVERT_ANTPCV','convert_antpcv'
     .        ,' ','Requested antenna found on input PCV file',0)
            finished = .true.
          endif
        else
          if( .not.found_ant ) then
            call report_stat('STATUS','CONVERT_ANTPCV'
     .             ,'convert_antpcv',' '
     .            ,'Complete copy, end reading input PCV file',0)
            finished = .true.
          endif
        endif
      
      enddo
      stop
      end

                  
c*************************************************************************

      Subroutine read_ngs( luin,ngs_code
     .                   , found,zen1,zen2,dzen
     .                   , offsetl1,offsetl2,elvtabl1,elvtabl2 )  

c     Input:
c       luin        logical unit number
c       ngs_code    20-character code of requested antenna; if blank read the next entry

c     Output: 
c       found       indicates requested antenna found
c       zen1        minimum zenith angle (max elev), usually = 0.
c       zen2        maxium zenith angle (min elev), usually 90.
c       dzen        increment of zenith angle of table, 5 deg for NGS files 
c       offsetl1(3) U N E constants for L1 phase center
c       offsetl2(3) U N E constants for L2 phase center
c       elvtabl1    array of elevation-dependent PCV corrections for L1 
c       elvtabl2    array of elevation-dependent PCV corrections for L2 
     

      implicit none
            
      include '../includes/dimpar.h'

      integer*4 luin,ioerr,i  
                              
      real*8 zen1,zen2,dzen,offsetl1(3),offsetl2(3)
     .     , elvtabl1(maxel),elvtabl2(maxel)
                    
      character*20 ngs_code,antenna
      character*256 line

      logical found,eof

c**   set default parameters for NGS ant_info files
  
      zen1 = 0.d0
      zen2 = 90.d0
      dzen = 5.d0        

c**   read the antenna entries for L1 and L2
                 
      found = .false.                    
      eof = .false.                               
      line = ' '            
c*      print *,'READ_NGS ngs_code found ',ngs_code,found
      do while (.not.found .and. .not.eof) 
        read(luin,'(a20)',iostat=ioerr) line(1:20)        
        if( ioerr.eq.-1 ) then 
           eof = .true. 
c           call report_stat('WARNING','CONVERT_ANTPCV'
c     .          ,'read_ngs',' ','EOF on NGS antenna file',ioerr)
        elseif (ioerr.ne.0 ) then 
           call report_stat('FATAL','CONVERT_ANTPCV'
     .          ,'read_ngs',' ','Error reading NGS antenna file',ioerr)
        else                    
          antenna = ' '                          
          read(line,'(a20)',iostat=ioerr) antenna      
c*          print *,'READ_NGS antenna found ',antenna,found
c         NGS has some non-standard names: convert them to IGS standards
          if( antenna(1:16).eq.'ASH701975.01A+GP') 
     .        antenna(1:16)='ASH701975.01Agp ' 
          if( antenna(1:16).eq.'ASH701975.01B+GP') 
     .        antenna(1:16)='ASH701975.01Bgp '  
          if( antenna(1:16).eq.'JPSMARANT_GGD   ') 
     .        antenna(1:16)='JNSMARANT_GGD   ' 
          if( antenna(1:16).eq.'NAVRT3010S      ') 
     .        antenna(1:16)='NAVCOM RT-3010S ' 
          if( antenna(1:16).eq.'NAVSF2040G      ') 
     .        antenna(1:16)='NAVCOM SF-2040G ' 
          if( antenna(1:16).eq.'SOK_RADIAN_IS   ') 
     .        antenna(1:16)='SOK RADIAN_IS   ' 
          if( antenna(1:16).eq.'SPP571908273+CR ')  
     .        antenna(1:16)='SPP571908273    '
                 
c          print *,'ioerr antenna(1:2) found ',ioerr,antenna(1:2),found  
c          print *,'      ngs_code(1:2) found ',ngs_code(1:2),found        

c         check for blank line at end
          if( ioerr.eq.-1 .or. antenna(1:10).eq.'          ' ) then
            eof = .true.
          elseif(antenna.eq.ngs_code .or. ngs_code(1:1).eq.' ') then
            found = .true.
            read(luin,*,iostat=ioerr) (offsetl1(i),i=1,3)   
            read(luin,*,iostat=ioerr) (elvtabl1(i),i=1,10) 
            read(luin,*,iostat=ioerr) (elvtabl1(i),i=11,19) 
            read(luin,*,iostat=ioerr) (offsetl2(i),i=1,3)  
            read(luin,*,iostat=ioerr) (elvtabl2(i),i=1,10)    
            read(luin,*,iostat=ioerr) (elvtabl2(i),i=11,19)    
            if( ioerr.ne.0 ) then   
               call report_stat('WARNING','CONVERT_ANTPCV'
     .          ,'read_ngs',' ','Error reading NGS table values',ioerr)
            endif
            ngs_code = antenna                
c **           print *,'ngs_code'
          endif
        endif 
      enddo

      return
      end

         
c*************************************************************************

      Subroutine write_antex( luout,anttype,sinex_code
     .                      , dazi,zen1,zen2,dzen
     .                      , offsetl1,offsetl2,elvtabl1,elvtabl2
     .                      , tablel1,tablel2 )

 
c     Write ANTEX records for a single antenna
     
      implicit none
            
      include '../includes/dimpar.h'
    
      integer*4 luout,jdstart,jdstop,nfrq,iyr,idoy,imo,iday,ihr,imin
     .        , nzen,naz,ifrq,i,j
                                                             
      real*4 sec,az

      real*8 dazi,zen1,zen2,dzen,offsetl1(3),offsetl2(3)
     .     , elvtabl1(maxel),elvtabl2(maxel)
     .     , tablel1(maxel,maxaz),tablel2(maxel,maxaz)                       

      character*1 svsys  
      character*10 sinex_code                       
      character*15 format
      character*20 anttype,antsn
                           
c** Serial number usually not available  
      antsn = ' ' 

c** Compute the elevation and azimuth limits

      if( dzen.gt.0.d0 ) then
        nzen = int((zen2-zen1)/dzen) + 1
      else       
        dzen = 90.d0
        zen1 = 0.d0
        zen2 = 90.0
        nzen = 2 
      endif

c** Write the antenna header information
                                                 
c      print *,'DEBUG anttype ',anttype
c      print *,'DEBUg antsn   ',antsn
      write(luout,'(60x,a)') 'START OF ANTENNA' 
      antsn = ' ' 
      write(luout,'(a20,a20,20x,a20)') anttype,antsn
     .         ,'TYPE / SERIAL NO    ' 
      write(luout,'(a10,50x,a20)') 'CONVERTED  '
     .   ,'METH / BY / # / DATE'
      write(luout,'(2x,f6.1,52x,a4,16x)') dazi,'DAZI'
      write(luout,'(2x,3f6.1,40x,a20)') zen1,zen2,dzen
     .          ,'ZEN1 / ZEN2 / DZEN  ' 
      nfrq = 2
      write(luout,'(i6,54x,a)') nfrq,'# OF FREQUENCIES'
c     no other format offers a date range, so these records can't be written
      jdstart = 0
      jdstop = 0
      if( jdstart.ne.0 ) then                                     
        call dayjul(jdstart,iyr,idoy)
        call monday( idoy,imo,iday,iyr )  
        ihr = 0                                  
        imin = 0
        sec = 0.
        write(luout,'(5i6,f13.7,17x,a20)') iyr,imo,iday,ihr,imin,sec
     .       ,'VALID FROM         '
      endif
      if( jdstop.ne.0 ) then                                     
        call dayjul(jdstop,iyr,idoy)
        call monday( idoy,imo,iday,iyr )  
        ihr = 0                                  
        imin = 0
        sec = 0.
        write(luout,'(5i6,f13.7,17x,a20)') iyr,imo,iday,ihr,imin,sec
     .       ,'VALID UNTIL        '    
      endif      
      write(luout,'(a10,50x,a)') sinex_code,'SINEX CODE'

c** Write the L1 values
                            
      svsys = 'G'
      ifrq = 1
      write(luout,'(3x,a1,i2,54x,a20)') svsys,ifrq
     .     ,'START OF FREQUENCY  ' 
      write(luout,'(3f10.2,30x,a)') offsetl1,'NORTH / EAST / UP'
      format = ' '            
      write(format,'(a7,i2,a5)') '(3x,a5,',nzen,'f8.2)'  
      write(luout,format) 'NOAZI',(elvtabl1(i),i=1,nzen)
      if( dazi.ne.0.0d0 ) then
        format = ' '
        write(format,'(a6,i2,a5)') '(f8.1,',nzen,'f8.2)'  
        naz = 360/int(dazi) + 1 
        az = 0.
        do j=1,naz
          write(luout,format) az,(tablel1(i,j),i=1,nzen)   
          az = az + dazi
        enddo
      endif
      write(luout,'(3x,a1,i2,54x,a)') svsys,ifrq
     .      ,'END OF FREQUENCY' 

c** Write the L2 values

      ifrq = 2
      write(luout,'(3x,a1,i2,54x,a)') svsys,ifrq
     .     ,'START OF FREQUENCY' 
      write(luout,'(3f10.2,30x,a)') offsetl2,'NORTH / EAST / UP' 
      format = ' '            
      write(format,'(a7,i2,a6)') '(3x,a5,',nzen,'f8.2)'
      write(luout,format) 'NOAZI',(elvtabl2(i),i=1,nzen)
      if( dazi.ne.0.0d0 ) then
        format = ' '
        write(format,'(a6,i2,a6)') '(f8.1,',nzen,'f8.2)'  
        az = 0.
        do j=1,naz
          write(luout,format) az,(tablel2(i,j),i=1,nzen)   
          az = az + dazi
        enddo
      endif
      write(luout,'(3x,a1,i2,54x,a)') svsys,ifrq
     .      ,'END OF FREQUENCY' 
 
c** Terminate antenna blocks

      write(luout,'(60x,a)') 'END OF ANTENNA'

      return
      end
         
c************************************************************************************

      Subroutine write_antmod( luout,gamit_code,amodel
     .                       , dazi,zen1,zen2,dzen
     .                       , offsetl1,offsetl2,elvtabl1,elvtabl2
     .                       , tablel1,tablel2 )


c     Write GAMIT antmod.dat records for a single antenna
     
      implicit none
            
      include '../includes/dimpar.h'
    
      integer*4 luout,nel,naz,incazi,incelv,i,j
                                                             
      real*8 dazi,zen1,zen2,dzen,offsetl1(3),offsetl2(3)
     .     , elvtabl1(maxel),elvtabl2(maxel)
     .     , tablel1(maxel,maxaz),tablel2(maxel,maxaz)   
     .     , fact                    
                      
      character*6 gamit_code        
      character*15 format     
      character*16 amodel
           
c**  Note:
c      GAMIT files have two identical entries for elevation-only; order
c      of entries is from horizon to zenith, opposite of ANTEX order
    
c** Set the constants and limits
             
      fact = 1.0
      if( dzen.eq.90.d0 ) then
c       no PCV, set az/el increments to 0
        nel = 0
        incelv = 0
        naz = 0
        incazi = 0
      else
        nel = int((zen2-zen1)/dzen) + 1   
        incelv = int(dzen)    
        if( dazi.eq.0.d0 ) then   
c         no az-dependent terms, GAMIT uses two entries with inc=360
          naz = 0
          incazi = 360
        else
          naz = 360/int(dazi) + 1   
          incazi = int(dazi)    
        endif
      endif   
      format = ' '
      write(format,'(a4,i2,a5)') '(6x,',nel,'f6.1)'       

c** write the L1 header and values  

c debug      print *,'amodel incazi incelv fact ',amodel,incazi,incelv,fact         
      write(luout,'(1x,a6,2x,a2,3f7.1,2x,a16,1x,i3,3x,i3,f7.1)')
     .   gamit_code,'L1',(offsetl1(i),i=1,3),amodel,incazi,incelv,fact
c     GAMIT files have two identical entries for elevation-only; order
c     of entries is from horizon to zenith, opposite of ANTEX order
c     If only two elevation values (0, 90., from ANTEX, assume constant)
      if( nel.gt.2 ) then
        if( incazi.eq.360 ) then 
          do j=1,2  
            write(luout,format) (elvtabl1(i),i=nel,1,-1) 
          enddo  
        else
          do j=1,naz
            write(luout,format) (tablel1(i,j),i=nel,1,-1)
          enddo
        endif   
      endif
 
c** write the L2 header and values

      write(luout,'(1x,a6,2x,a2,3f7.1,2x,a16,1x,i3,3x,i3,f7.1,2x,a80)')
     .   gamit_code,'L2',(offsetl2(i),i=1,3),amodel,incazi,incelv,fact   
      if( nel.gt.2 ) then
        if( incazi.eq.360 ) then 
          do j=1,2
           write(luout,format) (elvtabl2(i),i=nel,1,-1) 
          enddo  
        else
          do j=1,naz
            write(luout,format) (tablel2(i,j),i=nel,1,-1)
          enddo
        endif
      endif
      
      return
      end



