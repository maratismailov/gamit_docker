Copyright (c) Massachusetts Institute of Technology and the University of
California, San Diego. 1994.   All rights reserved.

      Subroutine read_antmod( luin,antcod,found,amodel
     .                      , dazi,zen1,zen2,dzen
     .                      , offsetl1,offsetl2,elvtabl1,elvtabl2
     .                      , tablel1,tablel2)
c
c Purpose    :  To read a GAMIT-style table of antenna phase-center variations
c
c Parameters :
c         in:  luin      : unit number of antmod.dat file                       i*4
c              antcod    : requested GAMIT antenna code name                    ch*6
c                          (if blank, read the next entry and return the antcod)
c                       
c        out : found                : .true. if requested antenna in table        logical
c              amodel               : full name of model from table               ch*16 
c              dazi                 : azimuth increment (=0 if elev-only)         r*8
c              zen1                 : minimum zenith angle (usually 0)            r*8
c              zen2                 : maxium zenith angle (usually 90)            r*8
c              dzen                 : zenith angle increment in deg               r*8
c              offsetl1(3)          : antenna ARP to L1 phase centre (U N E)      r*8
c              offsetl2(3)          : antenna ARP to L2 phase centre (U N E)      r*8
c              elvtabl1(maxel)      : elevation-dependent phase center offset L1  r*8
c              elvtabl2(maxel)      : elevation-dependent phase center offset L2  r*8
c              tablel1(maxel,maxaz) : elev & azimuth-dependent pc offset L1       r*8  
c              tablel2(maxel,maxaz) : elev & azimuth-dependent pc offset L1       r*8    
c                                     1st element is zenith angle, 0-90, 2d is azimuth, 0-360
c                                     Note that the order of vertical angles on output is
c                                     reversed from antmod.dat and the original read_antmod
c                                     in order to conform with ANTEX conventions.

c  Format of GAMIT antmod.dat (columns are 0-90 elevation angle, opposite of output)
c--------
c
c ASHL12  L1   79.9    0.5    0.3  I1_IGS01         360    05    1.0  ASHTECH 700228 A
c         0.0   0.0   0.9   0.7   1.5   2.3   2.6   2.8   3.0   3.3   3.0   2.4   2.1   2.1   1.8   1.2   0.5   0.1   0.0
c         0.0   0.0   0.9   0.7   1.5   2.3   2.6   2.8   3.0   3.3   3.0   2.4   2.1   2.1   1.8   1.2   0.5   0.1   0.0
c ASHL12  L2   79.2   -1.2    0.8  I1_IGS01         360    05    1.0  ASHTECH 700228 A
c         0.0   0.0  -1.8   1.6   2.4   1.8   1.5   1.7   1.9   2.0   2.1   2.3   2.2   1.8   1.6   1.5   1.1   0.4   0.0
c         0.0   0.0  -1.8   1.6   2.4   1.8   1.5   1.7   1.9   2.0   2.1   2.3   2.2   1.8   1.6   1.5   1.1   0.4   0.0         
c         0.0   0.0  -1.8   1.6   2.4   1.8   1.5   1.7   1.9   2.0   2.1   2.3   2.2   1.8   1.6   1.5   1.1   0.4   0.0         
c 
c where, in this case, the azimuth increment = 360 and identical lines mean that there is no azimuth 
c dependence (cf. ANTEX, which has azimuth increment = 0 and no repeated lines)


c created    :  2003/4/14 by R. King, from subset of code in read_antmod.f written
c               by S. McClusky (1995/10/17 (rest of code now in get_antpcv.f)

      Implicit none

      include '../includes/dimpar.h'
      

      character*1   cvalue
      character*2   afrq
      character*6   antcod,antcod1,antcod2
      character*16  amodel,amodel2
      character*80  prog_name
      character*256 line,message   

      integer*4 luin,len,rcpar,incazi,incelv,naz,nel
     .         ,ioerr,lerr,indx,i,j

      real*8  sign,rvalue,offsetl1(3),offsetl2(3)
     .      , elvtabl1(maxel),elvtabl2(maxel)
     .      , dazi,zen1,zen2,dzen
     .      , tablel1(maxel,maxaz),tablel2(maxel,maxaz)

      logical eof,found 
                

c     Get the calling module name for report_stat
      len = rcpar(0,prog_name)
          
c**   Initialize the output variables
                   
      found = .false.
      amodel = ' ' 
      amodel2 = ' '
      sign = 1.d0
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

c**   Start reading the table looking for the correct antenna type
                   
      eof = .false.
      do while ( .not.eof .and. .not.found ) 

        read(luin,'(a256)',iostat=ioerr) line 
c        print *,'READ file ',line 
        if(ioerr.eq.-1 .or. line(2:7).eq.'THEEND') then 
          eof = .true.   
        elseif( ioerr.ne.0 ) then  
          call report_stat('FATAL',prog_name,'lib/read_antmod',' ' 
     .     ,'Error reading antmod.dat ',ioerr)
        elseif( line.eq.' ' .or. line(1:1).ne.' ' ) then
c         skip blank line or comment (non-blank first column) 
          continue
        else
           
c         see if it's the requested antenna   
          if( line(2:7).eq.antcod .or. antcod.eq.'      ') then    
            found = .true. 
c            print *,'found true' 
c           need to assign antcod if it's blank
c           if so, decode the header line for L1
            indx = 1  
            call read_line(line,indx,'CH',lerr,rvalue,antcod1)
c           need to assign antcod for return if it's blank
            antcod = antcod1
            call read_line(line,indx,'CH',lerr,rvalue,afrq) 
            if( afrq.ne.'L1')
     .         call report_stat('FATAL',prog_name,'lib/read_antmod',' '
     .       ,'Cannot read old-style antmod.dat table--get a new one',0)  
c           read(line,'(50x,i4,2x,i4,f7.0)',iostat=ioerr) incazi,incelv
c     .         ,sign          
           if( ioerr.ne.0 ) call report_stat('FATAL',prog_name
     .       ,'lib/read_antmod',' '
     .       ,'Error decoding Az/El increments for L1',ioerr)
           do i=1,3
             call read_line(line,indx,'R8',lerr,offsetl1(i),cvalue)
           enddo                                                   
           call read_line(line,indx,'CH',lerr,rvalue,amodel)
           call read_line(line,indx,'I4',lerr,incazi,cvalue)
           call read_line(line,indx,'I4',lerr,incelv,cvalue)
           call read_line(line,indx,'R8',lerr,sign,cvalue)   
c         ("incazi" is the azimuth interval of the table in degrees 
c          360 degrees must evenly divisible by this value).    
           if( lerr.gt.0 ) call report_stat('FATAL',prog_name
     .      ,'lib/read_antmod',' ','Error decoding header for L1 ',lerr)
           if( incazi.gt.0 ) then
             naz=(360/incazi)+1   
           else
             naz = 0
           endif
c          if only two azimuth entries, they are 0 and 360, so elev-dependent model
           if( naz.eq.2 ) dazi = 0                
           if (incelv.eq.0) then
c            incelv =0 means constants-only
             nel = 0 
             zen1 = 0.d0
             zen2 = 90.d0
             dzen = 0.d0  
           else
             nel=(90/incelv)+1
             if( nel.gt.maxel ) then  
               write(message,'(a,i3,a,i3,a)') 
     .              'Number of elevation values (',nel
     .              ,') exceed maxel (',maxel,')' 
               call report_stat('FATAL',prog_name,'lib/read_antmod'
     .                         ,' ',message,0)
             endif
             zen1 = 0.
             zen2 = 90. 
             dzen = dfloat(incelv)  
           endif     
c          read and store the L1 values 
           if( naz.gt.0 ) then
             do j=1,naz
               read(luin,*,iostat=ioerr) (tablel1(i,j),i=nel,1,-1)  
               if( ioerr.ne.0 ) call report_stat('FATAL',prog_name
     .          ,'lib/read_antmod',' '
     .          ,'Error reading values from antmod.dat',ioerr)   
             enddo   
c            compute azimuth-averaged values     
             do i=1,nel
               do j=1,naz
                 elvtabl1(i) = elvtabl1(i) + tablel1(i,j)
               enddo
               elvtabl1(i) = elvtabl1(i)/dfloat(naz) 
             enddo
           endif
c          now read and decode the L2 header line
           read(luin,'(a256)',iostat=ioerr) line 
           if( ioerr.ne.0 ) then 
            call report_stat('FATAL',prog_name,'lib/read_antmod',' ' 
     .       ,'Error reading L2 corrections from antmod.dat ',ioerr)
           elseif( line.eq.' ' .or. line(1:1).ne.' ' ) then
c            skip blank line or comment (non-blank first column) 
             continue
           else      
             indx = 1  
             call read_line(line,indx,'CH',lerr,rvalue,antcod2)  
             call read_line(line,indx,'CH',lerr,rvalue,afrq)
             do i=1,3
              call read_line(line,indx,'R8',lerr,offsetl2(i),cvalue)
             enddo            
             call read_line(line,indx,'CH',lerr,rvalue,amodel2)
             call read_line(line,indx,'I4',lerr,incazi,cvalue)
             call read_line(line,indx,'I4',lerr,incelv,cvalue)
             call read_line(line,indx,'R8',lerr,sign,cvalue)
c            write(*,123) anttyp_table,ifreq(1),offsetl1(1)   
             if( lerr.gt.0 ) call report_stat('FATAL',prog_name
     .      ,'lib/read_antmod',' ','Error decoding header for L2 ',lerr)
c            make sure we have L2 and that the antenna code and model match the L1 entries    
             if( antcod2.ne.antcod1 .or. afrq.ne.'L2' .or.
     .           amodel2.ne.amodel ) then 
                write(message,'(a,a6 )')  
     . 'L2 line in antmod.dat does not match L1 line for antcod ',antcod
                call report_stat( 'FATAL',prog_name,'lib/read_antmod'
     ,                          ,' ',message,0 ) 
             endif
c            assume that the azimuth and elevation limits are the same
c            read and store the L1 values 
             if( naz.gt.0 ) then
               do j=1,naz
                 read(luin,*,iostat=ioerr) (tablel2(i,j),i=nel,1,-1)   
                 if( ioerr.ne.0 ) call report_stat('FATAL',prog_name
     .          ,'lib/read_antmod',' ','Error reading antmod.dat',ioerr)   
               enddo  
c              compute azimuth-averaged values     
               do i=1,nel
                 do j=1,naz
                   elvtabl2(i) = elvtabl2(i) + tablel2(i,j)
                 enddo
                 elvtabl2(i) = elvtabl2(i)/dfloat(naz) 
               enddo
             endif      
c          end if for valid L2 header record
           endif      

c         end if for correct antenna type
          endif                            

c       end if for valid antenna header record
        endif

c     end do over records until eof or found    
      enddo                                 
       
      return
      end

