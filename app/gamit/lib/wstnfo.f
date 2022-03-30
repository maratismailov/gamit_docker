      Subroutine wstnfo( lu, nlist, clist
     .                 , first_call, maxhdcmt, nhdcmt, head_comments
     .                 , ccom, sitcod, sname, anth, antn, ante, antdaz
     .                 , rcvcod, antcod, htcod, radome, swver
     .                 , rctype, rcvrsn, rcvers, anttyp, antsn  
     .                 , sessn, start, stop, line_comment )   


c     Write the records of station.info.  If first call, write the
c     header and and comments; on subsequent calls, write one record
c     R. King 2 May 2001
c      --header labels modified 12 December 2002   

c MOD TAH 200203: Added AntDAZ token for antenna  Alignment from True N

c     Currently coded with 19 values defined, and with the largest
c     field 19 columns wide

      implicit none

c        unit number
      integer*4 lu
             
c        number, token list, and values of items for each record  
      integer*4 nlist     
      character*6 clist(20)
      character*20 value(20)  

c       station.info variables   
      character*1 ccom
      character*4 sitcod
      character*5 htcod,radome  
      character*6 rcvcod,antcod
* MOD TAH 200213: Changed from  character*10 to character*(*)
      character*(*) rcvrsn
      character*15 anttyp
      character*16 sname
      character*20 rctype,rcvers,antsn 
      character*(*) line_comment     
      real*4 swver
      real*8 anth,antn,ante 
      real*8 antdaz  ! Antenna aligment from True N (deg).
      integer*4 sessn,start(5),stop(5)
            
c       flag indicating that the header is to be written
      logical first_call

c       comments and labels to be written into the header  
      integer*4 maxhdcmt,nhdcmt
      character*132 head_comments(maxhdcmt) 
      character*256 label_line 
    
c        other variables
      integer*4 trimlen,len,rcpar,icol,ioerr,i,j   
      character*20 format(20)
      character*80 prog_name
      character*256 line_fmt

      save format                

c       get program name for report_stat calls
      len = rcpar(0,prog_name)  

           
c On first call, write the header information and construct the format for values

      if( first_call ) then  

c     Write comments at the top of the file

         if( nhdcmt.gt.maxhdcmt ) call report_stat('FATAL',prog_name,
     .    'lib/wrtnfo',' ','Too many comments in station.info',0)
         do i=1,nhdcmt
           write(lu,'(a)') head_comments(i)(1:trimlen(head_comments(i)))   
         enddo   
         write(lu,'(a)') '*'

c     Construct and write the labels line and construct the format for writing values

c       assumes the order sitcod sname(optional), start, stop, then others
        label_line = ' '  
        if( clist(1).eq.'sitcod' ) then
            write(label_line(1:5),'(a5)') '*SITE' 
            format(1) = 'a1,a4' 
          else      
             call report_stat('FATAL',prog_name,'lib/wrstnfo',' '
     .        ,'sitcod not first item in station.info header list',0)
        endif 
        icol = 6
        do i= 2,nlist  
          if( clist(i).eq.'sname '  ) then    
            write(label_line(icol:icol+17),'(a)') '  Station Name    ' 
            icol = icol + 18     
            format(i)=  ',2x,a16'
          elseif( clist(i).eq.'start ' ) then
            write(label_line(icol:icol+18),'(a)') '  Session Start    '
            icol = icol + 19        
            format(i)= ',2x,a17'
          elseif( clist(i).eq.'stop  ' ) then
            write(label_line(icol:icol+18),'(a)') '  Session Stop     '
            icol = icol + 19              
            format(i)= ',2x,a17'   
          elseif( clist(i).eq.'sessn' ) then
            write(label_line(icol:icol+6),'(a)') ' SN' 
            icol = icol + 3   
            format(i) = ',2x,i1'
          elseif( clist(i).eq.'anth  ' ) then
            write(label_line(icol:icol+8),'(a)') '  Ant Ht '
            icol = icol + 9        
            format(i) = ',2x,a7'
          elseif( clist(i).eq.'htcod ' ) then
            write(label_line(icol:icol+6),'(a)') '  HtCod'
            icol = icol + 7  
            format(i) = ',2x,a5'
          elseif( clist(i).eq.'antn  ' ) then
            write(label_line(icol:icol+8),'(a)') '  Ant N '
            icol = icol + 9     
            format(i) = ',2x,a7'
          elseif( clist(i).eq.'ante  ' ) then
            write(label_line(icol:icol+8),'(a)') '  Ant E '
            icol = icol + 9     
            format(i) = ',2x,a7'                         
          elseif( clist(i).eq.'antdaz' ) then
            write(label_line(icol:icol+8),'(a)') '  AntDAZ'
            icol = icol + 9     
            format(i) = ',2x,a7'                         
          elseif( clist(i).eq.'rcvcod' ) then
            write(label_line(icol:icol+7),'(a)') '  RcvCod' 
            icol = icol + 8   
            format(i) = ',2x,a6'
          elseif( clist(i).eq.'rctype' ) then
            write(label_line(icol:icol+21),'(a)') 
     .                              '  Receiver Type       '
            icol = icol + 22   
            format(i) = ',2x,a20'       
          elseif( clist(i).eq.'rcvrsn' ) then
            write(label_line(icol:icol+21),'(a)') 
     .                                  '  Receiver SN        '
            icol = icol + 22   
            format(i) = ',2x,a20'
          elseif( clist(i).eq.'swver' ) then
            write(label_line(icol:icol+6),'(a)')    '  SwVer'
            icol = icol + 7 
            format(i) = ',2x,a5'      
          elseif( clist(i).eq.'rcvers' ) then
            write(label_line(icol:icol+21),'(a)') 
     .                                     '  Vers                '
            icol = icol + 22      
            format(i) = ',2x,a20' 
          elseif( clist(i).eq.'antcod' ) then
            write(label_line(icol:icol+7),'(a)') '  AntCod'
            icol = icol + 8   
            format(i) = ',2x,a6'
          elseif( clist(i).eq.'anttyp' ) then
            write(label_line(icol:icol+16),'(a)') '  Antenna Type   ' 
            icol = icol + 17  
            format(i) = ',2x,a15'     
          elseif( clist(i).eq.'antsn' ) then
            write(label_line(icol:icol+21),'(a)') 
     .                                       '  Antenna SN         '
            icol = icol + 22   
            format(i) = ',2x,a20'
          elseif( clist(i).eq.'radome' ) then
            write(label_line(icol:icol+6),'(a)') '  Dome '
            icol = icol + 7 
            format(i) = ',2x,a5'
          endif
cd          print *,'icol header ',icol,label_line
        enddo  
        write(lu,'(a)') label_line(1:icol)
c     end of first-call assignments
      endif  
 

c Write the entry in the new station.info    
                 
cd      print *,'in WSTNFO nlist clist ',nlist,(clist(i),i=1,nlist)
      do i= 1,nlist
        if( clist(i).eq.'sitcod' ) then
          value(i) = sitcod
        elseif( clist(i).eq.'sname'  ) then
          value(i)= sname
        elseif( clist(i).eq.'sessn' ) then 
          write(value(i),'(i1)') sessn
        elseif( clist(i).eq.'start'  ) then
          write(value(i),'(2i4,3i3)') (start(j),j=1,5)
        elseif( clist(i).eq.'stop'  ) then    
          write(value(i),'(2i4,3i3)') (stop(j),j=1,5)
        elseif( clist(i).eq.'anth' ) then 
          write(value(i),'(f7.4)') anth 
        elseif( clist(i).eq.'htcod') then
          value(i) = htcod
        elseif( clist(i).eq.'antn' ) then
          write(value(i),'(f7.4)') antn
        elseif( clist(i).eq.'ante' ) then
          write(value(i),'(f7.4)') ante
        elseif( clist(i).eq.'antdaz' ) then
          write(value(i),'(f5.0)') antdaz
        elseif( clist(i).eq.'rcvcod' ) then
          value(i) = rcvcod
        elseif( clist(i).eq.'rctype' ) then
          value(i) = rctype
        elseif( clist(i).eq.'rcvrsn' ) then
          value(i) = rcvrsn
        elseif( clist(i).eq.'swver' ) then
          write(value(i),'(f5.2)') swver
        elseif( clist(i).eq.'rcvers' ) then
          value(i) = rcvers
        elseif( clist(i).eq.'antcod' ) then
          value(i) = antcod
        elseif( clist(i).eq.'anttyp' ) then
          value(i) = anttyp
        elseif( clist(i).eq.'antsn' ) then
          value(i) = antsn
        elseif( clist(i).eq.'radome' ) then
          value(i) = radome
        endif
      enddo                   
cd      print *,'nlist value ',nlist,(value(i),i=1,nlist)
c     construct format from item list formats plus one extra field for comments
      write(line_fmt,'(a,(20a13),a)',iostat=ioerr) 
     .    '(',(format(i),i=1,nlist),',a)' 
cd      print *,'line_fmt ',line_fmt
      if( ioerr.ne.0 )  then
c         print *,'Internal I/O error writing line_fmt '
c         print *,'format list: ',(format(i),i=1,nlist)
        call report_stat('FATAL',prog_name,'lib/wstnfo',' '
     .      , 'Internal I/O error writing line_fmt',ioerr)
c         print *,'Internal I/O error writing line_fmt '
c         print *,'format list: ',(format(i),i=1,nlist) 
      endif                       
      write(lu,line_fmt,iostat=ioerr) ccom,(value(i),i=1,nlist)
     .                    ,line_comment(1:trimlen(line_comment))
cd      write(*,line_fmt,iostat=ioerr) ccom,(value(i),i=1,nlist)
cd     .                    ,line_comment(1:trimlen(line_comment))
      if( ioerr.ne.0 ) then
cd           print *,'Error writing values to station.info '
cd            print *,'values:',(value(i),i=1,nlist) 
         call report_stat('FATAL',prog_name,'lib/wrstnfo',' '
     .        ,'Error writing values to station.info',ioerr)
      endif

      return
      end



