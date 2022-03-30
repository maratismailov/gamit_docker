* Program to read two eq_files and merge the renames, giving extent priority 
* to a reference file, usually the official ITRF2008 or IGS list.
*
* R. King February 2013
*
* The format of entries is assumed to be 
*                             `163
*   rename ALIC ALIC_1PS 2010 3  6 15  0  0 2100 1 1  0  0 ! comment
*
* Non-blank first column is treated as a comment and not copied.
* 'break' entries are not handled, so you must first run 'break_to_rename'
* to get an extent assigned to the site at the break epoch. 
* 
* The run-string is
*
*   merge_renames [ref-file] [file2] [out-file]
*
*          where ref-file takes precedence.

      implicit none                          
                                 
      integer*4 maxren,maxsit       
*     change these in all subroutines if needed
      parameter(maxren=30000,maxsit=100)

      integer*4 list1_num,list1_sdates(5,maxren),list1_edates(5,maxren)
     .        , list2_num,list2_sdates(5,maxren),list2_edates(5,maxren)
     .        , numsit1,sdates1(5,maxsit),numsit2,sdates2(5,maxsit)
     .        , numren,sdates(5,maxsit),edates(5,maxsit)
     .        , iyear,imonth,iday,ihr,imn,isec,ihnsec
     .        , ilist1,iarg,ioerr,i,j

      character*4 sitcod,list1_sites(maxren),list2_sites(maxren)
     .          , list1_exts(maxren),list2_exts(maxren),exts1(maxsit)
      character*8 sitesrn(maxsit),site2rn      
      character*16 uname
      character*80 arg,infile1,infile2,outfile,cdum
      character*256 line,list1_comments(maxren),list2_comments(maxren)
     .           ,  comments1(maxsit),comments2(maxsit),comments(maxsit)


      logical samesite,end1,list2_used(maxren)
      logical debug/.false./

* FUNCTION
      integer*4 nblen,iclarg

* Get the run-string arguments and open the input files

      iarg = iclarg(1,infile1)
      if( iarg.gt.0 ) then 
        open(unit=1,file=infile1,status='old',iostat=ioerr) 
        if( ioerr.ne.0 ) then
          write(*,'(a)') 'Error opening reference file',infile1,ioerr
         stop
        else
         write(*,'(a)') 'Opened reference file ',infile1
         endif  
        iarg = iclarg(2,infile2)     
        open(unit=2,file=infile2,status='old',iostat=ioerr) 
        if( ioerr.ne.0 ) then
          write(*,'(a)') 'Error opening 2nd input file',infile2,ioerr
          stop
        else
          write(*,'(a)') 'Opened 2nd input file ',infile2
        endif  
      else
        write(*,'(3(/,a))') 'Missing runstring: '
     .        , ' merge_renames [file 1] [file 2] [out-file]'
     .        , ' where file 1 is the reference, taking priority'
        stop
      endif

* Open the output file and write a header
          
      iarg = iclarg(3,outfile)     
      open(unit=3,file=outfile,status='unknown',iostat=ioerr) 
      if( ioerr.ne.0 ) then
        write(*,'(a)') 'Error opening output file',outfile,ioerr
        stop
      else
        write(*,'(a)') 'Opened output file ',outfile
      endif                     
      call getdat(iyear,imonth,iday)
      call gettim(ihr,imn,isec,ihnsec)
      call getusr(uname) 
      write(3,'(a,i4,a,i3,a,i3,2x,i2,a,i2,a,i2,a,a16)')
     .  '* Written by MERGE_RENAMES   Run on ',iyear,'/',imonth,'/'
     .  ,iday,ihr,':',imn,':',isec,' by ',uname
      write(3,'(4a)') '* Reference file: ',infile1(1:nblen(infile1))
     .              ,'   File 2: ',infile2(1:nblen(infile2))
      write(3,'(a)') '*'
      write(3,'(2a)') '* Part I:  Reference file entries with original '
     .   ,'extents merged with File 2 entries with extents APS, BPS. ..'
      write(3,'(2a)') '* Part II (appended): File 2 sites not in '
     .                ,'Reference File (original extents)'          
      write(3,'(a)') '*'
      write(3,'(a)') '*   Part I'
     

        
* Read and store the entries from each file
                                
      call read_entries( 1,list1_num,list1_sites,list1_exts
     .                 , list1_sdates,list1_edates,list1_comments )
      call read_entries( 2,list2_num,list2_sites,list2_exts
     .                 , list2_sdates,list2_edates,list2_comments ) 
                      
      if(debug) then
         write(*,'(a,i4)') 'List 1 # ',list1_num
         do i=1,list1_num
           write(*,'(1x,a4,1x,a4,5i5)') 
     .       list1_sites(i),list1_exts(i),(list1_sdates(j,i),j=1,5)
     .      ,(list1_edates(j,i),j=1,5)
         enddo
        write(*,'(a,i4)') 'List 2 # ',list2_num
         do i=1,list2_num
           write(*,'(1x,a4,1x,a4,5i5)') 
     .       list2_sites(i),list2_exts(i),(list2_sdates(j,i),j=1,5) 
     .      ,(list2_edates(j,i),j=1,5)
         enddo    
      endif


                    
* Loop through list 1, checking for same-site in list 2, then append list 2 entries
               
      do i=1,list2_num
        list2_used(i) = .false.
      enddo        
      end1 = .false.
      ilist1 = 1 
      numsit1 = 1   
      if(debug ) print *,'ilist1 list1_num ',ilist1,list1_num
      do while( ilist1.le.list1_num ) 
        do j=1,5
          sdates1(j,numsit1) = list1_sdates(j,ilist1)
        enddo          
        exts1(numsit1) = list1_exts(ilist1)
        comments1(numsit1) = list1_comments(ilist1)        
        samesite = .true.   
        do while( samesite )         
          sitcod = list1_sites(ilist1)
          if(debug) print *,'numsit1+1 list1_sites sitcod '
     .          ,numsit1+1,list1_sites(numsit1+1),sitcod
          if(debug) print *,'samesite ',samesite
          if( list1_sites(numsit1+1).eq.sitcod ) then
            numsit1 = numsit1 + 1    
            ilist1 = ilist1 + 1
            exts1(numsit1) = list1_exts(ilist1)
            do j=1,5
              sdates1(j,numsit1) = list1_sdates(j,ilist1)
            enddo
            comments1(numsit1) = list1_comments(ilist1)
          else
            samesite = .false.        
          endif
        enddo  
        numsit2 = 0 
        do i=1,list2_num
          if( list2_sites(i).eq.sitcod ) then
            numsit2 = numsit2 + 1
            do j=1,5
              sdates2(j,numsit2) = list2_sdates(j,i)
            enddo
            comments2(numsit2) = list2_comments(i)
            list2_used(i) = .true.
          endif
        enddo  

        if( debug ) then
           write(*,'(a,i4,1x,a4)') 'numsit1 sitcod exts1 sdates1 '
     .          ,numsit1,sitcod
           do i=1,numsit1
             write(*,'(1x,a4,5i5)') exts1(i), (sdates1(j,i),j=1,5)
           enddo
           write(*,'(a,i4)') 'numsit2 sdates2 ',numsit2
           do i=1,numsit2
             write(*,'(5i5)')  (sdates2(j,i),j=1,5)
           enddo                       
         endif

        call merge_entries( sitcod,numsit1,exts1,sdates1,comments1
     .                    , numsit2,sdates2,comments2
     .                    , numren,sitesrn,sdates,edates,comments)
        do i=1,numren
          write(3,'(a,a4,1x,a8,2(i6,4i3),a)') 
     .             ' rename ',sitcod,sitesrn(i)
     .            , (sdates(j,i),j=1,5),(edates(j,i),j=1,5)
     .            , comments(i)(1:nblen(comments(i)))
        enddo
        ilist1 = ilist1 + 1
        numsit1 = 1
        if( debug ) write(*,'(a,i4)') 'End loop ilist1',ilist1
      enddo
                                                                   
      if( debug ) then 
        write(*,'(a)') 'list 2 flags '
        do i=1,list2_num
          write(*,'(a4,1x,l)') list2_sites(i),list2_used(i)
        enddo
      endif
     

c     Now add the list2 entries not merged with list1
                           
      write(3,'(a)') '*'
      write(3,'(2a)') '*  Part II: File 2 sites not in Reference File '
      do i=1,list2_num                       
        if( .not.list2_used(i) ) then
          site2rn = list2_sites(i)//list2_exts(i)
          write(3,'(a,a4,1x,a8,2(i6,4i3),a)') 
     .      ' rename ',list2_sites(i),site2rn,(list2_sdates(j,i),j=1,5)
     .      ,(list2_edates(j,i),j=1,5)
     .      ,list2_comments(i)(1:nblen(list2_comments(i)))
         endif
      enddo                                                                        

      stop
      end
                 
c-------------------------------------------------------------

      Subroutine READ_ENTRIES( lu,numsit,sites,exts,sdates,edates
     .                       , comments )
       

      integer*4 maxren       
      parameter(maxren=30000)

      integer*4 ioerr,lu,indx,idum,sdates(5,maxren),edates(5,maxren)
     .        , numsit,is,nblen,i,trimlen,indx_save
                  
      character*4 sites(maxren),exts(maxren)
      character*6 keywrd                      
      character*8 siteold,sitenew
      character*80 cdum
      character*256 line,comments(maxren)

      logical eof,newsite,debug/.false./

* Initialize the end dates in case they are not present (allowed)
               
      do i=1,maxren
        edates(1,i) = 2100
        edates(2,i) = 1
        edates(3,i) = 1
        edates(4,i) = 0
        edates(5,i) = 0
      enddo                                                           
      
* Loop through the files, storing all the inputs
      
      eof = .false.                             
      is = 0
      do while( .not.eof )     
        read(lu,'(a)',iostat=ioerr ) line
        indx = 2
        if( ioerr.eq.-1.or.line(1:3).eq.'   ' ) then
          eof = .true.
        elseif( ioerr.ne.0 ) then
          write(*,'(a,2i5)') 'Error reading file ',lu,ioerr
          stop
        elseif( line(1:1).ne.' ') then
          continue
        else
          call read_line(line,indx,'CH',ioerr,idum,keywrd)
          call lowers(keywrd)
          if( keywrd(1:3).ne.'ren' ) then
            write(*,'(2a)') 'Invalid keywrd: ',keywrd
            stop
          endif
          is = is + 1              
          if( is.gt.maxren ) then
            write(*,'(a,i5,a,i5 )') 'Number of entries ',is
     .         ,' > maxnum ',maxnum
            stop
          endif
          siteold = ' '
          sitenew = ' '                   
          call read_line(line,indx,'CH',ioerr,idum,siteold)   
          if( ioerr.eq.-1 ) then
            eof = .true. 
          elseif( ioerr.ne.0 ) then
            write(*,'(a,2i5)') 'Error reading site on lu ',lu,ioerr
            stop
          else
c          always use just a 4-character code for the initial name 
           sites(is) = siteold(1:4)
           call read_line(line,indx,'CH',ioerr,idum,sitenew)
           exts(is) = sitenew(5:8)
           call read_line(line,indx,'I4',ioerr,sdates(1,is),cdum)      
           if( ioerr.ne.0 ) then
c          if entry not numerical, assume an h-file code and skip
             call read_line(line,indx,'CH',ioerr,idum,cdum)      
             call read_line(line,indx,'I4',ioerr,sdates(1,is),cdum)      
           endif                                                 
           if( debug ) write(*,'(a,2i5,1x,a8,1x,a8,i5)') 
     .      'lu ioerr read sdate ',lu,ioerr,siteold,sitenew,sdates(1,is)
           call read_line(line,indx,'I4',ioerr,sdates(2,is),cdum)
           call read_line(line,indx,'I4',ioerr,sdates(3,is),cdum)
           call read_line(line,indx,'I4',ioerr,sdates(4,is),cdum)
           call read_line(line,indx,'I4',ioerr,sdates(5,is),cdum)
           if( ioerr.ne.0 ) then                                   
             write(*,'(a,i2,2x,1x,a8,1x,a8,1x,5i6)') 
     .          'Error reading site or start date from file ',lu
     .           ,  siteold,sitenew,(sdates(i,is),i=1,5)
             stop
           endif 
c          we don't know whether the next string is the end-date or a comment or blank   
           indx_save = indx
           call read_line(line,indx,'I4',ioerr,edates(1,is),cdum)
           if( ioerr.eq.0 ) then     
             if( debug ) write(*,'(a,2i5,i5)') 
     ,          'read ioerr edate ',ioerr,edates(1,is)
             call read_line(line,indx,'I4',ioerr,edates(2,is),cdum)
             call read_line(line,indx,'I4',ioerr,edates(3,is),cdum)
             call read_line(line,indx,'I4',ioerr,edates(4,is),cdum)
             call read_line(line,indx,'I4',ioerr,edates(5,is),cdum) 
           else
             indx = indx_save
           endif
           call read_line(line,indx,'CH',ioerr,idum,cdum) 
           if( debug ) write(*,'(2a)') 'read cdum ',cdum  
           if( ioerr.eq.0 ) then
              if ( cdum(1:1) .ne. '!' ) then
                write(*,'(a)') 'Apparent comment begins with ',cdum
                write(*,'(a)') '--must use !,  stop'
                stop
              else
                comments(is) = line(indx-1:)
              endif
	        else
              comments(is) = ' ' 
	        endif    
         endif
        endif 
      enddo  
      numsit = is   
      return
      end        

c---------------------------------------------

      Subroutine MERGE_ENTRIES( sitcod,numsit1,exts1,sdates1,comments1
     .                        , numsit2,sdates2,comments2
     .                        , numren,sitesrn,sdates,edates
     .                        , comments )

c     For a single site, merge the rename entries in two lists and insert
c     them into the full output list.  Retain the extents for list (reference, 
c     e.g. IGS), and assign appropriate new extents for list 2.

c   Input:
c       sitcod            c*4   4-character code for this site    
c       numsit1           i*4   number of entries for this site in list 1
c       exts1(maxsit)     c*4   extents associated with this site in list 1
c       sdates1(5,maxsit) i*4   start dates YYYY MM DD HH MM for entries in list 1
c       comments1(maxsit) c*256 comments for entries in list 1
c       numsit2           i*4   number of entries for this site in list 2
c       sdates2(5,maxsit) i*4   start dates YYYY MM DD HH MM for entries in list 2
c       comments2(maxsit) c*256 comments for entries in list 2 
c   Output
c       numren            i*4   number of merged entries
c       sitesrn(maxsit)   c*8   renamed sites in merged entries
c       sdates(5,maxsit)  i*4   start dates in merged entries
c       edates(5,maxsit)  i*4   end dates in merged entries
c       comments          c*256 comments in merged entries

      implicit none
        
                
      integer*4 maxsit
      parameter(maxsit=100)

      integer*4 numsit1,numsit2,sdates1(5,maxsit),sdates2(5,maxsit)
     .        , numren,sdates(5,maxsit),edates(5,maxsit),iorder
     .        , ns1,ns2,i,j
                          
      character*1 ext2char
      character*4 sitcod,exts1(maxsit),sites(maxsit)
      character*8 sitesrn(maxsit)
      character*256 comments1(maxsit),comments2(maxsit),comments(maxsit)

      real*8 tk,tkp,decyrs
          
      logical end1,end2, debug/.false./

c  FUNCTION
      integer*4 itimdif,nblen
                                
c  Initialize the extent character for the secondary lists (first new one will be 'A')
      ext2char = '@'

      if( debug ) then
        write(*,'(a,1x,a4,i4)') 'In MERGE_ENTRIES sitcod numsit1 '
     .       ,sitcod,numsit1
        do i=1,numsit1
          write(*,'(1x,a4,5i6,a)') exts1(i),(sdates1(j,i),j=1,5)
     .        ,comments1(i)(1:nblen(comments1(i)))
        enddo                      
        print *,'numsit2 ',numsit2
        do i=1,numsit2
          print *,(sdates2(j,i),j=1,5)
     .      ,comments2(i)(1:nblen(comments2(i)))
        enddo
      endif

c  Merge based just on start dates, then set the end dates to be consistent.
                            
c  See if one list is empty
                    
      numren = 0 
      if( numsit2.eq.0 ) then 
        do i=1,numsit1        
          numren = numren + 1
          sitesrn(numren) = sitcod//exts1(i)
          do j=1,5
            sdates(j,i) = sdates1(j,i)
          enddo               
          comments(i) = comments1(i)
        enddo
      elseif( numsit1.eq.0 ) then
        do i=1,numsit2
          numren = numren + 1                    
          if(i .eq.1 ) then 
            ext2char = 'A'
          else
            ext2char = char(ichar(ext2char)+1)
          endif
          ext2char = char(ichar(ext2char)+1)
          sitesrn(numren) = sitcod//"_"//ext2char//"PS"
          do j=1,5
            sdates(j,i) = sdates1(j,i)
          enddo             
          comments(i) = comments2(i)
        enddo
      else

c Go through the lists, interleaving entries
              
        ns1 = 0
        ns2 = 0 
        numren = 0 
        end1 = .false.
        end2 = .false.
        do while( .not.end1 .and. .not.end2 )
                           
          if( ns1.ge.maxsit.or.ns2.ge.maxsit.or.numren.ge.maxsit ) then
             write(*,'(a,1x,a4)') 'Too many entries for site ',sitcod
             stop
          endif 
          call next_entry( numsit1,ns1,sdates1,numsit2,ns2,sdates2
     .                   , iorder,end1,end2 )      
          if(debug) then
            print *,'aft next_entry ns1 end1 ns2 end2 iorder numren '
     .        ,ns1,end1,ns2,end2,iorder,numren
          endif
c          iorder = 0  next entries at the same time
c                   1  next entry is list 1
c                   2  next entry is list 2
                               
          if( .not.end1 .and.(iorder.eq.0 .or.iorder.eq.1) ) then
            numren = numren + 1 
            ns1 = ns1 + 1                
            if( iorder.eq.0 ) ns2 = ns2 + 1
            sitesrn(numren) = sitcod//exts1(ns1)
            do j=1,5
              sdates(j,numren) = sdates1(j,ns1)    
            enddo                            
            comments(numren) = comments1(ns1)
            if( debug ) write(*,'(a,2i4,1x,a8,5i4)') 
     .          'numren iorder sitesrn sdates'
     .          , numren,iorder,sitesrn(numren),(sdates(j,numren),j=1,5)
          elseif( .not.end2 .and.iorder.eq.2 ) then
            numren = numren +1
            ns2 = ns2 + 1
            sitesrn(numren) = sitcod
            ext2char = char(ichar(ext2char)+1)
            sitesrn(numren) = sitcod//"_"//ext2char//"PS" 
            do j=1,5
              sdates(j,numren) = sdates2(j,ns2)
            enddo                     
            comments(numren) = comments2(ns2)  
            if( debug ) write(*,'(a,i4,1x,a8,5i4)') 
     .        'numren iorder 2 sitesrn',numren,sitesrn(numren)
          else
c           assume both lists ended 
            if( debug ) write(*,'(a,2i4)') 'end of lists ns1 ns2 numren'
     .         ,ns1,ns2,numren
            continue
          endif
        enddo
      endif


c Now fix the end times to match the next start time

      do i = 1, numren-1
        do j=1,5
          edates(j,i) = sdates(j,i+1)
        enddo   
      enddo
      edates(1,numren) = 2100
      edates(2,numren) = 1 
      edates(3,numren) = 1
      edates(4,numren) = 0 
      edates(5,numren) = 0 
      if( debug ) then
         write(*,'(a,i4)') 'setting end times, numren ',numren
         do i=1,numren
           write(*,'(5i4)') (edates(j,i),j=1,5)
         enddo
      endif

      return
      end

c----------------------------------------------------------------------=      
      Subroutine NEXT_ENTRY( numsit1,ns1,sdates1,numsit2,ns2,sdates2
     .                     ,iorder,end1,end2 )

      implicit none 
                   
      integer*4 maxsit
      parameter(maxsit=100)

      integer*4 numsit1,ns1,sdates1(5,maxsit)
     .        , numsit2,ns2,sdates2(5,maxsit)
     .        , iorder
                             
      real*8 next1,next2,tol
                                                        
      logical end1,end2,debug/.false./ 

      real*8 timarg


c     Compare the next times for list 1 and list 2 and return a code:
c        iorder = 0   next entries are the same in both lists
c        iorder = 1   list 1 entry is next                   
c        iorder = 2   list 2 entry is next

c  In all of the time checks, use a tolerance of 1.1 days, on the assumption that 
c  a 1-day difference is probably an entry error, not a real difference in the time
c  of the change.    

      tol = 1.1d0/365.25d0

                                     
      if( debug ) write(*,'(a,4i4)') 
     .    'In NEXT_ENTRY ns1 ns2 numsit1 numsit2 '
     .     ,ns1,ns2,numsit1,numsit2
      if( ns1+1.le.numsit1 ) then
        next1 = timarg(sdates1(1,ns1+1))
      else
         end1 = .true.
      endif
      if( ns2+1.le.numsit2 ) then
        next2 = timarg(sdates2(1,ns2+1))
      else
        end2 = .true.
      endif      
      if( debug ) write(*,'(a,2f15.5,2L1)') ' next1 next2 end1 end2 '
     .          , next1,next2,end1,end2
      if( .not.end1 .and. .not.end2 ) then                            
        if( dabs(next1-next2) - tol.le.0 ) then
c         times are the same
          iorder = 0 
        elseif( (next1-next2).lt.0 ) then
c         list 1 next
          iorder = 1
        else
c         list 2 next
          iorder = 2
        endif
      elseif( end2 ) then
        iorder = 1
      elseif( end1 ) then
        iorder = 2
      endif    
      if( debug ) write(*,'(a,2f15.5,2L1,i3)') 
     .         ' next1 next2 end1 end2 iorder '
     .         , next1,next2,end1,end2,iorder


      return
      end

c-----------------------------------------------------------------------

      Function TIMARG( date )
               
c     Get days since 1900 to use as a time argument
    
      integer*4 date(5)
      real*8 timarg        
                                                                                              
      timarg =  dfloat(date(1)-1900) + 
     .      ( dfloat(idoy(date(1),date(2),date(3))) +
     .        dfloat(date(4))/24.d0 + dfloat(date(5))/1440.d0 )/365.25d0
      return
      end     


c---------------------------------------------------------------------------------------------
