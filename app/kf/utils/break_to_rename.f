* Program to convert a file of site discontinuities (dates only) to 
* GLOBK eq_file format with specified renames
* R. King December 2008 
* Modified January 2011 to first sort the break file    
*
* The format of the 'break' file is
*
*   break   ALIC  2003  6 15  0  0 
*   break   ALIC  2003 10 23  0  0 
*   break   ALIC  2008  4 16  0  0 RPS       
* 
* which result in the following renames
*   rename ALIC ALIC_1PS 2003  6 15  0  0 2003 10 23  0  0
*   rename ALIC ALIC_2PS 2003 10 23  0  0 2008  4 16  0  0
*   rename ALIC ALIC_RPS 2008  4 16  0  0 2100  1  1  0  0
*
* where the last entry includes an optional specified rename.
* Non-blank first column is treated as a comment and not copied.
* Break entries in the input file need not be in order and will
* sorted by site and time before assigning renames.
*
* The program can also be run in reverse to 'undo' a current file of
* non-sequential renames (but no sort is performed)
*
* The program takes two arguments, the input file name and an optional
* control to specify break-to_rename ('break', default) or rename-to-break
* ('rename').  A input file may have both breaks and renames, but 
* only the keyword specified will be read.
*
* Examples:
*
*   break_to_rename igs.breaks 
*   break_to_rename igs.renames rename
* 
* The output file uses the input file name with the extent '.out'

      implicit none                          
                                 
      integer*4 maxsit       
*     change this also in sorteq below
      parameter(maxsit=30000)

      integer*4 ioerr,indx,idum,idate1(5,maxsit),idate2(5,maxsit)
     .        , numsit,is,nblen,iclarg,iarg,i,trimlen
        
      character*1 extchar                        
      character*4 ext
      character*6 keywrd,filetype                      
      character*8 siteold(maxsit),sitenew(maxsit)
      character*80 arg,infile,outfile,cdum
      character*256 line,comment(maxsit)

      logical eof,newsite,nextnew

* Print help if no arguments

      iarg = iclarg(1,arg)
      if( iarg.eq.0 ) then
        write(*,'(a)') 
     .  'break_to_rn needs at least the file name as argument'      
        stop
      endif                                                  

* Get the run-string arguments and open the input file

      iarg = iclarg(1,infile)
      open(unit=1,file=infile,status='old',iostat=ioerr) 
      if( ioerr.ne.0 ) then
        write(*,'(a)') 'Error opening input file',infile,ioerr
        stop
      else
        write(*,'(a)') 'Opened input file ',infile
      endif  
      filetype = 'break'
      iarg = iclarg(2,arg)
      if( iarg.gt.0 ) then
        read(arg,*,iostat=ioerr) filetype 
        if( filetype(1:3).eq.'ren'.or. filetype(1:5).eq.'REN') then
          filetype = 'rename'
        else
          filetype = 'break '
        endif
      endif   

* Open the output file and write a header
                          
      outfile = infile(1:nblen(infile))//'.out'                     
      open(unit=2,file=outfile,status='unknown',iostat=ioerr)
      if( ioerr.ne.0 ) then
        write(*,'(a)') 'Error opening output file ',outfile,ioerr
        stop
      else
        write(*,'(a)') 'Opened output file ',outfile
      endif            
      if( filetype.eq.'break ' ) then
       write(2,'(2a)') '* GLOBK eq_file written from break file ',infile
       write(2,'(a)') '*'
      elseif( filetype.eq.'rename' ) then
       write(2,'(2a)') '* Break file written from eq_file file ',infile
       write(2,'(a)') '*'         
      else
       write(*,'(3a)') 'Filetype ',filetype,' not recognized'
      endif                                                                        

        
* Loop through the file, storing all the inputs
                          
      eof = .false.                             
      is = 0
      do while( .not.eof )     
        read(1,'(a)',iostat=ioerr ) line
        indx = 2
        if( ioerr.eq.-1 ) then
          eof = .true.
        elseif( ioerr.ne.0 ) then
          write(*,'(a)') 'Error reading the input file ',ioerr   
          stop
        elseif( line(1:1).ne.' ') then
          continue
        else
          call read_line(line,indx,'CH',ioerr,idum,keywrd)
          call lowers(keywrd)
cd        print *, 'keword, line ',keywrd(1:3),line
          if( filetype.eq.'break '.and.keywrd(1:3).eq.'bre' ) then  
            is = is + 1              
            if( is.gt.maxsit ) then
              write(*,'(a,i7,a)') 
     .          'Number of breaks greater than maxsit (',maxsit,')'
            endif
            siteold(is) = ' '
            call read_line(line,indx,'CH',ioerr,idum,siteold(is))
            call read_line(line,indx,'I4',ioerr,idate1(1,is),cdum)
            call read_line(line,indx,'I4',ioerr,idate1(2,is),cdum)
            call read_line(line,indx,'I4',ioerr,idate1(3,is),cdum)
            call read_line(line,indx,'I4',ioerr,idate1(4,is),cdum)
            call read_line(line,indx,'I4',ioerr,idate1(5,is),cdum)
            if( ioerr.ne.0 ) then    
              write(*,'(a,a8,1x,5i5)') 
     .           'Error reading site or date from break file: '
     .            ,  siteold(is),(idate1(i,is),i=1,5)
              stop  
            endif                    
            if( siteold(is)(5:8).eq.'    ') siteold(is)(5:8) = '_GPS'
            sitenew(is) = siteold(is)
c read optional new extent
            call read_line(line,indx,'CH',ioerr,idum,ext)   
            if( ioerr.eq.0 ) then
               call uppers(ext)
c reading ! means the rest of the line is a comment not ext name
	       if ( ext(1:1) .ne. '!' ) then
                 sitenew(is)(5:8) = '_'//ext
	       else
            comment(is) = line(indx-1:)
	       endif
	       
        endif
      elseif( filetype.eq.'rename'.and.keywrd(1:3).eq.'ren' ) then
            is = is + 1 
            siteold(is) = ' '
            call read_line(line,indx,'CH',ioerr,idum,siteold(is))   
            sitenew(is) = ' ' 
            call read_line(line,indx,'CH',ioerr,idum,sitenew(is))   
            call read_line(line,indx,'I4',ioerr,idate1(1,is),cdum)      
            if( ioerr.ne.0 ) then
c if entry not numerical, assume an h-file code and skip
              call read_line(line,indx,'CH',ioerr,idum,cdum)      
              call read_line(line,indx,'I4',ioerr,idate1(1,is),cdum)      
            endif   
            call read_line(line,indx,'I4',ioerr,idate1(2,is),cdum)
            call read_line(line,indx,'I4',ioerr,idate1(3,is),cdum)
            call read_line(line,indx,'I4',ioerr,idate1(4,is),cdum)
            call read_line(line,indx,'I4',ioerr,idate1(5,is),cdum)
            if( ioerr.ne.0 ) then                                   
              write(*,'(a,1x,a8,1x,a8,1x,5i6)') 
     .           'Error reading site or start date from eq file: '
     .            ,  siteold(is),sitenew(is),(idate1(i,is),i=1,5)
              stop
            endif    
            call read_line(line,indx,'I4',ioerr,idate2(1,is),cdum)
            call read_line(line,indx,'I4',ioerr,idate2(2,is),cdum)
            call read_line(line,indx,'I4',ioerr,idate2(3,is),cdum)
            call read_line(line,indx,'I4',ioerr,idate2(4,is),cdum)
            call read_line(line,indx,'I4',ioerr,idate2(5,is),cdum)
            if( ioerr.ne.0 ) then                                   
              write(*,'(a,5i5)') 'Error reading stop date from eq file:'
     .            ,  (idate1(i,is),i=1,5)
              stop
            endif    
c no check on illegal keyword since many possibilities foe eq_file
          endif 
        endif
      enddo  
      numsit = is


         
* Convert BREAK entries into RENAME entries

      if( filetype.eq.'break ' ) then

*       Sort the entries by station and time
        call sorteq( numsit,siteold,sitenew,idate1,idate2,comment )
        do is = 1, numsit
          newsite = .true.       
          nextnew = .true. 
          if( is.ne.1.and.siteold(is)(1:4).eq.siteold(is-1)(1:4) )
     .      newsite = .false.                                     
          if( newsite) then
c Start at _2PS to be consistant with itrf naming conventions.	  
            extchar = '2' 
          else          
            extchar = char(ichar(extchar)+1)
          endif
          sitenew(is)(6:6) = extchar
          if( is.ne.numsit.and.siteold(is+1)(1:4).eq.siteold(is)(1:4))
     .       nextnew = .false.
          if( nextnew ) then
            idate2(1,is) = 2100 
            idate2(2,is) = 1
            idate2(3,is) = 1
            idate2(4,is) = 0
            idate2(5,is) = 0 
          else
            do i=1,5
              idate2(i,is) = idate1(i,is+1)
            enddo
          endif 
        enddo       
      endif

      
* Write out the new file
         
      do is = 1,numsit
        if( filetype.eq.'break ' ) then
          write(2,'(a,2(1x,a8),2(i5,4i3),1x,a)') 
     .      ' rename' ,siteold(is)(1:4),sitenew(is)
     .      , (idate1(i,is),i=1,5),(idate2(i,is),i=1,5)
     .      , comment(is)(:trimlen(comment(is)))
     
        elseif( filetype.eq.'rename' ) then      
          write(2,'(a,a8,i5,4i3)') 
     .      'break' ,siteold(is),(idate1(i,is),i=1,5)
        endif
     
      enddo

      write(*,'(a)') 'Normal end'
      stop
      end

c------------------------------------------------------------

      Subroutine sorteq( numsit,siteold,sitenew,idate1,idate2,comment )
 
      implicit none

      integer*4 maxsit       
*     change this to match main routine
      parameter(maxsit=30000)

      integer*4 numsit,idate1(5,maxsit),idate2(5,maxsit)
     .        , index(maxsit),smallest_one,i,j,k,ik,ikp,is,ns
     .        , site_start(maxsit),last,itmp1(5,maxsit),iswap
     .        , itmp2(5,maxsit),jndex(maxsit)
                 
      character*4 sit4(maxsit),swap
      character*8 siteold(maxsit),sitenew(maxsit)
     .        , tmpold(maxsit),tmpnew(maxsit)
      character*256 comment(maxsit)
                                                       
      real*8 tk,tkp,dup_tol

      logical newsite
                              

*     Sort by name

      do i=1,numsit
        index(i) = i
        sit4(i) = siteold(i)(1:4)
      enddo
      call sort_snames(sit4,index,numsit)
      do i=1,numsit
        tmpold(i) = siteold(i)
        tmpnew(i) = sitenew(i)
        do j=1,5
          itmp1(j,i) = idate1(j,i)
          itmp2(j,i) = idate2(j,i)
        enddo
      enddo
      do i = 1,numsit
        siteold(i) = tmpold(index(i))
        sitenew(i) = tmpnew(index(i))
        do j=1,5
          idate1(j,i) = itmp1(j,index(i))
          idate2(j,i) = itmp2(j,index(i))
        enddo
      enddo           
cd      print *,'List after name sort  i index(i) numsit ',numsit
cd      do i=1,numsit   
cd        print *,i,sit4(i),siteold(i)
cd        if( i.gt.10 ) stop 
cd      enddo


* Within each site, sort by date
             
      is = 1
      do while( is.le.numsit)       
        ns = 1
        newsite = .false.
        do while(.not.newsite)        
cd          print *,'is sit4 sit4+1 newsite '
cd     .             ,is,sit4(is),sit4(is+1),newsite
          if( sit4(is+1).eq.sit4(is) ) then
             ns = ns + 1
             is = is + 1
cd             print *,'incremented is ns ',is,ns
          else
            if( ns.gt.1 ) then  
              newsite = .true.
cd              print *,'newsite caling sort_stimes is-ns+1',is-ns+1
              call sort_stimes(idate1(1,is-ns+1),jndex,ns)
cd             do i=1,ns
cd                print*,' ,i,jndex(i) ',i,jndex(i)
cd              enddo
cd              print *,'after sort_stimes is-ns+1,is',is-ns+1,is
              do i=1,ns
                do j=1,5
                  itmp1(j,i) =  idate1(j,is-ns+i)
                  itmp2(j,i) =  idate2(j,is-ns+i)
                enddo 
cd                print *,'i is-ns+i itmp1 ',i,is-ns+i,(itmp1(j,i),j=1,5)
              enddo
              do i=1,ns
                do j=1,5
                  idate1(j,is-ns+i) = itmp1(j,jndex(i))
                  idate2(j,is-ns+i) = itmp2(j,jndex(i))
                enddo    
cd               print *,'i is-ns+i,idate1 ',i,is-ns+i,(idate1(j,i),j=1,5)
              enddo           
            endif        
cd            print *,'after t-sort is-ns+1 is',is-ns+1,is
cd            do i= is-ns+1,is
cd              print *,sit4(i),(idate1(j,i),j=1,5),(idate2(j,i),j=1,5)
cd            enddo
            is = is + 1
          endif 
        enddo  
      enddo

cd      print *,'List after time sort  i index(i) '
cd      do i=1,numsit
cd        print *,i,index,i
cd      enddo 
cd      print *,'site start '
cd      do i=1,ns
cd        print *,site_start(i)
cd      enddo 

      return
      end

*----------------  

CTITLE SORT_SNAMES
 
      subroutine sort_snames( names, list, num)

      implicit none
 
*     This routine uses an exchange sort algormithm to sort
*     the list names into ascending alphabetical order.  
*     There are num values in list and names
 
*   num     - Number of values to be sorted
*   list(num)  - List to be sorted in to ascending order.
 
      integer*4 num, list(num)

*   names(num) -- Names to be sorted
      character*4 names(num)
 
* LOCAL VARIABLES
 
*   i,j     - Loop counters
*   smallest_one    - pointer to least value

      integer*4 i,j, smallest_one

*   swap   - swap string
      character*4  swap
 


****  Generate index list     
* This now done in mstinf
c      do i = 1, num
c         list(i) = i
c      end do
 
****  Start loop using exchange sort
      do i = 1, num
          smallest_one = i
          do j = i+1, num
              if( names(j).lt. names(smallest_one) ) then
                  smallest_one = j
              end if
          end do 
 
*****     See if we should swap
          if( smallest_one.gt. i ) then
              swap = names(smallest_one)
              names(smallest_one) = names(i)
              names(i) = swap
*             Now swap the index
              j = list(smallest_one)
              list(smallest_one) = list(i)
              list(i) = j

          end if
      end do
 
***** Thats all.  Now sorted in ascending order
      return
      end

CTITLE SORT_STIMES

      Subroutine sort_stimes( idate,list,num )

      implicit none
 
*     This routine uses an exchange sort algormithm to sort
*     the times into ascending .  
*     There are num values in list and names
 
*   num     - Number of values to be sorted
*   list(num)  - List to be sorted in to ascending order.
 
      integer*4 num, list(num), idate(5,num)

      real*8 t(num)
 
* LOCAL VARIABLES
 
*   i,j     - Loop counters
*   smallest_one    - pointer to least value

      integer*4 i,j, smallest_one

*   swap   - swap time
      real*8 swap
                     

cd      print *,'SORT_TIMES num idate ',num
      do i=1,num      
        t(i) = idate(1,i)+ (idate(2,i) + (idate(3,i) + 
     .         (idate(4,i) + 
     .         idate(5,i) /60.d0) /60.d0) /24.d0 )/365.25d0 
cd        print *,t(i),(idate(j,i),j=1,5)
      enddo

*  Generate index list     
      do i = 1, num
         list(i) = i
      end do
 
****  Start loop using exchange sort
      do i = 1, num
          smallest_one = i
          do j = i+1, num
            if( t(j).lt.t(smallest_one) ) then
               smallest_one = j
              end if
          end do 
*****     See if we should swap
          if( smallest_one.gt. i ) then
              swap = t(smallest_one)
              t(smallest_one) = t(i)
              t(i) = swap
*             Now swap the index
              j = list(smallest_one)
              list(smallest_one) = list(i)
              list(i) = j
          end if
      end do  

cd      print *,'SORT list '
cd      do i=1,num
cd        print *,list(i),t(i)
cd      enddo
 
***** Thats all.  Now sorted in ascending order
      return
      end

