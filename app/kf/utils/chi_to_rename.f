*     Program to read either an RNRP list from GLOBK or a rename file from TSFIT 
*     and use the chi increment from position or log estimates for renamed sites
*     to determine what renames, equates, and unequates are needed.  Write these
*     out in files to be used with GLOBK.   Optionally, an eq_file be read as input 
*     and written out with unneeded renames commented out (not yet implemented).

*     R.  King 21 September 2009

*     The program is command-line driven with the form:

*        chi2rename -org <org-file> -tsrn <ts_file> -chi <chi-value> -xps

*     Requires either
*         org_file   Glorg print file written with the RNRP option on;  OR 
*         ts_file    Rename file from tsfit

*     Optional:
*         eq_file    If given will read the file and write a new file with unneeded renames commented
*         chi-value  Maximum sqrt(chi^2) allowed for equate (default 3.0)
*         -xps       Include XPS renames from tsfit (default no)

      implicit none

      integer*4 maxsit
      parameter(maxsit=5000)

      logical xps,eof

      integer*4  luin,lueqin,lueqout,luequate,luunequate  
     .        ,  nr,len_run,ioerr,i
     .        ,  idate1(5),idate2(5)

      real*4 chimax,chi

* Function
      integer*4 nblen,rcpar,indx,idum

      character*3 intype            
      character*8 word,siteold,sitenew,cdum
      character*32 org_file,ts_file,infile,eqfin
     .          , eqfout,equatef,unequatef
      character*120 line,runstring

* Initialize
      luin = 1
      lueqin  = 2
      lueqout = 10
      luequate= 11
      luunequate = 12
      chimax = 3.0
      intype = '   '             
      org_file = ' ' 
      eqfin = ' '  
      eqfout = 'eq_file.out'
      equatef = 'equates.out'
      unequatef = 'unequates.out'
      xps = .false.

      
* Get the run-string arguments and open the files
    
      len_run = rcpar(1,runstring)
      if( len_run.le.0 ) then
         call proper_runstring('chi_to_rename.hlp','chi2rename',1)
      else
        nr = 0
        do while ( len_run.gt.0 )
          nr = nr + 1
          len_run = rcpar(nr, runstring )
          if( len_run.gt.0 ) then

*           See if we have new option
            word = runstring           
            call casefold(word)     
*           See if an org file is input          
            if( word(1:4).eq.'-ORG'  ) then
              nr = nr + 1 
              len_run = rcpar(nr,org_file)
              intype = 'org'
              print *,'Org file input not yet coded'
              stop
                          
*           See if a tsfit rename file is input
            else if( word(1:3).eq.'-TS' ) then
              nr = nr +1  
              len_run = rcpar(nr,ts_file)
              intype = 'tsf'

*           See if an input eq_file is to be used 
            else if( word(1:4).eq.'-EQF' ) then
              nr = nr + 1 
              len_run = rcpar(nr,eqfin) 
              print *,'Code for input eq_file not yet implemented'
              stop

*           See if the chi value is given 
            else if( word(1:4).eq.'-CHI ' ) then  
               nr = nr + 1
               len_run = rcpar(nr,runstring)
               read(runstring,*) chimax 

*           See if XPS commands from a tsfit rename file are to be included
            elseif( word(1:4).eq.'-XPS' ) then
               xps = .true.     

            endif 
          endif
        enddo 
      endif

       
*  Check for valid input

      if( intype.eq.'   ' ) then
        write(*,'(a)') 
     .   'Neither an org-file nor a tsfit rename file specified'
        stop
      elseif( intype.eq.'tsf'.and.org_file(1:1).ne.' ' ) then
        write(*,'(a)') 'Cannot use both org and tsfit files'
        stop
      endif

* Open the input file

      if( intype.eq.'org' ) then
        infile = org_file
      else
        infile = ts_file
      endif             
      open(unit=luin,file=infile,status='old',iostat=ioerr)
      if( ioerr.ne.0 ) then
         write(*,'(2a,i5)') 'Error opening input file ',infile,ioerr
           stop
      endif
      if( eqfin(1:1).ne.' ' ) then
        open(unit=lueqin,file='eqfin',status='old',iostat=ioerr)
        if( ioerr.ne.0 ) then
          write(*,'(2a,i5)') 'Error opening input eq_file ',eqfin,ioerr
          stop 
        endif
      endif
 
     
* Open the output files and write the headers

      open(unit=lueqout,file=eqfout,status='unknown',iostat=ioerr)
      if( ioerr.ne.0 ) then
        write(*,'(2a,i5)') "Error opening output eq_file ",eqfout,ioerr
         stop 
      endif          
      write(lueqout,'(2a,2x,a,f5.2 )')  
     .  "* Potential renames from CHI_TO_RENAME; input: "
     .    ,infile(1:nblen(infile)),"  chi max = ",chimax
             
      if( intype.eq.'org' ) then
        open(unit=luequate,file=equatef,status='unknown',iostat=ioerr)
        if( ioerr.ne.0 ) then
          write(*,'(2a,i5)') "Error opening output equate file "
     .      ,equatef,ioerr
          stop 
        endif 
        write(luequate,'(2a,2x,a,f5.2 )')  
     .    "* Potential equates from CHI_TO_RENAME; input: "
     .    , infile(1:nblen(infile)),"  chi max = ",chi
      endif
         
      open(unit=luunequate,file=unequatef,status='unknown',iostat=ioerr)
      if( ioerr.ne.0 ) then
        write(*,'(2a,i5)') "Error opening output unequate file "
     .     ,unequatef,ioerr
         stop 
      endif          
      write(luunequate,'(2a,2x,a,f5.2 )')  
     .   "* Potential unequates from CHI_TO_RENAME; input: "
     . , infile(1:nblen(infile)),"  chi max = ",chi
 

* Read the input file line-by-line and write the outputs
                        
      eof = .false.
      do while (.not.eof ) 
                     
        if( intype.eq.'org' ) then
c         this not yet coded
          continue
        
        else
          read( luin,'(a)',iostat=ioerr )  line
          if( ioerr.eq.-1 ) then
            eof = .true.
          elseif( line(1:1).ne.' ' ) then
            write(lueqout,'(a)') line
          else
            indx = 2
            call read_line(line,indx,'CH',ioerr,idum,word)
            call uppers(word)
            if( word(1:3).ne.'REN' ) then
              write(*,'(a)') 'TSFIT input file keyword not RENAME'
              write(*,'(a)') line
              stop
            endif                                         
            siteold = ' ' 
            sitenew = ' ' 
            call read_line(line,indx,'CH',ioerr,idum,siteold) 
            call read_line(line,indx,'CH',ioerr,idum,sitenew)
            call read_line(line,indx,'I4',ioerr,idate1(1),cdum)
            call read_line(line,indx,'I4',ioerr,idate1(2),cdum)
            call read_line(line,indx,'I4',ioerr,idate1(3),cdum)
            call read_line(line,indx,'I4',ioerr,idate1(4),cdum)
            call read_line(line,indx,'I4',ioerr,idate1(5),cdum)
            call read_line(line,indx,'I4',ioerr,idate2(1),cdum)
            call read_line(line,indx,'I4',ioerr,idate2(2),cdum)
            call read_line(line,indx,'I4',ioerr,idate2(3),cdum)
            call read_line(line,indx,'I4',ioerr,idate2(4),cdum)
            call read_line(line,indx,'I4',ioerr,idate2(5),cdum)
            call read_line(line,indx,'CH',ioerr,idum,cdum)  
            call read_line(line,indx,'CH',ioerr,idum,cdum)  
            call read_line(line,indx,'R4',ioerr,chi,cdum) 
            if( sitenew(6:8).eq.'XPS' ) then
              if( xps )  write(lueqout,'(a)') line  
            else
              if( chi.gt.chimax ) then
                write(lueqout,'(a,2(1x,a8,1x),2(i5,4i3,3x))') ' RENAME'
     .             ,siteold,sitenew,(idate1(i),i=1,5),(idate2(i),i=1,5)
               write(luunequate,'(a,a8,a)') ' UNEQUATE ',sitenew,' NDOT' 
               write(luunequate,'(a,a8,a)') ' UNEQUATE ',sitenew,' EDOT'
               write(luunequate,'(a,a8,a)') ' UNEQUATE ',sitenew,' UDOT'
              endif
            endif
          endif
        endif
      enddo
                                           
      if( intype.eq.'org' ) then
        write(*,'(2a,1x,a,1xa)') "Successfully wrote "
     .    ,eqfout(1:nblen(eqfout)),equatef(1:nblen(equatef))
     .    ,unequatef(1:nblen(eqfout))
      else
        write(*,'(2a,1xa)') "Successfully wrote "
     .     ,eqfout(1:nblen(eqfout)),unequatef(1:nblen(unequatef))
      endif

      stop
      end

                                                    


                                                                     

            
              
