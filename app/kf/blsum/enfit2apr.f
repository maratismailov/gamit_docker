      Program enfit2apr 

*      Routine to generate extended apr-file lines from the enfit summary file
*      for sites affected by earthquakes  
        
*      R. King 25 June 2002

       implicit none

                    
       character*3 type 
       character*8 site,arg 
       character*112 line(3)
       character*256 enfit_file,apr_file
          
       integer*4 len_run,rcpar,iymdhm1(5),iymdhm2(5),ierr,trimlen,i

       real*4 sigcut, snrcut, tau1, tau2, expo1(3), expo1_sig(3)
     .      , expo2(3), expo2_sig(3)

      logical found_expo1, found_expo2, eof

****  Get the runstring 

      len_run = rcpar(1,enfit_file)
      if( len_run.le.0 ) then  
         write(*,'(2a,/)') 'Write extended-format apr-file lines from'
     .            ,' an enfit summary file'
         write(*,'(a,/)') 'Runstring: enfit2apr file type <sig> <snr> ' 
         write(*,'(2a,/)') ' only type EXP is currently coded; ' 
     .                  , ' can extract up to 2'
         write(*,'(a)') '   sig is maximum sig (mm) (default 100 mm) ' 
         write(*,'(a)') '   snr is minimum value/sig (default 0)'

         stop
      end if
      len_run = rcpar(2,type(1:3))
      call casefold(type)
      if( len_run.le.0 ) then
         write(*,'(a)') 'Missing file name from runstring'
         stop
      endif

      len_run = rcpar(3,arg) 
      if( len_run.le.0 ) then
         sigcut = 100.   
      else
         read(arg,'(f8.0)') sigcut
      endif  

      len_run = rcpar(4,arg) 
      if( len_run.le.0 ) then
         snrcut = 0.    
      else
         read(arg,'(f8.0)') snrcut  
      endif


****  Open the files   
 
      open(101, file = enfit_file, iostat=ierr, status='old') 
      if( ierr.eq.0 ) then
         write(*,'(2a)') 'Opened input file ',enfit_file
      else
         write(*,'(2a)') 'Error opening input file ',enfit_file   
         stop
      endif
      apr_file = enfit_file(1:trimlen(enfit_file))//'.apr'
      open(102,file=apr_file,iostat=ierr,status='unknown')
      if( ierr.eq.0 ) then
         write(*,'(2a)') 'Opened output file ',apr_file
      else
         write(*,'(2a)') 'Error opening output file ',apr_file
         stop
      endif

      
****  Check for exponentials in headers, which can have up to 5 lines and then 
***** exactly form (with 'Site in column 1') 
*
*Site    ENU #         Offset at               Rate       Exponential 1  Exponential 2  wrms  nrms  line #
*                    1992 06 28 00 00 +-     mm/yr    +-  Tau   10.0 daysTau  200.0 days
*                                                         1992 06 28 00 01992 06 28 00 0 mm
      if( type.eq.'EXP') then    
        found_expo1 = .false. 
        found_expo2 = .false.
        do while (.not.found_expo1 )
          read(101,'(a)',iostat=ierr) line(1)   
          if( ierr.ne.0 ) then
             write(*,'(a)') 'Error reading first header line'
             stop
          endif
          if( line(1)(1:4).eq.'Site' ) then 
            read(101,'(a)') line(2) 
            read(101,'(a)') line(3)  
            if( index(line(1),'Exponential 1').ne.0 ) found_expo1=.true.
            if( index(line(1),'Exponential 2').ne.0 ) found_expo2=.true.   
            if( .not.found_expo1 ) then   
               write(*,*) 'No exponential found in file'
               stop
            else
               read(line(2),'(60x,f7.0)',iostat=ierr) tau1   
               if( ierr.ne.0 ) then
                  write(*,'(a)') 'Error reading tau1 '
                  stop
               endif 
*              note odd formating of HHMM 
               read(line(3),'(57x,i4,1x,i2,1x,i2,1x,i2,i2)',iostat=ierr) 
     .            (iymdhm1(i),i=1,5) 
               if( ierr.ne.0 ) then
                  write(*,'(a)') 'Error reading date for expo 1'
                  stop  
               endif
               if( found_expo2 ) then
                 read(line(2),'(75x,f7.0)',iostat=ierr) tau2 
                 if( ierr.ne.0 ) then
                    write(*,'(a)') 'Error reading tau2'
                    stop 
                 endif
                 read(line(3),'(72x,i4,1x,i2,1x,i2,1x,i2,i2)'
     .                ,iostat=ierr) (iymdhm2(i),i=1,5)     
                 if( ierr.ne.0 ) then
                    write(*,'(a)') 'Error reading date for expo 2'   
                    stop
                 endif
               endif 
            endif  
          endif
        enddo   
      else
        write(*,'(a,a3,a)') 'Input type ',type,' not yet allowed '
        stop
      endif
                          
****  Write some comment lines at the top of the output file

      write(102,'(2a)') '* Extended apr-file output from enfit file '
     .                ,enfit_file 
      write(102,'(a,f5.0,a,f3.1)')  
     .      '*   Maximum sigma (mm) ',sigcut,'   Minimum SNR ',snrcut 
      write(102,'(a,23x,a,13x,a,20x,a)') '*','Date'
     .                           ,'Tau (days)    NEU (m)','NEW sig (m)'  
      write(102,'(a)') '*'
        
****  Now read over the full file and extract the values for each station, assuming
*     that N, E, and U appear on successive lines
           
      eof = .false.
      do while (.not.eof )  
        do i=1,3
          read(101,'(a)',iostat=ierr) line(i) 
          if( ierr.eq.-1 ) then
            eof = .true. 
            goto 99
          elseif ( ierr.ne.0 ) then
            write(*,'(a)') 'Error reading station line'
            stop
          else
            read(line(i),'(a8,46x,1x,2f7.0)',iostat=ierr) 
     .          site,expo1(i),expo1_sig(i) 
            if( ierr.ne.0 ) then
                write(*,'(a,a8)') 'Error reading expo1 for site ',site 
                stop
            endif
            if( found_expo2 ) 
     .        read(line(i),'(70x,2f7.0)',iostat=ierr) 
     .           expo2(i),expo2_sig(i)   
              if( ierr.ne.0 ) then
                 write(*,'(a,a8)') 'Error reading expo2 for site ',site
                 stop
              endif
          endif 
        enddo 
*       write the values if at least one of the components has a significant exponential  
*       but zero-out the vertical if it fails to meet the SNR test    
              
        if( ( expo1_sig(1).le.sigcut .or. 
     .        expo1_sig(2).le.sigcut .or.   
     .        expo1_sig(3).le.sigcut ) .and.
     .       ( abs(expo1(1))/expo1_sig(1).ge.snrcut .or.
     .         abs(expo1(2))/expo1_sig(2).ge.snrcut .or.
     .         abs(expo1(3))/expo1_sig(3).ge.snrcut ) ) then  
           if( abs(expo1(3))/expo1_sig(3).le.snrcut ) expo1(3) = 0.
           write(102,'(a,a8,a,i4,4i3,f8.1,2x,3f8.4,4x,3f8.4)') 
     .          ' EXTENDED ',site,' EXP ',(iymdhm1(i),i=1,5),tau1
     .         ,(.001*expo1(i),i=1,3),(.001*expo1_sig(i),i=1,3)  
        endif
        if( found_expo2 ) then 
          if( ( expo2_sig(1).le.sigcut .or. 
     .          expo2_sig(2).le.sigcut .or.   
     .          expo2_sig(3).le.sigcut ) .and.
     .        ( abs(expo2(1))/expo2_sig(1).ge.snrcut .or.
     .          abs(expo2(2))/expo2_sig(2).ge.snrcut .or.
     .          abs(expo2(3))/expo2_sig(3).ge.snrcut ) ) then  
             if( abs(expo2(3))/expo2_sig(3).le.snrcut ) expo2(3) = 0.
             write(102,'(a,a8,a,i4,4i3,f8.1,2x,3f8.4,4x,3f8.4)') 
     .            ' EXTENDED ',site,' EXP ',(iymdhm2(i),i=1,5),tau2
     .           ,(.001*expo2(i),i=1,3),(.001*expo2_sig(i),i=1,3)
  
          endif
        endif

      enddo
           

   99 write(*,'(a)') 'Normal end of ENFIT2APR'
      stop
      end


