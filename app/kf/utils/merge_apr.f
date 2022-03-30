      program merge_apr

      implicit none 
 
*     This program is used to merge two .apr file together, where all entries from the 
*     primary file (first) file are included in the outpout along with unique entries 
*     from the second file.
*  
*     Runstring:
*     % merge_apr [reference apr] [update apr] [output apr] <yr> <doy>
*
*     where [reference apr ] contains the primary or reference apr file name.
*           [secondary apr] is the apr file to be checked for unique entries wrt the reference file.
*           [output apr] is the output apr file name.
* 
* MOD TAH 161001: Added <yr> and <doy> to runstring so that EXTENDED entries can be computed
*      based on values in reference apr file.
*
*       max_sites   - Total maxiumim number of sites (this is so all
*                   - duplicates than be removed.
 
      integer*4  max_sites

* MOD TAH 140310: Updated to 16384 sites from 5000. 
      parameter ( max_sites         = 16384 )

* MOD TAH 161001: Include the NON-sec inclused

      include '../gen_util/nonsec_util.h'
 
* PROGRAM VARIABLES
 
*   i,j,k       - Loop counters
*   trimlen     - Length of string
*   ierr        - IOSTAT error
*   rcpar       - Get runstring
*   len_run     - Length of runstring
*   indx        - Pointers in string
      integer*4 i,j,is,k, trimlen, ierr, rcpar, len_run, indx, jndx

      integer*4 num_site  ! Number of sites found

      integer*4 yr, doy   ! Year and day-of-year for non-sec evaluation
      integer*4 uns       ! Unit number of secondary file.  101 or -1 if not
                          ! specified. 
      
*   found   - Used to indicate unique site has been found in secondary file
      logical found

      real*8 ejd, eyr, sec  ! jd for extended entries and sec tag
      real*8 ref_xyz(3,2,max_sites)   ! Reference apr file XYZ position and rate
      real*8 ref_ep(max_sites)        ! Reference Epoch (jd)
      real*8 decyrs                    ! Deciminal year
      real*8 coords(3), dsol(3)       ! XYZ mapped to date, adjustments of non-secular terms

*   line        - Line read from unput files
*   ref_apr     - Reference input apr file
*   sec_apr     - Secondary input apr file
*   out_apr     - Output apr_file
      character*256 line, ref_apr, sec_apr, out_apr
      character*4 cdum   ! read_line, multiread

*   name        - List of reference site names read from reference file
*   site        - secondary site site name read from line.  
      character*8 site, name(max_sites)
   
** rwk 090929:  comment out this header 
cd      write(*,100)
cd 100  format(/' MERGE_APR: Merges two GLOBK apr station',
cd     .        ' coordinate files',/)
 
      len_run = rcpar(1, ref_apr)
      if( len_run.le.0 ) then
          call proper_runstring('merge_apr.hlp','merge_apr',1)
      end if
 
*     Open the name files
      open(100,file=ref_apr, iostat=ierr, status='old')
      call report_error('IOSTAT',ierr,'open',ref_apr,1,'REF APR FILE')
 
*     Get the secondary apriori file name
* MOD TAH 161011: Allow this to be optional to allow reference apr to
*     be mapped to specific date.
      len_run = rcpar(2,sec_apr)
      if( len_run.gt.0 ) then
          uns = 101
          open(uns,file=sec_apr, iostat=ierr, status='old')
          call report_error('IOSTAT',ierr,'open',sec_apr,1,
     .                      'SEC APR FILE')
      else
          uns = -1
      end if
 
*     Get name of output
      len_run = rcpar(3,out_apr)
      if( len_run.le.0 ) then
          call proper_runstring('merge_apr.hlp','merge_apr',1)
      end if

*     See of YR and DOY passsed
      sec = 0
      ejd = 0
      len_run = rcpar(4,line)
      if( len_run.gt.0 ) then
          read(line,*) yr
      else
          yr = 0
      endif
     
      len_run = rcpar(5,line)
      if( len_run.gt.0 ) then
          read(line,*) doy
      else
          doy = 0
      endif
      if( yr.gt.0 .and. doy.gt.0 ) then
          call yds_to_jd( yr, doy, sec, ejd)
          call jd_to_decyrs( ejd, eyr) 
      endif 
 
*     Check to make sure output does not equal input
      if( ref_apr.eq.out_apr ) then
          write(*,120)
 120      format(' ** ERROR ** Output apr file can be same name as',
     .            ' primary input file')
          stop 'MERGE_APR: Output and Input Files same name'
      end if
*     Check to make sure output does not equal input
      if( sec_apr.eq.out_apr ) then
          write(*,130)
 130      format(' ** ERROR ** Output apr file can be same name as',
     .            ' secondary input file')
          stop 'MERGE_APR: Output and Input Files same name'
      end if      
      
      open(200,file=out_apr, iostat=ierr, status='unknown')
      call report_error('IOSTAT',ierr,'open',out_apr,1,'OUT APR FILE')
  
****  Read the reference apr file: All valid entries including comments are 
*     directly written to the output file, site names are saved.
      ierr = 0
      i = 0
      num_site = 0
      num_nonsut = 0
      name(:) = ' '

      do while ( ierr.eq.0 )
          read( 100, '(a)', iostat=ierr ) line

* MOD TAH 070402: removed use of exit (not suupported on Solaris).
*         Also changed the counter i to be incremented when names
*         added.
*	  if ( ierr .ne. 0 ) exit
*         Save the name of reference file site and write to output file. 
          if( ierr.eq.0 .and. trimlen(line).gt.0 ) then

* MOD TAH 161001: See if these are extended entries and do not write the 
*             coordinates yet since we may need to modified to get to values
*             of date.
              if( line(1:1).ne.' ' )  
     .                         write(200,'(a)') line(1:trimlen(line))
              if ( ejd.eq.0 ) then 
*                 Just pass the line though as was originally done.
                  if( line(1:1).eq.' ' ) 
     .                         write(200,'(a)') line(1:trimlen(line))
              endif

              if( line(1:1).eq.' ' ) then
*                Get the names of all the sites in the reference file 
                 indx = 1
                 call getword( line, site, indx)
                 call casefold(site)
                 if( site.ne.'EXTENDED' ) then

                     call get_cmd(site, name, num_site, is,
     .                            indx )

                     if( is.le.0 ) then
                         num_site = num_site + 1 
                         if( i.gt.max_sites ) then
                             write(*,220) max_sites
 220                         format('**ERROR** Max sites ',i5,
     .                              ' exceeded in reference file')
                             stop 'Too many sites'
                         endif
                         name(num_site) = site
                         is = num_site
                     endif

*                    Decode the rest of the line
                     call multiread(line, indx, 'R8',ierr, 
     .                                       ref_xyz(1,1,is),cdum, 6)
                      
                     call read_line(line, indx, 'R8',ierr, 
     .                                               decyrs,cdum)
                     call decyrs_to_jd(decyrs, ref_ep(is))
                  else
* MOD TAH 161001: EXTENDED ENTRY: If yr and doy passed evaluate, other wise
*                    just write
                     call decode_nonsut( line, indx, name,
     .                                   num_site, ref_xyz)
                  endif 
              endif
	  endif
      enddo
       
      write(*,320) trim(ref_apr),num_site, num_nonsut, yr, doy, ejd
 320  format('In ',a,' there are ',i5,' sites, and ',i5,
     .          ' EXTENDED entries',/,
     .       'Mapping to year ',i4,' DOY ',I3.3,' JD ',F12.1,/,
     .       'Not mapped if JD is 0.00')

****  If ejd > 0, map values
      if( ejd.gt.0 ) then
*        Map to data and evaluate extended terms.
         do i = 1, num_site      
            do j = 1,3
               coords(j) = ref_xyz(j,1,i) + 
     .               ref_xyz(j,2,i)*(ejd-ref_ep(i))/365.25d0
            end do
****        Now evaluate the non-secular terms (0 sets that logs should be 
*           computed)
            call eval_nonsut(i, ejd, num_nonsut, param_nonsut,
     .                 apr_val_nonsut, dsol, 0 )

            do j = 1,3
               coords(j) = coords(j) + dsol(j)
            end do

****        Now write new values
            call jd_to_decyrs( ref_ep(i), decyrs ) 
            write(200,360) name(i), coords, ref_xyz(:,2,i), eyr, 
     .               trim(ref_apr), decyrs
 360        format(1x,a8,1x,3(F16.5,1x),1x,3(F10.5),1x,F12.4,
     .             ' ! From ',a,' Mapped from ',F12.4)

         enddo
      endif

      write(200,410) trim(sec_apr)
 410  format(/,'* Coordinates copied from ',a)

****  Read the secondary apr file: All entries unique to this file are 
*     written directly to the output file     
      ierr = 0
      j = 0
      do while ( ierr.eq.0 .and. uns.gt.0 )
          j = j + 1
          read( 101, '(a)', iostat=ierr ) line
* MOD TAH 070402: removed exit (not suppurted Solaris).
*	  if ( ierr .ne. 0 ) exit
* SKIP comments in secondary file (they will grow it we don't
*      do this).
          if( ierr.eq.0 .and. trimlen(line).gt.0 .and.
     .        line(1:1).eq.' ' ) then
*            Get the nae of this site
	     indx = 1
             call getword( line, site, indx)
             call casefold(site)
*            Check if this site unique to the secondary file. 
             found = .false. 
             do k = 1,num_site
                if ( site .eq. name(k) ) then
                  found = .true.
                endif
 	     enddo 
*            Skip extended lines and RFERENCE_FRAME line
             if( ejd.eq.0 ) site = 'X'  ! Allows the extedned lines to 
*                                         be printed.
	     if ( .not.found .and. site.ne.'EXTENDED' .and.
     .           line(1:4).ne.'+REF' ) then
                write(200,'(a)') line(1:trimlen(line))
	     endif
          endif
      enddo	  	        
	      
****  Thats all
      close(100)
      if( uns.gt.0 ) close(uns)
      close(200)
      end
 
