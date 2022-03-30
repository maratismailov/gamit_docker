      program tssum

      implicit none

*     Program to generate PBO time series files
*
*     Usages:
*     tssum <dir> <prod_id> <-R/A/C><K> <list of .org files>
*     
*     where <dir>  -- directory to put the time series in.  This directory
*      will be checked for existing files and these will be appended to
*           <prod_id> is product id with form: pbo.final_frame.  Characters
*                  5:9 will used for time series type (normally rapid or final)
*           <-R> causes exisiting time series files to be ignored
*           <-C> will keep exising entries but recompute the reference position
*           <-A> appends (any string will do this)
*           Adding K will keep lines that read with errors.
*           <list of .org files> glred/globk org files with output option PBOP
*      output option set.
*     
*     PROD_ID types valid in PBO:
*     Format <cen>.<series>_<frame>+<type>
*     <cen>     - Center 3 characters (bsl cwu pbo mit) or special product
*                 code e.g., aug for Augustine Volcano, aku for Akutan volcano
*     <series>  - Orbit series: final or rapid
*     <frame>   - Frame type: loose or frame or refernce name eg snf01
*     <type>    - Optional type.  If not given series name is used.  Additional
*                 type put in the final series is supplemental run (suppl).
*     
*     Stsndard PROD_ID
*     pbo.rapid_frame
*     pbo.final_frame
*     pbo.final_frame+suppl

 

      include 'tssum.h'

      integer*4 len_run, nr, ierr, ns, rcpar, trimlen
      character*4  institute ! Institute to generate product ID
      character*64 runstring 

****  OK, Read the runsting (for the moment just by position)
      len_run = rcpar(1,tsdir)
      ln_tsdir = len_run
      if( len_run.eq.0 ) then
        call proper_runstring('tssum.hlp','tssum',1)
      endif

*     Get prodct ID
      len_run = rcpar(2,prod_id)
* MOD TAH 101038: If prod_id is empty; then generate nama
      if( len_run.eq. 0 ) then
          call getenv('INSTITUTE',institute)
          call caseunfold(institute)
          prod_id = institute(1:3) // '.unkwn_unkwn' 
      end if

* MOD TAH 051129: Extract type from prod_id: Format pbo.final_frame+suppl
      ts_ref_type = prod_id(5:9)
      if( len_run.gt.15 ) then
          ts_ref_type = prod_id(17:21)
          prod_id = prod_id(1:15)
      endif

      ln_prod_id = trimlen(prod_id)
      ts_refresh = .false.
      ts_recomp = .false.
      ts_keep = .false.

      len_run = rcpar(3,runstring)
      if( runstring(1:2).eq.'-R' ) then
          ts_refresh = .true.   ! Creates new ts_files
      elseif( runstring(1:2).eq.'-C' ) then    
          ts_recomp = .true.    ! Retains ts_files but re-comps refpos.
      endif
* MOD TAH 140808: see if K is appended which will keep bad entries
      if( index(runstring,'K').gt.0 ) ts_keep = .true.
      
***** OK Loop over input files
      num_ent = 0
      num_site = 0
      num_code = 0 
      nr = 3
      len_run = 1
      tsprog = 'tssum'

      reference_frame = ' '

      do while ( len_run.gt.0 )
         
         nr = nr + 1
         len_run = rcpar(nr,in_file)
         if( len_run.gt.0 ) then
             call read_in
         endif
      end do

      call systime( date_rel, sec_rel) 

****  OK, now loop over all the codes
      do ns = 1, num_code

*        Generate TS file name
         ts_file = tsdir(1:ln_tsdir) // '/' // in_code(ns) // '.' // 
     .             prod_id(1:ln_prod_id) // '.pos'

*****    Try to open file unless we are refreshing
         if( .not.ts_refresh ) then
             open(200,file=ts_file,iostat=ierr,status='old')
             if( ierr.ne.0 ) then
                 new_ts = .true.
             else
                 new_ts = .false.
             end if
         else
             new_ts = .true.
         end if

*****    See if we should read old series
         num_ts = 0
         if( .not.new_ts ) then
             call read_tspos(200)
             if( ts_recomp ) then
                 call gen_refdata(ns)
             endif
         else
             call gen_refdata(ns)
         endif

* MOD TAH 080325: Before merging correct input series for any
*        East jumps due to latitude shifts

         call remove_ejmp(ns)

*****    OK, merge new entries with old
         call merge_ts(ns)

*****    Now close the infile and re-open to write
         close(200)
         open(200,file=ts_file,iostat=ierr,status='unknown')
         call write_ts(200, ns)
      end do

***** Thats all
      end

CITTLE READ_IN

      subroutine read_in

      implicit none

      include 'tssum.h'


      character*8 gsite_name 
      character*16 gsite_full

      integer*4 date(5), ierr, jerr, j, ns, cs, ne
      integer*4 trimlen, indx
      real*8  gmjd, pos_xyz_fin(3),xyz_std(6), unc_geod(3),llu_std(3),
     .        pos_neu_fin(3), neu_std(6), sec 

      character*512 line

****  Open the input file
      open(100,file=in_file, iostat=ierr, status='old')
      call report_error('IOSTAT',ierr,'open',in_file,0,'tssum')

      sec = 0.d0

      do while ( ierr.eq.0 )
         read(100,'(a)', iostat=ierr) line
* MOD TAH 080724: See if reference frame line
          if( line(1:26).eq.' Reference Frame        : ' .and.
     .        ierr.eq.0                            ) then
              reference_frame = line(27:)
c             print *,'RefFrame ',reference_frame, num_site, num_code,
c     .                ne, gmjd
          endif

* MOD TAH 060921: Replaced XPS with X since eq-renames may change the name
         if( ierr.eq.0 .and. line(1:4).eq.'pbo.' .and.
     .       line(11:11).ne.'X' ) then
            read(line,100,iostat=jerr) gsite_name, gsite_full,
     .                  date, gmjd,
     .                 (pos_xyz_fin(j), j=1,3), (xyz_std(j),j=1,6),
     .                 (unc_geod(j),j=1,3), (llu_std(j),j=1,3),
     .                 (pos_neu_fin(j),j=1,3), (neu_std(j),j=1,6)
 100        format(5x,a8,1x,a16,1x,i5,4(1x,i2.2),1x,F10.4,
     .                  1x,3F15.5,3F8.5,3F7.3,3x,
     .                  2F16.10,1x,F10.5,1x,2F8.1,1x,F10.5,3x,
     .                  2F15.5,1x,F10.5,3F8.5,3F7.3)
            call report_error('IOSTAT',jerr,'read',line,0,'read_in')
c100        format('pbo. ',a8,1x,a16,1x,i5,4(1x,i2.2),1x,F10.4,
c    .                  1x,3F15.5,3F8.5,3F7.3,' | ',
c    .                  2F16.10,1x,F10.5,1x,2F8.1,1x,F10.5,' | ',
c    .                  2F15.5,1x,F10.5,3F8.5,3F7.3)
            call ymdhms_to_mjd(date, sec, gmjd)

            if( jerr.eq.0 ) then

*             See if we can match site name
              indx = 0
              call get_cmd(gsite_name,in_site,num_site,ns,indx)
              if( ns.le.0 ) then
                  num_site = num_site + 1
                  ns = num_site
                  in_site(ns) = gsite_name
                  if( num_site.gt.max_site ) then 
                      call report_stat('FATAL','TSSUM','read_in','',
     .                     'Too many sites',max_site)
                  endif
              end if
              indx = 0
              call get_cmd(gsite_name(1:4),in_code,num_code,cs,indx)
              if( cs.le.0 ) then
                  num_code = num_code + 1
                  cs = num_code
                  in_full(cs) = gsite_full
                  in_code(cs) = gsite_name(1:4)
              end if

*****         OK save this entry
              num_ent = num_ent + 1
              if( num_ent.gt.max_ent ) then
                  call report_stat('FATAL','tssum','read_in','',
     .                'Too many enties Max ',max_ent)
              endif
              ne = num_ent
              in_ns(ne)  = ns
              in_cs(ne)  = cs
              in_mjd(ne) = gmjd
              do j = 1,3
                 in_xyz(j,ne) = pos_xyz_fin(j)
                 in_neu(j,ne) = pos_neu_fin(j)
                 in_llu(j,ne) = unc_geod(j)
              end do
              do j = 1,6
                 in_xyz_std(j,ne) = xyz_std(j)
                 in_neu_std(j,ne) = neu_std(j)
              end do
            end if
         end if
      end do

****  Tell user were we are
      write(*,110) in_file(1:trimlen(in_file)), num_site, num_code, 
     .             num_ent
 110  format('File: ',a,' Sites ',i5,' Codes ',i5,' Entries ',i10)
      return
      end 

