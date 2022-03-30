      program getrel

      implicit none 
 
*     This program will get the relative velocities of a group of sites
*     and print the results in a format comparable to the standard NE
*     velocity output from GLOBK.
*
*     Runstring is:
*     % getrel <site> <input file> <max_sigma>
*
*     where <site>       is the name of the site the velocities are
*                        to be given relative to.
*           <input file> in the name of the file output from GLOUT or
*                        GLORG.  Output options with bits 3 and 4 on
*                        should be set (i.e., 12 + origin options)
*           <max_sigma>  is the size of the maximum sigma to be output.
*                        (Default is to output all)
 
      include '../includes/kalman_param.h'
 
* VARIABLES
 
*   ierr, jerr  - IOSTAT Error
*   rcpar       - Function to read runstring
*   lenrun      - Length of runstring
*   trimlen     - Length of string
*   iel, indx   - Pointers to positions in strings
*   ns          - Number of sites found
*   nref        - Site number for reference site
*   is          - Sign to be applied to difference (to reverse
*               - baseline)
*   lref        - length of Reference site string
 
      integer*4 ierr, jerr, rcpar, lenrun, trimlen, iel, indx, ns, 
     .          is, i, nref, lref
 
*   found       - General logical to indicate we have found
*               - what we are searching for
      logical found
 
*   maxsigma        - Maxiumum sigma to allow on output
*   values(20)  - Values read form station position string
*   posll(max_glb_sites)    - Long and Lat of each site.
 
      real*8 maxsigma, values(20), posll(2,max_glb_sites)
 
*   refsite     - Name of the reference site to be used
*   diffsite        - Site that refsite is differenced from
*   snames(max_glb_sites)   - Name of sites read from list
*   cd          - Dummy string for readline
*   sity1 and site2 - Names pulled from rate line
 
      character*8 refsite, diffsite, snames(max_glb_sites), cd, 
     .            site1, site2
 
*   inputfile   - Name of the input file
 
      character*128 inputfile
 
*   line        - Line read from input
 
 
      character*256 line
 
****  Read the runstring
      lenrun = rcpar(1,refsite)
      if( lenrun.le.0 ) then
          call proper_runstring('getrel.hlp', 'GETREL', 1)
      end if
      call casefold(refsite)
 
      lenrun = rcpar(2, inputfile)
      if( lenrun.le.0 ) then
          call proper_runstring('getrel.hlp', 'GETREL', 1)
      end if
 
****  See if sigma passed
      lenrun = rcpar(3, line)
      if( lenrun.gt.0 ) then
          read(line,*, iostat=ierr) maxsigma
          call report_error('IOSTAT',ierr,'decode', line, 0, 'getrel')
          if( ierr.ne.0 ) then
              call proper_runstring('getrel.hlp', 'GETREL', 1)
          end if
      else
          maxsigma = 1000.d0
      end if
 
****  Now we seem to be set, open up the input
      call open_lu(100, inputfile, ierr, 'old')
      call report_error('IOSTAT',ierr,'open',inputfile,1,'getrel')
 
***** Now start reading the file until we get to the 'Long.' line
      ierr = 0
      found = .false.
      do while( ierr.eq.0 .and. .not.found )
          read(100,'(a)', iostat=ierr) line
          indx = index(line,'Long.')
*         MOD TAH 110308: Allow a range of values here for different 
*         formats.
          if( indx.ge.3 .and.indx.le.5 ) then
              found = .true.
          end if
      end do

***** See if line not found
      if( .not.found ) then
         write(*,110) 
 110     format('GETREL: Velocity summary not found. Use VSUM org_opt ',
     .          'in output')
         call report_error('VSUM velocity summary not found',ierr,
     .        'search',inputfile,1,'getrel')
      end if
 
***** See if we reached EOF
      if ( ierr.ne.0 ) then
          call report_error('Unexpected IOSTAT',ierr,'read', inputfile,
     .        0,'getrel')
          call proper_runstring('getrel.hlp', 'GETREL', 1)
      end if
 
*     Skip a line and read the long and lats and names of the sites
      read(100,'(a)' ) line
 
      found = .false.
      ns = 0
      do while ( ierr.eq.0 .and. .not.found )
          read(100,'(a)', iostat=ierr ) line
* MOD TAH 050407: Globk Vel 5.11 Check VEL STATISTICS line
          if( trimlen(line).eq.0 .or.line(1:3).eq.'VEL' ) then
*                                 ! We have reached to end of
              found = .true.
*                                 ! station list
          end if
 
*                                 ! Decode
          if( .not.found ) then
 
*             Get rid of any NaN or * in line
              call sub_char(line,'NaN', '0.0')
              call sub_char(line,'*', '1')
              indx = 1
              ns = ns + 1
              call multiread(line, indx, 'R8',jerr, values, cd, 12)
 
*             Save the long and lat of site and get the name
              posll(1,ns) = values(1)
              posll(2,ns) = values(2)
              call getword(line, snames(ns), indx)
          end if
      end do

****  Now get out the reference site postion.
      indx = 1
      call get_cmd(refsite, snames, ns, nref , indx)

* MOD TAH 050407: Add explicit message if reference site not fould
      if( nref.le.0 ) then
         write(*,130) refsite, nref
  130    format('GETREL: Reference site ',a8,' Not found.',
     .          'GetCommand return ',i3)
         call report_error('Reference site not found',nref,'decod',
     .        refsite, 1,'getrel')
      end if
 
***** See if we reached EOF
      if ( ierr.ne.0 ) then
          call report_error('Unexpected IOSTAT',ierr,'read', inputfile,
     .        0,'getrel')
          call proper_runstring('getrel.hlp', 'GETREL', 1)
      end if
 
****  If we have not reach EOF then we should now be able to read
      found = .false.
      do while ( ierr.eq.0 .and. .not.found )
          read( 100, '(a)', iostat=ierr ) line
          indx = index(line,'GLOBK: BASELINE COMPONENT RATES')
          if( indx.gt.0 ) then
              found = .true.
          end if
      end do
 
***** See if we reached EOF
      if ( ierr.ne.0 ) then
          call report_error('Unexpected IOSTAT',ierr,'read', inputfile,
     .        0,'getrel')
          call proper_runstring('getrel.hlp', 'GETREL', 1)
      end if
 
***** Now skip the next three lines
      read(100,'(a)') line
      read(100,'(a)') line
      read(100,'(a)') line

*     Write out some header information.
      write(*,190) snames(nref), inputfile(1:trimlen(inputfile))
  190 format(/,'* GETREL: Velocites relative to ',a8,' from input ',a)
 
      write( * ,200)
  200 format('   Long.       Lat. ',7x,'E & N Rate ',3x,
     .       ' E & N Adj. ',2x,' E & N +-',1x,
     .       ' RHO ',5x,' H Rate  H adj.   +-',1x,'SITE',/,
     .       2x,' (deg)      (deg) ',2x,3(6x,'(mm/yr)'),15x,
     .          '(mm/yr)' )

*     Now write out the values for the reference site 
      write(*,300) (posll(i,nref),i=1,2),
     .             0.00,  0.00,  0.00, 0.00, 0.00, 0.00, 0.00,
     .             0.00,  0.00,  0.00, snames(nref),'*'
 
***** Now start scaning the list of baseline rates
      found = .false.
      lref = trimlen(refsite)
      do while ( .not.found .and. ierr.eq.0 )
          read(100,'(a)', iostat=ierr ) line
          if( trimlen(line).eq.0 .or. ierr.ne.0 ) then
              found = .true.
*                             ! Decode the line
          else
              site1 = line(2:9)
              site2 = line(11:18)
 
              is = 0
              if( site1(1:lref).eq.refsite(1:lref) ) then
                  is = -1
                  diffsite = site2
              end if
              if( site2(1:lref).eq.refsite(1:lref) ) then
                  is = 1
                  diffsite = site1
              end if
 
*                                     ! Get velocities
              if( is.ne.0 ) then
                  indx = 1
                  call get_cmd(diffsite, snames, ns, iel, indx)
                  if( iel.le.0 ) then 
                      is = 0
                  end if
              end if
              if( is.ne.0 ) then
                  indx = 19
                  call sub_char(line,'NaN', '0.0')
                  call sub_char(line,'*', '1')
                  call multiread(line, indx, 'R8',jerr, values,
     .                cd, 13)
 
*                 Fix the signs if need be
                  values(4) = is*values(4)
                  values(5) = is*values(5)
 
                  values(7) = is*values(7)
                  values(8) = is*values(8)
 
                  values(11) = is*values(11)
                  values(12) = is*values(12)
 
*                 Now write out the results
                  if( values( 9).lt. maxsigma .and. 
     .                values( 6).lt. maxsigma )      then
                      write(*,300) (posll(i,iel),i=1,2),
     .                    values(7), values(4),
     .                    values(8), values(5),
     .                    values(9), values(6),
     .                    values(10),
     .                    values(11), values(12), values(13),
     .                    diffsite
 300                  format(2(1x,f10.5),1x,6(1x,f7.2),1x,f6.3,2x,
     .                    3(1x,f7.2), 1x,a8,a1)
*                         ! Refererence site found
                  end if
              end if
*                         ! Line not blank
          end if
*                         ! Reading input
      end do
 
****  Thats all
      end
 
 
 
 
