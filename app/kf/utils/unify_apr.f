 
      program unify_apr
 
      implicit none 

*     This program is used to ensure that their are no duplicate station
*     coordinates and that sites which have different names but are the
*     same site have the same coordinates.
*
*     Runstring:
*     % unify_apr [Names file] [input apr] <output apr>
*
*     where [Names file] contains the primary site name followed by all
*                        site names that should have the same coordinate
*           [input apr] is a standard globk apriori site coordinate file
*           [output apr] is the output apr file name.
 
*     Special lines in the Name files changes mode between making position
*     and velocity equal to just velocity being the same (and visa versa).
*     The mode change is done using POS_MODE and VEL_MODE lines. (Must be
*     fully spelled out).  Program starts in POS_MODE.
 
*       max_prime   - Maximum number of primary sites
*       max_com_per_prime       - Max number of common site names with
*                   - with the prime sites.
*       max_sites   - Total maxiumim number of sites (this is so all
*                   - duplicates than be removed.
 
      integer*4 max_prime, max_com_per_prime, max_sites
 
      parameter ( max_prime         = 5000 )
      parameter ( max_com_per_prime =  100 )
      parameter ( max_sites         = 5000 )
 
* PROGRAM VARIABLES
 
*   i,j,k       - Loop counters
*   trimlen     - Length of string
*   ierr        - IOSTAT error
*   rcpar       - Get runstring
*   len_run     - Length of runstring
*   indx, jndx  - Pointers in string
*   iel         - Name number in list
*   np, no      - NUmber of prime and other sites
*   nc(max_prime)   - Number of common stations per prime
*   nv(max_prime)   - Number of common velocity sites per prime
*   mode        - Mode for reading names file.
*               - 1 -- Position and velocity same (default)
*               - 2 -- Velocity only same.
*               - Mode is selected with pos_mode or
*               - vel_mode line.
*   len_rem     - Length of the remainder line
 
      integer*4 i,j,k, trimlen, ierr, rcpar, len_run, indx, jndx,
     .     iel, jel, np, no, nc(max_prime), nv(max_prime), mode, 
     .     len_rem, jerr
 
*   mode_change    - Logical to indicate that this line changes
*           - the mode and should be be processed for site names
*   other_used(max_sites)       - Set true when an other site name
*           - have been output afiliated with a primary site.
*   found   - Used to indicate quatnity has been found
 
      logical mode_change, other_used(max_sites), found, use
 
*   pr_pos(3), ot_pos(3)    - Position for prime and other site when
*           - when only velocity is to be made the same
*   pr_vel(3), ot_vel(3)    - Velocity for prime and other site when
*           - when only velocity is to be made the same
*   pr_ep, ot_ep        - Epoch for primary and other sites
 
      real*8 pr_pos(3), ot_pos(3), pr_vel(3), ot_vel(3), pr_ep, ot_ep
      real*8 p2_pos(3), p2_vel(3), p2_ep

*   line        - Line read from unput files
*   name_file   - Names files
*   in_apr      - Input apr file
*   out_apr     - Output apr_file
*   prime_coord(max_prime)  - Line of apr file after the site name
*   other_coord(max_sites)  - Line of apr file after the site name
*   line_rem    - Remainder of line after the reference epoch.
 
      character*256 line, name_file, in_apr, out_apr,
     .    prime_coord(max_prime), other_coord(max_sites), line_rem

*   cd          - Dummy string for readline

      character*8 cd   
 
*   prime_names(max_prime)  - Names of prime sites
*   copos_names(max_com_per_prime, max_prime)   - Names of the common
*                           - sites associated with each prime site
*   covel_names(max_com_per_prime, max_prime)   - Names of sites which
*                           - have same velocity (position is taken
*                           - from the original data line)
*   other_names(max_sites)  - Names of the other sites.
*   name        - Site name read from line
 
 
      character*8 prime_names(max_prime),
     .    copos_names(max_com_per_prime, max_prime),
     .    covel_names(max_com_per_prime, max_prime),
     .    other_names(max_sites), name
 
****  START decode the runstring
      write(*,100)
 100  format(/' UNIFY_APR: Unifies GLOBK apr station',
     .        ' coordinate files',/)
 
      len_run = rcpar(1, name_file)
      if( len_run.le.0 ) then
          call proper_runstring('unify_apr.hlp','unify_apr',1)
      end if
 
*     Open the name files
      open(100,file=name_file, iostat=ierr, status='old')
      call report_error('IOSTAT',ierr,'open',name_file,1,'NAMES FILE')
 
*     Get the apriori file name
      len_run = rcpar(2,in_apr)
      if( len_run.le.0 ) then
          call proper_runstring('unify_apr.hlp','unify_apr',1)
      end if
      open(101,file=in_apr, iostat=ierr, status='old')
      call report_error('IOSTAT',ierr,'open',in_apr,1,'IN APR FILE')
 
*     Get name of output
      len_run = rcpar(3,out_apr)
      if( len_run.le.0 ) then
          call proper_runstring('unify_apr.hlp','unify_apr',1)
      end if
 
*     Check to make sure output does not equal input
      if( in_apr.eq.out_apr ) then
          write(*,120)
 120      format(' ** ERROR ** Output apr file can not be same name',
     .            ' as input file')
          stop 'UNIFY_APR: Output and Input Files same name'
      end if
      open(200,file=out_apr, iostat=ierr, status='unknown')
      call report_error('IOSTAT',ierr,'open',out_apr,1,'OUT APR FILE')
 
***** We have all the files open:  Now read the names file first to get
*     all the prime and common site names
      np = 0
      mode = 1
      do while ( ierr.eq.0 )
          read(100,'(a)', iostat=ierr) line
 
*         If this is not a comment then process.  Set mode change true
*         so that we will not process line if blank or a comment, in
*         addition to case when mode really changes.
          mode_change = .true.
          if( ierr.eq.0 .and. trimlen(line).gt.0 .and.
     .        line(1:1).eq.' ' ) then

* MOD TAH 060327  Remove any thing after ! or #
              iel = 0
              iel = index(line,'!')
              if( iel.eq.0 ) iel = index(line,'#' )
              if( iel.gt.0 ) line(iel:) = ' '
 
*             Get the prime site name
              indx = 1
              call getword(line, name, indx)
 
****          See if mode change line
              call casefold(name)
              mode_change = .false.
              if( name.eq.'POS_MODE' ) then
                  mode = 1
                  mode_change = .true.
              else if( name.eq.'VEL_MODE' ) then
                  mode = 2
                  mode_change = .true.
              end if
          end if
 
****      If we are not changing mode and this is valid line, continue
*         processing
          if( .not.mode_change ) then
 
*             See if duplicate
              jndx = 1
              call get_cmd(name, prime_names, -np, iel, jndx)
              if( iel.le.0 ) then
*                 New prime site so add to list
                  np = np + 1
                  if( np.gt. max_prime ) then
                      write(*,200) max_prime
 200                  format('** DISTASTER ** Limit of ',i4,' prime ',
     .                    'stations reached')
                      stop 'UNIFY_APR: Too many prime stations'
                  end if

                  prime_names(np) = name
                  iel = np
                  nc(np) = 0
                  nv(np) = 0
              end if
 
*             Now get the list of common site names
              do while ( trimlen(name).gt.0 )
                  call getword(line, name, indx)
                  call casefold(name)
                  if( trimlen(name).gt.0 .and. mode.eq.1 ) then
                      nc(iel) = nc(iel) + 1
                      if( nc(iel).gt.max_com_per_prime ) then
                          write(*,220) max_com_per_prime,
     .                                 prime_names(iel)
 220                      format(' ** DISTASTER ** Limit of ',i4,
     .                        ' common names for site ',a8,' exceeded')
                          stop 'UNIFY_APR: Too many common names'
                      end if
                      copos_names(nc(iel),iel) = name
                  else if( trimlen(name).gt.0 .and. mode.eq.2 ) then
* MOD TAH 060327: Make sure that a velocity site does not appear in the 
*                 position list (velocity will be set with the position
*                 and hence not needed)
                      use = .true.
                      do k = 1,nc(iel)
                         if( name.eq.copos_names(k,iel) ) then
                             write(*,230) name
 230                         format('** WARNING ** Velocity equate',
     .                              ' site ',a8,' appears in Position',
     .                              ' equate.  Ignoring')
                             use = .false.
                         end if
                      end do 
                      if( use ) then 
                         nv(iel) = nv(iel) + 1
                         if( nv(iel).gt.max_com_per_prime ) then
                             write(*,240) max_com_per_prime, 
     .                                    prime_names(iel)
 240                         format('** DISTASTER ** Limit of ',i4,
     .                        ' common names for site ',a8,' exceeded')
                             stop 'UNIFY_APR: Too many common names'
                         end if 
                         covel_names(nv(iel),iel) = name
                      endif
                  end if
              end do
          end if
      end do
 
****  Tell user what is happening
 
      write(*,300) np, name_file(1:trimlen(name_file))
 300  format(' There are ',i4,' prime site names in ',a,/,
     .        ' PRIME   : Common Position and velocity names')
      do i = 1, np
          prime_coord(i) = ' '
          if( nc(i).gt.0 ) then
              write(*,320) prime_names(i), (copos_names(j,i), j=1,nc(i))
 320          format(1x,a8,': ',50(a8,1x))
          end if
      end do
      write(*,340)
 340  format(/,' PRIME   : Common velocity site names')
      do i = 1,np
          if( nv(i).gt.0 ) then
              write(*,320) prime_names(i), (covel_names(j,i), j=1,nv(i))
          end if
      end do
 
***** Write headers to output file
 
      write(200,380) name_file(1:trimlen(name_file)),
     .            in_apr(1:trimlen(in_apr))
 380  format('*',/,
     .       '* Apriori coordinate file unified with UNIFY_APR from:',/,
     .       '* NAMES FILE : ',a,/,
     .       '* APR   FILE : ',a,/,
     .       '*',/,'* Comments from original APR file',/,
     .       '* -------------------------------')
 
****  Now read the input apr file:  All comments are directly written to
*     the output file, and the prime and other site coordinate lines are
*     saved.
      ierr = 0
      no = 0
      do while ( ierr.eq.0 )
          read( 101, '(a)', iostat=ierr ) line
          if( ierr.eq.0 .and. trimlen(line).gt.0 .and.
     .        line(1:1).eq.' ' ) then
              indx = 1
              call getword( line, name, indx)
              call casefold(name)
 
*             See if prime
              jndx = 1
*             Use negative number of sites to force look up on full
*             name
              call get_cmd(line, prime_names, -np, iel, jndx)
*                                 ! Prime station
              if( iel.gt.0 ) then
 
*                 See if we have entry already
                  if( trimlen(prime_coord(iel)).gt.0 ) then
                      write(*,400) prime_names(iel)
 400                  format(' Duplicate coordinates in apr file',
     .                       ' primary   site ',a8)

*                 Save only the first occurrence of the site.
                  else
                       prime_coord(iel) = line(jndx:)
                  end if
*                         ! Save as other site
              else
 
*                 See if we already have this site
                  jndx = 1
                  call get_cmd(line, other_names, -no, iel, jndx)
                  if( iel.gt.0 ) then
                      write(*,410) other_names(iel)
 410                  format(' Duplicate coordinates in apr file',
     .                       ' secondary site ',a8)
                  else
                      no = no + 1
                      iel = no
                      if( no.gt.max_sites ) then
                          write(*,420) max_sites
 420                      format('** DISASTER ** Limit of ',i4,' for',
     .                        ' other sites exceeded')
                          stop 'UNIFY_APR: Too many other sites'
                      end if

* MOD TAH 980210:  Only save the site information if this is the
*                     first occurence. 
                      other_names(iel) = name
                      other_coord(iel) = line(jndx:)
                  end if
              end if
          else if( ierr.eq.0 .and. trimlen(line).gt.0 ) then
              write(200,'(a)') line(1:trimlen(line))
          end if
      end do
 
*     Mark all other sites are un-used.  When they have been output they
*     will be marked as used.
      do i = 1, no
          other_used(i) = .false.
      end do
 
****  Now we finished reading the input.  Now ouput, in order,
*     the prime sites with their position and velocity equivalents,
*     followed by the orther sites that have not been used.
 
      write(200,500)
 500  format('*',/,'* PRIMARY SITES and their affiliates',/,
     .       '* ---------------------------------')
 
      do i = 1, np
          if( trimlen(prime_coord(i)).eq.0 ) then
 
*             Prime site was not in file, try the sites in the
*             common position and velocity names.  Use one of these if
*             necessary
              found = .false.
              j = 0
              do while ( .not.found .and. j.lt. nc(i) )
                  j = j + 1
                  jndx = 1
                  call get_cmd(copos_names(j,i), other_names, -no,
     .                iel, jndx )
                  if( iel.gt.0 ) then
                      found = .true.
                      prime_coord(i) = other_coord(iel)
                  end if
              end do
          end if
 
*         Try again now with first common site name used.
          if( trimlen(prime_coord(i)).gt.0 ) then
 
*             Primary site found, so start output of primary site and the
*             affiliated sites:
              jndx = 1
              call multiread(prime_coord(i), jndx, 'R8', jerr,
     .                       pr_pos, cd, 3 )
              call multiread(prime_coord(i), jndx, 'R8', jerr,
     .                       pr_vel, cd, 3 )
              call read_line(prime_coord(i), jndx, 'R8', jerr,
     .                       pr_ep , cd )
              line_rem = prime_coord(i)(jndx:trimlen(prime_coord(i)))
              len_rem = max(1,trimlen(line_rem))
              write(200,520) prime_names(i), pr_pos, pr_vel, 
     .                pr_ep, line_rem(1:len_rem)
 520          format(1x,a8,3(1x,f14.4),1x,3(f9.5,1x),1x,f9.4,1x,a) 
 
*****         Now write the affliated positon and velocity sites
              do j = 1, nc(i)
                  jndx = 1
                  call get_cmd(copos_names(j,i), other_names, -no,
     .                    iel, jndx)
                  if( iel.gt.0 ) other_used(iel) = .true.
                  write(200,520) copos_names(j,i), pr_pos, pr_vel,
     .                pr_ep, line_rem(1:len_rem)
* MOD TAH 060327: See if this site is a prime for another group of 
*                 sites
                  do k = i+1, np
                     if( prime_names(k).eq.copos_names(j,i) ) then
                        write(*,530) k, prime_names(k), prime_names(i)
 530                    format('Updating Primary station ',i4,1x, a8,
     .                         ' with position velocity from ',a8)
                        write(prime_coord(k),535) pr_pos, pr_vel,
     .                       pr_ep, line_rem(1:len_rem)
 535                    format(3(1x,f14.4),1x,3(f9.5,1x),1x,f9.4,1x,a)
                     endif
                  end do
              end do
 
*****         Now do the sites with common velocity but different
*             positions
              do j = 1, nv(i)

* MOD TAH 060327: See if this site is a primary one
                  jndx = 1
                  call get_cmd(covel_names(j,i), prime_names, -no,
     .                    iel, jndx)
                  if( iel.gt.i ) then  
*                     This site that we equating velocities at is also
*                     a primary site for another purpose.  Apply the
*                     position equates first
                      if( trimlen(prime_coord(iel)).gt.0 ) then
 
*                        Primary site found, so start output of
*                        primary site and the affiliated sites:
                         jndx = 1
                         call multiread(prime_coord(iel), jndx, 'R8', 
     .                       jerr, p2_pos, cd, 3 )
                         call multiread(prime_coord(iel), jndx, 'R8', 
     .                       jerr, p2_vel, cd, 3 )
                         call read_line(prime_coord(iel), jndx, 'R8',
     .                       jerr, p2_ep , cd )
                         line_rem = prime_coord(iel)
     .                           (jndx:trimlen(prime_coord(iel)))
                         len_rem = max(1,trimlen(line_rem))
*                        Update the Prime coordinate velocity.  Note
                         write(prime_coord(iel),535) p2_pos, pr_vel,
     .                        p2_ep, line_rem(1:len_rem)
 
*****                    Now update the affliated positon and velocity sites
                         jel = 0
                         do k = 1, nc(iel)
                            jndx = 1
                            call get_cmd(copos_names(k,jel), 
     .                           other_names, -no, iel, jndx)
                            if( jel.gt.0 ) then 
                                other_coord(jel) = prime_coord(iel)
                            endif
                         enddo
                      endif
                  endif
*****             Now see if velocity is another site

                  jndx = 1
                  call get_cmd(covel_names(j,i), other_names, -no,
     .                    iel, jndx)
                  if( iel.gt.0 ) then
 
*                     In these cases we must have found the name so that
*                     we can output position
                      other_used(iel) = .true.
 
*                     Read the position and velocity from the other site
                      jndx = 1
                      call multiread(other_coord(iel), jndx, 'R8', jerr,
     .                               ot_pos, cd, 3 )
                      call multiread(other_coord(iel), jndx, 'R8', jerr,
     .                               ot_vel, cd, 3 )
                      call read_line(other_coord(iel), jndx, 'R8', jerr,
     .                               ot_ep , cd )
                      line_rem = other_coord(iel)(jndx:
     .                                     trimlen(other_coord(iel)))
                      len_rem = max(1,trimlen(line_rem))
                      write(200,520) covel_names(j,i), ot_pos, pr_vel, 
     .                        ot_ep, line_rem(1:len_rem)
 540                  format(1x,a8,1x,3(F14.4,1x),1x,3(F8.4,1x),
     .                       1x,F10.4)
                  end if
              end do
          else
              write(*,560) prime_names(i)
 560          format('**WARNING** No entry in apr file for ',a8,
     .            ' or its affiliated sites')
          end if
      end do
 
***** Now write other sites that we not menetioned.
      write(200,610)
 610  format(/,'*',
     .       /,'* OTHER sites not associated with primary sites',/,
     .         '* ---------------------------------------------',/)
 
      do i = 1, no
          if( .not.other_used(i) ) then
              jndx = 1
              call multiread(other_coord(i), jndx, 'R8', jerr,
     .                       ot_pos, cd, 3 )
              call multiread(other_coord(i), jndx, 'R8', jerr,
     .                       ot_vel, cd, 3 )
              call read_line(other_coord(i), jndx, 'R8', jerr,
     .                       ot_ep , cd )
              line_rem = other_coord(i)(jndx:trimlen(other_coord(i)))
              len_rem = max(1,trimlen(line_rem))
              write(200,520 ) other_names(i), ot_pos, ot_vel, ot_ep,
     .            line_rem(1:len_rem)
 620          format(1x,a8,1x,a)
          end if
      end do
 
****  Thats all
      close(100)
      close(101)
      close(200)
      end
 
