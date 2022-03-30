      program global_sel 

      implicit none 

*      This program will read an .apr file and given a reference station
*      will select a network on n stations which maxiumizes the separation
*      of the sites.

* MOD TAH 061104: Added feature where number of sites is negative which will
*      create this number of networks for all sites in apr file.  Used in 
*      forming clock networks.
*
* MOD SCM 061212: added feature to specify network size, the number of networks 
*      and number of reference sites to be used in in each network.

* PARAMETERS
* max_site  -- Maximum number of sites allowed
* max_bl    -- maximum number of baselines
* max_nets  -- maximum number of networks 
 
      integer*4 max_site, max_bl, max_nets

      parameter ( max_site = 4048)
      parameter ( max_nets = 50 )
      parameter ( max_bl   = max_site*(max_site+1)/2 )

* MAIN DECLARATIONS FOR THE PROGRAM
* 
* use_list(max_site)    -- List of sites to use
* all_list(max_site)    -- List of all sites.  Value set to -1
*                          when already used.
* tot_sites             -- Total number of sites
* sel_num               -- Number of sites selected by user
* use_num               -- Number of sites in use list
* num_refs              -- Number of reference sites total
* net_refs              -- Number of user defined reference sites / network
* num_avf               -- Number of files given with lists of
*                          available files.
* num_nets              -- Number of networks to generate

      integer*4 use_list(max_nets,max_site), all_list(max_site), 
     .          tot_sites, use_num, num_refs, num_avf,
     .          num_nets, sel_num, net_refs(max_nets), nrefs

* site_xyz(3,max_site)  -- Site coordinates
* site_vel(3,max_site)  -- Site velocity (used only for output)
* site_ep(max_site)     -- Site epoch (output only).
* blens(max_bl)         -- Lengths between all the sites 

      real*8 site_xyz(3,max_site), site_vel(3,max_site), 
     .       site_ep(max_site), blens(max_bl)  

* site_names(max_site)          -- Names of sites
* MOD TAH 200421: Made use_names 4-char to make sure we match in apr file.
* use_names(max_site)           -- Names of sites to use
* ref_names(max_site)  -- Names of reference sites (at least one must be
*                         given)
* name                 -- Generic name for site

      character*8 site_names(max_site), ref_names(max_site)
      character*4 use_names(max_site)   ! 4-char ID
      character*8 name

* apr_file   -- Name of apriori coordinate file
* runstring  -- Command line arguments
* line       -- Line read from file
* avail_file -- File containing a list of available sites

      character*256 apr_file, line, avail_file, test, message
      character*1024 net
      character*4096 runstring

* RUNNING DECLARATIONS
* len_run      -- Length of runstring
* ierr, jerr   -- IOSTAT errors
* indx, jndx   -- Pointers in strings
* trimlen      -- returns length of string
* i,j,k,l,m,n  -- Loop counters
* nb           -- Function to return baseline number
* iref         -- Reference site number
* rcpar        -- Returns runstring
* iel          -- Return from get_cmd (station number if it exist) 
* max_indx     -- Site with largest separation
* mode         -- mode in witch program called (mode=1,2,3)

      integer*4 len_run, ierr, jerr, indx, jndx, trimlen, i,j,k,l,m,n,
     .          nb, iref, rcpar, iel, max_indx, jel, mode
      integer*4 ir ! Loop counter for 4-char ref site lookup
      integer*4 num_missed ! Count of missed sites when avail list is short.

* max_dist, av_dist -- Max of average distance and average distance
* sq_dist           -- Sum of squares to distance to maxiumuze variance

      real*8 max_dist, av_dist, sq_dist, mean_dist, ri, rj, idotj

* use  -- Logical which is set false if site with 100 km of an
*         already used site.
* use_OK(max_site) -- Logical set true when ever a sites appears
*         in the available list being read.
      logical use, use_OK(max_site)
      logical use_found(max_site)  ! Set true when use site has been found
                                   ! in apr file (multiple _xPS names)

* MOD TAH 200419: Added skip to skip over extended lines
      logical skip   ! Set true for extended lines

* cd   -- Dummy string 
      character*8 cd



****  Read the runstring and print help
      mode = 1
      j = 0 
      runstring = ' '
      num_refs = 0
      sq_dist = 0.d0
      len_run = rcpar(1, runstring)
      if ( len_run.eq.0 ) then
          call proper_runstring('global_sel.hlp','global_sel',1)
      else
*         Decode the reference list (separated by colons)
          if ( index(runstring,"#") .ne. 0 ) then
             mode = 3
             call sub_char(runstring,'#',' ')
             net = 'XXX'
             indx = 0
             do while ( trimlen(net).gt.0 )
                j = j + 1
                call GetWord(runstring, net, indx)

                call sub_char(net,':',' ')
                name = 'XXX'
                net_refs(j) = 0
                jndx = 0
                do while ( trimlen(name).gt.0 )
                   name = ' '
                   call GetWord(net, name, jndx)
                   call casefold(name)
                   if( trimlen(name).gt.0 ) then
                     num_refs = num_refs + 1
                     net_refs(j) = net_refs(j) + 1
                     ref_names(num_refs) = name
                   end if
                end do
             end do 
             num_nets = j - 1
c             print*,'num_nets',num_nets
c             print*,'net_refs ',net_refs
c             print*,'ref_names ',ref_names
          else
             call sub_char(runstring,':',' ')
             name = 'XXX'
             indx = 0
             do while ( trimlen(name).gt.0 )
                name = ' '
                call GetWord(runstring, name, indx)
                call casefold(name)
                if( trimlen(name).gt.0 ) then
                     num_refs = num_refs + 1
                     ref_names(num_refs) = name
                end if
             end do  
          end if          
      end if

      len_run = rcpar(2, runstring) 
      if( len_run.eq.0 ) then
          call proper_runstring('global_sel.hlp','global_sel',1)
      else
          read(runstring,*,iostat=ierr) sel_num
          call report_error('IOSTAT',ierr,'decod',runstring,1,
     .                      'global_sel')
         
* MOD TAH 0611004: See if value is negative, if so make number of 
*          networks
           if( sel_num.lt.0 ) then
               num_nets = -sel_num
               mode = 2
           endif

      end if

      len_run = rcpar(3, apr_file) 
      if( len_run.eq.0 ) then
          call proper_runstring('global_sel.hlp','global_sel',1)
      end if

      

c      print*,'sel_num,num_nets,net_refs: ',sel_num,num_nets,net_refs
c      print*,'ref_names: ',ref_names

*     If avail_file given, read now to get names
      num_avf = 0
      use_num = 0
      do while ( len_run.gt.0 )
          len_run = rcpar(4+num_avf,test)
          if( len_run.gt.0 ) then 
              num_avf = num_avf + 1
              avail_file = test
              open(100, file=avail_file, status='old', iostat=ierr)
              call report_error('IOSTAT',ierr,'open',avail_file,1,
     .                          'global_sel')
              if( num_avf.eq.1 ) then
                 use_num = 0
                 do while ( ierr.eq.0 ) 
                     read(100,'(a)',iostat=ierr) line
                     if( ierr.eq.0 .and. line(1:1).eq.' ' .and.
     .                   trimlen(line).gt.0 ) then 
                         indx = 0
                         call GetWord(line, name, indx) 
                         call casefold(name) 
                         use_num = use_num + 1
                         use_names(use_num) = name(1:4)
                     end if
                 end do
              else 
*                Read the next file and mark all names that
*                appear in current list.
                 do i = 1, use_num
                    use_OK(i) = .false.
                 end do
                 ierr = 0
                 do while ( ierr.eq.0 ) 
                     read(100,'(a)',iostat=ierr) line
c                     print *, line(1:10), use_num
                     if( ierr.eq.0 .and. line(1:1).eq.' ' .and.
     .                   trimlen(line).gt.0 ) then 
                         indx = 0
                         call GetWord(line, name, indx) 
                         indx = 0
                         call get_cmd(name(1:4), use_names, use_num,
     .                                    iel, indx)
                         if( iel.gt.0 ) then
                             use_OK(iel) = .true.
                         end if
                     end if
                 end do

*****            Now see how many names have been marked as not OK.
*                Remove these names from the list
                 i = 0
                 do while ( i.lt.use_num )
                    i = i + 1
                    if( .not.use_OK(i) ) then
*                       Remove name
                        do j = i, use_num -1
                           use_names(j) = use_names(j+1)
                           use_OK(j) = use_OK(j+1)
                        end do
                        i = i - 1
                        use_num = use_num - 1 
                    endif
                 end do
             end if
             close(100)
         end if
      end do
!     write(*,50) use_num, use_names(1:use_num)
! 50  format('USE_NAMES # ',i5,/,200(20(a,1x),:/))

****  Open and read apr_file
      open(100,file=apr_file, status='old', iostat=ierr)
      call report_error('IOSTAT',ierr,'open',apr_file,1,
     .                      'global_sel')

      i = 0
      use_found = .false.  ! Mark all use sites as not found
      do while ( ierr.eq.0 ) 
         read(100,'(a)',iostat=ierr) line
* MOD TAH 200419: Explicitly ingore extended entries
         indx = 0
         skip = .false.
         if( ierr.eq.0 .and. line(1:1).eq.' ' .and.
     .        trimlen(line).gt.0 ) then 
            call GetWord(line, name, indx) 
            call casefold(name) 
            if( name(1:8).eq.'EXTENDED' ) then
              skip = .true.
            endif
         endif

         if ( ierr.eq.0 .and. line(1:1).eq.' ' .and.
     .        trimlen(line).gt.0 .and. .not.skip ) then 
            indx = 0
            call GetWord(line, name, indx) 
            call casefold(name) 
            jndx = 0
*           If avail site passed, make sure site is in list
            jel = 1
            iel = 0
            if( use_num.gt. 0 ) then
* MOD TAH 200421: Directly scan for just the 4-char code
!              call get_cmd(name, use_names, use_num, jel, jndx)
               jel = -1
               do ir = 1, use_num
                  if( name(1:4).eq.use_names(ir) ) then
                     if( .not.use_found(ir) ) then
                        jel = ir
                        use_found(ir) = .true.
                     else
                       jel = -1  ! Mark as not found
                     endif
                     exit
                  endif
               end do
            end if
            if( jel.gt.0 ) then
                jndx = 0
                call get_cmd(name, site_names, i, iel, jndx)
            end if
            if( iel.eq. -1 .and. jel.gt.0 ) then
*               This is a new name so add to list
                i = i + 1
                if( i.ge.max_site ) then
                    write(*,100) max_site
 100                format('Ending site reads at ',i4,' sites')
                    ierr = -1
                end if

                site_names(i) = name
                call multiread(line, indx, 'R8', jerr, site_xyz(1,i),
     .                          cd, 3)
                call multiread(line, indx, 'R8', jerr, site_vel(1,i),
     .                          cd, 3)
                call read_line(line, indx, 'R8', jerr, site_ep(i),cd)
                all_list(i) = 1
                if( jerr.ne.0 ) i = i - 1
            end if
         end if
      end do 

      tot_sites = i 
      write(*,120) tot_sites, apr_file(1:trimlen(apr_file)),
     .             use_num, avail_file(1:trimlen(avail_file))
 120  format('* GLOBAL_SEL: ',i5,' sites found in ',a,/,
     .       '*             ',i5,' sites are available in ',a)      
      k = 0
      do j = 1,num_nets
        if ( mode .eq.3 ) then
           write(*,130) j, (ref_names(i), i = k+1,k+net_refs(j))
           k = k + net_refs(j)
        else if ( mode .eq. 2 ) then
           write(*,130) j, (ref_names(i), i = 1,num_refs)
        endif
 130    format('* Network ',i2,' : Requested reference sites ',
     .        100(a8,1x))
      enddo
      if ( mode .eq. 1 ) then
           write(*,130) j, (ref_names(i), i = 1,num_refs)
      endif

*     Get the reference site from list
*     Mark all the reference sites
      j = 0
      l = 1
      m = 0
      k = net_refs(1)
      do i = 1, num_refs
         jndx = 0
         if ( mode . eq. 3. and. i .gt. k ) then
            k = k + net_refs(l+1)
            l = l+1
            n = 0
         endif 
         call get_cmd(ref_names(i), site_names, tot_sites, 
     .                    iref, jndx)
* MOD TAH 200421: If reference site not found; use just 4-char code
*        (case where no _GPS version of name appears in apriori file)
         if( iref.le.0 ) then
*           Scan for 4-character name and take first found
            do ir = 1, tot_sites
               if( ref_names(i)(1:4).eq.site_names(ir)(1:4) ) then
                  iref = ir
                  exit
               endif
            enddo
         endif
*        See if we are still missing site
         if( iref.le.0 ) then
             if ( mode .eq.3 ) net_refs(l) = net_refs(l) - 1
             write(*,140) ref_names(i)
 140         format('Reference site ',a,' not found. ')
         else
             j = j + 1
             n = n + 1
             if ( mode .eq. 2 ) then 
                 do m = 1,num_nets
                     use_list(m,j) = iref
                 enddo
             else if ( mode .eq. 3 ) then
                 use_list(l,n) = iref
             else
                 use_list(1,j) = iref
             endif
             all_list(iref) = -1
         endif
      end do 
      num_refs = j

* Check that all networks have at least 1 reference sties (necessary to start search)
      if ( mode .eq. 3 ) then
          do i = 1, num_nets
              if ( net_refs(i) .le. 0 ) then
                  write(*,145) i
 145              format('* Network ',i2,
     .                   ': Zero valid reference sites found - Stop')
                  write(message,145) i
                  call report_stat('FATAL','GLOBAL_SEL','main',
     .                             ' ',message,0) 
                  stop
              endif
          end do
      endif

* MOD SCM 061212: Three modes of running. Orginal: (1) One network with sel_num
*     sites, (2) multiple networks (#sites/net defined by user) simulataneously 
*     selected using all available sites in .apr file. All reference sites (num_ref)
*     common to all networks. (3) multiple networks (#sites/net defined by user, #networks 
*     defined by user). Unique reference sites for each network (net_refs) common to all
*     networks. Mode 2 and 3 have networks simulataneously selected from all available sites 
*     in .apr file. 

* MOD TAH 061104: See how many sites per network are needed
      if( num_nets.eq.0 .and. net_refs(1) .eq. 0 ) then
          write(*,150) sel_num, num_refs
 150      format('* Network Mode 1: Selecting a network with ',i3,
     .           ' sites, ',i3,' Reference sites')
      else if ( num_nets.gt.0 .and. net_refs(1) .eq. 0 ) then
          sel_num = (tot_sites-num_refs)/num_nets + num_refs
          write(*,160) num_nets, sel_num, num_refs
 160      format('* Network Mode 2: ',i3,' nets each with ',i3,
     .           ' sites, ',i3,' Common reference sites per network')
      else if ( num_nets.gt.0 .and. net_refs(1) .gt. 0  ) then
          write(*,170) num_nets, sel_num
 170      format('* Network Mode 3: Forming ',i2,' networks containing '
     .           ,i3,' sites each')
      endif

*     Now compute all the baseline lengths
      do i = 1, tot_sites-1
         ri = sqrt(site_xyz(1,i)**2+site_xyz(2,i)**2+
     .             site_xyz(3,i)**2)
         do j = i+1, tot_sites
            k = nb(i,j)
            rj = sqrt(site_xyz(1,j)**2+site_xyz(2,j)**2+
     .                site_xyz(3,j)**2)
            idotj = (site_xyz(1,i)*site_xyz(1,j)+
     .               site_xyz(2,i)*site_xyz(2,j)+
     .               site_xyz(3,i)*site_xyz(3,j))/(ri*rj)
            blens(k) = acos(idotj)*6378153.d0

*            blens(k) =  sqrt((site_xyz(1,i)-site_xyz(1,j))**2+
*     .                      (site_xyz(2,i)-site_xyz(2,j))**2+
*     .                      (site_xyz(3,i)-site_xyz(3,j))**2)
* MOD TAH 031125: Replace with arc length
            
c            write(*,190) i,j,k, site_names(i), site_names(j), 
c     .                   blens(k)/1000.d0
c 190        format(3i4,1x,2(a8,1x),F10.4)

         end do
      end do

****  Now start forming the list of stations
* MOD TAH 061104: Added zero show part of all networks

* If mode 1 or 2 (common reference sites)
      if( net_refs(1) .eq. 0 ) then 
          do i = 1, num_refs
              write(*,200) site_names(use_list(1,i)), 
     .             (site_xyz(j,use_list(1,i)),j=1,3),
     .             (site_vel(j,use_list(1,i)),j=1,3),
     .             site_ep(use_list(1,i)),1
  200              format(1x,a8,1x,3(F15.4,1x),1x,3(F8.4,1x),
     .                    1x,F10.4,2x,I2,'Ref     1N')
          end do
* If mode 3 (net_refs unique reference sites)
      else
         i = 0 
         do n  = 1,num_nets
            do k = 1, net_refs(n)
               i = i + 1
               write(*,205) site_names(use_list(n,k)), 
     .              (site_xyz(j,use_list(n,k)),j=1,3),
     .              (site_vel(j,use_list(n,k)),j=1,3),
     .               site_ep(use_list(n,k)),k, n, n
 205               format(1x,a8,1x,3(F15.4,1x),1x,3(F8.4,1x),
     .                  1x,F10.4,3x,I2,2x,I2,'Ref',5x,I2,'N')
            enddo
         enddo
      endif
  
* MOD TAH 061104: Two modes of running. Orginal: One network with sel_num
*     sites or new mode, mulitple networks simulataneosly selected.
* MOD 200421: Add some booking keeping for less out when small numbers
*     of sites are used for core-selection.
      num_missed = 0
      if( num_nets.eq. 0 ) then 
         do i = num_refs+1, sel_num
            max_dist = 0.d0 
            max_indx = 0

*           Get the mean distance between each site and the list of
*           sites already obtained
            do j = 1, tot_sites

*              Compute mean distance from this site to ones in use_list.
*              (provided site not in use_list)
               if( all_list(j).eq.1 ) then
*                  Site not in list
                   av_dist = 0
*    MOD TAH 050407: Initialise sq_dist summation variable.
                   sq_dist = 0
                   use = .true.
                   do k = 1, i-1
c                      if( blens(nb(j,use_list(k))).lt.
c        .                               20000.d3/sel_num ) then
c                          use = .false.
c                      end if

                      av_dist = av_dist + blens(nb(j,use_list(1,k)))
                      sq_dist = sq_dist + blens(nb(j,use_list(1,k)))**2
                   end do
                   mean_dist = av_dist/(i-1)
C                  av_dist = sqrt(sq_dist-mean_dist**2*(i-1))/(i-1)
C                  if( i-1.lt.1000 ) av_dist = mean_dist
                   av_dist = mean_dist
                   if( av_dist.gt. max_dist .and. use ) then
                       max_dist = av_dist
                       max_indx = j
                   end if
               end if
C              write(*,*) i,j, mean_dist/1000, av_dist/1000, max_indx, 
C        .                 max_dist
           end do 

****    See who is max
           if( max_indx.gt.0 ) then
               use_list(1,i) = max_indx
               all_list(max_indx) = -1
               write(*,210) site_names(use_list(1,i)), 
     .                (site_xyz(j,use_list(1,i)),j=1,3),
     .                (site_vel(j,use_list(1,i)),j=1,3),
     .                site_ep(use_list(1,i)), i, (max_dist)/1000.d0
 210            format(1x,a8,1x,3(F15.4,1x),1x,3(F8.4,1x),
     .               1x,F10.4,1x,i4,F10.2,'   1N')
           else
               num_missed = num_missed + 1
!              write(*,*) 'No more sites'
           endif
         end do
         if( num_missed.gt.0 ) write(*,215) num_missed 
 215     format('Short Fixed list.  Unable to assign ',i4,' sites')
      else
* MOD TAH 061104: New mode simular code but simulteneous network selection.
* Get the correct number of refs included in each network.
         if ( net_refs(1) .eq. 0 ) then
              nrefs = num_refs + 1
         else
              nrefs = 1
         endif
         do i = nrefs, sel_num

*           MOD TAH 061104: New loop over the networks
            do n  = 1,num_nets
               max_dist = 0.d0 
               max_indx = 0

* Check if this network should be skipped this iteration
               if ( net_refs(n) .lt. i ) then

*              Get the mean distance between each site and the list of
*              sites already obtained
                  do j = 1, tot_sites

*                 Compute mean distance from this site to ones in use_list.
*                 (provided site not in use_list)
                     if( all_list(j).eq.1 ) then
*                       Site not in list
                        av_dist = 0
*       MOD TAH 050407: Initialise sq_dist summation variable.
                        sq_dist = 0
                        use = .true.
                        do k = 1, i-1
c                          if( blens(nb(j,use_list(k))).lt.
c                            20000.d3/sel_num ) then
c                              use = .false.
c                          end if

                         av_dist = av_dist + 
     .                             blens(nb(j,use_list(n,k)))
                         sq_dist = sq_dist + 
     .                             blens(nb(j,use_list(n,k)))**2
                        end do
                        mean_dist = av_dist/(i-1)
C                       av_dist = sqrt(sq_dist-mean_dist**2*(i-1))/(i-1)
C                       if( i-1.lt.1000 ) av_dist = mean_dist
                        av_dist = mean_dist
                        if( av_dist .gt. max_dist .and. use ) then
                           max_dist = av_dist
                           max_indx = j
                        end if
                     end if
c                     write(*,*) i,j, mean_dist/1000, av_dist/1000, 
c     .                          max_indx, max_dist
                  end do 

****          See who is max
                  if( max_indx.gt.0 ) then
                     use_list(n,i) = max_indx
                     all_list(max_indx) = -1
                     write(*,220) site_names(use_list(n,i)), 
     .                     (site_xyz(j,use_list(n,i)),j=1,3),
     .                     (site_vel(j,use_list(n,i)),j=1,3),
     .                     site_ep(use_list(n,i)), i, 
     .                     (max_dist)/1000.d0, n
 220                 format(1x,a8,1x,3(F15.4,1x),1x,3(F8.4,1x),
     .                    1x,F10.4,1x,i4,F10.2,2x,I2,'N')
                  else
                     num_missed = num_missed + 1
!                    write(*,*) 'No more sites'
                  endif
* end of skip network this round
               endif
            end do   ! Looping over networks
         end do      ! Loop over number of sites needed
         if( num_missed.gt.0 ) write(*,225) num_missed 
 225     format('Short Available list.  Unable to assign ',i4,' sites')
      end if

****  Thats all
      end


CTITLE NB
      integer*4 function nb( i, j  )

      implicit none 

*     Function to return baseline number

      integer*4 i,j

      if( j.gt.i ) then
          nb = (j-1)*(j-2)/2 + i
      else if( i.gt.j ) then
          nb = (i-1)*(i-2)/2 + j
      else
          write(*,*) 'Invalid baseline in NB ',i,j
          stop
      end if

      return
      end




