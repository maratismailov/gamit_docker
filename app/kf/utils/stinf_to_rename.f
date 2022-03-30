      program stinf_to_rename

      implicit none

*     Program to read a new format station.info file and a list of GPS
*     site names and generate rename lines for each antenna and radome
*     change
*
*     MOD 081008 McClusky - Added receiver swap renames and more 
*     output documentation.
*     MOD 090128 Herring - Added option to ignore receiver changes (third 
*                runstring parameter -norec)
*     MOD 121012 Herring - Introduced new scheme: 2-9 and then ASCII (J-Z)
*               for 10-26 renames. Above 26 there are problems (left at Z).
*               UNKN==NONE in new definitions for radomes.

      integer*4 max_stinf   ! Maximum number of station.info lines
      integer*4 max_site    ! Maximum number of sites to check
      integer*4 max_confc   ! Maximum number of configuration changes

*     parameter ( max_stinf = 20000 )
C MOD TAH 180315: Increased number of records and sites llowed in station.info
      parameter ( max_stinf = 30000 )
      parameter ( max_site  = 10000 )
      parameter ( max_confc =   128 )

****  Main program variables
      integer*4 num_site ! Number of sites
     .,         num_stinf  ! Nomber of station.info entries
     .,         num_confc  ! Number of configuration changes
     .,         sf_start(5), sf_stop(5)
     .,         sf_sdate(5,max_stinf), sf_edate(5,max_stinf)
     .,         ierr

      real*8 sf_smjd(max_stinf), sf_emjd(max_stinf) ! Start and stop MJD

      integer*4 rcpar, len_run, trimlen, i, j,k, date(5)
     .,         nent, old_j, end_j,  il
     .,         jc(max_confc)   ! Indices to records with configuration changes

      real*8 sectag, jd

      real*8 post_mjd   ! MJD after which entries should be output (set with
                ! -post option 

      character*4 site_names(max_site)
     .,           sf_names(max_stinf)

      character*1 typec   ! Type of change: A -- antenna, R -- Radome

      character*16 sf_antstr(max_stinf), ref_ant

      character*20 sf_recstr(max_stinf), ref_rec

      character*5  sf_dome(max_stinf), ref_dome

      character*10 sf_antsn(max_stinf), ref_antsn

      character*128 in_site, in_stinf
      character*256 line

      integer*4 nc, na, nr, ns   ! Counts for recv, ant, radome and single counter)
      character*1 antc(15), radc(10), recc(8)

      character*1 cmapns   ! Function to map ns to character

      logical norec   ! Set true if receiver changes are to be ignored.
      logical noend   ! Set true if stop times to be omitted throughout
      logical oldform ! Set true to use old scheme with different letters
                      ! for different types (A-> for antenna, R-> for radome
                      ! 2-> receiver changes).
      data antc / 'A','B','C','D','E','F','H','I','J','K',
     .            'L','M','N','O','1' /
      data recc / '2','3','4','5','6','7','8','9'/
      data radc / 'R','S','T','U','V','W','Y','Z','Q','P' / 


****  Get the file with list of sites
      len_run = rcpar(1, in_site)
      if( len_run.gt.0 ) then
        open(100,file=in_site, status='old',iostat=ierr)
        call report_error('IOSTAT',ierr,'open',in_site,1,
     .       'stinf_to_rename')
      else
        call proper_runstring('stinf_to_rename.hlp',
     .       'stinf_to_rename/in_site',1)
      end if
      len_run = rcpar(2, in_stinf)
      if( len_run.gt.0 ) then
        open(101,file=in_stinf, status='old',iostat=ierr)
        call report_error('IOSTAT',ierr,'open',in_stinf,1,
     .       'stinf_to_rename')
      else
        call proper_runstring('stinf_to_rename.hlp',
     .       'stinf_to_rename/in_stinf',1)
      end if     
      norec   = .false.
      noend   = .false.
      oldform = .false.

      post_mjd = 15020.0d0

      nr = 2
      do while ( len_run.gt.0 ) 
         nr = nr + 1
         len_run = rcpar(nr,line)
         if( len_run.gt.0 ) then 
*           See which option passed
            if( line(1:4).eq.'-nor' ) then
               norec = .true.    
            elseif (line(1:4).eq.'-noe' ) then
               noend = .true.
            elseif (line(1:4).eq.'-old' ) then
               oldform = .true.
            elseif ( line(1:2).eq.'-p' ) then
*              -post option given, read the date
               do j = 1, 5
                  len_run = rcpar(nr+j,line)
                  if( len_run.gt.0 ) then
                      read(line,*) date(j)
                  else
                      date(j) = 0
                  end if
               end do
               sectag = 0
               call ymdhms_to_mjd(date,sectag,post_mjd)
               nr = nr + 5
            end if
         endif
      enddo

****  Read in the list of sites
      num_site = 0
      do while ( ierr.eq.0 )
         read(100,'(a)',iostat=ierr) line
         call trimlead(line)
         call casefold(line)
         if( ierr.eq.0 ) then
             num_site = num_site + 1
             site_names(num_site) = line(1:4)
         end if
      enddo
      write(*,110) num_site, in_site(1:trimlen(in_site))
 110  format('* STINF_to_RENAME: ',i5,' Sites found in ',a)
      if( norec ) write(*,120)
 120  format('* Receiver changes not included')
      if( noend ) write(*,130) 
 130  format('* No end times will be written')
      if( post_mjd.gt. 15021.d0 ) then
          call mjd_to_ymdhms(post_mjd, date,sectag)
          write(*,140) date
 140      format('* Only new entries after ',i4,4i3,' will be output')
      endif
      if( oldform ) then
          write(*,'(a)') '* OLD A->, R->, 2-> form used'
      else
          write(*,'(a)')
     .             '* NEW STYLE: Sequential 2-9,J-Z (A-I for unknown)'
      end if


****  OK Read station.info
      num_stinf = 0
      ierr = 0
      sectag = 0
      do while ( ierr.eq.0 )
         read(101,'(a)',iostat=ierr) line
         if( ierr.eq.0 .and. line(1:1).eq.' ' .and.
     .       trimlen(line).gt.0 ) then
*            Decode line
             num_stinf = num_stinf + 1
             if( num_stinf.gt.max_stinf ) then
                write(*,145) in_stinf(1:trimlen(in_stinf)), max_stinf
 145            format('**DISASTER** Too many entries in ',a,/,
     .                 '             Max allowed (max_stinf) is ',I8)
                stop 'Station.info too large'
             end if
             i = num_stinf
             read(line,150) sf_names(i), (sf_start(j),j=1,5),
     .           (sf_stop(j),j=1,5)
 150         format(1x,a4,20x,i4,1x,i3,3i3,2x,i4,1x,i3,3i3)
             call casefold(sf_names(i))
             sf_recstr(i) = line(98:117)
             sf_antstr(i) = line(171:185)
             sf_dome(i)   = line(188:192)
             sf_antsn(i)  = line(195:204)
cd            write(*,160) sf_names(i), (sf_start(j),j=1,5),
cd     .          (sf_stop(j),j=1,5),sf_recstr(i),sf_antstr(i),sf_dome(i)
cd160         format(1x,a4,1x,i4,1x,i3,3i3,2x,i4,1x,i3,3i3,
cd     .             2x,a16,2x,a16,1x,a5)
             if( sf_stop(1).eq.9999 ) sf_stop(1) = 2100
             if( sf_stop(2).eq. 999 ) sf_stop(2) = 1
             if( sf_dome(i).eq.'-----' ) sf_dome(i) = 'NONE'
             if( sf_dome(i).eq.'UNKN ' ) sf_dome(i) = 'NONE'
             if( sf_dome(i).eq.'null ' ) sf_dome(i) = 'NONE'

*            Convert the dates
             date(1) = sf_start(1)
             date(2) = 1
             date(3) = sf_start(2)
             date(4) = sf_start(3)
             date(5) = sf_start(4)
             call ymdhms_to_mjd(date,sectag,jd)
             call mjd_to_ymdhms(jd,sf_sdate(1,i),sectag)
             sf_smjd(i) = jd
             date(1) = sf_stop(1)
             date(2) = 1
             date(3) = sf_stop(2)
             date(4) = sf_stop(3)
             date(5) = sf_stop(4)
             call ymdhms_to_mjd(date,sectag,jd)
             call mjd_to_ymdhms(jd,sf_edate(1,i),sectag)
             sf_emjd(i) = jd           
         endif
      end do
      write(*,180) num_stinf, in_stinf(1:trimlen(in_stinf))
 180  format('* STINF_to_RENAME: ',i5,' entries found in ',a)

****  OK: Now generate the renames that are needed.
      il = 0
      do i = 1, num_site

*****    Scan accross station.info
         nent = 0
         num_confc = 0
         nc = 0
         na = 0
         nr = 0
         ns = 1  ! First used value is 2PS
         do j = 1, num_stinf
            if( sf_names(j).eq.site_names(i) )  then
               nent = nent + 1
               il = j
               if( nent.eq.1 ) then
                  ref_ant = sf_antstr(j)
                  ref_rec = sf_recstr(j)
                  ref_dome = sf_dome(j)
                  ref_antsn = sf_antsn(j)
                  num_confc = 1
                  jc(1) = j
               else
****              See if change
* MOD TAH 090128: Check to see if receivers are to be considered
                  if( norec ) then  ! Test only antenna and receiver
                     if( ref_ant  .ne.sf_antstr(j) .or.
     .                   ref_dome .ne.sf_dome(j)   .or.
     .                   ref_antsn.ne.sf_antsn(j)  ) then
*                        OK: There is a change in the configuration
*                        Save index of change
                         num_confc = num_confc + 1
                         jc(num_confc) = j
*                        Save the new configuration
                         ref_rec = sf_recstr(j)
                         ref_ant = sf_antstr(j)
                         ref_dome = sf_dome(j)
                         ref_antsn = sf_antsn(j)
                     end if

                  else   ! Original code, consider receiver and antenna
* MOD TAH 090910: Added antenna serial number check to standard case.
                     if( ref_rec.ne.sf_recstr(j) .or.
     .                   ref_ant.ne.sf_antstr(j) .or.
     .                   ref_dome.ne.sf_dome(j)  .or. 
     .                   ref_antsn.ne.sf_antsn(j) ) then
*                        OK: There is a change in the configuration
*                        Save index of change
                         num_confc = num_confc + 1
                         jc(num_confc) = j
*                        Save the new configuration
                         ref_rec = sf_recstr(j)
                         ref_ant = sf_antstr(j)
                         ref_dome = sf_dome(j)
                         ref_antsn = sf_antsn(j)
                     end if
                  endif
               endif
            end if
         end do

*        OK: Now see how many changes
*        If there is more than the initial configuration print the changes
         if( num_confc.gt.1 ) then
            typec = ' '
            do j = 2,num_confc
               if( sf_antstr(jc(j-1)).ne.sf_antstr(jc(j)) .or.
     .             sf_antsn(jc(j-1)) .ne.sf_antsn(jc(j))  ) then
                   na = na + 1
                   ns = ns + 1
                   if( na.gt.15 .and. oldform ) then
                      write(*,'(a)') 
     .                     '** WARNING: Max antenna renames reached: 16'
                      na = 15
                   endif
                   if( oldform ) then
                      typec = antc(na)
                   else
                      typec = cmapns(ns) 
                   end if
                   write(line,220) sf_antstr(jc(j-1)),sf_antstr(jc(j)),
     .                             sf_recstr(jc(j-1)),sf_recstr(jc(j)),
     .                             sf_dome(jc(j-1)),sf_dome(jc(j))
 220               format('Antenna  swap ',a20,' to ',a20, 
     .                    ' Receiver ',a20,' to ',a20, 
     .                    ' Dome ',a5,' to ',a5)
               elseif ( sf_dome(jc(j-1)).ne.sf_dome(jc(j)) ) then
                   nr = nr + 1
                   ns = ns + 1
                   if( nr.gt.10 ) then
                      write(*,'(a)') 
     .                      '** WARNING: Max dome renames reached: 11'
                      nr = 10
                   endif
                   if( oldform .and. oldform ) then
                      typec = radc(nr)
                   else
                      typec = cmapns(ns) 
                   end if
                   write(line,240) sf_dome(jc(j-1)),sf_dome(jc(j)),
     .                             sf_antstr(jc(j-1)),sf_antstr(jc(j)),
     .                             sf_recstr(jc(j-1)),sf_recstr(jc(j))
 240               format('Dome swap     ',a5,' to ',a5, 
     .                    ' Antenna  ',a20,' to ',a20, 
     .                    ' Receiver ',a20,' to ',a20)
               elseif( sf_recstr(jc(j-1)).ne.sf_recstr(jc(j)) ) then
                   nc = nc + 1
                   ns = ns + 1
                   if( nc.gt.8 .and. oldform ) then 
                      write(*,'(a)') 
     .                   '** WARNING: Max receiver renames reached: 8'
                      nc = 8
                   endif
                   if( oldform ) then
                      typec = recc(nc)
                   else
                      typec = cmapns(ns) 
                   end if
                   write(line,260) sf_recstr(jc(j-1)),sf_recstr(jc(j)),
     .                             sf_antstr(jc(j-1)),sf_antstr(jc(j)),
     .                             sf_dome(jc(j-1)),sf_dome(jc(j))
 260               format('Receiver swap ',a20,' to ',a20, 
     .                    ' Antenna  ',a20,' to ',a20, 
     .                    ' Dome ',a5,' to ',a5) 
               endif
               old_j = jc(j)
* MOD TAH 140317: The reason to treat the last entry differently is
*              the end of the rename is given by the start of the next
*              one.  The last rename will be open ended.
               if( j.lt.num_confc ) then
                  end_j = jc(j+1)   
* MOD RWK 091001: Add open-ended option
* MOD TAH 100618: Check date to see if just new entries.
                  if( noend ) then
* MOD TAH 110719: Removed _GPS in name 1 so all names will match
!                     if( sf_smjd(old_j).ge.post_mjd ) 
* MOD TAH 140317: Here old_j points to the end of the j'th rename and
*                 we should test if the end of this renames falls after the
*                 post_mjd we are interested in.  The end_j entry here is
*                 last entry for this site which is 2100/01/01 in general.
*                 Changed end_j to old_j in tests below.  
!                    if( sf_emjd(end_j).ge.post_mjd ) 
!                    if( sf_emjd(old_j).ge.post_mjd ) 
* MOD TAH 140317: Use the start time of the next entry (end times may not
*                    be correct wheh multiple station.info entries with the
*                    same antenna. 
                     if( sf_smjd(end_j).ge.post_mjd )  
     .               write(*,280) site_names(i), site_names(i), typec,
     .                  (sf_sdate(k,old_j),k=1,5), line(1:trimlen(line)) 	     
 280                 format(1x,'rename ',a4,'     ',a4,'_',a1,'PS ',
     .                  (I4,4i3,1x),1x,'!',1x,a)
*                     Check the end date when closed because this 
*                     will change.
                  else
* MOD TAH 110719: Removed _GPS in name 1 so all names will match
!                    if( sf_emjd(end_j).ge.post_mjd )  
!                    if( sf_emjd(old_j).ge.post_mjd )  
                     if( sf_smjd(end_j).ge.post_mjd )  
     .               write(*,281) site_names(i), site_names(i), typec,
     .                  (sf_sdate(k,old_j),k=1,5), 
     .                  (sf_sdate(k,end_j),k=1,5), line(1:trimlen(line)) 	     
 281                 format(1x,'rename ',a4,'     ',a4,'_',a1,'PS ',
     .                  2(I4,4i3,1x),1x,'!',1x,a)
                  endif
               else
c                  write(*,280) site_names(i), site_names(i), typec,
c     .                (sf_sdate(k,old_j),k=1,5), 
c     .                (sf_edate(k,il),k=1,5), line(1:trimlen(line)) 
* MOD TAH 090910: Make the last entry run to 2100 1 1 0 0 
* MOD RWK 091001: Add open ended option
                  if( noend ) then
* MOD TAH 140317: Test on the end time of the rename to see if should be
*                 output (was tested on start time).
!                    if( sf_smjd(old_j).ge.post_mjd ) 
!                    if( sf_emjd(old_j).ge.post_mjd ) 
* MOD TAH 161119: Always output the last entry if typec is not blank
*                 (old code had problems with receiver change when -norec
*                  used if the change was after the last antenna change
*                  and before the -post date.
                     if( typec.ne.' ' ) 
     .               write(*,280) site_names(i), site_names(i), typec,
     .                 (sf_sdate(k,old_j),k=1,5),line(1:trimlen(line)) 	     
                  else
!                    if( sf_smjd(old_j).ge.post_mjd ) 
!                    if( sf_emjd(old_j).ge.post_mjd ) 
* MOD TAH 161119: Always output the last entry if typec is not blank
                     if( typec .ne.' ') 
     .               write(*,281) site_names(i), site_names(i), typec,
     .                (sf_sdate(k,old_j),k=1,5), 
     .                2100,1,1,0,0, line(1:trimlen(line)) 	     
                  endif
              end if
           end do
         end if
 
      end do

****  Thats all
      end


      character*(1) function cmapns(ns)
      implicit none

*     Function to map sequenctial change number to a character:
*     2-9 Numerical value; 10-26 J-Z

      integer*4 ns

***** See if less than 10
      if( ns.le.9 ) then
         cmapns = char( 48 + ns ) ! ns = 2; ascii 50 = '2'
      elseif( ns.le.26 ) then
         cmapns = char( 64 + ns ) ! ns = 10; ascii 74 = 'J' 
      else
         write(*,120) ns
 120     format('** WARNING ** too many renames (',i3,'). Code ',
     .          'kept at Z')
         cmapns = 'Z'
      endif

****  Thats all
      return
      end

 

