CTITLE PREREAD_SNX
 
      subroutine preread_snx(unit, line, np, ierr, sepoch)
 
*     This routine will pre-read the sinex file to find
*     out how many parameters were estimated in this
*     solution.

* MOD TAh 981020: Count all EOP parameters since with 1.02 global
*     files we will be be able to save all values.

      include '../includes/kalman_param.h'
      include '../includes/htoglb_comm.h'

 
* PASSED VARIABLES
 
*   unit        - Unit number for input file
*   np          - Number of parameters
*   ierr        - IOSTAT error
 
      integer*4 unit, np, ierr

*   speoch      - Reference time for solution

      real*8 sepoch
 
*   line        - Line read from file
 
      character*(*) line
 
* LOCAL VARIABLES
 
*   indx        - Position in string
 
      integer*4 indx, cyr, cdoy, csec, ne

*   snx_vers        - Version number for sinex file
*   epvals(max_glb_parn) -- Epoch numbers of parameters which
*     we take the median as the value.

      real*8 snx_vers, epvals(max_glb_parn) 
 
*   start_found     - Set true when end of SOLUTION/ESTIMATE
*               - block is found
*   end_found   - End of block found
 
      logical start_found, end_found
 
****  Start reading down file looking for the SOLUTION/ESTIMATE
*     line
 
      start_found = .false.
      end_found = .false.
      ierr  = 0
      
* MOD TAH 960820: Moved getting of site list to the pre-read stage
*     so that only those sites used are moved forward.    
      qnum_sites = 0
      
* MOD TAH 960430: Get the version number
      read(line,'(5x,F5.2,1x,a3)', iostat=ierr) snx_vers, qowner
      call report_error('IOSTAT',ierr,'read',line,0,'Preread_snx')
      
      do while ( .not.end_found .and. ierr.eq.0 )
 
          read(unit,'(a)',iostat=ierr ) line

          if( ierr.eq.0 ) then
 
*             See if we have already found the start
              if( .not.start_found ) then
                  indx = index(line,'+SOLUTION/ESTIMATE')
                  if( indx.eq.1 ) then
                      start_found = .true.
                      np = 0
                      ne = 0
                  end if
              else
 
*                 See if end has been found
                  indx = index(line,'-SOLUTION/ESTIMATE')
                  if( indx.eq.1 ) then
                      end_found = .true.
                  else if( line(1:1).eq.' ' ) then
 
*                     Increment number of parameters.  See if
*                     will use this parameter
                      if( line(8:10).eq.'STA' ) then
                      

*                         Pull off the station name (we do this to make
*                         sure that only those sites which are 
*                         estimated make it through to the site name
*                         list.
                          if( line(8:11).eq.'STAX' ) then
                              call add_site_name( line, snx_vers )
                          end if                      
                      
                          np = np+1
*                         Get the time from the line.
                          if( snx_vers.eq.0.05d0 ) then
                              read(line,210,iostat=ierr) cyr, cdoy, csec
 210                          format(25x,i2,1x,i3,1x,i5)
                          else
                              read(line,215,iostat=ierr) cyr, cdoy, csec
 215                          format(27x,i2,1x,i3,1x,i5)
                          end if
                          ne = ne + 1
                          call yds_to_jd( cyr, cdoy, csec, epvals(ne))
                      end if

                      if( line(8:10).eq.'VEL' ) np = np+1

*                     Check translation and scale
                      if( line(8:11).eq.'TRAN' ) np = np +1
                      if( line(8:12).eq.'SCALE') np = np +1

*                     Now check polar motion/ut1
                      if( line(8:10).eq.'LOD' .or.
     .                    line(8:11).eq.'LODR') then
                          np = np+1
                      end if
                      if( line(8:10).eq.'UT ' .or.
     .                    line(8:10).eq.'UTR' ) then
                          np = np+1
                      end if
                      if( line(8:11).eq.'XPO ' ) then
                          np = np+1
                      end if
                      if( line(8:11).eq.'XPOR'  ) then
                          np = np+1
                      end if
                      if( line(8:11).eq.'YPO ' ) then
                          np = np+1
                      end if
                      if( line(8:11).eq.'YPOR'  ) then
                          np = np+1
                      end if
                     
                  end if
              end if
          end if
      end do
      
****  Now get the median epoch of the station positions
      call mdian2(epvals, ne, sepoch)
 
****  Finished: Report any errors and rewind file
      if ( ierr.ne.0 ) then
          call report_error('IOSAT',ierr,'SOLUTION find',line,0,
     .        'preread_snx')
      end if
 
      write(*,300) np
 300  format(' PreRead_snx: ',i5,' Parameters found')
      rewind(unit)

*     Read the first line of file again.
      read(unit,'(a)' ) line
      return
      end
 
CTITLE ADD_SITE_NAME

      subroutine add_site_name( line, snx_vers )
      
*     Routine to add site name to list based on estimated parameters.

      include '../includes/kalman_param.h'
      include '../includes/htoglb_comm.h'

* PASSED VARIABLES

* line      - Line read from file
* snx_vers  - Version of sinex file

      character*(*) line
      real*8 snx_vers      

* LOCAL VARIABLES

* ierr  - IOSTAT Error
* jp    - Parameter numbers
* stype - Parameter type
* code  - Site code
* pt    - Point number for site
* occ   - occupation code
* chocc - Character occupation code

      integer*4 ierr, jp, occ, indx, sn
      character*8 stype, code, pt, test_name 
      character*4 chocc
            
            
***** See what version of sinex we have
      if( snx_vers.eq.0.05d0 ) then
           read(line,210,iostat=ierr) 
     .         jp, stype, code,pt,occ
 210       format(i6,1x,a4,1x,a4,1x,a2,1x,i4)
      else
           read(line,215,iostat=ierr) 
     .          jp, stype, code,pt,chocc
 215              format(i6,1x,a6,1x,a4,1x,a2,1x,a4)
       
*          Check the occupancy value
           if( chocc.eq.'----' ) chocc = '   1'
           read(chocc,*,iostat=ierr) occ
      end if

* MOD TAH 970716: Check for blank pt names, replace with ' A'
      if( pt(2:2).eq.' ' .or. pt(2:2).eq.'-' ) pt = ' A'
      
****  Now construct the site name
      write(test_name,310) code(1:4),occ,pt(2:2)
 310  format(a4,i3.3,a1)
 
****  Increment the number of sites and save the name
      call casefold(test_name)               
      indx = 1
      call get_command(test_name, qsite_names, qnum_sites,
     .                 sn, indx )
      if( sn.gt.0 ) then
          write(*,410) line(1:80), test_name, sn
 410      format('FATAL ERROR: site in SOLUTION/ESTIMATE line',/,
     .           a,/,' Matches site name ',a8,' as site number ',i4)
          call report_stat('FATAL','htoglb','Preread_snx',test_name,
     .           'Duplicate site name',0)
      end if
 
***** All seems OK, so add site to list.      
      qnum_sites = qnum_sites + 1
      qsite_names(qnum_sites) = test_name
      
***** Thats all
      return
      end
      
