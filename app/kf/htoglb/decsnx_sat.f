CTITLE DECSNX_SAT

      subroutine decsnx_sat( unit, line)

      implicit none

*     Routine to decode the satellite information blocks

      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
      include '../includes/glb_hdr_def.h'
      include '../includes/svinf_def.h'

* PASSED VARIABLES
 
*   unit        - Unit number
 
      integer*4 unit

*   line        - Line read from input
 
      character*(*) line
 
* LOCAL VARIABLES

      integer*4 indx, trimlen
 
*   block_found - True is block found
 
      logical block_found

****  Start decoding the types of site blocks

      block_found = .false.
 
      indx = index(line,'SATELLITE/ID' )
      if( indx.gt.0 ) then
          block_found = .true.
          call dec_sat_id( unit, line)
      endif 

      indx = index(line,'SATELLITE/PHASE_CENTER' )
      if( indx.gt.0 ) then
          block_found = .true.
          call dec_sat_phs( unit, line)
      endif 

****  Thats all the blocks so far
      if( .not.block_found ) then
          write(*,500) line(1:trimlen(line))
500       format('DECSNX_SAT: Unknown block type',/a)
          call snx_finbl(unit)
      end if

****  Thats all
      return
      end

CTITLE DEC_SAT_ID

      subroutine dec_sat_id ( unit, line)

      implicit none

*     Routine to write the satellite ID block.

      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/svinf_def.h'
      include '../includes/sln_def.h'

* PASSED VARIABLES
 
*   unit        - Unit number
 
      integer*4 unit

*   line        - Line read from input
 
      character*(*) line
 
* LOCAL VARIABLES

c      real*8 dtl   ! Total Days since launch

      integer*4 dyrl  ! Years since launch
     .,         ddayl ! Days since launch   
     .,         dsl   ! Seconds part of day since launch
     .,         dyre  ! Years until end
     .,         ddaye ! Days until end  
     .,         dse   ! Seconds part of day until end
     .,         i,j   ! Loop counter
     .,         blk   ! Block number
     .,         ierr  ! IOSTAT Error
     .,         jerr  ! IOSTAT Error (decoding lines)
     .,         trimlen  ! Length of string

*   end_block   - Set true at end of block
      integer*4 conoff  ! Function to return constellation offset value

      logical end_block

      character*2 pt   ! Point code
* ssys -- Satelite system.  Fo the moment we only handle GPS (G)
      character*1 tab, ssys 

      character*1 offcon    ! Function to return constellation letter based on
                            ! offset PRN (original PRN for mod(prn,100))

      character*20 in_ant_name

****  Start decoding
 
      end_block = .false.
      tab = char(9)
      ierr = 0
      i = 0
      do while ( ierr.eq.0 .and. .not.end_block )
 
          read(unit,'(a)', iostat=ierr) line
          if( ierr.eq.0 .and. line(1:1).eq.' ' ) then
              i = i + 1
              if( i.gt.max_glb_svs ) then
                  call report_stat('WARNING','htoglb','dec_sat_id',' ',
     .                'Too Many satellites',max_glb_svs)
                  i = i - 1
              endif
              read(line, 120,iostat=jerr) ssys, qsvi_svn(i), 
     .             qsvi_prn(i),  pt, dyrl, ddayl, dsl,
     .             dyre, ddaye, dse,in_ant_name
 120          format(1x,a1,I3.3,1x,i2.2,1x,9x,1x,a1,1x,
     .             I2.2,1x,i3.3,1x,I5.5,1x,I2.2,1x,i3.3,1x,I5.5, 
     .             1x, a20)
              call report_error('IOSTAT',jerr,'decod',line,1,
     .             'DEC_SAT_ID')

* MOD TAH 180401: Update for the constellation type (multi-gnss)
              qsvi_prn(i) = qsvi_prn(i) + conoff(ssys)
              call name_to_blk(+1, in_ant_name, blk )

              if( blk.eq.0 ) then
                  write(*,220) line(1:trimlen(line)), 
     .                         in_ant_name(1:trimlen(in_ant_name))
 220             format('*** WARNING *** No antenna type match in: ',
     .                  a,':',a,':')
              endif
              qsvi_block(i) = blk
              call yds_to_jd(dyrl, ddayl, dsl, qsvi_launch(i))

          end if
          if( line(1:1).eq.'-' .or. ierr.ne.0 ) then
              end_block = .true.
          end if

      end do

*     Save the number of satellites read
      qnum_svs = i

*     Save the SVS and SVS and PRN numbers as a name
      do i = 1, qnum_svs
          ssys = offcon(qsvi_prn(i))
          write(qsvs_names(i),320) ssys, qsvi_svn(i), 
     .          mod(qsvi_prn(i),100)
 320      format(a1,I3.3,'_',I2.2)
      end do

****  Thats all
      return
      end

CTITLE DEC_SAT_PHS

      subroutine dec_sat_phs( unit, line)

      implicit none

*     Routine to write the satellite Phase model block.

      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/svinf_def.h'

* PASSED VARIABLES
 
*   unit        - Unit number of sinex file
 
      integer*4 unit
 
*   line        - Line read from file
 
      character*(*) line
 
* LOCAL VARIABLES
      integer*4 i, j
     .,         tsvs  ! Test of SVS number
     .,         ierr  ! IOSTAT Error

*   end_block   - Set true at end of block
 
      logical end_block

* ssys  -- Satellite system (GPS only for the moment)
      character*1 ssys

      end_block = .false.
      ierr = 0
      i = 0
      if( qnum_svs.eq.0 ) then
         call report_stat('WARNING','htoglb','dec_sat_phse',
     .                   ' ','No SATELLITE/ID Block found',0)
      endif
         
      do while ( ierr.eq.0 .and. .not.end_block )

          read(unit,'(a)', iostat=ierr) line
* MOD TAH 200131: Changed to work for all systesm
          if( ierr.eq.0 .and. line(1:1).eq.' ' ) then
*             Code below works for both GPS and GLONASS
              read(line, 120) ssys, tsvs
*             See if we have satellites records
              i = 0
              if( qnum_svs.gt.0 ) then
                  do j = 1, qnum_svs
                     if( qsvi_svn(j).eq.tsvs ) i = j
                  end do
              else
                  i = i + 1
                  if( i.gt.max_glb_svs ) then
                     call report_stat('WARNING','htoglb','dec_sat_phse',
     .                    ' ','Too Many satellites',max_glb_svs)
                     i = i - 1
                  endif
              endif
*             Read if satellite found
              if( i.gt.0 ) then 
                 read(line, 120) ssys, qsvi_svn(i), qsvi_ocode(i)(1:1),
     .              qsvi_antpos(3,1,i), qsvi_antpos(1,1,i),
     .              qsvi_antpos(2,1,i), qsvi_ocode(i)(2:2),
     .              qsvi_antpos(3,2,i), qsvi_antpos(1,2,i),
     .              qsvi_antpos(2,2,i), qsvi_antmod(i)(1:10),
     .              qsvi_ocode(i)(3:3), qsvi_ocode(i)(4:4)
 
 120             format(1x,a1,I3.3,2(1x,a1,3(1x,F6.4)),1x,a10,1x,
     .               a1,1x,a1)
              end if 
          end if
          if( line(1:1).eq.'-' .or. ierr.ne.0 ) then
              end_block = .true.
          end if

      end do

****  If the number of satellites passed is zero, save as SVS numbers
      if( qnum_svs.eq.0 ) then
          qnum_svs = i
          do j = 1, qnum_svs
             qsvi_prn(j) = qsvi_svn(j)
             write(qsvs_names(i),320) qsvi_prn(j)
 320         format('SVS_',i3.3)
          end do
      end if

            

****  Thats all
      return
      end


         

