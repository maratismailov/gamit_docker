CTITLE RX_HEAD
 
      subroutine read_xf_head

      implicit none
 
*     Routine to read the header from the X data files.
*     Returns will be:
*     Approximate site position
*     station name
 
      include '../includes/xfile_def.h'
      
* Local variables

*    swver  - Rinex software version
      character*20 flg1*1, flg2*1
            
*   ierr    - IOSTAT error
*  trimlen  - Length of string
*   indx    - Position of substring in a string
 
      integer*4 ierr, trimlen, indx, i
      integer*4 hr, min, dlat,mlat, dlon, mlon
      real*8 slat,slon,radius
 
*   eoh     - End of header flags
*   readx   - True if x file
       logical eoh, readx, site, firsteph, antenna
 
*   line    - Line read from file
 
 
      character*256 line
      character*1 cr
      character*2 datatypes(xf_maxdat)
 
****  Open the data file
      cr = char(13) 
      call openfile(obs_lu, xf_rfile, 'old',' ','ASCII', ierr)
      call report_error('IOSTAT',ierr,'open',xf_rfile, 1,
     .        'read_rinex_file')

*     Start looping over the header

*     See if we find X file
      read(obs_lu,'(a)', iostat=ierr) line
      indx = index(line,'Phase and Pseudorange')
      if( indx.ge.5 )  readx = .true.

*     Escape the first four lines
      do i =1, 4
        read(obs_lu,'(a)', iostat=ierr) line
      end do
      
      eoh = .false.
      site = .false.
      firsteph = .false.
      antenna = .false.
      xf_ntext = 0
      
      do while ( .not.eoh )
          read(obs_lu,'(a)', iostat=ierr) line
          write(*,*) line(1:84)
c          call sub_char(line,cr,' ')
c          call casefold(line)

          if( index(line,'EPOCH # FLG').gt.0 ) eoh = .true.

*     See if we need to change position format

          indx = index(line,'STATIC')
          if( indx.gt.0 )   
     .       write(line,'(a)') ' EXPANDED DATA FORMAT KINEMATIC'

*         Save the header into xf_test

          xf_ntext = xf_ntext + 1
          write(xf_text(xf_ntext),'(a)') line(1:trimlen(line))
 
*         Get RINEX file name
          indx = index(line,'FROM FILE:')
          if( indx.gt.5 ) read(line(16:), '(a)') xf_rfile
          
         
*         See if marker name
          indx = index(line,'site code  :')
          if( indx.gt.5 ) xf_makernam = line(21:24)

*         See if operator name
          indx = index(line,'Operator name    :')
          if( indx.gt.0 ) xf_openam = line(21:30)
                    
*         See if receiver type and serial # , software version
          indx = index(line,'Receiver serial #')
          if( indx.gt.0 ) then
            read(line(21:), 100 ) xf_rcvnum
          end if
 100      format(a20)

          indx = index(line,'Receiver         :')
          if( indx.gt.0 ) then   
            read(line(21:), 100 )  xf_rctype
          end if
          
*         See if antenna type and serial #.
c          indx = index(line,'Antenna          : ')
c          if( indx.gt.0 ) then
c            read(line, 110 ) xf_antnum, xf_anttyp
c          end if
   
                     
*         See if position
          if (site) then
             read(line, 120) xf_sitnam, flg1, dlat, mlat, slat,
     .          flg2, dlon, mlon, slon, radius, xf_rctype, xf_swver
             write(xf_posflg, '(a1,a1,"U")') flg1, flg2
             
c             if (flg1.eq.'N'.or.flg1.eq.'S') xf_posflg(3:3)='U'

             if (xf_pos(1).eq.0)
     .         call geoxyz(dlat,mlat,slat,flg1,dlon,mlon,slon,
     .             flg2,radius,xf_pos(1), xf_pos(2), xf_pos(3))
             site = .false.
           end if
 120       format(1x,a16,1x,a1,i2,1x,i2,1x,f8.5,1x,a1,i3,1x,i2,
     .         1x,f8.5,f13.4,2x,a3,1x,f5.2 )
           
          if(index(line,'SS.SSSSS').gt.0) site = .true.
            

*         See if antenna offset here.
          if(antenna) then
              read(line,140 )xf_anttyp, xf_offsL1, xf_offsL2
              antenna = .false.
           endif
  140      format(1x,a6,9x,2(3x,3f8.4))
          if (index(line,'ANT OFFSETS (M)').gt.0) antenna = .true.

  
*         See if data types, and satellite number
          indx = index(line,'DATA TYPES:')
          if( indx .gt.9 ) then
* MOD TAH 090114: Allow for more than 9 entries on a line
              read(line,160) xf_nsat, xf_ndat,
     .                    (datatypes(i),i=1, max(xf_ndat,9)) 
 160          format(i3, 23x, i2, 12x, 8(1x,a2))


              do i = 1, xf_ndat 
                xf_dattyp(i) = i
              end do
              
              do i=1, xf_nsat
                read(obs_lu,'(a)', iostat=ierr) line
                xf_ntext = xf_ntext + 1
                write(xf_text(xf_ntext),'(a)') line(1:trimlen(line))
                 read(line, 200) xf_ischan(i), xf_prn(i)
              enddo
                                      
          end if
 200      format(8x, i3,6x, i3)
          
*         See if  last obs. epoch.

          if (firsteph) then
             read(line, *) xf_srtdoy, hr,min, xf_srtsec, xf_inter,
     .             xf_ircint, xf_isessn
             xf_srtsect = hr*3600+min*60 + xf_srtsec
             firsteph = .false.
           end if
          indx = index(line,'INTERVAL(SECS)')
          if( indx.gt.0 )   then
             firsteph = .true.
          end if
*         See if  # of epochs.
          indx = index(line,'EPOCHS')
          if( indx.eq.5 )  read(line,*) xf_nepoch

      end do

 
*     See if we found xf file
      if( readx ) then
          write(*,150) xf_rfile(1:trimlen(xf_rfile)),
     .        xf_sitnam, xf_pos
 150      format(/'* For X data file ',a,/,
     .        '* Site ',a,' Aprrox. position ',3F15.3)
      else
          write(*,170) xf_rfile(1:trimlen(xf_rfile))
 170      format(a,' Does not appear to be X  file')
          stop 'rx_head: Wrong type data file'
      end if

      
****  Thats all
      return
      end
 
