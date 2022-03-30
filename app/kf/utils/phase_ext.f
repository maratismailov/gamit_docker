      program phase_ext

      implicit none 

*     Rouitne to extract the phase data from DPH files the same as the
*     shell script below
C   set cnt = 32
C   \rm tmp.prn?? tmp.edt?? >& /dev/null

C   while ( $cnt )
C     if ( $cnt < 10 ) then
C       set prn = '0'$cnt
C     else
C       set prn = $cnt
C     endif
C     if ( -e ${ext}.${base}.PRN${prn} ) then
C       grep '^ ' ${ext}.${base}.PRN${prn}| awk '{ if( $12==0 ) print $0}' | awk -v ot=0 '{if( $1-ot > 10 ) print ">"} { if ( $4^2 < 1 ) {print $10,90-$11,$4, $1}} {ot=$1}' >! tmp.prn${prn}
C       grep '^ ' ${ext}.${base}.PRN${prn}| awk '{ if( $12==1 ) print $0}' | awk -v ot=0 '{if( $1-ot > 10 ) print ">"} {print $10,90-$11," 0.001", $1} {ot=$1}' >! tmp.edt${prn}
C     endif
C     @ cnt = $cnt - 1
C   end

      character*16 ext  ! The extend of the file name (DPH normally)
      character*4  base  ! Site name
      character*256 infile  ! Name of inoput file
      character*256 outprn ! Name of output file (tmp.prnXX)
      character*256 outedt ! Name of edit file (tmp.edtXX)
      character*256 line

      integer*4 prn   ! PRN number
      integer*4 lastep, lasted   ! Epoch numbers for last epochs
      integer*4 ierr, jerr, kerr  ! IOSTAT errors
      integer*4 lene, lenb, rcpar   ! Length of runstring entries
      real*8 values(14)  ! Values on line up to the edit flag

      logical openprn, openedt  ! Set true once files are opened

* MOD TAH 051101: Updated for P1 and P2 residuals being added to files.
      logical newformat
      integer*4 voff     ! Offset in values for addition of P1 and P2 (either 0 or 2)


****  Get the ext base
      lene = rcpar(1,ext)
      lenb = rcpar(2,base)
      if( lene.eq.0 .or. lenb.eq.0 ) then
         write(*,50) 
 50      format('PHASE_EXT runstring:',/,
     .          '% phase_ext <code> <site>',/,
     .          'where <code> is lead characters on phase files ',
     .          '(e.g., DPH)',/,
     .          '      <site> is 4-character site name',/
     .          'Program generates tmp.* files for sh_oneway')
         stop 'PHASE_EXT: Incomplete runstring'
      end if

      do prn = 1, 32
         write(infile,120) ext(1:lene), base(1:lenb), prn 
 120     format(a,'.',a,'.PRN',i2.2)
         lastep = 0
         lasted = 0
         open(50,file=infile,status='old',iostat=ierr)
         openprn = .false.
         write(outprn,140) prn 
 140     format('tmp.prn',i2.2)
         openedt = .false.
         write(outedt,160) prn 
 160     format('tmp.edt',i2.2)

         newformat = .false.
         voff = 0 

         if( ierr.eq.0 ) then    ! OK File exists
              do while ( ierr.eq.0 )
                 read(50,'(a)',iostat=ierr ) line
                 if( ierr.eq.0 .and. line(1:1).eq.' ' ) then
*                    Decode the line
                     call sub_char(line,'*','0')
                     read(line,*,iostat=jerr) values
*                    See which file this should goto
                     if( values(12+voff).eq.0 ) then   ! OK, good data

*                       Only output residuals less than 1 cycle and 
*                       see if we need a line separator
                        if( lastep.ne.0 .and. 
     .                      values(1)-lastep.gt.10 ) then
                            if( .not.openprn ) then
                                open(60,file=outprn,status='unknown',
     .                               iostat=kerr)
                                openprn = .true.
                            end if
                            write(60,'(a)',iostat=kerr) '>'
                            write(60,'(a)',iostat=kerr) '>'
                       endif
*                       Now write out line
                        if( .not.openprn) then
                            open(60,file=outprn,status='unknown',
     .                  	 iostat=kerr)
                            openprn = .true.
                        end if
*  MOD TAH 051101: Include offser value
                        if( abs(values(4+voff)).lt.1 ) 
     .                  write(60,220) values(10+voff), 
     .                                90-values(11+voff),
     .                                values(4+voff),int(values(1)) ! $10,90-$11,$4, $1}
 220                    format(2F12.4,1x,f12.3,1x,i8)
                        lastep = values(1)
                     elseif( values(12+voff).eq.1 ) then   ! bad data point

*                       Only output residuals less than 1 cycle and 
*                       see if we need a line separator
                        if( lasted.ne.0 .and. 
     .                      values(1)-lasted.gt.10 ) then
                            if( .not.openedt ) then
                                open(61,file=outedt,status='unknown',
     .                               iostat=kerr)
                                openedt = .true.
                            end if
                            write(61,'(a)',iostat=kerr) '>'
                        endif
*                       Now write out line
                        if( .not.openedt) then
                            open(61,file=outedt,status='unknown',
     .                  	 iostat=kerr)
                            openedt = .true.
                        end if
                        write(61,220) values(10+voff), 
     .                                90-values(11+voff),
     .                                0.001 ,int(values(1)) ! $10,90-$11,$4, $1}
                        lasted = values(1)
                     end if
                 else        ! Check out the format
                     if( index(line,'* Epoch').eq.1 .and. 
     .                   index(line,' P1 ').gt.0 ) then
                         newformat = .true.
                         voff = 2
                     endif 
                 end if
              enddo
              if( openprn ) close(60,iostat=kerr)
              if( openedt ) close(61,iostat=kerr)
          end if
      end do

****  Thats all
      end



                        
                           
                   



     
