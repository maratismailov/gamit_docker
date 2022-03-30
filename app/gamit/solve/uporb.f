Copyright 1995 Massachusetts Institute of Technology and The Regents of the
c              University of California, San Diego.  All Rights Reserved.

      Subroutine UPORB
C
C     Update a G-file from an orbit solution.
C     Y.Bock 1986; R. King June 1987
C     Clean up to handle case of non-existant G-file: Kurt & Dong July 91
c     Add ability to handle additional orbital parameters: King June 95       

c     In this routine, isat, nsat, and isprn refer to satellites on the 
c     observation files, igsat and ngsat(not actually used) to the g-file

      include '../includes/dimpar.h'
      include 'solve.h'        
      include 'parameters.h' 
      include 'models.h'

      integer*4 igsat,ind,indx,ioerr,iparm,isat,i

      CHARACTER*3 END,UPPERC
      character*4 parm(maxorb)
      CHARACTER*16 SATNAM,BUF16
      CHARACTER*16 TMPNAM,ORBNAM,LOWERC
      CHARACTER*80 LINE,message

      logical gerror,fcheck

      character*4 ecom1_par(15),ecom2_par(15),ecomc_par(19)

      data end/'END'/

      data ecom1_par/'X   ','Y   ','Z   ','XDOT','YDOT','ZDOT'
     . ,'DRAD','YRAD','BRAD','DCOS','DSIN','YCOS','YSIN','BCOS','BSIN'/ 
      data ecom2_par/'X   ','Y   ','Z   ','XDOT','YDOT','ZDOT'
     . ,'DRAD','YRAD','BRAD','BCOS','BSIN','DCS2','DSN2','DCS4','DSN4'/
      data ecomc_par/'X   ','Y   ','Z   ','XDOT','YDOT','ZDOT'
     . ,'DRAD','YRAD','BRAD','DCOS','DSIN','YCOS','YSIN','BCOS','BSIN'
     . ,'DCS2','DSN2','DCS4','DSN4'/

c       Open the original G-file

      orbnam = ginf
      if (fcheck (orbnam)) then
         OPEN(UNIT=16,
     .        FILE=LOWERC(ORBNAM),
     .        STATUS='OLD',
     .        iostat=ioerr)
         gerror = .false.
         if (ioerr .ne. 0) then
            call report_stat('WARNING','SOLVE','update',tmpnam
     .          , 'Cannot open G-file--no update',ioerr)
            gerror = .true.
         endif
      else
           call report_stat('WARNING','SOLVE','update',tmpnam
     .          , 'Cannot find G-file--no update',0)
         gerror = .true.
      endif


c         Open the output G-file

      if (.not. gerror) then
         tmpnam=goutf
        open (unit   = 17,
     .        file   = lowerc(tmpnam),
     .        status = 'unknown',
     .        iostat = ioerr,
     .        err=920)


c           Write the G-file names to the log and Q-file

         if( logprt )write(6,'(/,2a,/,2a)' ) 'Updating G-file : ',ginf
     .                      , 'New G-file      : ',goutf
         if( iqflag.eq.1 )
     .   write(10,'(/,2a,/,2a)' ) 'Updating G-file : ',ginf
     .                       , 'New G-file      : ',goutf


c            Read and write the first line of the G-file

c        format:
c       (i2,1x,i3,1x,i2,1x,i2,1x,i2    20x,  a4,1x,a5,1x,a5,1x,a5  )
c        yy ddd hh mm ss                    GPST J2000 IAU76 BERNE
         read (16,fmt='(a80)') line
c        make explicit the model values
         if( line(36:39).eq.'    ' ) line(36:39)='UTC '
         line(41:45) = frame
         line(47:51) = precmod
         line(53:57) = srpmod    
         line(59:63) = nutmod
         line(65:69) = gravmod
         line(71:75) = eradmod
         line(77:80) = antradmod
         write (17,fmt='(a80)') line


c            Read the second line of the G-file and write it out with current parameters

         read(16,'(i2,1x,19(a4,1x))') iparm,(parm(i),i=1,iparm)
c         print *,'UPORB iparm parm ',iparm,(parm(i),i=1,iparm)
c         print *,' norb ',norb
         if( iparm.lt.norb ) then
           write(message,'(2a)') 'Input G-file has fewer parameters'
     .             ,' than solution: increase on output'
           write(10,'(a)') message
           call report_stat('WARNING','SOLVE','uporb',' ',message,0)
         elseif( iparm.gt.norb) then
           write(message,'(2a)') 'Input G-file has more parameters'
     .             ,' than solution: decrease on output'
           write(10,'(a)') message
           call report_stat('WARNING','SOLVE','uporb',' ',message,0)
         endif
c         print *,'UPORB srpmod ',srpmod
         if ( srpmod.eq.'BERNE' .or. srpmod.eq.'ECOM1'.or.
     .        srpmod.eq.'UCLR1' .or. srpmod.eq.'UCLR2'  ) then
              do i=1,15
                parm(i) = ecom1_par(i)
              enddo
         elseif ( srpmod.eq.'ECOMC' ) then 
           do i=1,19
             parm(i) = ecomc_par(i)
           enddo
         else
           write(message,'(2a)') 'Solar radiation pressure model ('
     .          ,srpmod,') not recognized '
           write(10,'(a)') message
           call report_stat('FATAL','SOLVE','uporb',' ',message,0)
         endif
         write(17,'(i2,1x,19(a4,1x))') norb,(parm(i),i=1,norb)


c             Read and write the comment lines of the G-file

150      read (16,fmt='(a80)') line
         if (upperc(line(1:3)).eq.upperc('END')) goto 160
         if( logprt ) write (6,'(1x,a80)') line
         write (17,'(a80)') line
         go to 150
160      if (gerror) goto 950
         write(17,fmt='(A80)') gline
         write(17,fmt='(a3)') end


c             Determine the index of orbital parameters in the adjust array

         indx = 0
         do i=1,ntpart
           if( islot1(i).ge.501.and.islot1(i).le.2400 ) goto 200
           indx = indx + 1
         enddo
  200    continue


C        Loop over the satellites on the input G-file
C        updating the parameters for any that were observed

         do 500 ind=1,maxsat

           read (16,'(a16)',end=915,iostat=ioerr,err=915) satnam
           write(buf16,'(a16)') satnam

c          read satellite name - convert to satellite number
           if (upperc(satnam(1:3)).eq.upperc('END')) then
              write(17,'(a3)') 'END'
              goto 950
           else if( upperc(satnam(1:1)) .eq. upperc('P') ) then
C             Label is 'PRN nn' -pre-GNSS naming convention 
              read(buf16,'(4x,i2)') igsat  
           else 
c             Label is GNSS code followed by PRN number
              read(buf16,'(1x,i2)') igsat
           endif

           write (17,'(a16)') satnam

C          Loop over the observed satellites
C
           do isat=1,nsat
C            If we've updated the sat in SOLVE, rewrite the params
             if( isprn(isat).eq.igsat ) then
               do i=1,norb
                  if( i.le.iparm ) then
c                   advance one line in G-file
                    read(16,fmt='(a80)') line
                  endif
                  write(17,'(d20.13,60x)')
     .               postvl(indx+norb*(isat-1)+i)
               enddo
c              --break out of the channel loop once the satellite is found
               goto 500
             endif
           enddo
c           --if satellite not adjusted, copy the old parameters
           do i=1,norb
              if( i.le.iparm ) then
                  read(16,fmt='(a80)') line
               else
                 line(1:20) = '                0.d0'
               endif
                  write(17,fmt='(a80)') line
           enddo

C        --End loop over G-file satellites
  500    continue


c        Need to write 'END' to output file if NSAT=MAXSAT
         write(17,'(a3)') 'END'
      else
c        Could not find G-file
c        still need to read 2 lines from input
c        READ(5,FMT='(A80)') LINE
c        READ(5,FMT='(A80)') LINE
         call report_stat('WARNING','SOLVE','uporb',orbnam
     .          , 'Cannot find G-file--no upodate',0)
      endif


c     come here on decoding error of SATNAM
 915  continue
      if (ioerr .ne. 0) then
        call report_stat('WARNING','SOLVE','uporb',orbnam
     .          , 'Cannot open G-file--no update',ioerr)

        write(message,'(a,a16)')
     .     'Missing END or bogus SATNAM on input G-file: ',satnam
           call report_stat('WARNING','SOLVE','uporb',' ',message,ioerr)
      endif


c     come here on end of file or error
 920  continue
      if (ioerr .ne. 0) then
         call report_stat('WARNING','SOLVE','update',orbnam
     .          , 'Error or end in G-file--no update',ioerr)
      endif

950   CLOSE(UNIT=16)
      CLOSE(UNIT=17)

      RETURN
      END

