Copyright (c) Massachusetts Institute of Technology and the University of
California at San Diego, 1995.  All rights reserved.
       subroutine rot_gfile(gfname,arc_frame,gframe,arc_prec,gfprec
     .                      ,moveold,iugout,iuin,program
     .                      ,te,jde,tdtoff,inut,iut1,ipole )

c
c  subroutine to rotate a gfile in "gframe" with "gfprec" model into "arc_frame"
c  with "arc_prec" model
c
c  input gfile will be moved to a different filename, and a new one created in the
c  requested reference frame and precession model
c
c
c  P. Tregoning 18th April, 1995
c
c  Modified by PT 9th May, 1995 to allow rotations between any combinations of
c  B1950, J2000, IAU68 and IAU76
c
c  Modified by PT 27th July, 1995 to remove common statements, to be called by GTOG
c     as well as ARC, and moved to the library.
c
c INPUT:  gfname               : input gfile name
c         arc_frame, arc_prec  : requested inertial frame and precession model
c         gframe,gfprec        : current gfile frame and precession model
c         moveold              : logical - T to replace old gfile, F to create new one
c         iugout               : unit number of output gfile (if moveold = .false.)
c         iuin                 : unit number of input gfile
c         program              : name of calling program
c
c         jde                  : julian day number of ics epoch
c         te                   : seconds of day of ics epoch
c
c         inut,iut1,ipole      : nutabl. ut1. and pole. unit numbers

       implicit none

      include '../includes/dimpar.h'

       character*4 program
       character*5 arc_frame,arc_prec,gframe,gfprec
       character*10 gfname
       character*19 dattim
       character*22 gfold
       character*80 junk(6),line,prog_name
       character*256 message
       real*8 satic(3),satic_v(3),satic_e(3),satic_v_e(3),rot(3,3)
       real*8 rotdot(3,3),sidtm,xpole,ypole
       real*8 rotsav(3,3),rotdif(3,3)
c      real*8 tempve(3)
       integer k,blank,iugout,iout,iuin
       integer ierr,j,i,iut1pol
       logical moveold,svflag/.false./

       real*8 te,tdtoff
       integer*4 jde

       integer*4 inut,iut1,ipole,len,rcpar

c       print *,'In ROT_GFILE gfname,arc_frame,gframe,arc_prec,gfprec '
c     .                     , gfname,arc_frame,gframe,arc_prec,gfprec

c      get name of calling program for report_stat
       len = rcpar(0,prog_name)


       gfold = ' '
       gfold(1:10) = gfname


c rewind input gfile
       rewind(iuin)

c set output unit numbers, depending on whether call to this subroutine
c requires a gfile to be replaced (ARC) or a new one created (GTOG)
       if(moveold)then
         iout = 51
       else
         iout = iugout
       endif

       if(moveold)then
c  determine whether frame, or precession model or both need changing
         if(arc_frame.ne.gframe)then
           gfold(11:16) = '.'//gframe
         endif
         if(arc_prec.ne.gfprec)then
           blank = index(gfold,' ')
           gfold(blank:blank+5) = '.'//gfprec
         endif

c make a copy of gfile
         open(unit=50,file=gfold,status='unknown',form='formatted')
c can't have named scratch files when compiling for the HP!!! McClusky 950701
c        open(unit=51,file='gfile.tmp',status='scratch')
         open(unit=iout,status='scratch',form='formatted')

c blank out buffer arrays
         do i=1,6
          junk(i)=' '
         enddo
         line = ' '

c  rewind gfname and write it to gfold_name. 
         ierr = 0
         do while (ierr.eq.0)
           read(iuin,'(a80)',iostat=ierr)junk(1) 
           if( ierr.ne.-1 ) then 
             if( .not.svflag .and.junk(1)(1:3).eq.'END' ) then
               line(1:30)='Gfile moved by arc by ARC  '
               call runtim(dattim)
               line(31:51)='  '//dattim
               write(50,'(a80)') line
               write(50,'(a3)')'END'
               svflag = .true.
               line = ' '
             else
               write(50,'(a80)')junk(1)
             endif  
          else
            continue
          endif
        enddo
        call report_stat('WARNING',prog_name,'lib/rot_gfile',' '
     .     , 'Mismatch between requested frame/precession and G-file',0)
        write(message,'(a,a5,2x,a5,a,a5,2x,a5)')
     .  '(cont)  G-file: ',gframe,gfprec,'   Input: ',arc_frame,arc_prec
        call report_stat('WARNING',prog_name,'lib/rot_gfile'
     .                  ,' ',message,0)
        write(message,'(a,a10,a,a22)') '(cont) ',gfname,' moved to '
     .                 , gfold
        call report_stat('WARNING',prog_name,'lib/rot_gfile'
     .                  ,' ',message,0)
        rewind(iuin)
        svflag = .false.

      else
        write(message,'(2a,a5,a,a5)')'Creating gfile in inertial  '
     .           ,'frame - ',arc_frame,' Precession - ',arc_prec
        call report_stat('STATUS',prog_name,'lib/rot_gfile'
     .                  ,' ',message,0)
      endif

c  now read header line and add new frame and precession to line
      read(iuin,100)junk(1)
100   format(a80)
      junk(1)(41:45) = arc_frame
      junk(1)(47:51) = arc_prec

c  write out first line of gfile - now indicating ARC frame and precession
      write(iout,100)junk(1)

c  now read gfile until the END of comments is found
      do 1000 j=1,(15+10*maxsat)
        read(iuin,'(a80)',iostat=ierr)junk(1)

c if the end of file is reached, return gfile to original position in file and return
        if(ierr.ne.0) then

c if moveold = false then we are finished
           if(.not.moveold) return

c          copy temporary gfile to actual gfile
c          rewind gfname and write it to gfold_name.
           rewind(iuin)
           rewind(iout)
           ierr = 0
           do while (ierr.eq.0)
             read(iout,'(a80)',iostat=ierr)junk(1)
             write(iuin,'(a80)')junk(1)
           enddo
           write(message,120) gfname,arc_frame,arc_prec
120        format('Satellite IC''s in new ',a10
     .       ,' are now consistent with ',a5,' and ',a5,'.')
           call report_stat('STATUS',prog_name,'lib/rot_gfile',' '
     .                     , message,0)
           rewind(iuin)
           read(iuin,100)junk(1)
c          stop
           close(50)
           close(iout)
           return

        else

          if( .not.svflag .and.junk(1)(1:3).eq.'END' ) then
            write(line,125) program,arc_frame,arc_prec
125         format('ICs converted by ',a4,' to ',a5,' and ',a5
     .             ,' precession')
            call runtim(dattim)
            line(52:72)='  '//dattim
            write(iout,'(a80)') line
            svflag = .true.
          endif

          if( junk(1)(1:2).eq.'PR' .or. junk(1)(1:2).eq.'NS') then

c           write out the SV name line
            if( junk(1)(1:2).eq.'PR' ) then
              junk(1)(7:80) = ' '
            elseif ( junk(1)(1:2).eq.'NS' ) then
              junk(1)(5:80) = ' '
            endif
            write(iout,'(a80)') junk(1)

c           then read in the IC's of the satellite and convert to e-fixed using
c           gfile inertial frame and precession, and then to inertial using ARC
c           inertial frame and precession
            do i=1,6
              read(iuin,100)junk(i)
              if(i.le.3)then
                read(junk(i)(1:20),'(d20.0)')satic(i)
c               print*,'satic ',i,satic(i)
              else
                read(junk(i)(1:20),'(d20.0)')satic_v(i-3)
c               print*,'satic ',i,satic_v(i-3)
              endif
            enddo
c           rotate the coords
c           rotate from wrong inertial back to date of epoch
c           original coding has iut1pol set to zero (ie ignore) so we must
c           undo the same
* MOD TAH 200505: Read the sestbl. to get the correct value   
C           iut1pol = 0   ! Original line.
            call get_iut1pol( iut1pol ) 

c       print*,' calling rotsnp'
            call rotsnp( -1,jde,te,tdtoff,iut1pol    
     .                 , gframe,gfprec,iut1,ipole,inut
     .                 , rot,rotdot,sidtm,xpole,ypole ) 
c      print*,' called rotsnp'
c            print *,'1 rot ',rot
c            print *,'rotdot ',rotdot
            do i=1,3
              do k=1,3
                rotsav(i,k) = rot(i,k)
              enddo
            enddo
c            print *,'satic_v ',satic_v
            call matmpy( rot,satic,satic_e,3,3,1)
            call matmpy( rot,satic_v,satic_v_e,3,3,1)
c            call matmpy( rotdot,satic,tempve,3,3,1)
c            do i=1,3
c              satic_v_e(i) = satic_v_e(i) + tempve(i)
c            enddo
c           now rotate the Earth-fixed at date to the new inertial frame
c            print *,'2 rot ',rot
c            print *,'rotdot ',rotdot
            call rotsnp( 1,jde,te,tdtoff,iut1pol    
     .                 , arc_frame,arc_prec,iut1,ipole,inut
     .                 , rot,rotdot,sidtm,xpole,ypole ) 
            call matmpy( rot,satic_e,satic,3,3,1)
            call matmpy( rot,satic_v_e,satic_v,3,3,1)
c            print *,'3 rot ',rot
c            print *,'rotdot ',rotdot
c            print *,'ROT_GFILE stop '
c            stop
            do i=1,3
               do k=1,3
                rotdif(i,k) = rot(i,k)-rotsav(k,i)
c                print *,'i j rotdif ',i,k,rotdif(i,k)
               enddo
            enddo
c            print *,'satic_v ',satic_v
c            call matmpy( rotdot,satic_e,tempve,3,3,1)
c            do i=1,3
c               satic_v(i) = satic_v(i) + tempve(i)
c            enddo
c           write the rotated IC's back into the character buffer
            do i=1,3
              write(junk(i)(1:20),'(d20.13,60x)') satic(i)
              write(junk(i+3)(1:20),'(d20.13,60x)') satic_v(i)
            enddo
c           write the information back to the gfile
            do i=1,6
c              write(*,100)junk(i)
               write(iout,100)junk(i)
            enddo

          else
            write(iout,'(a80)') junk(1)
          endif

        endif

1000  continue
c     end loop on all records of input g-file

      end
