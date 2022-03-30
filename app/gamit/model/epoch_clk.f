         Subroutine EPOCH_CLK( ipass,iter,jdsend,tsend,ichan,none )
c
c PURPOSE: subroutine to compute the receiver clock correction from
c         satellite clocks, pseudo-range, and the theoretical delay
c
c PARAMETERS:
c         IN: ipass   : clock iteration pass counter                    I*4
c             iepoch  : epoch number                                    I*4 model.h
c             iprnt   : unit number for pfile output file               I*4 model.h
c             iter    : iteration number                                I*4
c             iuj     : unit number of jfile                            I*4 model.h
c             nchan   : number of satellites                            I*4 model.h
c             ischan  : vector of satellite numbers                     I*4(maxsat) model.h
c             jdsend  : julian day of the current epoch send time       I*4
c             tsend   : time of day in secs of current epoch send time  R*8
c             ishan   : satellite array index                           I*4
c             klock   : type of receiver clock correction requested     I*4 model.h
c             jfiln   : jfile name                                      C*16 model.h
c             none    : debug flag                                      C*4
c             obs     : observed psuedo ranges                          R*8(maxdat,maxsat) model.h
c             vel_light: velocity of light in m/s                       R*8  constat_param.h
c             svantbody: SV type, needed only for GPS BLOCK I           C*20 model.h
c             rclock0 : mean rclock correction from iteration one       R*8  model.h
c             delay   : theoretical delay from coords and sat position  R*8(maxdat,maxsat) model.h
c
c        OUT: ier     : flagged epoch bad data                          I*4(maxsat) model.h
c             rclock  : receiver clock correction at epoch.             R*8 model.h
c
c SUBROUTINES CALLED: readj,
c
c CREATED: 27th DEC 1993               LAST MODIFIED: 27th DEC 1993
c
c AUTHOR: rwk, kf, put in SR by Simon McClusky.
c
c COPYRIGHT: DEPARTMENT OF EARTH AND PLANETRY SCIENCES
c            M.I.T. 1993
c
         implicit none
c
      include '../includes/dimpar.h'
      include '../includes/units.h'
      include '../../libraries/includes/const_param.h'
      include '../includes/errflg.h'
      include '../includes/model.h'
c
         character*4 none
         character*256 message
c
         integer*4 ipass,iter,jdsend,icall,jdc,ichan
c
         real*8 tsend,svdt,tc,svcepc,svcrat,svcacc,valid,prangt
c
cd         print*,' EPOCH_CLK ichan  tsend  klock rclock0 delay '
cd     .         ,            ichan,tsend,klock,rclock0,delay(1,ichan)
         if(klock .eq. 3) then
c
c           get the 3 SV clock poly coeffs from J-file
c
c           secret back door for debugging:
c           Can avoid use of J-file by specifying NONE for name
                if(index(jfiln,none) .eq. 0) then
            icall= 1
            call readj ( iuj,ischan,nchan,ischan(ichan),
     .      jdsend,tsend,svdt,icall,jdc,tc,svcepc,svcrat,svcacc,valid)
                endif
            if (ipass .eq. 1) then
               prangt= obs(3,ichan)/(vel_light)

c       print*,' in epec_clk observed pseudorange is ',prangt
c              svdt may be contaminated by selective availability.
c              if we have a block II sat, then we will have a set of sv clock
c              parameters every epoch, and first order term is best.
               if (svantbody(ichan)(1:8).eq.'BLOCK I ' ) then
                  rclock = prangt - delay(1,ichan) + svdt
               else
c PT 950619: use satellite clock polynomial for all other SVs
c                  rclock = prangt - delay(1,ichan) + svcepc
                  rclock = prangt - delay(1,ichan) + svdt
               endif
            else
c              second pass, use averaged estimate
               rclock = rclock0
            endif
         endif

c         print*,' in epoch_clk prangt, delay,svdt ' ,prangt
c     .              ,delay(1,ichan), svdt

c        check receiver clock at this point
         if(dabs(rclock).gt.2.0d0 ) then
            write(iprnt,805) ischan(ichan),iepoch,rclock
805         format(1x,'MODEL warning in EPOCH_CLK: '
     .                 ,'receiver clock correction from PRN'
     .                 ,i2.2,' at epoch',i5,' is enormous: ',1pe12.4
     .                 ,' seconds!')
            write(message,810) ischan(ichan),iepoch,rclock
810         format('Receiver clock correction from PRN',i2.2,
     .            ' at epoch',i5,' is enormous: ',1pe12.4,' seconds!')
            call report_stat('WARNING','MODEL','epoch_clk',' '
     .      ,message,0)
c           if the clock correction is really big, set it back to
c           zero to avoid infinite loops in the table lookup routines.
            if(dabs(rclock).gt.1000.0d0) then
               write(iprnt,815) iepoch,ischan(ichan),iter
815            format(1x,'MODEL warning: '
     .                 ,'unweighting data at epoch'
     .                 ,i5,' from PRN',i2.2,' iter = ',i4)
               write(message,820) iepoch,ischan(ichan),iter
820            format('Outragous clock corrn, unweighting data at epoch'
     .        ,i5,' from PRN',i2.2,' iter = ',i4)
               call report_stat('WARNING','MODEL','epoch_clk',' '
     .         ,message,0)
               rclock = 0.0d0
               ier(ichan) = igunwt
            endif
         endif
         return
         end
