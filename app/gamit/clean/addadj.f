Copyright (c) Massachusetts Institute of Technology and the University of
California at San Diego, 1995. All rights reserved.

      real*8 function addadj( isat,dt,iepoch,norb,nparts,tmpart
     .                      , coord_index,atm_index1,atm_index2
     .                      , clock_index,orb_index,svant_index
     .                      , eop_index,grad_index,grad_parts )

C     Return the correction to O-C to get post-fit residuals

c     Rewritten by R. King 28 Aug 95 using pre-computed pointers to avoid
c     looping through the adjust array for every data point.

c  Input:
c     isat                  =  the sat index
c     dt                    =  time from the initial epoch in seconds
c     iepoch                =  epoch number   
c     norb                  =  number of orbit + SV ant partials (3,6,9,12,15,or 18)
c     nparts                =  number of partials on C-file
c     tmpart(maxlab)        =  partials from C-file
c     Index in adjust array for the first parameter in each group for input site and sat
c       coord_index         - station coordinates (3)
c       atm_index1          - average zenith delay parameter (1)
c       atm_index2          - time-dependent zenith delay parameters (numzen)
c       clock_index         - station clock parameters (1)
c       orb_index           - orbital parameters (norb = 9, 12 or 15)
c       svant_index         - satellite antenna offsets (3)
c       eop_index           - Earth rotation parameters (6)
c       grad_index(2)       - N/S and E/W atmospheric gradient parameters (2)
c       grad_parts          - computed atmospheric gradient partials (2)

c   In common /zendel/
c     numzen =  number of zenith-delay parameters per site
c     zenmod =  model for zenith-delay parameters
c     idtzen =  break epochs for zenith-delay model

c C-file slots for partials
c
c      tmpart(1)  ! latitude
c      tmpart(2)  ! longitude
c      tmpart(3)  ! radius(height)
c      tmpart(4)  ! atmosphere
c      tmpart(5)  ! clock epoch
c      tmpart(6)  ! x
c      tmpart(7)  ! y
c      tmpart(8)  ! z
c      tmpart(9)  ! xdot
c      tmpart(10) ! ydot
c      tmpart(11) ! zdot
c      tmpart(12) ! radiation 1  (direct)
c      tmpart(13) ! radiation 2  (y-bias)
c      tmpart(14) ! radiation 3  (x, z, or b-vec bias--see arc/ertorb.f)
c      tmpart(15) ! radiation 4  (cos 1/rev direct or sin 1/rev x-axis)
c      tmpart(16) ! radiation 5  (sin 1/rev direct or sin 3/rev x-axis)
c      tmpart(17) ! radiation 6  (cos 1/rev y-axis or sin 1/rev z-axis)
c      tmpart(18) ! radiation 7  (sin 1/rev y-axis)
c      tmpart(19) ! radiation 8  (cos 1/rev b-axis)
c      tmpart(20) ! radiation 9  (sin 1/rev b-axis)
c      tmpart(21) | sv antenna offset x-axis 
c      tmpart(22) | sv antenna offset y-axis
c      tmpart(23) | sv antenna offset z-axis
c      tmpart(24) ! ut1
c      tmpart(25) ! ut1dot
c      tmpart(26) ! xpole
c      tmpart(27) ! xpole dot
c      tmpart(28) ! ypole
c      tmpart(29) ! ypole dot

      implicit none

      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
c
      integer*4 isat,iepoch,coord_index,atm_index1,atm_index2
     .        , clock_index,nparts,orb_index,svant_index,eop_index
     .        , grad_index(2),eoppart,i,len,rcpar,norb
      real*8 update,dt,ff,tmpart(maxlab),grad_parts(2)
      character*3 upperc
      character*80 prog_name

      logical debug
      data debug/.false./

                         
      if( debug ) then
        print *,'ADDADJ isat dt iepoch nparts tmpart '
     .        ,       isat,dt,iepoch,nparts,(tmpart(i),i=1,nparts)
        print *,' indices: crd atm clk svant eop orb ',
     .  coord_index,atm_index1,atm_index2,clock_index,svant_index
     .   ,eop_index,orb_index
      endif

c-----Initialize the increment

      update = 0.d0

c-----Update for site coordinates

      if( coord_index.gt.0 ) then
c         latitude
        update = update + tmpart(1) * adjust(coord_index)
c         longitude
        update = update + tmpart(2) * adjust(coord_index+1)
c         radius
        update = update + tmpart(3) * adjust(coord_index+2)
c          if( iepoch.eq.78 ) then
c           print *,'coord_index tmpart adjust update '
c     .            ,coord_index,tmpart(3),adjust(coord_index+2),update
c          endif
      endif

c------Update for zenith delays

      if( atm_index1.gt.0 ) then
c       average zenith delay
        update = update + tmpart(4) * adjust(atm_index1)   
        if( atm_index2.gt.0 ) then 
c         multiple zenith delays
          do i=1,numzen-1
c          --constant (step) model
           if( upperc(zenmod).eq.'CON' ) then
             if( iepoch.ge.idtzen(i).and.iepoch.lt.idtzen(i+1) ) then
c                for all epochs except the last, check if GE tabular pt 'i' and LT tabular pt i+1
c                (this logic avoids double adding at tabular points)
                 update = update + tmpart(4) * adjust(atm_index2-1+i)
             else if( i.eq.numzen-1 .and. iepoch.eq.idtzen(i+1) ) then
c                special logic for last point, skipped in other test
                 update = update + tmpart(4) * adjust(atm_index2-1+1)
             endif
           elseif (upperc(zenmod).eq.'PWL' ) then
c             see comment above for span logic
              if( iepoch.ge.idtzen(i).and.iepoch.lt.idtzen(i+1) ) then
                  ff = dfloat(iepoch-idtzen(i)) /
     .                   dfloat(idtzen(i+1)-idtzen(i))

c          if( iepoch.ge.74 .and. iepoch.le.78 .and.
c    .         ((atm_index2.ge.22 .and. atm_index2.le.25)
c    .         .or. (atm_index2.ge.42.and.atm_index2.le.45)) .and.
c    .         (isat.eq.5 .or. isat.eq.2) ) then
c                  print *,'atm_index2 i ff tmpart update '
c    .                 ,atm_index2,i,ff,tmpart(4),update
c           endif

                  update = update
     .                + (1.d0-ff) * tmpart(4) * adjust(atm_index2-1+i)
     .                       + ff * tmpart(4) * adjust(atm_index2-1+i+1)

              elseif( i.eq.numzen-1 .and. iepoch.eq.idtzen(i+1) ) then
                  update = update * tmpart(4) * adjust(atm_index2-1+i+1)

c
c          if( iepoch.ge.74 .and. iepoch.le.78 .and.
c    .         ((atm_index2.ge.22 .and. atm_index2.le.25)
c    .         .or. (atm_index2.ge.42.and.atm_index2.le.45)) .and.
c    .         (isat.eq.5 .or. isat.eq.2) ) then
c                print *,'  adjust1 adjust2 update',
c    .             update,adjust(atm_index2-1+i),adjust(atm_index2-1+i+1)
c           endif

               endif
           else
               len = rcpar(0,prog_name)
               call report_stat('STATUS',prog_name,'addadj',' ',
     .                         'Error, ZENMOD not defined in ADDADJ ',0)
           endif
          enddo
        endif
      endif

c------Update for atmospheric gradient parameters

      if( grad_index(1).gt.0 ) then  
c           N/S gradient
        if( numgrad.eq.1 ) then
          update = update + grad_parts(1) * adjust(grad_index(1))
c           print *,'grad_index(1) update ',grad_index(1),update  
        else 
c         multiple gradient parameters
          do i=1,numgrad-1
c          --constant (step) model
           if( upperc(gradmod).eq.'CON' ) then
             if( iepoch.ge.idtgrad(i).and.iepoch.lt.idtgrad(i+1) ) then
c                for all epochs except the last, check if GE tabular pt 'i' and LT tabular pt i+1
c                (this logic avoids double adding at tabular points)
                 update = update + tmpart(4) * adjust(grad_index(1)-1+i)
             else if( i.eq.numgrad-1 .and. iepoch.eq.idtgrad(i+1) ) then
c                special logic for last point, skipped in other test
                 update = update + tmpart(4) * adjust(grad_index(1)-1+1)
             endif
           elseif (upperc(gradmod).eq.'PWL' ) then
c            see comment above for span logic
             if( iepoch.ge.idtgrad(i).and.iepoch.lt.idtgrad(i+1) ) then
                  ff = dfloat(iepoch-idtgrad(i)) /
     .                   dfloat(idtgrad(i+1)-idtgrad(i))
c          if( iepoch.ge.74 .and. iepoch.le.78 .and.
c    .         ((grad_index(1).ge.22 .and. grad_index(1).le.25)
c    .         .or. (grad_index(1).ge.42.and.grad_index(1).le.45)) .and.
c    .         (isat.eq.5 .or. isat.eq.2) ) then
c                  print *,'grad_index1 i ff tmpart update '
c    .                 ,grad_index(1),i,ff,tmpart(4),update
c          endif
             endif
           endif
          enddo
        endif
      endif
      if( grad_index(2).gt.0 ) then
c         E/W gradient
        if( numgrad.eq.1 ) then
          update = update + grad_parts(2) * adjust(grad_index(2))
c          print *,'grad_index(2) update ',grad_index(2),update
        else
c         multiple gradient parameters
          do i=1,numgrad-1
c          --constant (step) model
           if( upperc(gradmod).eq.'CON' ) then
             if( iepoch.ge.idtgrad(i).and.iepoch.lt.idtgrad(i+1) ) then
c                for all epochs except the last, check if GE tabular pt 'i' and LT tabular pt i+1
c                (this logic avoids double adding at tabular points)
                 update = update + tmpart(4) * adjust(grad_index(2)-1+i)
             else if( i.eq.numgrad-1 .and. iepoch.eq.idtgrad(i+1) ) then
c                special logic for last point, skipped in other test
                 update = update + tmpart(4) * adjust(grad_index(2)-1+1)
             endif
           elseif (upperc(gradmod).eq.'PWL' ) then
c            see comment above for span logic
             if( iepoch.ge.idtgrad(i).and.iepoch.lt.idtgrad(i+1) ) then
                 ff = dfloat(iepoch-idtgrad(i)) /
     .                  dfloat(idtgrad(i+1)-idtgrad(i))
c          if( iepoch.ge.74 .and. iepoch.le.78 .and.
c    .         ((grad_index(2).ge.22 .and. grad_index(2).le.25)
c    .         .or. (grad_index(2).ge.42.and.grad_index(2).le.45)) .and.
c    .         (isat.eq.5 .or. isat.eq.2) ) then
c                  print *,'grad_index2 i ff tmpart update '
c    .                 ,grad_index(2),i,ff,tmpart(4),update
c          endif
             endif
           endif
          enddo
        endif
      endif

c------Update for site clock

      if( clock_index.gt.0 ) then
c         epoch
        update = update + tmpart(5) * adjust(clock_index)
c        if( iepoch.eq.78 ) print *,'after clock ',update
c** adjustments for rate and acceleration no longer supported
cc         rate
cc        update = update + tmpart(5) * dt * adjust(clock_index+1)
cc         acceleration
cc        update = update + 0.5d0*tmpart(5)*dt*dt/86400.d0 *
cc     .               adjust(clock_index)
c        print *,'clock_index update ',clock_index,update
      endif

c-----Update for orbital parameters

      if( orb_index.gt.0 ) then
        do i=1,norb
          update = update + tmpart(5+i) * adjust(orb_index-1+i)  
c          if( iepoch.eq.78 ) print *,'norb update ',i,update
        enddo
c        update = update + tmpart(6) * adjust(orb_index)
c     .                  + tmpart(7) * adjust(orb_index+1)
c     .                  + tmpart(8) * adjust(orb_index+2)
c     .                  + tmpart(9) * adjust(orb_index+3)
c     .                  + tmpart(10)* adjust(orb_index+4)
c     .                  + tmpart(11)* adjust(orb_index+5)
c        update = update + tmpart(12)* adjust(orb_index+6)
c     .                  + tmpart(13)* adjust(orb_index+7)
c     .                  + tmpart(14)* adjust(orb_index+8)
c        if( nparts.gt.20 ) then
c          update = update + tmpart(15) * adjust(orb_index+9)
c     .                    + tmpart(16) * adjust(orb_index+10)
c     .                    + tmpart(17) * adjust(orb_index+11)
c     .                    + tmpart(18) * adjust(orb_index+12)
c     .                    + tmpart(19) * adjust(orb_index+13)
c     .                    + tmpart(20) * adjust(orb_index+14)
c        print *,'orb_index update ',orb_index,update
c        endif
      endif

c-----Update for SV antenna offsets

      if( svant_index.gt.0 ) then 
        do i=1,3
          update = update + tmpart(5+norb+i) * adjust(svant_index-1+i)   
c         if( iepoch.eq.78 ) then
c          print *,'svant_index tmpart ',svant_index-1+i,tmpart(5+norb+i)  
c         endif
        enddo
      endif

c-----Update for Earth orientation parameters

      if( eop_index.gt.0 ) then
         eoppart = nparts - 5
         update = update + tmpart(eoppart)   * adjust(eop_index)
     .                   + tmpart(eoppart+1) * adjust(eop_index+1)
     .                   + tmpart(eoppart+2) * adjust(eop_index+2)
     .                   + tmpart(eoppart+3) * adjust(eop_index+3)
     .                   + tmpart(eoppart+4) * adjust(eop_index+4)
     .                   + tmpart(eoppart+5) * adjust(eop_index+5)
c      if (iepoch.eq.78)  print *,'eop_index nparts update '
c     .      ,eop_index,nparts,update
      endif

      addadj = update

      return
      end
