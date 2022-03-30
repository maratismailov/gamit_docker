      Subroutine SET_PARA(wcmd,job)
c
c     setup parameter mode
c
c     job = 1:   estimate parameters
c     job =-1:   fix parameters
c
      include '../includes/dimpar.h'
      include 'solve.h'    
      include 'parameters.h'
c
      character*(*) wcmd
      character*16 code
      character*16 code_lower
      character*16 item
      character*20 label1,label2
      character*4 snam,upperc
      character*5 lowerc
      integer job,job2,job3,i0,i2,ic,type
      integer ib,i,match_name,id
      integer count_arg,lift_arg

      logical debug/.false./

c     zero out the snam variable
      snam='    '
c
c set job2 and job3 depending whether fixing of freeing parameters
c job2 = fixing/freeing a group of parameters on all sites/satellites
c job3 = fixing/freeing a group of parameters at a specific site/satellite
c -ve job indicates fix parameters, +ve job indicates estimate parameters
      if(job.gt.0) then
        job2 = 2
        job3 = 3
      else
        job2 = -2
        job3 = -3
      endif                                 
c
c     decompose the command line
      ic = count_arg(wcmd)
c     not enough arguments
      if (ic.le.1) goto 200
c     pointer to
      ib = lift_arg(wcmd,code,1)
      if (ib.le.0) goto 200
      code_lower = code
      call lowers(code_lower)
c
c     determine the type of parameter to be estimated/fixed
c     type=1 global, type=2 all_sit, type=3 all_sat, type=4 spacific parameter choice,
c     ie site name or satellite number.
      type = 0
      if (code_lower(1:6).eq.'global') then
        type = 1
      elseif (code_lower(1:6).eq.'all_si') then
        type = 2
      elseif (code_lower(1:6).eq.'all_sa') then
        type = 3
      else
        type = 4
        snam = upperc(code(1:4))
      endif
c
c     read site name, or sat number and match its index
      do 150 i2 = 2,ic
         ib = lift_arg(wcmd,item,i2)
         if (ib.le.0) goto 150
c        all parameters
         if (lowerc(item(1:5)).eq.'all_p') then
            do i = 1,ntpart
            if(debug.and.i.eq.615) print *,'i 615 job2 islot1 free '
     .            ,job2,islot1(i),free(i)
              if (type.eq.1.and.islot1(i).ge.80001.and.islot1(i)
     .            .le.80006) then
                if(abs(free(i)).le.abs(job)) free(i) = job 
                print *,'Setting EOP ',i,free(i),job
                
              elseif (type.eq.2.and.
     .            ((islot1(i).ge.1.and.islot1(i).le.500).or.
     .            (islot1(i).ge.2501.and.islot1(i).le.2700).or.
     .            (islot1(i).ge.2901.and.islot1(i).le.29000))) then
                  if(abs(free(i)).le.abs(job)) free(i) = job

              elseif (type.eq.3.and.(islot1(i).ge.501.and.
     .                               islot1(i).le.2700) ) then
                  if(abs(free(i)).le.abs(job)) free(i) = job

              else
                  label1 = rlabel(i)
                  i0 = match_name(1,4,label1,snam)
                 if (i0.gt.0.and.abs(free(i)).le.abs(job3)) free(i)=job3
               endif
            if(debug.and.i.eq.615) print *,'i 615 job2 islot1 free '
     .            ,job2,islot1(i),free(i)
            enddo
         endif
c
c        coordinates
         if (lowerc(item(1:5)).eq.'coord')then
            do 20 i = 1,ntpart
               if((islot1(i).gt.300)) go to 20
               if (type.eq.2) then
                  if(abs(free(i)).le.abs(job2)) free(i) = job2
               else
                  label1 = rlabel(i)
                  i0 = match_name(1,4,label1,snam)
                 if (i0.gt.0.and.abs(free(i)).le.abs(job3)) free(i)=job3
               endif
 20         continue
         endif
c        zenith delay
         if (lowerc(item(1:5)).eq.'zenit' .or. 
     .       lowerc(item(1:5)).eq.'atmos') then
            do 30 i = 1,ntpart
* MOD TAH 190615: Based on other construct changes here to replace
*     .gt.<start_slop>-1 with .ge.<start_slot> where in this case
*     start_islot is 300, replaced .gt. with .ge. but with same 301 value.
              if( (islot1(i).ge.301.and.islot1(i).le.400) .or.
     .             islot1(i).gt.11501.and.islot1(i).le.24000) then
                if (type.eq.2) then
                  if(abs(free(i)).le.abs(job2)) free(i) = job2
                else
                  label1 = rlabel(i)
                  i0 = match_name(1,4,label1,snam)
                  if(i0.gt.0.and.abs(free(i)).le.abs(job3)) free(i)=job3
                endif
              endif
 30         continue
         endif
c        atmospheric gradients   
         if (lowerc(item(1:5)).eq.'grad ') then
            do 35 i = 1,ntpart
            if (islot1(i).ge.24001.and.islot1(i).le.29000) then
               if (type.eq.2) then
                  if(abs(free(i)).le.abs(job2)) free(i) = job2
               else
                  label1 = rlabel(i)
                  i0 = match_name(1,4,label1,snam)
                 if (i0.gt.0.and.abs(free(i)).le.abs(job3)) free(i)=job3
               endif
            endif
 35         continue
         endif
c        receiver clock
         if (lowerc(item(1:5)).eq.'clock') then
            do 40 i = 1,ntpart
               if((islot1(i).lt.401).or.(islot1(i).gt.500)) go to 40
               if (type.eq.2) then
                  if(abs(free(i)).le.abs(job2)) free(i) = job2
               else
                  label1 = rlabel(i)
                  i0 = match_name(1,4,label1,snam)
                 if (i0.gt.0.and.abs(free(i)).le.abs(job3)) free(i)=job3
               endif
 40         continue
         endif
c        satellite orbital initial conditions
         if (lowerc(item(1:5)).eq.'orbit') then
            do 50 i = 1,ntpart
               if((islot1(i).lt.501).or.(islot1(i).gt.1100)) go to 50
               if (type.eq.3) then
                  if(abs(free(i)).le.abs(job2)) free(i) = job2
               else
                  label1 = rlabel(i)
                  id = lift_arg(label1,label2,4)
                  i0 = match_name(1,4,label2,snam)
                 if (i0.gt.0.and.abs(free(i)).le.abs(job3)) free(i)=job3
               endif
 50         continue
         endif
c        satellite orbital non-gravitational force (radiation) parameters
         if (lowerc(item(1:5)).eq.'radia') then
            do 60 i = 1,ntpart
               if((islot1(i).lt.1101).or.(islot1(i).gt.2400)) go to 60
               if (type.eq.3) then
                  if(abs(free(i)).le.abs(job2)) free(i) = job2
               else
                  label1 = rlabel(i)
                  id = lift_arg(label1,label2,4)
                  i0 = match_name(1,4,label2,snam)
                 if (i0.gt.0.and.abs(free(i)).le.abs(job3)) free(i)=job3
               endif
 60         continue
         endif
c        satellite antenna offsets      
         if (lowerc(item(1:5)).eq.'svant') then  
            do 80 i = 1,ntpart     
               if((islot1(i).lt.2401).or.(islot1(i).gt.2800)) go to 80
               if (type.eq.3) then
                  if(abs(free(i)).le.abs(job2)) free(i) = job2
                  if(debug) print *,'SET_PARA svant i islot1 free job2 '
     .                     ,  i,islot1(i),free(i),job2  
               else
                  label1 = rlabel(i)
                  id = lift_arg(label1,label2,4)
                  i0 = match_name(1,4,label2,snam)
                 if (i0.gt.0.and.abs(free(i)).le.abs(job3)) free(i)=job3
               endif
 80         continue
         endif
c        L1-bias parameters      
         if (lowerc(item(1:5)).eq.'l1-bi') then
            do 90 i = 1,ntpart         
               if((islot1(i).lt.2901).or.(islot1(i).gt.7400)) go to 90
               if (type.eq.2) then
                  if(abs(free(i)).le.abs(job2)) free(i) = job2
               else
                  label1 = rlabel(i)
                  i0 = match_name(1,4,label1,snam)
                 if (i0.gt.0.and.abs(free(i)).le.abs(job3)) free(i)=job3
               endif    
 90         continue
         endif                 
c        L2-bias parameters
         if (lowerc(item(1:5)).eq.'l2-bi') then
            do 100 i = 1,ntpart
               if((islot1(i).lt.7401).or.(islot1(i).gt.11900)) go to 100
               if (type.eq.2) then
                  if(abs(free(i)).le.abs(job2)) free(i) = job2
               else
                  label1 = rlabel(i)
                  i0 = match_name(1,4,label1,snam)
                 if (i0.gt.0.and.abs(free(i)).le.abs(job3)) free(i)=job3
               endif
 100        continue
         endif
c        X, Y, pole position
         if (lowerc(item(1:5)).eq.'wob  '.or.
     .       lowerc(item(1:5)).eq.'wobbl') then
           do 110 i = 1,ntpart
             if((islot1(i).ne.80001).and.(islot1(i).ne.80003)) go to 110
             if (type.eq.1) then
                 if(abs(free(i)).le.abs(job2)) free(i) = job2
             endif
 110       continue
         endif
c        X, Y, pole rates
         if (lowerc(item(1:5)).eq.'wob_r') then
           do 120 i = 1,ntpart
             if((islot1(i).ne.80002).and.(islot1(i).ne.80004)) go to 120
             if (type.eq.1) then
                if(abs(free(i)).le.abs(job2)) free(i) = job2
             endif
 120       continue
         endif
c        UT1 correction
         if (lowerc(item(1:5)).eq.'ut1  ') then
            do 130 i = 1,ntpart
               if((islot1(i).ne.80005)) go to 130
               if (type.eq.1) then
                  if(abs(free(i)).le.abs(job2)) free(i) = job2
               endif
 130         continue
         endif
c        UT1 rate
         if (lowerc(item(1:5)).eq.'ut1_r') then
            do 140 i = 1,ntpart
               if((islot1(i).ne.80006)) go to 140
               if (type.eq.1) then
                  if(abs(free(i)).le.abs(job2)) free(i) = job2
               endif
 140         continue
         endif
 150  continue
c
 200  continue  
      return
      end

