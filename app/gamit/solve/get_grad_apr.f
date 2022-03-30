Copyright 1993 Massachusetts Institute of Technology and The Regents of the
c              University of California, San Diego.  All Rights Reserved.

      Subroutine GET_GRAD_APR( indgrad )

      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'
      include 'parameters.h'

      character*256 message
      character*120 wcmd
      character*16 code
      character*4 snam,upperc,buf4

      integer i,j,i1,i2,i0,ib,ic
      integer lcmd,count_arg,lift_arg
      integer indgrad

      real*8 m2cyc,at10deg,temp(3)    


c      Input:  rlabel(maxlab)  in common /params/  - parameter labels
c              indgrad          index of first gradient parameter in parameter array

c      Output: grad_apr, grad_apr2, grad_mar, grad_mar2, grad_tau, grad_tau2
c              -- a priori and Markov constraints (and time constrant) for gradient delays

        
c      Constant scale factors  

      m2cyc = 1/(299792458.d0/1.57542d9)
c Changed as the G77 compiler does not support intrinsic trig functions working with degrees 
c      at10deg = 1/(sind(10.d0)*tand(10.d0)+0.003)
      at10deg = 1/(sin(10.d0/180.d0*pi)*tan(10.d0/180.d0*pi)+0.003)
 
             
c      Get number of gradient-delay parameters and type of model

      ngrad=1
      call getcmd(5,'aprior',wcmd,lcmd,2)
      call getcmd(5,'gradient',wcmd,lcmd,3)
      if (lcmd.le.0 ) then
        write(message,'(2a)')
     .   'No gradient-delay controls in batch file, set defaults:'
     .  ,'  gradmod=con  nungrad=1 '
        call report_stat('WARNING','SOLVE','get_grad_apr',' ',message,0)
        gradmod = 'con'
        goto 100
      else    
       ic = count_arg(wcmd)
       if (ic.le.0) goto 100
c      decompose command line
       ic = count_arg(wcmd)
c      not enough arguments
       if (ic.le.1) goto 100
       ib = lift_arg(wcmd,code,1)
       if (ib.le.0) goto 100
c      all sites                 
       if (upperc(code(1:4)).eq.'ALL_') then
             read (wcmd(ib+1:lcmd),*) ngrad
             if ( ic.ge.3) then
               ib = lift_arg(wcmd,code,3)
               if (ib.le.0 )
     .           call report_stat('FATAL','SOLVE','get_grad_apr',' '
     .                      ,'gradmod not found ',0)
               gradmod = code(1:3)
               call uppers(gradmod)
               if( gradmod.ne.'CON'.and.gradmod.ne.'PWL') then
                  write(message,'(2a)')
     .                 'Gradient model not recognized--'
     .               , 'piecewise linear (PWL) assumed'
                  call report_stat('WARNING','SOLVE','get_grad_apr',' '
     .                            , message,0)
                endif
              endif
        else
           write(message,'(2a)') 'Temporary restriction: '
     .                 , 'all sites must have equal nungrad'
           call report_stat('WARNING','SOLVE','get_grad_apr',' '
     .                     , message,0)
           read (wcmd(ib+1:lcmd),*) ngrad
        endif
        if( gradmod.eq.'PWL'.and.ngrad.eq.1) then
          write(message,'(2a)')
     .          'Piecewise-linear model not allowed'
     .         ,'  with NGRAD=1; constant used'
          call report_stat('WARNING','SOLVE','get_grad_apr',' '
     .                            , message,0)
          gradmod = 'CON'
        endif
        if( ngrad.gt.maxgrad/2 ) then
           write(message,'(a,i2,a,i2,a)')
     .       'Number of gradient-delay parameters (',ngrad
     .      ,') greater than MAXGRAD/2 (',maxgrad/2,')'
          call report_stat('FATAL','SOLVE','get_grad_apr',' ',message,0)
        endif
        ngrad = ngrad
      endif
 100  continue

                   
c      Set default constraints
     
       do i=1,nsite 
         do j=1,2
           grad_apr(i,j) = 0.03d0      
           grad_apr2(i,j) = 0.03d0 
           if( ngrad.eq.1 ) then
             grad_mar(i,j) = 0.d0   
             grad_mar2(i,j) = 0.d0
             grad_tau(i,j) = 0.d0  
             grad_tau2(i,j) =0.d0 
           else   
             grad_mar(i,j) = 0.01d0   
             grad_mar2(i,j) = 0.01d0
             grad_tau(i,j) = 100.d0  
             grad_tau2(i,j) = 100.d0 
           endif
         enddo
       enddo

c      Read in constraints for first ('tight') solution

      i2 = 0
      call getcmd(5,'aprior',wcmd,lcmd,2)
      do 150 i0 = 1,1000
        call getcmd(5,'tight_apr_gr',wcmd,lcmd,3)
        if (lcmd.gt.0) then
          i2 = i2+1
c         decompose command line
          ic = count_arg(wcmd) 
c         not enough arguments
          if (ic.lt.4) goto 145
c         get the first argument and count the number of characters 
          ib = lift_arg(wcmd,code,1)   
          if (ib.le.0) goto 145  
          if (upperc(code(1:4)).eq.'ALL_') then
            read (wcmd(ib+1:lcmd),*) (temp(j),j=1,ic-1)
            if ( temp(1).eq.0.d0 .or. temp(2).eq.0 ) 
     .      call report_stat('FATAL','SOLVE','get_grad_apr',' ',
     .            'Gradient constraints for ALL are 0. in batch file',0)
            do i=1,nsite
              do j=1,2
                grad_apr(i,j) = temp(1) 
                if( ngrad.gt.1 ) then
                  grad_mar(i,j) = temp(2)   
                  grad_tau(i,j) = temp(3) 
                endif
              enddo 
            enddo
            goto 155
          else
             snam = upperc(code(1:4))         
             read (wcmd(ib+1:lcmd),*) (temp(j),j=1,ic-1) 
             if ( temp(1).eq.0.d0 .or. temp(2).eq.0 ) then
               write(message,'(2a,a4,a)') 'Gradient tight constraints '
     .           ,'are 0. for site ',snam,' in batch file'
               call report_stat('FATAL','SOLVE','get_grad_apr',' '
     .                          , message,0)  
             endif   
             i1 = indgrad
             if( j.eq.2 ) i1 = indgrad + nsite*ngrad
             do i=1,nsite
                buf4=upperc(rlabel(i1)(1:4))
                if (buf4.eq.snam) then  
                  do j=1,2
                    grad_apr(i,j) = temp(1)
                    if( ngrad.gt.1 ) then
                      grad_mar(i,j) = temp(2)
                      grad_tau(i,j) = temp(3)  
                    endif
                  enddo
                  goto 150
                else               
                  i1= i1 + ngrad
                endif
             enddo  
c            ignore constraint in batch file not in station list
             goto 150  
c         end if on ALL vs station constraint
          endif 
          goto 150  
c       end if on valid command 
        endif
        goto 150
 145    write(message,'(3a)') 
     .   'Missing or incomplete a priori tight constraints for gradient'
     .    ,' delays in batch file--setting defaults:  0.03 m '
     .    ,' 0.01 m/sqrt(hr)  100 hr'
        call report_stat('WARNING','SOLVE','get_grad_apr',' ',message,0)
c     150: found one or none, continue searching  
 150  continue    
c     155: found all
 155  continue
c     use different formats for single and multiple gradient parameters
      if ( ngrad.gt.1 ) then
        if( logprt ) write(6,160)
        write(10,160)  
      else     
        if( logprt ) write(6,162)
        write(10,162)
      endif                       
      do i=1,nsite  
        if( ngrad.gt.1 ) then  
          if( logprt ) write( 6,161) i,rlabel((i-1)*3+1)(1:4),sitnam(i)
     ,       ,ngrad,(grad_apr(i,j),grad_mar(i,j),grad_tau(i,j),j=1,2)
          write(10,161) i,rlabel((i-1)*3+1)(1:4),sitnam(i),ngrad
     .                ,(grad_apr(i,j),grad_mar(i,j),grad_tau(i,j),j=1,2)  
        else    
          if( logprt ) write( 6,163) i,rlabel((i-1)*3+1)(1:4),sitnam(i)
     .                  ,(grad_apr(i,j),j=1,2)
          write( 10,163) i,rlabel((i-1)*3+1)(1:4),sitnam(i)
     .                  ,(grad_apr(i,j),j=1,2)
        endif
c       convert units 
        grad_apr(i,1)=grad_apr(i,1)/at10deg*m2cyc  
        grad_apr(i,2)=grad_apr(i,2)/at10deg*m2cyc 
        if( ngrad.gt.1 ) then
          grad_mar(i,1)=grad_mar(i,1)/at10deg*m2cyc  
          grad_mar(i,2)=grad_mar(i,2)/at10deg*m2cyc 
        endif
      enddo         
  160 format(/,
     .    3x,'A priori atmospheric gradient error at 10 degrees'
     .   ,' elevation angle',/,
     .   'Station                 # N/S: A priori(m) '
     .   ,' Mar(m/sqrt(hr)) Correl(hr)   E/W: A priori(m) '
     .  ,' Mar(m/sqrt(hr)) Correl(hr)',/)   
  161     format(i3,1x,a4,2x,a12,1x,i3,2(6x,f8.5,f16.5,f12.1))
  162 format(/,
     .    3x,'A priori atmospheric gradient error in meters at ',
     .       '10 degrees elevation angle',/,
     .       'Station                  North-South  East-West',/)
  163     format(i3,1x,a4,2x,a12,3x,2(3x,f8.5)) 

    
c      Get a priori for loosely constrained solution

      i2 = 0
      call getcmd(5,'aprior',wcmd,lcmd,2)
      do 250 i0 = 1,1000
         call getcmd(5,'loose_apr_gr',wcmd,lcmd,3)
           if ( lcmd.gt.0 ) then
             i2 = i2+1
c            decompose command line
             ic = count_arg(wcmd)
c            not enough arguments
             if (ic.lt.4) goto 245  
c            get the first argument and count the number of characters 
             ib = lift_arg(wcmd,code,1)
             if (ib.le.0) goto 245 
             if (upperc(code(1:4)).eq.'ALL_') then
                read (wcmd(ib+1:lcmd),*) (temp(j),j=1,ic-1)
               if ( temp(1).eq.0.d0 .or. temp(2).eq.0 ) 
     .           call report_stat('FATAL','SOLVE','get_grad_apr',' ',
     .     'Gradient loose constraints for ALL are 0. in batch file',0)
               do i=1,nsite
                  do j=1,2
                    grad_apr2(i,j) = temp(1)/at10deg*m2cyc 
                    if( ngrad.gt.1 ) then
                      grad_mar2(i,j) = temp(2)/at10deg*m2cyc   
                      grad_tau2(i,j) = temp(3)  
                    endif 
                  enddo
               enddo
               goto 300
             else
               snam = upperc(code(1:4))   
               read (wcmd(ib+1:lcmd),*) (temp(j),j=1,ic-1)
               if ( temp(1).eq.0.d0 .or. temp(2).eq.0 ) then
                 write(message,'(a,a4,a)') 'Gradient loose constraints '
     .           ,'are 0. for site ',snam,' in batch file'
                call report_stat('FATAL','SOLVE','get_grad_apr',' '
     .                          , message,0)  
               endif
               i1 = indgrad
               if( j.eq.2 ) i1 = indgrad + nsite*ngrad
               do i=1,nsite
                 buf4=upperc(rlabel(i1)(1:4))
                 if (buf4.eq.snam) then
                   do j=1,2
                     grad_apr2(i,j) = temp(1)/at10deg*m2cyc
                     if( ngrad.gt.1 ) then
                       grad_mar2(i,j) = temp(2)/at10deg*m2cyc
                       grad_tau2(i,j) = temp(3)  
                     endif
                   enddo  
                   goto 250
                 else
                   i1= i1 + ngrad
                 endif
               enddo  
c              ignore constraint in batch file not in station list  
               goto 250 
c            end if on ALL vs station constraint
             endif     
c          end if on valid command
           endif
           goto 250
 245       write(message,'(3a)') 
     .  'Missing or incomplete a priori loose constraints for gradient'
     .    ,' delays in batch file--setting defaults:  0.03 m '
     .    ,' 0.01 m/sqrt(hr)  100 hr'
        call report_stat('WARNING','SOLVE','get_grad_apr',' ',message,0)   
c     250: found one or none, continue searching
 250  continue       
c     300: found 'all', no more to do
 300  continue
      return
      end
