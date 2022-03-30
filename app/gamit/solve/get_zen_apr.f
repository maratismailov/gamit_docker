Copyright 1993 Massachusetts Institute of Technology and The Regents of the
c              University of California, San Diego.  All Rights Reserved.

      Subroutine GET_ZEN_APR ( indzen )

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
      integer indzen

      real*8 temp(3)    


c      Input:  rlabel(maxlab)  in common /parameter/  - parameter labels
c              indzen          index of first zenith-delay  in parameter array

c      Output: zen_apr, zen_apr2, zen_mar, zen_mar2, zen_tau, zen_tau2
c              -- a priori and Markov constraints (and time constrant) for zenith delays


c      Get number of zenith-delay parameters and type of model

      call getcmd(5,'aprior',wcmd,lcmd,2)
      call getcmd(5,'zenith',wcmd,lcmd,3)
      if (lcmd.le.0 ) then
        write(message,'(2a)')
     .   'No zenith-delay controls in batch file, set defaults:'
     .  ,'  zenmod=con  nzen=1 '
        call report_stat('WARNING','SOLVE','get_zen_apr',' ',message,0)
        zenmod = 'con'
        nzen = 1
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
             read (wcmd(ib+1:lcmd),*) nzen
             if ( ic.ge.3) then
               ib = lift_arg(wcmd,code,3)
               if (ib.le.0 )
     .           call report_stat('FATAL','SOLVE','get_zen_apr',' '
     .                      ,'zenmod not found ',0)
               zenmod = code(1:3)
               call uppers(zenmod)
               if( zenmod.ne.'CON'.and.zenmod.ne.'PWL') then
                  write(message,'(2a)')
     .                 'Zenith-delay model not recognized--'
     .               , 'piecewise linear (PWL) assumed'
                  call report_stat('WARNING','SOLVE','get_zen_apr',' '
     .                            , message,0)
                endif
              endif
        else
           write(message,'(2a)') 'Temporary restriction: '
     .                 , 'all sites must have equal nzen'
           call report_stat('WARNING','SOLVE','get_zen_apr',' '
     .                     , message,0)
           read (wcmd(ib+1:lcmd),*) nzen
        endif
        if( zenmod.eq.'PWL'.and.nzen.eq.1) then
          write(message,'(2a)')
     .          'Piecewise-linear model not allowed'
     .         ,'  with NZEN=1; constant used'
          call report_stat('WARNING','SOLVE','get_zen_apr',' '
     .                            , message,0)
          zenmod = 'CON'
        endif
        if( nzen.gt.maxatm ) then
           write(message,'(a,i2,a,i2,a)')
     .       'Number of zenith-delay parameters (',nzen
     .      ,') greater than MAXATM (',maxatm,')'
           call report_stat('FATAL','SOLVE','get_zen_apr',' ',message,0)
        endif
      endif
 100  continue

                   
c      Set default constraints
     
       do i=1,nsite
          zen_apr(i) = 0.5d0   
          zen_apr2(i) = 0.5d0
          zen_mar(i) = 0.02d0   
          zen_mar2(i) = 0.02d0
          zen_tau(i) = 100.d0  
          zen_tau2(i) = 100.d0
       enddo

c      Read in constraints for first ('tight') solution

      i2 = 0
      call getcmd(5,'aprior',wcmd,lcmd,2)
      do 150 i0 = 1,1000
         call getcmd(5,'tight_apr_ze',wcmd,lcmd,3)
         if (lcmd.gt.0) then
           i2 = i2+1
c          decompose command line
           ic = count_arg(wcmd) 
c          not enough arguments
           if (ic.lt.4) goto 145
c          get the first argument and count the number of characters 
           ib = lift_arg(wcmd,code,1)   
           if (ib.le.0) goto 145
           if (upperc(code(1:4)).eq.'ALL_') then
              read (wcmd(ib+1:lcmd),*) (temp(j),j=1,ic-1)
              if ( temp(1).eq.0.d0 .or. temp(2).eq.0 ) 
     .         call report_stat('FATAL','SOLVE','get_zen_apr',' ',
     .     'Zen delay tight constraints for ALL are 0. in batch file',0)
              do i=1,nsite
                zen_apr(i) = temp(1)
                zen_mar(i) = temp(2)   
                zen_tau(i)=temp(3)  
              enddo
              goto 155
           else
              snam = upperc(code(1:4))
              read (wcmd(ib+1:lcmd),*) (temp(j),j=1,ic-1)
              if ( temp(1).eq.0.d0 .or. temp(2).eq.0 ) then
               write(message,'(2a,a4,a)') 'Zen delay tight constraints '
     .           ,'are 0. for site ',snam,' in batch file'
                call report_stat('FATAL','SOLVE','get_zen_apr',' '
     .                          , message,0)  
              endif
              i1 = indzen  
              do i=1,nsite      
                 buf4=upperc(rlabel(i1)(1:4))    
                 if (buf4.eq.snam) then
                    zen_apr(i) = temp(1)
                    zen_mar(i) = temp(2)
                    zen_tau(i)=temp(3)   
                    goto 150
                  else              
                       i1= i1 + 1
                  endif
              enddo 
c             ignore constraint in batch file not in station list
           endif  
c          end if on ALL vs station constraint
           goto 150
         endif 
c        end if on valid command
         goto 150
 145    write(message,'(3a)') 
     .    'Missing or incomplete a priori tight constraints for zenith '
     .    ,'delays in batch file--setting defaults:  0.5 m '
     .    ,' 0.02 m/sqrt(hr)  100 hr'
         call report_stat('WARNING','SOLVE','get_zen_apr',' ',message,0)
c     150: found one or none, continue searching
 150  continue            
c     155: found all
 155  if( logprt ) write( 6,160) zenmod
      write(10,160) zenmod
 160  format(/
     1  ,'   A priori zenith delay     Model = ',a3,/
     2  ,'Station                  # A priori (m) '
     2  ,'Markov (m/sqrt(hr))  Correlation time (hrs)',/)
      do i=1,nsite
          if( logprt) write( 6,165) i,rlabel((i-1)*3+1)(1:4),sitnam(i)
     .              ,nzen,zen_apr(i),zen_mar(i),zen_tau(i)
          write(10,165) i,rlabel((i-1)*3+1)(1:4),sitnam(i),nzen
     .                  ,zen_apr(i),zen_mar(i),zen_tau(i)
 165      format(i3,2x,a4,1x,a12,1x,i3,3x,f8.3,8x,f8.3,12x,f8.3)
      enddo

c      Get a priori for loosely constrained solution

      i2 = 0
      call getcmd(5,'aprior',wcmd,lcmd,2)
      do 250 i0 = 1,1000
         call getcmd(5,'loose_apr_ze',wcmd,lcmd,3)
           if ( lcmd.gt.0 ) then
           i2 = i2+1
c          decompose command line
           ic = count_arg(wcmd)
c          not enough arguments
           if (ic.lt.4) goto 245  
c          get the first argument and count the number of characters 
           ib = lift_arg(wcmd,code,1)
           if (ib.le.0) goto 245 
           if (upperc(code(1:4)).eq.'ALL_') then
              read (wcmd(ib+1:lcmd),*) (temp(j),j=1,ic-1)
             if ( temp(1).eq.0.d0 .or. temp(2).eq.0 ) 
     .         call report_stat('FATAL','SOLVE','get_zen_apr',' ',
     .     'Zen delay loose constraints for ALL are 0. in batch file',0)
              do i=1,nsite
                zen_apr2(i) = temp(1)
                zen_mar2(i) = temp(2)   
                zen_tau2(i)=temp(3)  
              enddo
              goto 300
           else
              snam = upperc(code(1:4))
              read (wcmd(ib+1:lcmd),*) (temp(j),j=1,ic-1)
              if ( temp(1).eq.0.d0 .or. temp(2).eq.0 ) then
                write(message,'(a,a4,a)') 'Zen delay tight constraints '
     .           ,'are 0. for site ',snam,' in batch file'
                call report_stat('FATAL','SOLVE','get_zen_apr',' '
     .                          , message,0)  
              endif
              i1 = indzen
              do i=1,nsite
                 buf4=upperc(rlabel(i1)(1:4))
                 if (buf4.eq.snam) then
                    zen_apr2(i) = temp(1)
                    zen_mar2(i) = temp(2)
                    zen_tau2(i)=temp(3)  
                    goto 250
                  else
                    i1= i1 + 1
                  endif
              enddo 
c             ignore constraint in batch file not in station list
           endif 
           goto 250
        endif 
        goto 250
 245    write(message,'(3a)') 
     .    'Missing or incomplete a priori loose constraints for zenith '
     .    ,'delays in batch file--setting defaults:  0.5 m '
     .    ,' 0.02 m/sqrt(hr)  100 hr'
        call report_stat('WARNING','SOLVE','get_zen_apr',' ',message,0)   
c     250: found one or none, continue searching
 250  continue       
c     300: found 'all', no more to do
 300  continue
      return
      end
