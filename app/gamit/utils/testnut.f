      Program TESTNUT

c     Temporary program to test new coding for the nutation
c      rwk 180327

      implicit none
     
      integer*4 inut,ioerr
      logical fcheck             
                          
c     primary variables
      real*8 pjd,tjd,dpsi,deps
c     dummy for MHB_2000 call
      real*8 dpsi_ls,deps_ls,dpsi_plan,deps_plan,dpsi_fcn,deps_fcn
     .     , dpsi_prec,deps_prec

      inut = 1
     
      pjd = 2458181.25d0 
      pjd = 2458181.25d0
      tjd = 2458180.750d0 

      if( fcheck('nbody') ) then
        print *,'nbody exists, call MHB_2000 '
        call MHB_2000(tjd,dpsi_ls,deps_ls
     .                   ,dpsi_plan, deps_plan
     .                   ,dpsi_fcn,  deps_fcn
     .                   ,dpsi_prec, deps_prec
     .                   ,dpsi      ,deps )
        print *,'MHB_2000 tjd dpsi deps ',tjd,dpsi,deps        
        print *,' deps_plan deps_fcn deps_prec '
     .         ,  deps_plan,deps_fcn,deps_prec 
      else 
        print *,'nbody link empty, call nutred  '
        OPEN (UNIT=INUT,FILE='nutabl.',STATUS='OLD',iostat=ioerr)
        if (ioerr .ne. 0) then
           call report_stat('FATAL','UTILS','testnut','nutabl.',
     .    'Error opening nutation table',ioerr)
        endif
        call nutred(inut,pjd,dpsi,deps)
        print *,'nutabl pjd dpsi,deps ',pjd,dpsi,deps 
      endif

      stop
      end 

     



