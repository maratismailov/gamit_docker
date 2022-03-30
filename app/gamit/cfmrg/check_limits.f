      Subroutine check_limits( nprm,nlab,nbin,maxbin )

c     Check that the number of total parameters, C-file labels, 
c     and parameter types,do not exceed the maximum dimensions
                      
      implicit none

      include '../includes/dimpar.h'
   
      integer*4 nprm,nlab,nbin,maxbin
      character*256 message


c      Check the parameter index vs maxprm (in dimpar.h)

      if( nprm.gt.maxprm ) then
         write(message,'(a,i4,a,i4)') 'Too many parameters: n=',nprm
     .                            ,' maxprm=',maxprm
         call report_stat('FATAL','CFMRG','check_limits',' ',message,0)
      endif

c        Check the label index vs maxlab (in dimpar.h)

      if( nlab.gt.maxlab ) then
         write(message,'(a,i4,a,i4)') 'Too many labels: n=',nlab
     .                            ,' maxlab=',maxlab
         call report_stat('FATAL','CFMRG','check_limits',' ',message,0)
      endif
            
c        Check the bin index vs maxbin 

      if( nbin.gt.maxbin ) then
         write(message,'(a,i4,a,i4)') 'Too many bins: nbin=',nbin
     .                            ,' maxbin=',maxbin
         call report_stat('FATAL','CFMRG','check_limits',' ',message,0)
      endif


      return
      end

