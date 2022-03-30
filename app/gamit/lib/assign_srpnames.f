      Subroutine assign_srpnames(srpmod,nics,srpnam)

c     Assign the 4-character parameter names based on the satellite radiation 
c     pressure model requested

      implicit none

      
      integer*4 nics,len,rcpar,i

      character*4  ecom1_par(9),ecom2_par(9),ecomc_par(13),srpnam(13)
      character*5 srpmod
      character*80 prog_name
      character*80 message 


      data ecom1_par/'DRAD','YRAD','BRAD','DCOS','DSIN','YCOS','YSIN'
     .             , 'BCOS','BSIN'/                      
                               
      data ecom2_par/'DRAD','YRAD','BRAD','BCOS','BSIN'
     .             ,'D2CS','D2SN','D4CS','D4SN'/                      

      data ecomc_par/'DRAD','YRAD','BRAD','DCOS','DSIN','YCOS','YSIN'
     .             , 'BCOS','BSIN','D2CS','D2SN','D4CS','D4SN'/                      


c  Get calling program name for report_stat

      len =  rcpar(0,prog_name)

                           
c  Set the parameter names according to the model

c   ---at present all models have the constant and one-per-rev terms avaiable 
      if( srpmod.eq.'BERNE'.or.srpmod.eq.'ECOM1'.or.srpmod.eq.'UCLR1' 
     .    .or.srpmod.eq.'UCLR1'.or.srpmod.eq.'ECOM2'
     .    .or.srpmod.eq.'ECOMC') then
        nics = 15
        do i=1,9 
          srpnam(i) = ecom1_par(i)
        enddo
        if( srpmod.eq.'ECOM2' ) then
          do i=1,9
            srpnam(i) = ecom2_par(i)
          enddo
        endif 
        if( srpmod.eq.'ECOMC' ) then
          nics = 19
          do i=10,13
            srpnam(i) = ecomc_par(i)
          enddo
        endif
      else                      
       write(message,'(a,a5,a)') 'SRP model (',srpmod,') not recognized'
           call report_stat('FATAL',prog_name,'lib/assign_srpnames'
     .                     ,' ',message,0)
      endif

      return
      end 
