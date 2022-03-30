      Subroutine set_dcb_flag( sitcod,rcvcod, pcncod, rxver, illi
     .                       , dcb_override, idcb) 

c     Determine from the receiver type and whether A/S is on whether corrections need to
c     be made to the C1 and P2 observables to account for biases between correlating 
c     and non-correlating receivers.  This routine sets a flag for the data line of the
c     X-file, which MODEL then converts into a setting for bits 28 and 29 of the AUTCLN 
c     data_flag on the C-file.  The actual corrections are determined and made in MODEL
c     by reading file p1c1bias, which records the values estimated by AIUB monthly.

c     R. King 30 September 2004

      implicit none                                              

c  Input
c  -----                            

c    Site code (for report_stat message)
       character*4 sitcod
           
c    Receiver code
       character*6 rcvcod

c    Code from rcvant.dat indicating the type of receiver  
       character*1 pcncod     
c         'P' : receiver is cross-correlating and requires correction of P2' and C1 
c                 Rogue SNR, Trimble 4000, etc.
c         'C' : receiver is non-cross-correlating but reports C1 instead of P1
c                 Trimble 4700, 5700, Leica RS500, CRS1000, SR9600, etc. 
c                  unless AS is off
c         'N' : receiver is non-cross-correlating and reports true P1, P2

                         
c    RINEX file version (needed to check loss-of-lock indicator)
       real*4 rxver

c    Loss-of-lock indicator from RINEX file
       integer*4 illi 
c        Bit 2 ('4') set if AS is on; some receivers are cross-correlating under AS
c        but non-cross-correlating with AS off.  (Trimble 4000, 4700, 5700)    

c    Flag to override rcvant.dat and illi if CC2NONCC was run on RINEX file
       logical dcb_override
c         false by default, true if CC2NONCC run
               

c  Output
c  ------

c    Pseudo-range code for X-file
       integer*4 idcb
c         0 : No corrections needed (non-cross-correlating receiver)
c         1 : Correct C1 only; set bit 28 of data_flag
c         2 : Correct C1 and P2; set bit 29 of data_flag
                  
c   External function
c   -----------------
      logical kbit
                     
c   Local
c   -----
      logical as           
      character*256 message  


c   See if AS is on 

      if( rxver.ge.2.0 ) then
         if( kbit(illi,3) ) then
             as = .true.
         else
             as = .false.
         endif
      else
c       if RINEX version 1, can't tell: assume AS on since for early data, 
c       we won't  have corrections anyway
         as = .true.
      endif
                 
      if( pcncod.eq.'N' ) then
          
        idcb = 0

      elseif( pcncod.eq.'C' ) then 

        idcb = 1

      elseif( pcncod.eq.'P' ) then
               
c        make correction unless AS off and certain receivers

         if( .not.as ) then   
           if( rcvcod.eq.'TRMSST' .or. rcvcod.eq.'TRMSSE' .or.
     .         rcvcod.eq.'TRMSSI' .or. rcvcod.eq.'TR4SIS' .or.
     .         rcvcod.eq.'TR4700' .or.
     .         rcvcod.eq.'TR5700' ) then  
             idcb = 0  
           else
            idcb = 2
           endif
         else   
           idcb = 2    
         endif  
       
      elseif( dcb_override ) then   
        idcb = 0

      else  

        write(message,'(a4,a,a6,a)') sitcod
     .   ,' PCN-code missing for receiver type '
     .   ,rcvcod,' in rcvant.dat, cannot set C1/P2 correction flag'
         call report_stat('FATAL','MAKEX','set_dcb_flag',' ',message,0)  

      endif   
                         
      return
      end


        

