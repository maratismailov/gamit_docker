      Subroutine settyp ( debug,rxver,nobtyp,rxobtyp,iobtypx
     .                  , nwave2,nsat,rcvrsw,swver,ndat,dattyp,lambda )

c     Set the x-file header values for type of observable and wavelength factor
c     Trap some special cases.  

c     R. King 17 July 1989, revised and moved to makex 1 October 2015

      implicit none

      include '../includes/dimpar.h'    
      include '../includes/makex.h'

c INPUT

c   RINEX version
      real*4 rxver
c   Number and observables on the RINEX file
c     RINEX versions 1 and 2 have 2-character obervables, e.g. 'L1 ','P2 ',
c     'C2 '; RINEX version 3 has 3-character observables which complete
c     describe the signals, e.g. 'L1C','C2P','C2W'. 
      integer*4 nobtyp       
      character*3 rxobtyp(maxobt)
c   Indices in the rxobtyp array for the two phase and two pseudorange obervables 
c   to be written on the x-file
      integer*4 iobtypx(6)
c   Wavenlength factor for L2 phase from the RINEX header (L1 always assumed to be '1')
c   If = 0, not specified yet: If RINEX 3, determine from observable, if RINEX 1 or 2, 
c   assume to be 1 but issue warning.
      integer*4 nwave2
c   Number of satellites
      integer*4 nsat                                                                                                  
c   GAMIT receiver and firmware codes
      character*3 rcvrsw
      real*4 swver

      logical debug 
      
c OUTPUT             
c ** RWK 160113: The code for these has been temporarily removed in favor
c                of keepiing the 1 2 3 4  scheme for phase (higher, lower freq)
c                and pseudorange (higher, lower freq) since the only way 
c                these codes were used was in xtorx  and in cfmrg to see
c                if dattyp(3) is zero to indicate single-frequency.  In their
c                place, the actual 3-character observable codes are printed
c                at the end of the line for documention, used only by xtrox.

c   Integer codes for the data types on the X-file
c       1 : L1 phase  (all systems except Beidou)
c       2   L2 phase  (GPS, Glonass, QZSS)
c       3   L1 pseudo-range (all systems except Beidou)
c       4   L2 pseudo-range (GPS, Glonass, QZSS) 
c     ( 5   L1 C/A code pseudorange no longer seperately indicated)
c       6   E5a/L5 phase (GPS, Galileo, SBAS, QZSS)
c       7   E5a/L5 pseudo-range (GPS, Galileo, SBAS, QZSs)
c       8   E5b/B2/L7 phase (Galileo, Beidou)
c       9   E5b/B2/L7 pseuro-range (Galileo, Beidou)
c      10   E6/L6  phase (Galileo, QZSS)
c      11   E6/L6  pseudo-range (Galileo, QZSS)

      integer*4 ndat,dattyp(maxdat)
c   Wavelength code indicating whether the observable is present for each
c   satellite (channel), whether it is ambiguous, and whether
c   the ambiguity spacing is full or half wavelength.  All observables
c   are converted to full wavelength.
      integer*4 lambda(maxsat,maxdat)      
c        0   observable is not present for this channel 
c       -1   ambiguity is one wavelength (e.g. 19 cm for L1 phase)
c       -2   ambiguity is one-half wavelength (i.e., frequency-doubled phase from
c                from a codeless receiver)
c        1  no ambiguity
c        2  no ambiguity (original observable was half wavelength)

c LOCAL

      character*256 message
      integer*4 i,j


c   Initialize arrays
      do i=1,maxdat
         dattyp(i) = 0
         do j=1,maxsat
           lambda(j,i) = 0
         enddo
      enddo
                   
      ndat = 4
      if(debug) print *,'in SETTYP iobtypx rxobtyp ',iobtypx,rxobtyp
      dattyp(1) = 1             
      if( iobtypx(2).eq.0 ) then    
        dattyp(2) = 0 
        call report_stat('WARNING','MAKEX','settyp'
     .        ,' ','Single-frequency observations',0)
* MOD TAH 200512: Added L8 for Galileo E8 observables.
      elseif( rxobtyp(iobtypx(2))(1:2).eq.'L2' .or.
     .        rxobtyp(iobtypx(2))(1:2).eq.'L5' .or.
     .        rxobtyp(iobtypx(2))(1:2).eq.'L7' .or.
     .        rxobtyp(iobtypx(2))(1:2).eq.'L8' .or.
     .        rxobtyp(iobtypx(2))(1:2).eq.'L6' .or.
     .        rxobtyp(iobtypx(2))(1:2).eq.'L9' ) then
          dattyp(2) = 2           
c*      elseif( rxobtyp(iobtypx(2))(1:2).eq.'L5') then  
c*        dattyp(2) = 6
c*      elseif( rxobtyp(iobtypx(2))(1:2).eq.'L7') then
c*        dattyp(2) = 8   
c*      elseif( rxobtyp(iobtypx(2))(1:2).eq.'L6' ) then
c*        dattyp(2) = 10
      else     
        write(message,'(a,f4.2,a,a2)') 'Invalid RINEX ',rxver
     .      ,' L2 phase observable ',rxobtyp(iobtypx(2))
        call report_stat('FATAL','MAKEX','settyp',' ',message,0)
      endif                     
      dattyp(3) = 3   
      if( iobtypx(4).eq.0 ) then 
        dattyp(4) = 0 
        call report_stat('WARNING','MAKEX','settyp'
     .         ,' ','No L2 pseudorange: codeless tracking',0)
* MOD TAH 200512: Added C8 and fixed reporting error (Lx reported)
      elseif( rxobtyp(iobtypx(4))(1:2).eq.'C2' .or.
     .        rxobtyp(iobtypx(4))(1:2).eq.'P2' .or.
     .        rxobtyp(iobtypx(4))(1:2).eq.'C5' .or.
     .        rxobtyp(iobtypx(4))(1:2).eq.'C7' .or.
     .        rxobtyp(iobtypx(4))(1:2).eq.'C8' .or.
     .        rxobtyp(iobtypx(4))(1:2).eq.'C6' .or.
     .        rxobtyp(iobtypx(4))(1:2).eq.'C9' ) then
        dattyp(4) = 4           
c*      elseif( rxobtyp(iobtypx(4))(1:2).eq.'C5') then
c*        dattyp(4) = 7
c*      elseif( rxobtyp(iobtypx(4))(1:2).eq.'C7') then
c*        dattyp(4) = 9   
c*      elseif( rxobtyp(iobtypx(4))(1:2).eq.'C6' ) then
c*        dattyp(4) = 11                   
      else  
        write(message,'(a,f4.2,a,a2)') 'Invalid RINEX ',rxver
     .      ,' L2 pseudorange observable ',rxobtyp(iobtypx(4))
        call report_stat('FATAL','MAKEX','settyp',' ',message,0)
      endif     

c  Set the wavelength factor for L2 phase (all others full wavelength)

      do i=1,nsat
        lambda(i,1) = -1 
        lambda(i,2) = -1
        lambda(i,3) =  1
        lambda(i,4) =  1   
        if( dattyp(2).eq.0 ) then
           lambda(i,2) = 0 
        elseif( nwave2.eq.2 ) then
           lambda(i,2) = -2
        elseif( rxver.ge.3.0 .and. rxobtyp(iobtypx(2)).eq.'L2N' ) then
           lambda(i,2) = -2
        endif  
      enddo
      if( dattyp(4).eq.0 ) lambda(i,4) = 0 


c  ** Trap some special cases **
c        TEQC now uses L2 wavelength factor = 1 in the header for all 
c        receivers, and invokes factor = 2 for codeless receivers via
c        the LLI bit.  Rather than having to scan the files in advance,
c        we'll use the firmware version and observation types to infer
c        the L2 wavelength factor.  An incorrect header wavelength
c        factor can also arise from the Berne translators if it was set
c        incorrectly in the input control file (common with TRMSST). 
c        An exception to this is the Trimble 'serial P-code receiver, from
c        which C2 and full wavelength L2 will be available for non-AS satellites
c        (Block I) satellites.  Since these are usually only a small subset and
c        C2 is not useful for editing or bias-fixing, force these to be treated 
c        as codeless if the input firmware version is set (artificially) to 4.00. 
cc RWK 180522: Check for dattyp=0 instead of =4 
cc      if( dattyp(2).ne.0.and.dattyp(4).ne.4 .and. nwave2.eq.1 ) then      
      if( dattyp(2).ne.0.and.dattyp(4).eq.0 .and. nwave2.eq.1 ) then
        write(message,'(a)') 
     .     'No P-code but nwave2=1; reset nwave2 = 2 '
        call report_stat('WARNING','MAKEX','settyp'
     .                   ,' ',message,0)
        write(uinfor,'(2a)') 'WARNING: ',message 
        if( rcvrsw.eq.'ASH' .and.swver.ne.1.0 ) then
          write(message,'(a)') 
     .        'Apparently codeless Ashtech but swver.ne.1.0'
          call report_stat('WARNING','MAKEX','settyp'
     .                     ,' ',message,0)
          write(uinfor,'(2a)') 'WARNING: ',message    
        elseif( (rcvrsw.eq.'TRM' .or. rcvrsw.eq.'SSP' 
     .             .or. rcvrsw.eq.'SST') .and. swver.eq.4.0 ) then
          write(message,'(a)')
     .      'P2 and full wavelength L2 overridden by input swver 4.0' 
          call report_stat('WARNING','MAKEX','settyp'
     .                     ,' ',message,0)
          write(uinfor,'(2a)') 'WARNING: ',message  
          dattyp(3) = 5
          dattyp(4) = 0
        endif
        do i=1,nsat
          lambda(i,2) = -2
        enddo
      endif

c        It is also possible (and too frequent!) for users to translate 
c        codeless (e.g Ashtech L-XII or Trimble SST) data with the RINEX 
c        opt file input set for a code-tracking receiver (e.g SSE).  In 
c        In this case there will be no P2 data (or it will be all 0.) and
c        AUTCLN will throw out all the data.  For the Ashtech, use the
c        GAMIT firmware code to fix this.  For Trimble, we cannot unambiguously 
c        distinguish this case from that of a serial P-code SST, so issue a
c        warning here and also in AUTCLN.  
       if ( dattyp(4).eq.4 ) then  
         if((rcvrsw.eq.'TRM'.and.swver.eq.4.0 ) .or.
     .      (rcvrsw.eq.'ASH'.and.(swver.ge.1.0.and.swver.lt.2.0)) ) 
     .     then
           write(message,'(a,a3,f5.1,a)') 
     .        'P2 in header but firmware (',rcvrsw,swver
     .         ,') implies codeless--fix header'
           call report_stat('WARNING','MAKEX','settyp'
     .                       ,' ',message,0)
           write(uinfor,'(2a)') 'WARNING: ',message 
           if( rcvrsw.eq.'ASH') then
             dattyp(4) = 0
             do i=1,nsat
               lambda(i,2) = -2   
               lambda(i,4) = 0 
             enddo  
           endif
         endif
       endif

c        It is also possible for users to translate L1-only data (e.g. TR4600LS)
c        as dual-frequency so that the header has L1 L2 C1 P2 P1 D1 (e.g.).
c        For the TR4600LS case, we can detect this from the artificially assigned
c        firwmare code 2.99.
      if( rcvrsw.eq.'TRM'.and.swver.eq.2.99 ) then
        write(message,'(2a)')
     .       'Dual-freq observable codes erroneously assigned '
     .       ,'for L1-only TR4600 receiver: override'
        call report_stat('WARNING','MAKEX','settyp'
     .                      ,' ',message,0)
        write(uinfor,'(2a)') 'WARNING: ',message  
        do i=1,nobtyp
            if(dattyp(i).eq.2.or.dattyp(i).eq.4 ) dattyp(i)=0 
        enddo                   
        do i=1,nsat
          lambda(i,2) = 0 
          lambda(i,4) = 0 
        enddo
      endif
      if(debug)  print *,'end of settyp dattyp lamba(1) '
     .        ,dattyp,(lambda(1,i),i=1,4) 
      return
      end
