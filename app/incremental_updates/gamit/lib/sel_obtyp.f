      Subroutine sel_obtyp( gnsslf,nobtyp,rxobtyp,iobtypx )

c     From the full suite of observables available from a RINEX
c     file, select two phase and two pseudorange observables to
c     write on the x-file.   The priorities currently set (Sept 2015)
c     reflect the observables available from all satellites in a
c     Wettzell RINEX file from a Leica GR25 from 2014 096. 
c
c     There are some RINEX files produced by teqc that have both
c     C1 and P1 listed (and possibly C2 and P2) but where P1
c     is absent at all epochs (perhaps also only at some epochs).
c     To handle this case, set two additional flags as iobtypx(5-6)
c     as backups.
                                                                   
      implicit none
             
      include '../includes/makex.h'

c       SV system
      character*1 gnss
      character*(*) gnsslf

c       Number and codes for observables from the RINEX header
      integer*4 nobtyp
      character*3 rxobtyp(maxobt)
c     If RINEX 2, the 3rnd character will be blank

c       Indices in rxobtyp for the preferred observables 
c         (phase1, phase2, PR1, PR2, backup PR1, backup9 PR2)
      integer*4 iobtypx(6)
           
c       Local
      integer*4 len,i  
c     function
      integer*4 rcpar
      character*80 prog_name
      character*256 message

       

c        Get the program name for report_stat calls
      len = rcpar(0,prog_name)
      gnss = gnsslf(1:1)
  
c       Initialize the index array and RINEX observation type
* MOD TAH 200526: Dimensioned to 6 so initialize to 6 not 4.
      do i=1,6
        iobtypx(i) = 0 
      enddo                                                  
           
c       Select the preferred observables for the system requested
c       Each call to fill_obtypx will overwrite any previous assignment, so
c       order the calls in inverse order of priority. Two-character codes
c       apply to RINEX 2, 3-character to RINEX 3

      if( gnss.eq.'G') then                                  

c         Higher-frequency carrier phase        
        call fill_obtypx(1,'LB ',nobtyp,rxobtyp,iobtypx)  
        call fill_obtypx(1,'LA ',nobtyp,rxobtyp,iobtypx)  
        call fill_obtypx(1,'L1 ',nobtyp,rxobtyp,iobtypx)  
        call fill_obtypx(1,'L1C',nobtyp,rxobtyp,iobtypx)

c         Lower-freqeuncy carrier phase 
* MOD TAH 200511: Added L5 frequency choice (default L2)
        if( gnsslf.eq.'G5' ) then         
           call fill_obtypx(2,'L5 ',nobtyp,rxobtyp,iobtypx) 
           call fill_obtypx(2,'L5X',nobtyp,rxobtyp,iobtypx) 
           call fill_obtypx(2,'L5I',nobtyp,rxobtyp,iobtypx)
           call fill_obtypx(2,'L5Q',nobtyp,rxobtyp,iobtypx)
        else
           call fill_obtypx(2,'LC ',nobtyp,rxobtyp,iobtypx) 
           call fill_obtypx(2,'L2 ',nobtyp,rxobtyp,iobtypx) 
           call fill_obtypx(2,'L2Y',nobtyp,rxobtyp,iobtypx)
           call fill_obtypx(2,'L2S',nobtyp,rxobtyp,iobtypx)
           call fill_obtypx(2,'L2W',nobtyp,rxobtyp,iobtypx)
        endif

c         Higher-frequency pseudorange       
        call fill_obtypx(3,'CB ',nobtyp,rxobtyp,iobtypx) 
        call fill_obtypx(3,'CA ',nobtyp,rxobtyp,iobtypx) 
        call fill_obtypx(3,'C1 ',nobtyp,rxobtyp,iobtypx)   
        call fill_obtypx(3,'P1 ',nobtyp,rxobtyp,iobtypx)
        call fill_obtypx(3,'C1W',nobtyp,rxobtyp,iobtypx)   
        call fill_obtypx(3,'C1C',nobtyp,rxobtyp,iobtypx)

c         backup higher-frequency pseudorange
        call fill_obtypx(5,'C1 ',nobtyp,rxobtyp,iobtypx)  
        call fill_obtypx(5,'C1C',nobtyp,rxobtyp,iobtypx) 

c         Lower-frequency pseudorange                    
* MOD TAH 200511: Added L5 frequency choice (default L2)
        if( gnsslf.eq.'G5' ) then   
           call fill_obtypx(4,'C5 ',nobtyp,rxobtyp,iobtypx)
           call fill_obtypx(4,'C5X',nobtyp,rxobtyp,iobtypx)
           call fill_obtypx(4,'C5I',nobtyp,rxobtyp,iobtypx)
           call fill_obtypx(4,'C5Q',nobtyp,rxobtyp,iobtypx)
        else   ! Default selection
           call fill_obtypx(4,'CC ',nobtyp,rxobtyp,iobtypx)
           call fill_obtypx(4,'C2 ',nobtyp,rxobtyp,iobtypx)
           call fill_obtypx(4,'P2 ',nobtyp,rxobtyp,iobtypx)
           call fill_obtypx(4,'C2Y',nobtyp,rxobtyp,iobtypx)
           call fill_obtypx(4,'C2S',nobtyp,rxobtyp,iobtypx)
           call fill_obtypx(4,'C2W',nobtyp,rxobtyp,iobtypx) 
c            backup lower-frequency pseudorange
           call fill_obtypx(6,'C2 ',nobtyp,rxobtyp,iobtypx)
           call fill_obtypx(6,'C2S',nobtyp,rxobtyp,iobtypx)
           call fill_obtypx(6,'C2X',nobtyp,rxobtyp,iobtypx)
         endif

      elseif( gnss.eq.'R' ) then  
c         Higher-frequency carrier phase 
        call fill_obtypx(1,'L1 ',nobtyp,rxobtyp,iobtypx)
        call fill_obtypx(1,'L1C',nobtyp,rxobtyp,iobtypx)

c         Lower-freqeuncy carrier phase
        call fill_obtypx(2,'L2 ',nobtyp,rxobtyp,iobtypx)  
        call fill_obtypx(2,'L2C',nobtyp,rxobtyp,iobtypx)
        call fill_obtypx(2,'L2P',nobtyp,rxobtyp,iobtypx)

c         Higher-frequency pseudorange       
        call fill_obtypx(3,'C1 ',nobtyp,rxobtyp,iobtypx)  
        call fill_obtypx(3,'P1 ',nobtyp,rxobtyp,iobtypx)
        call fill_obtypx(3,'C1C',nobtyp,rxobtyp,iobtypx)
c         backup higher-frequency pseudorange
        call fill_obtypx(5,'C1 ',nobtyp,rxobtyp,iobtypx)   

c         Lower-frequency pseudorange
        call fill_obtypx(4,'P2 ',nobtyp,rxobtyp,iobtypx)   
        call fill_obtypx(4,'C2 ',nobtyp,rxobtyp,iobtypx)  
        call fill_obtypx(4,'C2P',nobtyp,rxobtyp,iobtypx)
        call fill_obtypx(4,'C2C',nobtyp,rxobtyp,iobtypx) 
cd        print *,'R nobtyp rxobtyp iobtypx ',nobtyp,rxobtyp,iobtypx
c         backup lower-frquency pseudorange
* MOD TAH 200526: Changed C2 to P2 because C2 would be selected
*       above if it is available.
        call fill_obtypx(6,'P2 ',nobtyp,rxobtyp,iobtypx)     
        call fill_obtypx(6,'C2P',nobtyp,rxobtyp,iobtypx)          
cd        print *,'back R nobtyp rxobtyp iobtypx ',nobtyp,rxobtyp,iobtypx
cd        stop 
      elseif( gnss.eq.'C' ) then                

c         Higher-frequency carrier phase                  
        call fill_obtypx(1,'L1 ',nobtyp,rxobtyp,iobtypx)
        call fill_obtypx(1,'L1I',nobtyp,rxobtyp,iobtypx)  
        call fill_obtypx(1,'L1X',nobtyp,rxobtyp,iobtypx)  
c         RINEX 3.03 change
        call fill_obtypx(1,'L2I',nobtyp,rxobtyp,iobtypx)   
        call fill_obtypx(1,'L2X',nobtyp,rxobtyp,iobtypx)   

c         Lower-freqeuncy carrier phase  -- force to B2 ('7') for now
* MOD TAH/RWK 200223: Changed choice to L6 (B3) for lower frequency
*       (was L7 etc).  We keep L7 option since GAMUT should 
*       be able to handle mutltiple frequencies although ambiguity
*       resolution may be affected.
* MOD TAH 200511: Added L7 frequency choice (default L6). This
*       change undoes 200223 mod above so that frequencies are
*       no longer mixed.
        if( gnsslf.eq.'C7' ) then   
           call fill_obtypx(2,'L7 ',nobtyp,rxobtyp,iobtypx)
           call fill_obtypx(2,'L7X',nobtyp,rxobtyp,iobtypx)  
           call fill_obtypx(2,'L7I',nobtyp,rxobtyp,iobtypx)    
           call fill_obtypx(2,'L7Q',nobtyp,rxobtyp,iobtypx)    
* MOD MAF 210729: Added L5 lower frequency choice
        elseif( gnsslf.eq.'C5' ) then   
           call fill_obtypx(2,'L5 ',nobtyp,rxobtyp,iobtypx)
           call fill_obtypx(2,'L5X',nobtyp,rxobtyp,iobtypx)  
           call fill_obtypx(2,'L5P',nobtyp,rxobtyp,iobtypx)  
        else 
* Higher prority to L6 (highest priority at end list 
           call fill_obtypx(2,'L6 ',nobtyp,rxobtyp,iobtypx)
           call fill_obtypx(2,'L6I',nobtyp,rxobtyp,iobtypx)
           call fill_obtypx(2,'L6A',nobtyp,rxobtyp,iobtypx)    
           call fill_obtypx(2,'L6X',nobtyp,rxobtyp,iobtypx)   
        endif 
       
c         Higher-frequency pseudorange       
        call fill_obtypx(3,'C1 ',nobtyp,rxobtyp,iobtypx)
        call fill_obtypx(3,'C1I ',nobtyp,rxobtyp,iobtypx)  
        call fill_obtypx(3,'C1X ',nobtyp,rxobtyp,iobtypx)  
c         RINEX 3.03 change 
        call fill_obtypx(3,'C2I ',nobtyp,rxobtyp,iobtypx) 
        call fill_obtypx(3,'C2X ',nobtyp,rxobtyp,iobtypx)

c         Lower-frequency pseudorange
* MOD TAH 200511: Added L7 frequency choice (default L6). This
*       change undoes 200223 mod above so that frequencies are
*       no longer mixed.
        if( gnsslf.eq.'C7' ) then   
           call fill_obtypx(4,'C7 ',nobtyp,rxobtyp,iobtypx)
           call fill_obtypx(4,'C7X',nobtyp,rxobtyp,iobtypx)
           call fill_obtypx(4,'C7I',nobtyp,rxobtyp,iobtypx)
           call fill_obtypx(4,'C7Q',nobtyp,rxobtyp,iobtypx)
* MOD MAF 210729: Added C5 lower frequency choice
        elseif( gnsslf.eq.'C5' ) then   
           call fill_obtypx(4,'C5 ',nobtyp,rxobtyp,iobtypx)
           call fill_obtypx(4,'C5X',nobtyp,rxobtyp,iobtypx)
           call fill_obtypx(4,'C5P',nobtyp,rxobtyp,iobtypx)
         else
* NOD TAH/RWK Add C6 series of pseudo range
* Higher prority to L6 (highest priority at end list( 
           call fill_obtypx(4,'C6 ',nobtyp,rxobtyp,iobtypx)
           call fill_obtypx(4,'C6I',nobtyp,rxobtyp,iobtypx)
           call fill_obtypx(4,'C6A',nobtyp,rxobtyp,iobtypx)    
           call fill_obtypx(4,'C6X',nobtyp,rxobtyp,iobtypx) 
         endif   
 

      elseif( gnss.eq.'E' ) then
c         Higher-frequency carrier phase 
        call fill_obtypx(1,'L1 ',nobtyp,rxobtyp,iobtypx)  
        call fill_obtypx(1,'L1X',nobtyp,rxobtyp,iobtypx)
        call fill_obtypx(1,'L1C',nobtyp,rxobtyp,iobtypx)
    
c         Lower-freqeuncy carrier phase
* MOD TAH 200511: Added L6, L7,L8 frequency choice (default L5)
        if( gnsslf.eq.'E6' ) then   
           call fill_obtypx(2,'L6A',nobtyp,rxobtyp,iobtypx)   
           call fill_obtypx(2,'L6B',nobtyp,rxobtyp,iobtypx)
           call fill_obtypx(2,'L6C',nobtyp,rxobtyp,iobtypx)
           call fill_obtypx(2,'L6X',nobtyp,rxobtyp,iobtypx)
           call fill_obtypx(2,'L6Z',nobtyp,rxobtyp,iobtypx)
        elseif( gnsslf.eq.'E7' ) then   
           call fill_obtypx(2,'L7 ',nobtyp,rxobtyp,iobtypx)   
           call fill_obtypx(2,'L7X',nobtyp,rxobtyp,iobtypx)
           call fill_obtypx(2,'L7I',nobtyp,rxobtyp,iobtypx)
           call fill_obtypx(2,'L7Q',nobtyp,rxobtyp,iobtypx)
        elseif ( gnsslf.eq.'E8' ) then
           call fill_obtypx(2,'L8 ',nobtyp,rxobtyp,iobtypx)   
           call fill_obtypx(2,'L8X',nobtyp,rxobtyp,iobtypx)
           call fill_obtypx(2,'L8I',nobtyp,rxobtyp,iobtypx)
           call fill_obtypx(2,'L8Q',nobtyp,rxobtyp,iobtypx)
        else   ! Default L5
           call fill_obtypx(2,'L5 ',nobtyp,rxobtyp,iobtypx)   
           call fill_obtypx(2,'L5X',nobtyp,rxobtyp,iobtypx)
           call fill_obtypx(2,'L5I',nobtyp,rxobtyp,iobtypx)
           call fill_obtypx(2,'L5Q',nobtyp,rxobtyp,iobtypx)
        endif

c         Higher-frequency pseudorange       
        call fill_obtypx(3,'C1 ',nobtyp,rxobtyp,iobtypx)                
        call fill_obtypx(3,'C1X',nobtyp,rxobtyp,iobtypx)
        call fill_obtypx(3,'C1C',nobtyp,rxobtyp,iobtypx)
        
c         Lower-frequency pseudorange
* MOD TAH 200511: Added L7 or L8 frequency choice (default L5)
        if( gnsslf.eq.'E6' ) then   
           call fill_obtypx(4,'C6A',nobtyp,rxobtyp,iobtypx)    
           call fill_obtypx(4,'C6B',nobtyp,rxobtyp,iobtypx)      
           call fill_obtypx(4,'C6C',nobtyp,rxobtyp,iobtypx)      
           call fill_obtypx(4,'C6X',nobtyp,rxobtyp,iobtypx)      
           call fill_obtypx(4,'C6Z',nobtyp,rxobtyp,iobtypx)  
        elseif( gnsslf.eq.'E7' ) then   
           call fill_obtypx(4,'C7 ',nobtyp,rxobtyp,iobtypx)    
           call fill_obtypx(4,'C7X',nobtyp,rxobtyp,iobtypx)      
           call fill_obtypx(4,'C7I',nobtyp,rxobtyp,iobtypx)      
           call fill_obtypx(4,'C7Q',nobtyp,rxobtyp,iobtypx)  
        elseif ( gnsslf.eq.'E8' ) then
           call fill_obtypx(4,'C8 ',nobtyp,rxobtyp,iobtypx)    
           call fill_obtypx(4,'C8X',nobtyp,rxobtyp,iobtypx)      
           call fill_obtypx(4,'C8I',nobtyp,rxobtyp,iobtypx)      
           call fill_obtypx(4,'C8Q',nobtyp,rxobtyp,iobtypx)  
        else    ! Default
           call fill_obtypx(4,'C5 ',nobtyp,rxobtyp,iobtypx)    
           call fill_obtypx(4,'C5X',nobtyp,rxobtyp,iobtypx)      
           call fill_obtypx(4,'C5I',nobtyp,rxobtyp,iobtypx)      
           call fill_obtypx(4,'C5Q',nobtyp,rxobtyp,iobtypx)  
        endif    

      elseif( gnss.eq.'I' ) then

c         Higher-frequency carrier phase 
        call fill_obtypx(1,'L9 ',nobtyp,rxobtyp,iobtypx)
        call fill_obtypx(1,'L9C',nobtyp,rxobtyp,iobtypx)
    
c         Lower-freqeuncy carrier phase
        call fill_obtypx(2,'L5 ',nobtyp,rxobtyp,iobtypx)
        call fill_obtypx(2,'L5C',nobtyp,rxobtyp,iobtypx)

c         Higher-frequency pseudorange       
        call fill_obtypx(3,'C9 ',nobtyp,rxobtyp,iobtypx)
        call fill_obtypx(3,'C9C',nobtyp,rxobtyp,iobtypx)
        
c         Lower-frequency pseudorange
        call fill_obtypx(4,'C5 ',nobtyp,rxobtyp,iobtypx)
        call fill_obtypx(4,'C5C',nobtyp,rxobtyp,iobtypx)      

      else
        write(message,'(a,a1,a)') 'GNSS ',gnss,'not recognized '
        call report_stat('FATAL',prog_name,'lib/sel_obtyp',' '
     .                  , message,0 )
      endif

      return
      end
                              

      Subroutine fill_obtypx(islot,obtyp,nobtyp,rxobtyp,iobtypx )

c     Insert the obervable into the appropriate x-file slot

      implicit none                                

      include '../includes/makex.h'

      integer*4 islot,nobtyp,iobtypx(4),i

      character*3 obtyp,rxobtyp(maxobt)

cd      print *,'FILL islot obtyp nobtyp rxobtyp iobtypx '
cd     .            , islot,obtyp,nobtyp,rxobtyp,iobtypx
      do i=1,nobtyp        
cd      print *,'i obtyp rxobtyp(i) ',i,obtyp,rxobtyp(i)                           
        if( obtyp.eq.rxobtyp(i) ) iobtypx(islot) = i 
      enddo
      return
       end 
                                  



      



                     

