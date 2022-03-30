      Subroutine tdtrit(iungs,iut,spver,jds,ts,jdf,tf,delt,nepchs,nintrs
     .                  ,gnss_sel,numsp3sv,numsat,itsat    
     .                 , pvsigb,clksigb,otlmod )

c       Read an NGS SP file and write the data records of a T-file
c       R. King  March 1988
c       P. Fang  March 1993 mod. to allow SP#3 format input 
c       R. King  Dec   2002 mod. to allow reading SP3-c format
c       R> King  Nov   2015 mods for GNSS

c       SV clocks, orbit accuracy and event flags now available but not yet used

c     Input

c       iungs    : unit number for sp3 file
c       iut      : unit number for t-file
c       spver    : SP3 version 
c       jds,ts   : start epoch (PEP JD, sec-of-day)
c       jdf,tf   : stop epoch 
c       delt     : t-file interval (sec)
c       nepchs   : number of epochs (same on SP3 and t-file)
c       nintrs   : number of components to write (3 or 6)
c       gnss_sel : GNSS code for SVs selected
c       numsp3sv : number of SVs on SP3 file 
c       numsat   : number of SVs for t-file (all matching gnss_sel)
c       itsat    : PRNSs for SVs on t-file
c       psvsigb  : base for pos/vel values
c       clksigb  : base for clock values
c       otlmod   : OTL model
 
c     Ouput variables: none
c       
      implicit none

      include '../includes/dimpar.h'

      integer*4 iungs,iut,julday,jd,jds,jdf,nintrs
     .        , numsp3sv,numsat,itsat,nepchs,iepchs
     .        , iyr,imon,iday,ihr,imin,id,sigx
     .        , i,j,k,len,rcpar,ioerr,notl

      real*8 delt,delts,tf,ts,t,x,y,sec,clock,pvsigb,clksigb
     .     , docmc(3)

      character*1 spver,eflg,cflg,mflg,pflg,gnss_sel,asvid
      character*8 otlmod
      character*80 prog_name,line
      character*256 message                          

      logical debug/.false./
    
      dimension itsat(maxsat),x(6),y(3,maxsat),sigx(6)
    
c     Get the program name calling tdtrit
      len = rcpar(0,prog_name)
c             

c     Code assumes that the SP3 and t-file headers have been
c     read and written, thus positioned to read and write
c     the first epoch

      do  iepchs=1,nepchs

c       read, but don't use, the time record of each group
        if(spver.eq.'1') then
           read(iungs,'(4x,i4,4(1x,i2),1x,f10.7)',iostat=ioerr)
     .      iyr,imon,iday,ihr,imin,sec
        else
          if( debug ) then
            read(iungs,'(a)') line
            print *,'line of epoch ',line  
             backspace iungs
          endif
          read(iungs,'(3x,i4,4(1x,i2),1x,f10.7)',iostat=ioerr) 
     .           iyr,imon,iday,ihr,imin,sec
        endif  
        if( ioerr.ne.0 ) then 
          write(message,'(a,f5.2,a)') 'Error reading epoch ',iepchs
     .        ,' for SP ',spver,'file'  
          call report_stat('FATAL',prog_name,'orbits/tdtrit',' '
     .           ,message,ioerr)
        endif
        
        do j=1,numsp3sv
          if (spver.eq.'1') then
            read(iungs,'(3x,i2,3(1x,f12.5),3(1x,f11.8))',iostat=ioerr)
     .             id,(x(i),i=1,6)
          else
            read(iungs,'(a)',iostat=ioerr) line   
c           assume that sp3-c has the full 80 columns, though some can be blank
            read(line,'(1x,a1,i2,4f14.6,3i3,i4,1x,2a1,2x,2a1)'
     .          ,iostat=ioerr) asvid,id,(x(i),i=1,3),clock
     .            , (sigx(i),i=1,4),eflg,cflg,mflg,pflg     
            if( asvid.eq.' ' ) asvid = 'G'
            if(debug) print *,'read epoch asvid id x '
     .         ,iepchs,asvid,id,(x(i),i=1,3)
          endif
          if( ioerr.ne.0 ) then
            write(message,'(a,i5,a,f5.2,a)') 
     .         'Error reading coordinates at epoch ',iepchs
     .            ,'for SP ',spver,' file'
            call report_stat('FATAL',prog_name,'orbits/tdtrit',' '
     .                      ,message,ioerr)
          endif        
          if(debug) print *,'Looping over numsat itsat ',numsat
          do  k=1,numsat
            if(debug)  print *,' k asvid id  itsat)k)  '
     .          , k,asvid,id,itsat(k)
            if( asvid.eq.gnss_sel.and.id.eq.itsat(k) ) then
               if(debug) print *, 'Set y(k) = x'
               do i=1,3
                 y(i,k)= x(i) 
               enddo
             endif    
           enddo
c         end on SVs
          enddo      
          
c** RWK 160602: I had changed to above code to the following, which doesn't work for
c               cases of omitted SVs. I don't remember why I made the change, but
c               undo it for now.
c          if( asvid.eq.gnss_sel ) then
c            isatcnt = isatcnt + 1
c            do i=1,3
c              y(i,isatcnt)= x(i) 
c            enddo
c          endif        
c        enddo
c        if( isatcnt.ne.numsat ) then
c          write(message,'(2(a,i3))') 'SV count for t-file',isatcnt
c     .                 ,' not equal numsat ',numsat
c          call report_stat('FATAL',prog_name,'orbits/tdtrit',' '
c     .                      ,message,ioerr)
c        endif
 
c       apply the correction from CE to CM for ocean loading 
        if( otlmod(1:1).ne.' ' ) then   
          jd= julday( imon,iday,iyr )
          t= dble(ihr)*3600.d0 + dble(imin)*60.d0 + sec
c         hard-wire the # of components
          notl = 11
          call otlcmc( jd,t,otlmod,notl,2,docmc ) 
cd           write(6,'(a,i8,f8.1,3f7.4)') 'CMC db ',jd,t
cd            ,(docmc(i),i=1,3)
          do j=1,numsat
            do i=1,3
              y(i,j) = y(i,j)  + docmc(i)/1.d3
            enddo
          enddo
        endif
        write(iut) ((y(i,j),i=1,nintrs),j=1,numsat)
            
c     end loop on epochs
      enddo

      write(message,100) iepchs
  100 format('Successfully wrote Earth-fixed T-file ',i5
     .       ,' epochs written on T-file')
      call report_stat('STATUS',prog_name,'orbits/tdtrit',' ',message,0)

      return
      end
