      Subroutine read_sp3(lu)

*      Read an SP3-C file into storage

      implicit none

      include 'orbitx.h'

      integer*4 lu
                        
*   Values from SP3 header (local)
      character*1 spver,pvflag,gdum 
      character*2 linesym,ftype
      character*3 orbtyp,time_sys,clkmod,orbmod
      integer*4 igpswk,mjd,idumy(85),iaccsv(maxsv),i
      real*8    fmjd,pvsigb,clksigb

        
      logical debug/.true./
     
 read(lun,'(a2,1x,a2,1x,a3)',iostat=ioerr) linesym,ftype,time_sys 
c    
c Determine the file format using the first two characters of the first line

c   SP1    ' # '
c   SP3-a  '# '  or '# a'
c   SP3-b  '#b'
c   SP3-c  '#c'

c  For now, support only SP3-c

c  Read the first and second records of the SP3 file header
       
      read(lun,'(a2,a1,i4,4i3,1x,f11.8,3x,i5,7x,a5,1x,a3,1x,a4)'
     .    ,iostat=ioerr) linesym,pvflag
     .    ,fstr_year,fstr_mon,fstr_day,fstr_hr,fstr_min,fstr_sec
     .    ,nfepoch,crdsys,orbtyp,agency   
      if( ioerr.ne.0 ) call report_stat('FATAL','ORBITS','read_sp3'
     .                 ,' ','Error reading 1st line of sp3 file',ioerr)  
      if( linesym.ne.'# '.and.linesym.ne.'#a'.and.linesym.ne.'#b'.and.
     .    linesym.ne.'#c' ) then
        call report_stat('FATAL','ORBITS','read_sp3'
     .        ,' ','Invalid version type on line 1 of sp3 file',0)
      else
        spver = linesym(2:2)
      endif
      call check_y2k(iyear)  
      read(lun,'(3x,i4,17x,f14.8,1x,i5,1x,f15.13)',iostat=ioerr)
     .            igpswk,finterval,mjd,fmjd    
      if( ioerr.ne.0 ) call report_stat('FATAL','ORBITS','read_sp3'
     .                 ,' ','Error reading 2d line of sp3 file',ioerr)  

c  Read the PRN numbers and their accuracy codes
            
c     lines 3-7: PRNSs 
      read(lun,'(a60)',iostat=ioerr) line
      if( ioerr.ne.0 ) call report_stat('FATAL','ORBITS','read_sp3'
     .             ,' ','Error reading line 3 of SP3 file',ioerr)  
      read(line,'(4x,i2)') numsv    
      if( ioerr.ne.0 ) call report_stat('FATAL','ORBITS','read_sp3'
     .             ,' ','Error decoding line 3 of SP3 file',ioerr)  
      if( numsv.gt.maxsorbsat ) then
          write(message,'(a,i2,a,i2,a)') 
     .      'Number of satellites on SP3 file (',numsv,') > maxsv ('
     .        ,maxsv,')'
        call report_stat('FATAL','NGSTOT','orbits/read_sp3',' ',message,0)
      endif 
      backspace(lun)
      read(lun,'(4x,i2,3x,17(a1,i2),4(/,9x,17(a1,i2)))',iostat=ioerr)
     .     nsv,(gdum,isv(i),i=1,nsv),(gdum,idumy(i),i=nsv+1,85)
       if( ioerr.ne.0 ) call report_stat('FATAL','ORBITS','read_sp3'
     .                ,' ','Error reading lines 3-7 of sp3 file',ioerr)
c     lines 8-12: accuracy codes  
      read(lun,'(9x,17i3,4(/,9x,17i3))',iostat=ioerr)
c      read(lun,'(9x,17i3)',iostat=ioerr)   
     .    (iaccsat(i),i=1,numsv),(idumy(i),i=numsv+1,85)  
      if( ioerr.ne.0 ) call report_stat('FATAL','ORBITS','read_sp3'
     .                ,' ','Error reading lines 8-12 of sp3 file',ioerr)
c     convert the accuracy codes from 2**n mm to floating point m
      do i=1,nsv
        sigp(i) = 1.d-3*2.d0**iaccsv(i)
      enddo  

c Read the file and time types and exponent base for time-dependent position and velocity sigmas    
                                                    
c     line 13
      read(lun,'(a2,1x,a2,1x,a3)',iostat=ioerr) linesym,ftype,time_sys 
c     these will be dummy but readable for sp3-a
      if( ioerr.ne.0 ) call report_stat('FATAL','ORBITS','read_sp3'
     .                ,' ','Error reading line 13 of sp3 file',ioerr) 
c    line 14
      read(lun,'(1x)',iostat=ioerr) 
c    line 15
      read(lun,'(3x,f10.7,1x,f12.9)',iostat=ioerr) pvsigb,clksigb
      if( ioerr.ne.0 ) call report_stat('FATAL','ORBITS','read_sp3'
     .             ,' ','Error reading line 14 or 15 of sp3 file',ioerr)    

c Rewind the file and look for models in the comment line

      rewind(lun,iostat=ioerr)
      if( ioerr.ne.0 ) call report_stat('FATAL','ORBITS','read_sp3'
     .     ,' ','Error rewinding SP3 file after reading header',ioerr)
      eoh = .false.
      pcvmod = ' '
      otlmod = ' '
      atlmod = ' '
      clkmod = ' '
      orbmod = ' '
      do while (.not.eoh ) 
        read(lun,'(a)',iostat=ioerr) line   
        if( ioerr.ne.0 ) then
          call report_stat('FATAL','ORBITS','read_sp3'
     .                ,' ','Error reading header for models line',ioerr)
        elseif( line(1:1).eq.'*') then 
c         end of comments
          eoh = .true.
        elseif( line(14:15).eq.'TL' ) then 
c         read comment line for models   
c         this group for temporary (few weeks, Nov 2006) format
          pcvmod = line(8:12) 
          otlmod = line(17:24) 
          otlflg = line(28:28)
          atlmod = line(30:37)
          atlflg = line(41:41) 
          clkmod = line(47:49)
          orbmod = line(55:57)
c       elseif( line(18:19).eq.'OL') then
c         this group for revised format
          pcvmod = line(8:17) 
          otlmod = line(25:32) 
          otlflg = line(43:43)
          atlmod = line(34:41)
          atlflg = line(44:44) 
          clkmod = line(58:60)
          orbmod = line(50:52)
        endif                                                  
      enddo

c  File should now be positioned at the first epoch; backspace one
  
      backspace(lun,iostat=ioerr) 
      if( ioerr.ne.0 )  call report_stat('FATAL','ORBITS','read_sp3'
     .   ,' ','Error backspacing after reading comments',ioerr)

c  Echo the SP-3 file headers

      if( debug ) then 
        write(6,30) spver,orbtyp,crdsys,org,igpswk
     1          , iyear,imon,iday,ihr,imin,sec,mjd,fmjd,nepoch,delt
     2          , numsv,(issat(i),i=1,numsv)
* MOD TAH 200618: Updated 32I to 50I to allow for 35 Beidou satellites
   30 format(//,1x,'Header records from NGS Standard Product file:'
     1      , // 
     a      ,1x,'SP version: ',a1,//
     2      ,1x,'Orbit Type: ',a3,'   Coordinate System : ',a5
     3      , '  Organization : ',a4,'   GPS Week : ',i4,/
     4      ,1x,'Start epoch (GPST): ',i4,4i3,1x,f10.7,/
     5      ,1x,'            MJD   : ',i5,1x,f15.14,/
     6      ,1x,'Number epochs     : ',i6,'  Interval :',f7.2,' sec',//
     7      ,1x,i2,' satellites, PRN #s = ',50i3,/
     8      ,5x,     'accuracy (2**n mm)= ',50i3,/) 
        write(6,31) (iaccsat(i),i=1,numsv)  
   31   format(5x,'accuracy (2**n mm)= ',50i3,/)    
        if( pcvmod(1:1).eq.' ') then
          write(6,'(a)') 'No comment line for models'
        else
        write(6,32) pcvmod,otlmod,otlflg,atlmod,atlflg,clkmod,orbmod
   32     format(1x,'PCV model : ',a5,/
     .          ,1x,'OTL model : ',a8,'  applied = ',a1,/
     .          ,1x,'ATL model : ',a8,'  applied = ',a1,/
     .          ,1x,'CLK corrections : ',a3,/
     .          ,1x,'ORB corrections : ',a3)
        endif
      endif   

*   Now read all the values into storage

      do i = 1,nepoch

c       read, but don't use, the time record of each group
        read(iungs,12,end=110) iyr,imon,iday,ihr,imin,sec

        do j=1,numsv

          read(iungs,'(a)',iostat=ioerr) line   
c         assume that sp3-c has the full 80 columns, though some can be blank
          if( spver.eq.'c') then
            read(line,'(2x,i2,4f14.6,3i3,i4,1x,2a1,2x,2a1)'
     .            ,iostat=ioerr,end=130) id,(x(i),i=1,3),clock
     .            , (sigx(i),i=1,4),eflg,cflg,mflg,pflg
          else
            read(line,'(2x,i2,3(1x,f13.6))',iostat=ioerr,end=130)
     .                 id,(x(i),i=1,3) 
          endif
          if( ioerr.ne.0 ) 
     .      call report_stat('FATAL',prog_name,'orbits/read_sp3',' '
     .         ,'Error decoding coordinates line SP3 file',ioerr)
          endif   

c*          write(6,'(1x,a,2i3,3f11.3)') 'J,ID,X=',j,id,(x(i),i=1,3)
            do  k=1,nsat
                if( itsat(k).eq.issat(j) ) then
                  do i=1,3
                    y(i,k)= x(i) 
                  enddo
                endif    
              enddo
c           end loop on SVs
            enddo
                            
c           apply the correction from CE to CM for ocean loading 
            if( otlmod(1:1).ne.' ' ) then   
              jd= julday( imon,iday,iyr )
              t= dble(ihr)*3600.d0 + dble(imin)*60.d0 + sec
c             hard-wire the # of components
              notl = 11
              call otlcmc( jd,t,otlmod,notl,2,docmc ) 
c*              write(6,'(a,i8,f8.1,3f7.4)') 'CMC db ',jd,t
c*     .           ,(docmc(i),i=1,3)
              do j=1,nsat
                do i=1,3
                  y(i,j) = y(i,j)  + docmc(i)/1.d3
                enddo
              enddo
            endif
            write(iut) ((y(i,j),i=1,nintrs),j=1,nsat)
            iepcht= iepcht + 1
            if( iepcht.eq.nepcht ) goto 100
            
c        end loop on epochs
         enddo

      else

         do  i=1,numsv
           read(iungs,'(1x)',end=120)
         enddo
         goto 10

      endif

  100 write(message,101) iepcht
  101 format('Successfully wrote earth-fixed T-file ',i5
     1       ,' epochs written on T-file')
      call report_stat('STATUS',prog_name,'orbits/tdtrit',' ',message,0)

      return

  110 write(message,111) iyr,imon,iday,ihr,imin,sec
  111 format('End of file while looking for starting time on NGS file '
     1       ,'Last epoch read: ',i4,4(1x,i2),1x,f10.7 )
      call report_stat('FATAL',prog_name,'orbits/tdtrit',' ',message,0)

  120 write(message,121) jds,ts,jd,t
  121 format('Unexpected end of file in data records,'
     1      ,' Last computed: JDS,TS,JD,T ',2(i7,f8.5) )
      call report_stat('FATAL',prog_name,'orbits/tdtrit',' ',message,0)

  130 write(message,131) iepchs,iepcht
  131 format('End of NGS SP file  IEPOCH=',i4,' T-file IEPOCH=',I4)
      call report_stat('FATAL',prog_name,'orbits/tdtrit',' ',message,0)

      end





      return
      end



























