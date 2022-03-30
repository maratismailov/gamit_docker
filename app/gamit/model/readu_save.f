      Subroutine READU ( iuu,sitcod )

                                                         
       
      implicit none
                                 
      include '../includes/modkin.h' 

      logical eoh,eof,found_site,newblock,end_site
      integer*4 iuu,ioerr,i,j

      real*4 uvers
                       
      character*4 sitcod,upperc,table_type
      character*8 keyword,keywrd1
      character*80 message
      character*256 line
            
                          
       found_site = .false.
       eof = .false.
cc       do while (.not.found_site .and. .not.eof ) 
 
        read(iuu,'(a)',iostat=ioerr) line 
ccc        if( ioerr.eq.0 ) then    
          call report_stat('WARNING','MODEL','readu',' '
     .           ,'Error reading line of u-file',ioerr)  
ccc        elseif( line(1:1).ne.' ') then
          continue
ccc        elseif( line(2:8).eq.'ENDFILE' ) then
          eof = .true.
ccc        else 
          if( line(2:8).eq.'STATION') then 
            if( upperc(line(11:14)).eq.upperc(sitcod) ) then  
c             read blocks for site  
              found_site = .true.
              do while (.not.end_site) 
c               this assumption in this loop is that within each model block,
c               all of the available records will be read, so that the 
c               next non-comment line should be a new model, a new station,
c               or an ENDFILE. 
                read(iuu,'(a)',iostat=ioerr) line
                 if( ioerr.ne.0 ) then 
                    call report_stat('FATAL','MODEL','readu',' '
     .               ,'Error reading MODEL line of u-file',ioerr) 
                 elseif( line(1:1).ne.' ') then
                   continue      
                 elseif( line(2:8).eq.'STATION') then
                    end_site = .true.
                 else
                   keyword = line(2:8)     
                 endif
                 if( line(11:15).ne.'MODEL' )  
     .              call report_stat('FATAL','MODEL','readu',' '
     .                ,'Keyword in first block line not MODEL',0)
                 keyword = line(2:9)

                 if( keyword.eq.'ATMOSLOD') then
c                  read the model line and then an arbitrary number of epochs
                   read(line,'(16x,a8)',iostat=ioerr) atmlmod
                   if( ioerr.ne.0 )  call report_stat('FATAL','MODEL'
     .                ,'readu',' '
     .               ,'Error reading MODEL line for ATMOSLOD in u-file'
     .                 ,ioerr)
                   natml = 0 
                   do while( .not.newblock )   
                     read (iuu,'(a)',iostat=ioerr) line
                     if( ioerr.ne.0 ) then  
                        call report_stat('FATAL','MODEL','readu',' '
     .                ,'Error reading ATMOSLOD MODEL from u-file',ioerr) 
                     elseif( line(2:8).ne.keyword ) then
                       newblock = .true.
                     elseif( line(2:8).eq.'ENDFILE') then
                       newblock = .true.
                       eof = .true.  
                     else
                       natml = natml + 1
                       read(iuu,'(1x,a8,1x,f6.0,3f15.0)',iostat=ioerr ) 
     .                      atml_time(natml),(atml_val(j,natml),j=1,3)
                        if( ioerr.ne.0 ) call report_stat('FATAL'
     .                    ,'MODEL','readu',' '
     .                    ,'Error reading ATMOSLOD values from u-file'
     .                    ,ioerr)
                     endif
                   enddo
                                  

                 else
                   continue 
c                  endif for selecting block (MODEL line) 
                 endif 
                    
c              enddo on site block
               enddo
            
              else 
                continue
c               endif for search on site code
              endif

            else
              continue
c             endif for 'STATION'
            endif

ccc          else    
            continue
c           endif for reading next line
ccc          endif
            
c       enddo on search over all sites
cc        enddo
     

      if( .not.found_site ) then
         write(message,'(a,a4,a)') 'Site ',sitcod,' not found on u-file'
         call report_stat('WARNING','MODEL','readu',' ',message,0)
      endif
            
      return
      end

