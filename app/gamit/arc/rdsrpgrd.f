       Subroutine rdsrpgrd (filename,xllcorner,yllcorner
     .   , cellsize,accgrid,nrows,ncols)
c      Program to read UCL SRP grid files
C      format:
c       nrows 100
c       xllcorner -91.817374366161616
c       yllcorner -181.78825066616162
c       cellsize 3.6347487323232324
c       nodata_value 1.7014100000000001E+038
c       4.2304456171657069E-009 -6.8459395072139169E-010 1.2942214728790698E-009
        implicit none
C Declare the variables
       integer*4 iusrp,nrows,ncols,ioerr,test,i,j
       real*8 xllcorner,yllcorner,cellsize,nodata_value
     .    ,accgrid(100,120)
       character*80 filename, row,col,celltxt,xllctxt
     .    ,yllctxt,nanvaltxt

C      open the file
       iusrp=99
       print*,"debug: arc/dsrpgrd filename:",filename
       open(unit=iusrp,file=filename,status='old',
     .    action='read',iostat=ioerr)
      if (ioerr .ne. 0) then
          call report_stat('FATAL','ARC','rdsrpgrd',filename,
     .    'Error opening grid file',ioerr) 
        endif  
C Set up the array
       do i=1,100
         do j=1,120
           accgrid(i,j)=0.d0
         enddo
       enddo

Cd       print*,'preread'
c read the nrows and ncols
       read(iusrp,'(a5,1x,i3)',iostat=ioerr) col  ,ncols  
      if (ioerr .ne. 0) then
          call report_stat('FATAL','ARC','rdsrpgrd',filename,
     .    'Error reading grid file',ioerr) 
        endif  

       read(iusrp,'(a5,1x,i3)') row,nrows
       read(iusrp,'(a9,1x,f20.8)') xllctxt,xllcorner
       read(iusrp,'(a9,1x,f20.8)') yllctxt,yllcorner
       read(iusrp,'(a8,1x,f18.8)') celltxt,cellsize
       read(iusrp,'(a12,1x,E23.15)') nanvaltxt,nodata_value
       do i=1,nrows
c          read(99,'(50(E23.15,1x))') (accgrid(i,j)
c     .       ,j=1,ncols)
           read(iusrp,*) (accgrid(i,j),j=1,ncols)
c           print*,i
       enddo

Cd       print*, row,nrows,ioerr
Cd       print*, col,ncols
Cd       print*, xllctxt,xllcorner,yllctxt,yllcorner
Cd       print*,celltxt,cellsize,nodata_value,nanvaltxt
Cd       print*,accgrid(1,2)
Cd       print*,accgrid(100,51)
Cd       print*,accgrid(101,52)
       close(unit=iusrp)
       end
