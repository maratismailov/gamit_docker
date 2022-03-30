      program fixhr

      implicit none

*     Program to take SOPAC format station.info and reconcatenate lines that have
*     spit on day boundaries
*     Also replaces any 1900 001 entries with 1990 001
*     Run with:
*     fixhr < station.info > station.out
*

      character*256 line1, line2, lineo
      integer*4 ierr, cnt
      logical out

      read(*,'(a)',iostat=ierr) line1
      cnt = 0
      do while ( ierr.eq.0 )
         read(*,'(a)',iostat=ierr) line2
         out = .true.
         if( ierr.eq.0 ) then 
            if( line2(1:1).eq.' ' ) then
*               Replace 1900 lines
                call sub_char( line2,'1900 001','1990 001')

*               Fix problem with "DHTCR" "DHARP" at sites with  TRM41249.00
*               antennas
                if( index(line2,"DHTCR").gt.0 ) then
                    if( index(line2,"TRM41249.00").gt.0 ) then
                        call sub_char( line2,'DHTCR','DHARP')
                    endif
                endif
*               Fix GeoNet site height type.
                if( index(line2,"DHBGP").gt.0 ) then
                    if( index(line2,"TRM29659.00").gt.0 .or.
     .                  index(line2,"TPSCR4").gt.0  ) then
                        call sub_char( line2,'DHBGP','DHARP')
                    endif
                endif

*               Compare to previous end time
                if( line1(45:52).eq.line2(45:52) .and.
     .              line1( 1: 5).eq.line2( 1: 5) ) then
                    lineo = line1(1:44) // line2(45:61) // line1(62:)
                    line1 = lineo
                    out = .false.
                    write(*,'(a)') trim(line1)
                    read(*,'(a)',iostat=ierr) line2
                endif 

*               Skip sites with bad entrys
                if( line2(2:5).eq.'COT1' .or. line2(2:5).eq.'MAD2' .or.
     .              line2(2:5).eq.'MAJU' .or. line2(2:5).eq.'TSKB' .or.
     .              line2(2:5).eq.'TSK2' .or. line2(2:5).eq.'VBCA' .or.
     .              line2(2:5).eq.'ELRO' .or. line2(2:5).eq.'GMSD' .or.
     .              line2(2:5).eq.'MAJU' .or. line2(2:5).eq.'MTJO' .or.
     .              line2(2:5).eq.'SOLO' .or. line2(2:5).eq.'TONG' .or.
     .              line2(2:5).eq.'OUSD' .or. line2(2:5).eq.'ULAB' .or.
     .              line2(2:5).eq.'NETP' .or. line2(2:5).eq.'TUBI' .or.
     .              line2(2:5).eq.'P655' .or. line2(2:5).eq.'P659' .or. 
     .              line2(2:5).eq.'P664' .or. line2(2:5).eq.'P665' .or.
     .              line2(2:5).eq.'P667'                          )then
                      line2(1:1) = '#'
                end if
            end if
            if( out .and. line1(1:1).ne.'#' ) write(*,'(a)') trim(line1)
            line1 = line2
            cnt = cnt + 1
            if( cnt.eq.1 ) write(*,120)
 120        format('* Updated with program fixhr to re-concatentate',
     .             ' split day entries')
         else
            if( line1(1:1).ne.'#' ) write(*,'(a)') trim(line1)
         end if
      end do

      end
             
!Example
!ACOR  A Coruna          2008 112  9 30  0  2008 340 12 30  0   3.0460  DHPAB   0.0000   0.0000  LEICA GRX1200PRO      5.62                   5.62  459187                LEIAT504   
!ACOR  A Coruna          2008 340 12 30  0  2010  34  0  0  0   3.0460  DHPAB   0.0000   0.0000  LEICA GRX1200PRO      6.02                   6.02  459187                LEIAT504   
!ACOR  A Coruna          2010  34  0  0  0  2010  34 12  0  0   3.0460  DHPAB   0.0000   0.0000  LEICA GRX1200PRO      7.53                   7.53  459187                LEIAT504   
!ACOR  A Coruna          2010  34 12  0  0  2010  43  0  0  0   3.0460  DHPAB   0.0000   0.0000  LEICA GRX1200PRO      7.53                   7.53  459187                LEIAT504   
!ACOR  A Coruna          2010  43  0  0  0  9999 999  0  0  0   3.0460  DHPAB   0.0000   0.0000  LEICA GRX1200PRO      7.80                   7.80  459187                LEIAT504   

!ACOR  A Coruna          2008 112  9 30  0  2008 340 12 30  0   3.0460  DHPAB   0.0000   0.0000  LEICA GRX1200PRO      5.62                   5.62  459187                LEIAT504   
!ACOR  A Coruna          2008 340 12 30  0  2010  34 12  0  0   3.0460  DHPAB   0.0000   0.0000  LEICA GRX1200PRO      6.02                   6.02  459187                LEIAT504   
!ACOR  A Coruna          2010  34 12  0  0  2010  43  0  0  0   3.0460  DHPAB   0.0000   0.0000  LEICA GRX1200PRO      7.53                   7.53  459187                LEIAT504   
!ACOR  A Coruna          2010  43  0  0  0  9999 999  0  0  0   3.0460  DHPAB   0.0000   0.0000  LEICA GRX1200PRO      7.80                   7.80  459187                LEIAT504   

!SG01  SG01 (EF13, Lamo  1900 001 00 00 00  2005 302 00 00 00   0.0000  DHPAB   0.0000   0.0000  TRIMBLE 4700          1.35                   1.35  0220217852            TRM33429.20+GP   NONE   0220199216          
