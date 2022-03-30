      Subroutine block_name(iblk,ablk)

      implicit none

      integer*4 iblk,rcpar,len
      character*5 ablk
      character*80 prog_name,message
                   
      if ( iblk.eq.1 ) then
        ablk = 'I    ' 
      elseif ( iblk.eq.2 ) then
        ablk = 'II   '
      elseif ( iblk.eq.3 ) then
        ablk = 'IIA  ' 
      elseif ( iblk.eq.4 ) then
        ablk = 'IIR-A'   
      elseif ( iblk.eq.5 ) then
        ablk = 'IIR-B'
      elseif ( iblk.eq.6 ) then
        ablk = 'IIR-M'
      elseif ( iblk.eq.7 ) then
        ablk = 'IIF  '
      else
       len = rcpar(0,prog_name) 
       write(message,'(a,i3,a)') 'SV Block ',iblk,' not recognized'
       call report_stat('FATAL',prog_name,'lib/block_names',' '
     .                 ,message,0)
      endif 
                                    
      return
      end

