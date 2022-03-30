       SUBROUTINE  ACMAKE( istat,dfile,cfile,autcln_cmd,delaut
     .                  , typana,avedit_opt,eclopt,pre_post
     .                  , autcln_clk,mfile,ierr_sestbl )

C     subroutine to make AUTCLN/SINCLN batch files
C
C     R. King    91/7/1/; last modified 13/03/28
C
C     calling arguments
C        istat      : number of stations                   input
c        dfile      : D-file name for site list            input
c        cfile      : array of C-file names                input/output  
c        autcln_cmd : base autcln command file             input
c        delaut     : 'Y' to always delete input  C-files  input
c                     'N' to keep all input C-files        
c                     'I' to keep only the n-1 C-file for postfit AUTCLN 
c        typana     : type of analysis (from sestbl)       input  
c        avedit_opt : delete bad-clock data ('YES'/'NO ')  input 
c        eclopt     : full string for sh_autecl            input
c        pre_post   : this call ('PRE '/'POST'/'CLK')      input
c        autcln_clk : there is to be an extra autcln       input 
c        mfile      : m-file name for postfit edits        input
c        ierr_sestbl: 0 if not sestbl. error               output
                     

      include '../includes/dimpar.h'
                                
      integer*4 ierr_sestbl,istat,iper,nblen,i

      character* 1  delaut,autcln_clk,inchr,outchr,lowerc
      character* 3  outcmd,avedit_opt    
      character* 4  pre_post
      character* 5  typana,postopt
      character*7   avopt
      character*16  cfile(maxsit),cflout,dfile,mfile
     .              ,mfile_aut          
      character*30  eclopt
      character*80  autcln_cmd,autcmd         

      ierr_sestbl = 0 
 
c-------------------------------------------------------------------
c     Sequence for AUTCLN

c       create autcln.cmd.[prefit/postfit] by invoking script sh_autedit to
c       read the autcln base command file and 
c         1) uncomment the 'post' lines
c         2) add eclipse edit_sv_site command (call sh_autecl)
c         3) add bad-clock edit_sv_site commands (call sh_avedit)
                             
	if( pre_post.ne.'CLK '  ) then
          if( pre_post.eq.'post'.or. pre_post.eq.'POST' ) then  
             autcmd = 'autcln.cmd.postfit'
             postopt = '-post'
          else
             autcmd = 'autcln.cmd.prefit'
             postopt = ' '               
          endif     
          if( avedit_opt.eq.'yes' .or. avedit_opt.eq.'YES' ) then
             avopt = '-avedit'
          else
             avopt = ' '
          endif
          write(17,'(a,1x,a,1x,a,1x,a,1x,a)' ) 
     .     'sh_autedit -base ',autcln_cmd(1:nblen(autcln_cmd))
     .                      , postopt, avopt, eclopt
        endif
c       input C-file series
        iper = index(cfile(1),'.')
        inchr = cfile(1) (iper-1:iper-1)

C       write a comment to screen that AUTCLN is running
c** no longer needed with status files?  RWK: restored 100716
        write(17,'(A)') 
     .     'echo AUTCLN is running--see autcln.out for messages'

c       output C-file series - use the delete control to decide whether
c       to overwrite or increment except for the autcln run prior to
c       a requested extra run for clocks, where  we need to save the
c       original (MODEL-written) c-files to be used in the clock run
        if( autcln_clk.eq.'Y'.and.pre_post.eq.'POST' ) then
	  outcmd = ' + '
	else
          if( delaut.eq.'Y' ) 
     .      then
c           AUTCLN will overwrite
            outcmd = ' . '
          elseif (delaut.eq.'N' ) then
c           AUTCLN will increment
            outcmd = ' + ' 
          elseif (delaut.eq.'I' ) then
c           AUTCLN will increment for for prefit editing (this leaves excess
c           C-files in multiple iteration situations, but it is designed 
c           mainly for 0-ITER, where the first set of C-files are saved)
            if( pre_post.eq.'PRE ' ) then
              outcmd = ' + '
            else
              outcmd = ' . '
            endif
	  endif
        endif 
    
c       write the AUTCLN execution line into the batch file
        write(17,'(6a)') 'autcln ',autcmd(1:nblen(autcmd)),outcmd,dfile
     .                 ,  inchr,' >! autcln.out'
        write(17,'(a)') 'if( -e GAMIT.fatal ) exit'

c       if autcln executed in overwrite mode, need to update C-files
        if( autcln_clk.ne.'Y' .and. 
     .    ( delaut.eq.'Y'.or. (delaut.eq.'I'.and.pre_post.eq.'POST'))) 
     .     then
           outchr = inchr
           call newchr(outchr)
           inchr=lowerc(inchr)
           outchr=lowerc(outchr)
           if(typana.eq.'PREFI') return
           
           write(17,'(4a)') 'mvcf ',inchr,' ',outchr
        endif
c       update the FIXDRV C-file names in either case
        do i=1,istat
          call upnam1( cfile(i),cflout )
          cfile(i) = cflout
        enddo    

c       if postfit, copy the M-file used by AUTCLN to get phase clocks so that it's  
c       not overwritten by CFMRG and SOLVE, and can be used to display one-ways in CVIEW
        if( pre_post.eq.'POST' ) then
          mfile_aut = mfile(1:nblen(mfile))//'.autcl'  
          write(17,'(a,1x,a,1x,a)') 
     .    "/bin/cp ",mfile(1:nblen(mfile)),mfile_aut(1:nblen(mfile_aut))
        endif

         
      RETURN
      END
