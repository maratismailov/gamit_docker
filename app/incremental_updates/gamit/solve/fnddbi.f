c
      subroutine fnddbi(free,ia1,last_bias)
c
c     Find the dependent biases by detecting when the bias-part of the 
c     normal matrix goes singular; then remove the row and column
c     Written by D. Dong  Mar 1991
c     Modifications to catch problems  rwk Oct 93, Mar 96 
c     Minor modifications  rwk/scm Jan 97
c     Commented and rewritten structurally (no change in algorithms) in an attempt 
c      to avoid what appears to be a compiler bug with HP fort77 with optlevel 3.   rwk Oct 97 
c       ** failed
c     Float all calculations of the form k(k-1)/2 to avoid bug.  Works.  tah 12 Nov 97
c
      implicit none
      include '../includes/dimpar.h'
      include 'solve.h'
c
      logical ifirst
      integer*4 free,ia1,last_bias
     .        , ii,jj,kk,i2,i20,in22,nomal,ij2,ijk,ijkk,ij20
     .        , krow,kcol,ier,kbad,keig,jelf
c      these added for bias apr debug
       integer*4 i

c      real*8 i2_mid

      real*8 eig,rcond,rcond_saved,rcond_ratio
      character*256 message
                       
      logical debug/.false./

      dimension eig(maxobs),keig(maxobs),free(maxprm)
                
c      call report_stat('STATUS','SOLVE','fnddbi',' ','Called fnddbi',0)


c     initialization to avoid compiler warning
      nomal = 0  
      ij20 = 0
c     pointer to elements removed
      do ii=1,maxobs
        keig(ii) = 0
      enddo
      
c     get the size of the bias matrix
      in22=(last_bias-ia1)/iband   
                                 
c     copy the A matrix into the AN22 matrix`     
      if(debug) then 
        print *,'FNDDBI last_bias nlive ',last_bias,nlive
        print *,'FNDDBI maxwm1 in22 ,ia1 ',maxwm1,in22,ia1
      endif        
      call fillan22( in22,ia1 )       
      if(debug) print *,'FNDDBI a(139655) an22(2) ',a(139655),an22(2)
                   
c      Write the condition number criterion into the q-file
       write(message,'(a,f9.1)')
     .   'Condition-number ratio for removing dependent biases is '
     .   ,bias_rcond  
       write(10,'(a)') message
       if( logprt ) write(6,'(a)') message

c     invert it to see if it's singular (ier=130; dependent biases exist)
c      call invers(an22,eig,4,in22,ier)  
      call inver2(an22,eig,4,in22,rcond,ier) 
      if(debug) then 
         print *,'FNDDBI 1st inversion ier ',ier
         print *,(an22(i),i=1,5)     
      endif 

                              
c*** make this nonsensical to force checking biases and getting debug
      if( ier.eq.999 ) then
c**      if( ier.ne.130 ) then
c        no dependent biases--exit
         return

      else

c        sequentially remove rows and columns until there are no dependent biases

c        first refill the an22 matrix since it was overwritten in the inversion
         call fillan22( in22,ia1 )
       
c        initialize the number of dependent biases counter
         kbad=0
         ifirst = .true.

c        outer loop over all elements of the matrix
         do 90 ijk=1,in22

c          initialize the row counter for each pass through
           krow=0 

c          loop over all elements up to the current limit, skipping any rows and
c          columns corresponding to dependent biases already found 

c          loop over rows      
           krow=0
           do 80 ii=1,ijk
            if (kbad.gt.0) then
              do  kk=1,kbad
                if (ii.eq.keig(kk)) goto 80
              enddo
            endif
            krow=krow+1 
c           replace by function call
c            i20=nint(float(krow)*float((krow-1))/2)
            i20 = jelf(krow)
c*           Fix Bug
c            i2_mid = float(ii)*float((ii-1))
c            i2=nint(i2_mid/2)   
            i2 = jelf(ii)
                     
c           loop over columns
            kcol=0
            do 60 jj=1,ii
               if (kbad.gt.0) then
               do kk=1,kbad
                  if (jj.eq.keig(kk)) goto 60
               enddo
               endif
               kcol=kcol+1
               ij2=i2+jj
               ij20=i20+kcol
               bn22(ij20)=an22(ij2)
 60         continue

 80      continue
        
c        now invert the compressed matrix to see if it's singular         
         ij2=ijk-kbad   

         if(debug) then 
           print *,'FNNDBI bef inver2 ij2 bn22 ' ,ij2
           print *,(bn22(i),i=1,5)
         endif 
c         call invers(bn22,eig,4,ij2,ier)  
         call inver2(bn22,eig,4,ij2,rcond,ier)  
         if(debug) then 
           print *,'FNDDBI 2nd inversion ier bn22 ',ier
           print *,(bn22(i),i=1,5)     
          endif 



c If matrix is truely singular skip condition number ratio calc since RCOND must be 0.
         if ( ier .ne. 130 ) then
c Don't throw out the first bias.
           if (ifirst) then
              rcond_ratio = 1.0D0
              ifirst = .false.
           else
              rcond_ratio = rcond_saved/rcond
           endif
c
c Check the condition number ratio of this inversion VS. the last inversion.
c Throw out this bias if the ratio is below the set tolorance bias_rcond     
           if ( rcond_ratio .gt. bias_rcond ) then
              write(message,85) rcond, rcond_ratio
 85           format('Bias matrix ill conditioned - bias removed with ',
     .               'rcond: ',e12.6,' ratio: ',e12.3)
              call report_stat('STATUS','SOLVE','FNDDBI',' ',message,0)
              ier = 140
           else
c Save this condition number for comparison with the next inversion estimate.
              rcond_saved = rcond
           endif
         else
           rcond_ratio = 99999999.999999
         endif
c
c Debug         
      if( bias_debug ) then
       write(6,87) ier,ij20,bn22(ij20),rcond,rcond_saved,rcond_ratio
       write(10,87) ier,ij20,bn22(ij20),rcond,rcond_saved,rcond_ratio
87     format('ier, ij20, bn22(ij20), rcond, rcond_saved, rcond_ratio: '
     . ,i4,i8,f14.5,f16.9,f16.9,f10.2)     
      endif

c         print *,'IER: ',ier,dabs(bn22(ij20))
c         print *,'inverted eig ier',eig,ier
c         print *,(bn22(i),i=1,5)   
c
c Use the bad inversion flag and also a small-number limit to test for
c bias matrix inversion problem.
         if( dabs(bn22(ij20)).lt.1.0d-6 .or. ier.ge.130  ) then

cs            print *,'FNDDBI:ikj ij2 ij20 bn22',ijk,ij2,ij20,bn22(ij20)
            kbad=kbad+1
            keig(kbad)=ijk
            ij2=ia1+ijk
            free(ij2)=0
            ijkk=ij2-lpart
            idxb(ijkk)=-ij2
            nlive=nlive-1
            if(logprt) write(6,110) ij2,dabs(bn22(ij20)),rcond_ratio,ier
            write(10,110) ij2,dabs(bn22(ij20)),rcond_ratio,ier
            if (l2flag.ge.2) then
               free(ij2+in22)=0
               idxb(ijkk+in22)=-ij2-in22
               nlive=nlive-1
               if( logprt ) write(6,110) ij2+in22
               write(10,110) ij2+in22
            endif    
            if(debug) 
     .        print *,'FNDDBI ij20 bn22 ratio ier nlive ij2 in22 free '
     .      ,  ij20,bn22(ij20),rcond_ratio,nlive,ij2,in22,free(ij2+in22)
         else    
c ** uncomment these for debugging 
c            write(6,120) ia1+ijk, dabs(bn22(ij20)), ier 
c            write(10,120) ia1+ijk, dabs(bn22(ij20)), ier
         endif
 90   continue  
           
c     endif on initial existence of dependent biases
      endif      

 110  format(1x,' Fix dependent bias param. of index ',i4,
     .       ' Tol ',e14.3,' rcond_ratio ',e12.3,' ier ',I4)
c 120  format(1x,' OK  dependent bias param. of index ',i4,
c     .       ' Tol ',e14.3,' ier ',I4)
      
      return
      end

       
      Subroutine fillan22 ( in22,ia1 )

      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'

      integer*4 in22,ii,jj,ia1,ia2,i2,i4,ij2,ij4,jelf
      
      logical debug/.false./
                               
c      print *,'a'
c      print *,(a(ii),ii=1,100) 

      if(debug) print *,'FILLAN22 in22 ia1 A ',in22,ia1,(a(ii),ii=1,5)
      do ii=1,in22
         ia2=ia1+ii  
c        replace by call to function   tah/rwk 980121
c         i4=nint((float(ia2)*float((ia2-1))/2))
c         i2=nint((float(ii)*float((ii-1))/2))  
         i4 = jelf(ia2)
         i2 = jelf(ii)  
         if(debug)  print *,'  ia1 ii ia2 i4 i2 ',ia1,ii,ia2,i4,i2
         do jj=1,ii
            ij4=i4+ia1+jj
            ij2=i2+jj
            an22(ij2)=a(ij4) 
            if(debug)  print *,'  ij4 ij2 ',ij4,ij2
         enddo
       enddo    
c       print *,'a* ',a(22155),a(22365)
c       print *,'an22 ',an22(1),an22(2)           
       if(debug) print *,'FNDDBI a(22365) an22(2) ',a(22365),an22(2)
       return
  
       return
       end

