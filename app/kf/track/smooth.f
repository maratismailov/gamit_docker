CTITLE WR_COV_PARM

      subroutine wr_cov_parm(ep, option, data_type)

      implicit none

*     Routine to write/read the covariance matrix, solution vector,
*     and information about the ambiquity numbers into scratch
*     array and write it to disk or read from disk.

      include 'track_com.h'

* PASSED VARIABLES
* ep     -- Epoch being processed

      integer*4 ep

* option -- Set to W (write) or R (read).
* data_type -- Type of data being used. (L1,L2, or LC)

      character*(*) option, data_type

* LOCAL VARIABLES
* i, j, k  -- Loop variables
* ne  -- Number of bias parameters types (L1, L2, LC)
* np, nj -- Parameter numbers
* an  -- Ambiquity number
* npdest(max_parm) -- Destination parameter number of each of the
*        entries in the read buffer
* jel -- Element for computing position in lower-diagonally stored
*        matrix
* ierr -- IOSTAT error
* used_parm -- Number of used parameters at specific epoch.
 
      integer*4 i, j, k, ne, np, nj, an, npdest(max_parm), jel,
     .          ierr, used_parm

* idabuf(max_dabuf) -- Integer buffer for write cov_parm, sol_vec
*     and bias information to disk
* rdabuf(max_dabuf/2) -- Real*8 buffer for above information.

      integer*4 idabuf(max_dabuf)
      real*8    rdabuf(max_dabuf/2)

      equivalence ( idabuf, rdabuf )


***** See if write or read
      if( option(1:1).eq.'W' ) then
*         Clear the buffer before we start (Note: darecl is
*         in I*4 words).
          do i = 1, darecl
             idabuf(i) = 0
          end do

*         OK, Save the covariance, solution and bias numbers
          k = 0
          do i = 1, num_parm
             do j = 1, i
                k = k + 1
                rdabuf(k) = cov_parm(i,j)
             end do
          end do

*         Now do solution vector.  Here what we have depends
*         on whether or not the sites are considered static
*         or kinematic.  In the kinematic case, the position is
*         updated and this updated position will be used in the
*         back solution so here we set the sol_vec entry to
*         zero.  
*         Compute the start of the sol_vec entries based on total
*         number of parameters at any time because number of parameters
*         at any epoch in the forward and backward solutions may 
*         differ.
          k = (tac_parm+1)*tac_parm/2

          do i = 1, num_kine
             do j = 1,3
                k = k + 1
!               Zero sol_vec entry if curr_xyz has been updated.
!               debug_start = -3 stops this happening,
                if( .not.static(i) .and. debug_start.ne.-3 ) then
                    rdabuf(k) = 0.d0
                else
                    rdabuf(k) = sol_vec((i-1)*3+j)
                end if
             end do
          end do
*         Now save the rest of the solution vector
          do i = num_kine*3+1,num_parm
             k = k + 1
             rdabuf(k) = sol_vec(i)
          end do

*         Now save the ambiquity numbers.  These we save in the
*         integer array
          k = (tac_parm+1)*tac_parm/2 + tac_parm
          do i = 1, num_ambs
             if( amb_parn(1,i).ne.0 ) then
*                OK, this ambiquity is estimated and its parameter
*                number is given by amb_parn(1,i) for L1, (2,i) for
*                L2 (if used).
                 idabuf(2*k+abs(amb_parn(1,i))) = i
                 if( amb_parn(2,i).ne.0 ) 
     .               idabuf(2*k+abs(amb_parn(2,i))) = i
             endif
          end do
c     k = (tac_parm+1)*tac_parm/2 + tac_parm
c     write(*,160) ep, (idabuf(i),i=2*k+1,2*k+num_parm)
c160  format('AMBS ',I6,' A ', 100I5)

*         OK: Now write the buffer the to disk
          write(119, iostat=ierr, rec=ep) (idabuf(i),i=1,darecl)
          call report_error('IOSTAT',ierr,'writ',
     .         'Direct Access Smooth File',1,'wr_cov_parm')

          if( ep.ge.debug_start .and. ep.le.debug_end ) then
             call dump_cov(ep, rdabuf, idabuf, num_parm,
     .                     npdest,'N')
          end if

      else

****      OK: We need to read the direct access file.
          read(119,iostat=ierr, rec=ep) (idabuf(i),i=1,darecl)
          call report_error('IOSTAT',ierr,'read',
     .         'Direct Access Smooth File',1,'wr_cov_parm')

****      Initialize the saved covariance matrix and solution
*         vector.  This is to ensure that any parameters that
*         are not used will be correctly initialized.  Do the
*         whole matrix to make sure that the correlations 
*         between the non-ambiquity parameters and non-used
*         ambiquities is zero.
          do i = 1, tac_parm
             sol_vec_sav(i) = 0.d0
             do j = 1,i-1
                cov_parm_sav(i,j) = 0.d0
                cov_parm_sav(j,i) = 0.d0
             end do
             cov_parm_sav(i,i) = 10.d0
          end do

*         Move the buffer into place. Start with the covariance
*         matrix.   This gets tricky because we want to move the
*         entries for the ambiquity parameters into to their
*         correct positions.
          k = 0
          do i = 1, non_amb_parm
             do j = 1, i
                k = k + 1
                cov_parm_sav(i,j) = rdabuf(k)
                cov_parm_sav(j,i) = rdabuf(k)
             end do
          end do

*         Copy the static part of the solution vector.
          k = (tac_parm+1)*tac_parm/2
          do i = 1, non_amb_parm
             k = k + 1
             sol_vec_sav(i) = rdabuf(k)
          end do

****      Now move the ambiquity parameters.  Here we need to
*         be sure that we move them into the correct slots and
*         re-arrange the rows and columns
          k = (tac_parm+1)*tac_parm/2 + tac_parm
*         See how many ambiquities we need to estimate
          ne = 0 
          if( index(data_type,'L1').gt.0 .or.
     .         index(data_type,'LC').gt.0     ) ne = ne + 1
          if( index(data_type,'L2').gt.0     ) ne = ne + 1

*****     Set the mapping for the fixed parameters.  Here
*         the order is the same in the saved and current
*         matrix
          do i = 1, non_amb_parm
             npdest(i) = i
          end do

****      Check actual number of parameters used
          used_parm = non_amb_parm
          do i = non_amb_parm+1, tac_parm, ne
             if( idabuf(2*k+i).gt.0 ) used_parm = used_parm + ne
          end do

          do i = non_amb_parm+1, used_parm, ne

*            Clear the destination parameter numbers
             do j = 1, ne
                npdest(i+j-1) = 0
             end do

*            Get the ambiquity number for this parameter
             an = idabuf(2*k+i)

*            Now see what parameter number this is in this 
*            direction filter run
             if ( an.gt. 0 ) then
                np = abs(amb_parn(1,an)) 
                if( np.eq.0 ) then
*                   We have a problem.  This ambiquity does not seem
*                   to be estimated.
                    write(*,300) ep, i, an
 300                format('**ERROR alligning bias parameters at ',
     .                     ' epoch ',i6,' Entry ',i4,' Ambiquity ',
     .                     i4)
                    print *,'AmbParn', (amb_parn(1,j),j=1,num_ambs)
                    print *,'TAC ',k,tac_parm, num_parm, used_parm,
     .                      (npdest(j),j=1,used_parm)
                    if( ep.gt.10 ) stop

                else

*                   Save where the parameter needs to go.  Do this
*                   for the number of ambiquity types (ie. L1, L2, LC)
                    do j = 1, ne
                       npdest(i+j-1) = abs(amb_parn(j,an))
                    end do
                end if
             end if
          end do

*         OK, we now have the destination parameter numbers for
*         each of the input bias parameter numbers.
          do i = non_amb_parm+1, used_parm
             np = npdest(i)
             if( np.gt.0 ) then
                 sol_vec_sav(np) = 
     .             rdabuf((tac_parm+1)*tac_parm/2+i)

*               Now loop over the covariance matrix entries.
*               Compute the position of the diagonal element
*               of the parameter immediately before this one.
                jel = i*(i-1)/2
                do j = 1, non_amb_parm
                   cov_parm_sav(np,j) = rdabuf(jel+j)
                   cov_parm_sav(j,np) = rdabuf(jel+j)
                end do

*               Now do the rest of the ambiquity parameters
                jel = jel + non_amb_parm
                do j = non_amb_parm+1,i
                   nj = npdest(j)
* MOD TAH 080626: Only increment counter and save results if there is
*                 a destination parameter
                   if( nj.gt.0 ) then 
                      jel = jel + 1
                      cov_parm_sav(nj,np) = rdabuf(jel)
                      cov_parm_sav(np,nj) = rdabuf(jel)
                   else
                      jel = jel + 1
                   end if
                end do
             end if
          end do

****      See if debug on:
          if( ep.ge.debug_start .and. ep.le.debug_end ) then
             call dump_cov(ep, rdabuf, idabuf, num_parm,
     .                     npdest,'S')
          end if

      end if

***** Thats all
      return 
      end

CTITLE DUMP_COV

      subroutine dump_cov( ep, rdabuf, idabuf, np, npdest, option) 

      implicit none

*     Test routine to see if we are dumping things correctly

      include 'track_com.h'

      integer*4 ep, np, idabuf(*), npdest(*)
      real*8 rdabuf(*)
      character*(*) option

      integer*4 i,j,k 

*     Write out the buffer
      write(*,100) ep, np, num_parm, tac_parm
 100  format('DABUF DUMP: Ep ',i6,' NP, num, tot ',3i5,/,
     .       'Lower Diagonal Cov_parm ')
      k = 1
      do i = 1, np
         k = k + i - 1
         write(*,120) i, (rdabuf(j),j=k,k+i-1) 
  
 120     format(i5,(100D12.2))
      end do
      k = (tac_parm+1)*tac_parm/2
      write(*,140) (rdabuf(i),i=k+1,k+np)
 140  format('SOL: ',(100D12.2))

      k = (tac_parm+1)*tac_parm/2 + tac_parm
      write(*,160) (idabuf(i),i=2*k+1,2*k+np)
 160  format('AMBS ',100I5)

*     Now see if should write the reconstructed matrix
      if( option(1:1).eq.'S' ) then

          write(*,200) 
 200      format('COV_PARM_SAV')
          do i = 1, num_parm
             write(*,220) i,(cov_parm_sav(i,j),j=1,num_parm)
 220         format(i5,(100D12.2))
          end do
          write(*,240) (sol_vec_sav(i),i=1,num_parm)
 240      format('SOLS ',(100D12.2))
          write(*,260) (npdest(i),i=1,num_parm)
 260      format('NPDST',(100I5))
      end if

****  Thats all
      return
      end

CTITLE SMOOTH

      subroutine smooth(ep)

      implicit none

*     Routine to average the forward and backward running filters.
*     It is called after the predicition step of the forward
*     running filter since the backward running filter is run first.

* ep - epoch number

      integer*4 ep

      include 'track_com.h'

* LOCAL PARAMETERS
* i,j  -- Loop counters
* ipivot(max_parm) -- Pivot elements for inversion

      integer*4 i, j, ipivot(max_parm)

* jmat(max_parm, max_parm) -- Kalman gain matrix for doing the 
*     averaging
* scale(max_parm)          -- Scale for invert_vis plus general array
*     for saving intermediate results

      real*8 jmat(max_parm, max_parm), scale(max_parm) 


****  If this is the first epoch, just use the saved results
      if( ep.eq.1 ) RETURN

****  First add the forward propagated covariance matrix to the 
*     one from the backward run
      do i = 1, num_parm
         call dwadd(cov_parm_sav(1,i),1,cov_parm(1,i),1,
     .              cov_parm_sav(1,i),1, num_parm)
      end do

****  Now invert this matrix
      call invert_vis(cov_parm_sav, sol_vec_sav, scale, ipivot,
     .                num_parm, max_parm, 0)


****  Now form the kalman gain matrix.
      do i = 1, num_parm
         do j = 1, num_parm
            call dwdot(scale(j),cov_parm(1,i),1,
     .                          cov_parm_sav(1,j), 1, num_parm)
         end do
*        Move the column into jmat
         call dwmov(scale, 1, jmat(1,i),1, num_parm)
      end do

*     Form the difference between the two solution vectors
      call dwsub(sol_vec_sav,1, sol_vec, 1, sol_vec_sav,1, num_parm)

*     Multiply the difference by the Kalman gain
      do i = 1, num_parm
         call dwdot(scale(i), jmat(1,i),1,sol_vec_sav,1,num_parm)
      end do

*     Save the final answer in sol_vec_sav
      do i = 1, num_parm
         sol_vec_sav(i) = sol_vec(i) + scale(i)
      end do

*     Now complete the covariance matrix: Again save result in 
*     cov_parm_sav
      do i = 1, num_parm
         do j = 1, i
            call dwdot(scale(j),jmat(1,j),1, cov_parm(1,i),1, 
     .                 num_parm)
            cov_parm_sav(i,j) = cov_parm(i,j) - scale(j)
            cov_parm_sav(j,i) = cov_parm_sav(i,j)
         end do
      end do

****  Thats all:
      return 
      end


         





