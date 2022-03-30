CTITLE WRITE_GLB_COREL
      subroutine write_glb_corel( iout, options, cov_parm,
     .                            sol_parm )
 
      implicit none
 
*     Routine to write out the correlation matrix
*
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
 
*   iout    - Output LU
*   i,j,k   - Loop counters
*   options - OPtions parameter for output
 
      integer*4 iout, i,j, k, options

*   kbit    - Checks bit status:  1 - correlation matrix
*                                13 - Covariance matrix

      logical kbit
 
 
*   cov_parm(num_glb_parn,num_glb_parn) - Covariance matrix
*           - of estimated parameter
*   sol_parm(1j)    - Solution vector
*   covar(max_glb_parn)     - Correlations for one line
 
      real*8 cov_parm(num_glb_parn,num_glb_parn),
     .    sol_parm(num_glb_parn), covar(max_glb_parn)

* Format
      character*64 neucform   ! Format for NECU output specific
                              ! to number of paranmeters
 
 
      common covar
 
*     Write out header

      if( kbit(options,1) ) then
 
          write(iout,100)
  100     format(/,' CORRELATION MATRIX')
 
          do i = 1, num_glb_parn
              do j = 1,i
                  if(  cov_parm(i,i).gt.0.d0 .and. 
     .                 cov_parm(j,j).gt. 0.d0 ) then
                       covar(j) = cov_parm(i,j)/ sqrt( cov_parm(i,i)*
     .                                            cov_parm(j,j) )
                  else
                       covar(j) = 0.d0
                  end if
              end do
 
              write(iout,120) i,(covar(j),j=1,i)
  120         format(i5,'. ', 10(1x,f7.4):,/,200(7x,10(1x,f7.4),:/) )
          end do
      end if 

***** See if full covariances
      if( kbit(options,13) ) then
          write(iout,200)
  200     format(/,' COVARIANCE MATRIX (lower diagonal)',/,
     .           '   P1    P2    Covariance (internal units) ')

          do i = 1, num_glb_parn
             do j = 1, i
                  write(iout,220) i,j, cov_parm(i,j)
  220             format(2i5,1x,D25.16)
             end do
          end do
      end if

***   See if the print NEU correlation has been selected.
      if( kbit(options, 31) ) then

         write(iout,300)
  300    format(/,'N CORRELATION MATRIX WITH NEUC (lower diagonal)',/,
     .          'Label   Row Percent correlations by column ')


****     Rotate the covariance matrix into NEU form
         call cov_xyz_neu( +1, cov_parm, sol_parm)

*        Now write out values
C MOD TAH 161204: Create format lines
         write(neucform,410) num_glb_parn
 410     format('("NEUC ",I6,2x,',I5,'(I4,1x))')
         do i = 1, num_glb_parn

C            write(iout,320) (i, j*20+1,
C    .         (nint(100*cov_parm(i,j*20+k)/ sqrt( cov_parm(i,i)*
C    .                   cov_parm(j*20+k,j*20+k)+1.d-20)),
C    .                         k=1,min(20,i-j*20)),j=0,(i-1)/20)
C 320        format('N ',I5,1x,I5,20(1x,I4),:,/,
C    .                    1000(:,'N ',I5,1x,i5,20(1x,I4),/))
C MOD TAH 161204: Updated format so that file can be read directly
C        in matlab (after stripping NEUC from start of lines).
             write(iout,neucform) I,
     .            (nint(100*cov_parm(i,k)/ sqrt( cov_parm(i,i)*
     .                   cov_parm(k,k)+1.d-20)),k=1,num_glb_parn)
         end do
 
****     Now map the covarinace matrix back
         call cov_xyz_neu( -1, cov_parm, sol_parm)
      endif
 
****  Thats all
      return
      end
 
 
