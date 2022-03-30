CTITLE READ_APR_COV

      subroutine read_apr_cov(unit, line)
      
      implicit none

*     Routine to read the apriori covariance matrix from a
*     GAMIT hfile

      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'

* PASSED VARIABLES
 
*   unit        - INput hfile unit
 
      integer*4 unit

*   line        - Line read from input

      character*(*) line

* LOCAL VARIABLES
*   ierr        - IOSTAT error
*   trimlen     - Length of string
*   i,j         - Loop counter
*   nq, np      - Parmeter number from the original hfile
*   mq, mp      - Parameter number in the output solution
*   ns, nv      - Site and satellite numbers
*   jerr        - Error reading line
* MOD TAH 080416: Account for possible presence of atmospheric delays.
*   begsvs_parn - Parameter number before start of satellites

 
      integer*4 ierr, trimlen, i,j,k, nq, np, mp, mq,
     .          ns, nv, jerr, begsvs_parn

*   cov         - Covariance element read.

       real*8  cov

*   finished    - Indicates that we have finished reading these
*                 entries

      logical finished
       
****  Start looping over the entries
      finished = .false.

* MOD TAH 080416: Get the paramter number before start of SVS elements.
*     Code also accounts for possible non-estimated site at end of list.
      do i = 1, qnum_sites
         if( qparn_sites(3,i).gt.0 ) 
     .       begsvs_parn = qparn_sites(3,i)
      enddo
*     See if atmospheric delays are present
      do i = 1, qnum_sites
         if( qparn_atm(i).gt.0 ) 
     .       begsvs_parn = qparn_atm(i)
      enddo

      
      call init_qapr_cov
           
      do while ( .not.finished )
         read(unit,'(a)', iostat=ierr ) line
         if( ierr.ne.0 .or. trimlen(line).eq.0 ) then
            finished = .true.
            line = ' '
         else
            read(line,*, iostat=jerr ) np, nq, cov
            if( jerr.ne.0 ) finished = .true.
         end if
         
         if( .not.finished ) then
*           Now decide what type of parameter this is.  Check the
*           to see if the parameter is still in the solution
            mp = 0
            mq = 0
            do i = 1, qnum_parn
               if( np.eq.atoo(i) ) mp = i
               if( nq.eq.atoo(i) ) mq = i
            end do
            
            if( mp.gt.0 .and. mq.gt.0 ) then

*               These parameters are still in the solution.  Find
*               out what they are.
                if( mp.le. qparn_sites(3,qnum_sites)) then

*                   This is a site position.  Save in the site position
*                   array.
                    ns = (mp-1)/3 + 1
                    k  = mp - (ns-1)*3
                    qapr_cov(k,k,ns) = cov

****            Ok See if satellite orbits IC (ie.X Y Z Xdot Ydot Zdot)
                else if( mp.gt.begsvs_parn .and.
     .              mp.le.qparn_svs(6,qnum_svs)    ) then


*                   See if IC
                    nv = (mp-begsvs_parn-1)/(6+rad_num)+1
                    k = mp-begsvs_parn
     .                    - (nv-1)*(6+rad_num)
                    if( k.le.6 ) then
*                       Get the lower diagonal element number                     
                        j = mq-mp+k
                        qapr_svs(j,k,nv) = cov
                        qapr_svs(k,j,nv) = cov
                    else
                        qnum_apr_diag =  qnum_apr_diag + 1
                        qapr_diag(mp) = cov
                    end if
                else

****                Must be diagonal element for EOR                
                    qnum_apr_diag =  qnum_apr_diag + 1
                    qapr_diag(mp) = cov       
                end if                                                                
            end if
         end if
      end do
      
***** Save the numbers of values
      qnum_apr_cov = qnum_sites
      qnum_apr_svs = qnum_svs

***** Now do the conversion of units.  First scale all the entries
*     As read from the file the LLR entries are diagonal so thats
*     all we need to scale
   
      do i = 1, qnum_sites
         do j = 1,3
            np = qparn_sites(j,i)
            qapr_cov(j,j,i) = qapr_cov(j,j,i)*qscale(np)**2
         end do
      end do
      
      do i = 1, qnum_parn
         qapr_diag(i) =  qapr_diag(i)* qscale(i)**2
      end do

*     Now do the satellites.  These are non-diagonal and so we need
*     to do completely
      do i = 1, qnum_svs
          do j = 1, 6
             np = qparn_svs(j,i)             
             do k = 1,6
                nq = qparn_svs(k,i)
                qapr_svs(k,j,i) =  qapr_svs(k,j,i)* qscale(np)*
     .                                              qscale(nq) 
             end do
          end do
      end do           
      
      if( datum_type.eq.'LLR' ) then
      
*****    Now rotate the covariance matrx
         do j = 1, qnum_sites
             np = qparn_sites(1,j)
             if( np.gt.0 ) then
                 call rotate_cov(rot_mat(1,1,j),
     .                        qapr_cov(1,1,j), 3,
     .                        rot_mat(1,1,j), 3, 3 )
              end if
         end do
      end if
      
***** Thats all
      return
      end
      
