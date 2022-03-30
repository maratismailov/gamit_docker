CTITLE RD_COV

      subroutine rd_cov( unit, cov_parm, sol_parm,qnum_parn
     .                   ,qnum_sites,qnum_svs,apr,qsol,qparn_sites
     .                   ,qparn_svs,qparn_pmu,num_orb,scale_sv
     .                   ,scale_pmu,site_pos,llr_xyz,qscale )

*     This routine reads the covarinance matrix. Once it is read,
*     it and the solution vector are scaled and rotated (for the
*     site coordinates)
c
c NOTE: This is essentially the Tom Herring routine READ_COV.f from htoglb
c       but slightly modified to produce the information required for
c       creating a sinex file
c
c Modifications by P Tregoning
c August, 1995


* PASSED VARIABLES

*   unit        - Input hfile unit

      implicit none

      include '../includes/dimpar.h'

      integer*4 unit,qnum_parn,max_parn,qnum_sites,qnum_svs,num_orb
      parameter(max_parn = maxsit*3+maxsat*maxorb+6)

*   cov_parm(qnum_parn, qnum_parn)    - Covariance matrix
*   sol_parm(qnum_parn)          - Solution vector

      real*8 cov_parm(qnum_parn,qnum_parn),sol_parm(qnum_parn)
     .      ,qscale(qnum_parn),qsol(qnum_parn),apr(qnum_parn)

* LOCAL VARIABLES
*   ierr        - IOSTAT error
*   trimlen     - Length of string
*   indx        - Pointer in string
*   i,j         - Loop counter
*   nq, np      - Parmeter number for rotating solution and
*               - covariance matrix

      integer*4 ierr, trimlen, indx, i,j,k,l, nq, np,count,
     .    qparn_sites(3,qnum_sites),
     .    qparn_svs(maxorb,qnum_svs),
     .    qparn_pmu(2,3)


*   umat        - 3*3 unit matrix used for rotating the covariance
*                 elements between the station coordinates and the
*                 satellite orbits.  (For convenience we do the
*                 sv elements 3 at a time.  There are 9 for sv)
c   scaling factors to use for coords, orbits and eop

       real*8 umat(3,3),llr(3,maxsit),site_pos(3,maxsit)
     .        ,rot_mat(3,3,maxsit)
     .        ,svs_pos(maxorb,maxsat),pmu_pos(2,3),scale_sv(15)
     .        ,scale_pmu(6),llr_xyz(3*maxsit,3)

       character*3 datum_type
       parameter (datum_type = 'LLR')

*   line        - Line read from file

      character*256 line,message

      data umat / 1.d0, 0.d0, 0.d0,  0.d0, 1.d0, 0.d0,
     .            0.d0, 0.d0, 1.d0 /


c  first assign the coords to llr matrix
          do i =1,qnum_sites
             llr(1,i) = apr(i*3-2)
             llr(2,i) = apr(i*3-1)
             llr(3,i) = apr(i*3)
          enddo

c       call printmat(llr,3,qnum_sites,'sites ')

c  set up the qscale matrix
      do i=1,qnum_sites
           qscale(qparn_sites(1,i)) = llr(3,i)*1000.d0
           qscale(qparn_sites(2,i)) = llr(3,i)* cos( llr(1,i) )*1000.d0
           qscale(qparn_sites(3,i)) = 1000.d0
      enddo

      do i=1,qnum_svs
        do j=1,num_orb
          qscale(qparn_svs(j,i)) = scale_sv(j)
        enddo
      enddo

      do i=1,3
        do j=1,2
          if(qparn_pmu(j,i).ne.0)then
            qscale(qparn_pmu(j,i)) = scale_pmu(2*i-2+j)
          endif
        enddo
      enddo

c       call printmat(qscale,qnum_parn,1,'qscale ')


****  Skip the next (check to see that it says covariance)
      read(unit,'(a)', iostat=ierr ) line
      indx = index(line,'Covariance')
*                             ! We have a problem
      if( indx.eq.0 ) then
          write(*,100) line(1:max(1,trimlen(line)))
 100      format(' ** DISASTER ** Covariance matrix not where'
     .        ,' expected.  Line read is:',/,a)
          RETURN
      end if

***** Now loop over all of the parameters which have been estimated
*     Read the column of the matrx
      print*,' Reading the covariance matrix ....'
      do i = 1, qnum_parn
c       print*,' reading hfile',i
          read(unit,200, iostat=ierr) (cov_parm(j,i), j= 1,i)
 200      format(1x,5x,5000(5(1x,d23.16),:,/,6x))
          if( ierr.ne.0 ) then
              cov_parm(i,i) = 1.d0 
              write(message,'(a,i4,a)') 'Error reading row',i
     .          ,' of covariance matrix'
              call report_stat('WARNING','HTOSNX','utils/rd_cov'
     .          ,' ',message,ierr) 
          end if
      end do

****  Now fill in the symetric part of the matrix (This code actually
*     copies the diagonal as well but this should be no problem)
      do i = 1, qnum_parn
          call dwmov(cov_parm(1,i),1, cov_parm(i,1), qnum_parn,i)
      end do

****  Now rescale the matrix (meter everywhere except solar radiation
*     which will be ratio to direct effect)
      do i = 1, qnum_parn
        do j=1,qnum_parn
           cov_parm(i,j) = qscale(i)*cov_parm(i,j)
           cov_parm(j,i) = qscale(i)*cov_parm(j,i)
        enddo
        sol_parm(i) = qsol(i)*qscale(i)

      end do
****  Now rotate the station coordinates from local frame to cartesian
*     XYZ.  If needed

      if( datum_type.eq.'LLR' ) then

c  now assign the sat coords to svs_pos
          do i = 1,qnum_svs
             do j=1,num_orb
               svs_pos(j,i) = apr(qparn_svs(j,i))*scale_sv(j)
             enddo
          enddo
c       call printmat(svs_pos,maxorb,qnum_svs,'satellites ')

c  now assign eop values
          count = 0
          do i=1,3
            do j=1,2
            count = count + 1
              if(qparn_pmu(j,i).ne.0)then
                pmu_pos(j,i) = apr(qparn_pmu(j,i))*scale_pmu(count)
              else
                pmu_pos(j,i) = 0.d0
              endif
            enddo
          enddo

          do i = 1, qnum_sites

*             Rotate the site coordinates
              call llr_to_xyz( llr(1,i), site_pos(1,i), rot_mat(1,1,i))
c PT950818: store the rotation matrix for use later
              do k = 1,3
                do l = 1,3
                  llr_xyz(3*i-3+k,l) = rot_mat(k,l,i)
                enddo
              enddo
          end do

*****     Now rotate the solution vector and the covariance matrix
          do i = 1, qnum_sites

*             Solution vector
              nq = qparn_sites(1,i)
              if( nq.ne.0 ) then
                  call rotate_sol( rot_mat(1,1,i), sol_parm(nq))

*****             Now rotate the covariance matrx
                  do j = 1, qnum_sites
                      np = qparn_sites(1,j)
                      if( np.gt.0 ) then
                          call rotate_cov(rot_mat(1,1,i),
     .                        cov_parm(nq,np), qnum_parn,
     .                        rot_mat(1,1,j), 3, 3 )
                      end if
                  end do

*****             Now do the rotation of the satellite/site covariance
*                 elements.  Do in groups of three moving accross the
*                 rows.  When we are finished we will do the move the
*                 values to the columns.
                  do j = 1, qnum_svs
                     do k = 1,3
                        np = qparn_svs((k-1)*3+1,j)
                        if( np.gt.0 ) then
                            call rotate_cov(rot_mat(1,1,i),
     .                           cov_parm(nq,np), qnum_parn,
     .                           umat, 3, 3 )
                        end if
                     end do
                  end do

*****             Now do the polar motion/Ut1 estimates.  Do one at
*                 a time since we don't know how many there will be
                  count = 0
                  do j = 1, 3
                     do k = 1,2
                         count = count + 1
                        np = qparn_pmu(k,j)
                        if( np.gt.0 ) then
                            call rotate_cov(rot_mat(1,1,i),
     .                           cov_parm(nq,np), qnum_parn,
     .                           umat, 3, 1 )
                        end if
                     end do
                  end do

              end if
*                     ! Looping over sites
          end do

****      Now fill in the symetric part of the matrix (This code actually
*         copies the diagonal as well but this should be no problem)
          do i = 1, qnum_parn
              call dwmov(cov_parm(1,i),1, cov_parm(i,1), qnum_parn,i)
          end do
*                     ! coordinates need rotating
      end if

****  Now construct the total solution (add apriori and adjustment
      do i = 1,qnum_sites
          do j = 1,3
              nq = qparn_sites(j,i)
              if( nq.gt.0 ) then
                  sol_parm(nq) = site_pos(j,i) + sol_parm(nq)
              end if
c PT951002: currently cannot handle velocity estimates so comment this out
c              nq = qparn_vel(j,i)
c              if( nq.gt.0 ) then
c                  sol_parm(nq) = site_vel(j,i) + sol_parm(nq)
c              end if
          end do
      end do

***** Repeat for satellites
      do i = 1,qnum_svs
          do j = 1,num_orb
              nq = qparn_svs(j,i)
              if( nq.gt.0 ) then
                  sol_parm(nq) = svs_pos(j,i) + sol_parm(nq)
              end if
          end do
      end do

***** Now add adjusmtent to get total values
      do i = 1,3
         do j = 1,2
            nq = qparn_pmu(j,i)
            if( nq.gt.0 ) then
               sol_parm(nq) = pmu_pos(j,i) + sol_parm(nq)
            end if
         end do
      end do

****  Thats all.  We are now ready to generate the solution
      return
      end


