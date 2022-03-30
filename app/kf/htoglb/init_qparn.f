      subroutine init_qparn( nsites, nsvs )
 
      implicit none

*     This routine will clear the site and svs parameter number
*     pointer arrays so that they will be set freshly for each
*     solution in the hfiles and for subsquent hfiles during a
*     run

      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
 
*   nsites, nsvs    - Number of sites and satellites
 
      integer*4 nsites, nsvs
 
* LOCAL PARAMETERS
 
*   i,j             - Loop counters
 
 
      integer*4 i,j, k
 
****  Loop over each parn array and set to zero
      atm_num = 0
      do i = 1, nsites
          do j = 1,3
              qparn_sites(j,i) = 0
              qparn_vel(j,i) = 0
          end do
	  qparn_axo(i) = 0
          qparn_atm(i) = 0
      end do

*     PMU
      do i = 1,3
          do j = 1,2
              qparn_pmu(j,i) = 0
              do k = 1, max_mul_pmu
                 qparn_mul_pmu(j,i,k) = 0
              end do
              qnum_mul_pmu(j,i) = 0
          end do
      end do

 
 
*     Satellites
      do i = 1, nsvs
          do j = 1,max_svs_elem
              qparn_svs(j,i) = 0
          end do
      end do

      qnum_misc = 0
 
****  Thats all, any other parameter arrays should be added
*     later.
      return
      end
 
CTITLE INIT_SITE
 
      subroutine init_site( pos, vel, num)
 
      implicit none

*     This routine will clear the site positions and velocities
*     before a h-file is read.  (Here we want the velocity cleared
*     in case there are none in the h-file).
 
*   num     - Number of sites
 
      integer*4 num
 
*   pos(3,num)  - Site positions (XYZ)
*   vel(3,num)  - Site velocities (XYZ dot)
 
      real*8 pos(3,num), vel(3,num)
 
* LOCAL VARIABLES
 
*   i,j     - Loop counters
 
      integer*4 i,j
 
****  Loop over arrays clearing the values
 
      do i = 1, num
          do j = 1,3
              pos(j,i) = 0.d0
              vel(j,i) = 0.d0
          end do
      end do
 
***** Thats all
      return
      end

CTITLE INIT_SVS 
 
      subroutine init_svs( svs_pos, num_svs, num_elem )

      implicit none

*     Routine to initialze satellite parameters

      integer*4 num_svs, num_elem
      real*8 svs_pos(num_elem, num_svs )

      integer*4 i,j

      do i = 1, num_svs
         do j = 1, num_elem
            svs_pos(j,i) = 0.d0
         end do
      end do

      return
      end

CTITLE INIT_SVINF

      subroutine init_svinf

      implicit none

*     Rouitne to clear all the satellite infomration
      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'

      integer*4 i,j,k

* MOD TAH 050622: Clear all the sinex entrires
      do i = 1, max_glb_svs
         qsvi_prn(i) = 0
         qsvi_svn(i) = 0
         qsvi_block(i) = 0
         qsvi_antmod(i) = 'NONE'
         qsvi_ocode(i) = ' '
         do j = 1,3
            do k = 1,2 
               qsvi_antpos(j,k,i) = 0.d0
            end do
         end do
         qsvi_launch(i) = 0.d0
      end do


      return
      end
 
