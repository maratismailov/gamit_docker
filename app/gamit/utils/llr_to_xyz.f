CTITLE LLR_TO_XYZ

      subroutine llr_to_xyz ( llr, xyz, rot_mat )

*     Routine to convert from lat, long and radius to XYZ cartersian
*     and to returm the rotation matrix from llr to xyz.

* PASSED VARIABLES

*   llr(3)  - lat, long, radius (rad, rad, kms)
*   xyz(3)  - XYZ geocentric (m)
*   rot_mat(3,3)    - Rotation from NEU ro XYZ


      real*8 llr(3), xyz(3), rot_mat(3,3)

****  Get the XYZ coordinates
      xyz(1) = llr(3)*cos(llr(1))*cos(llr(2)) * 1000.d0
      xyz(2) = llr(3)*cos(llr(1))*sin(llr(2)) * 1000.d0
      xyz(3) = llr(3)*sin(llr(1)) * 1000.d0

****  Now compute the rotation matrix
      rot_mat(1,1) = -sin(llr(1)) * cos(llr(2))
      rot_mat(1,2) = -sin(llr(2))
      rot_mat(1,3) =  xyz(1)/llr(3) / 1000.d0

      rot_mat(2,1) = -sin(llr(1)) * sin(llr(2))
      rot_mat(2,2) =  cos(llr(2))
      rot_mat(2,3) =  xyz(2)/llr(3) / 1000.d0

      rot_mat(3,1) =  cos(llr(1))
      rot_mat(3,2) =  0.d0
      rot_mat(3,3) =  xyz(3)/llr(3) / 1000.d0

****  Thats all
      return
      end


