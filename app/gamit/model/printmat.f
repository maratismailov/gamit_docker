        subroutine printmat(mat,m,n)

c  prints out a m x n matrix
c  PT 9/6/92

        implicit none
        integer m,n,i,j
        real*8 mat(m,n)
        character*1 ent

        do 10 i=1,m
            write(*,*)(mat(i,j),j=1,n)
10      continue

        print*,' '
        write(*,'(a)')' press return to continue : '
        read(*,'(a)')ent

        return
        end

