!!$!===========================================================================
!!$!=======Inverse of a general matrix ========================
!!$!===========================================================================
MODULE invmod
  USE nrtype
contains
  SUBROUTINE  inverse(a,b)
    USE f95_lapack
    IMPLICIT NONE

    !	purpose: getting the inverse b of a general matrix a
    !	input: the general matrix a (nxn) to be inversed
    !	output: the inverse matrix b

    integer :: info, n
    real(DP) :: rcond
    real(DP), dimension(:,:) :: a, b
    integer, allocatable, dimension(:) :: indx

    n = size(a,1)

    if (size(a,2) /= n) then
       write(0,*) 'matrix a not square in inverse(a,b)!! quitting'
       STOP 77
    end if

    allocate(indx(n))

    b(:,:) = a(:,:)

    CALL LA_GETRF(b,indx,rcond,INFO=info)
    CALL LA_GETRI(b,indx,info)

    deallocate(indx)

    return
  end SUBROUTINE inverse
end MODULE invmod
