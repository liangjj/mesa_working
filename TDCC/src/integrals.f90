module integrals
  interface integral
     module procedure integral_r, &
          &           integral_c
  end interface
  interface integral2d
     module procedure integral2d_r, &
          &           integral2d_c
  end interface
  interface integral3d
     module procedure integral3d_r, &
          &           integral3d_c
  end interface
contains
  !======================================================================
  !======================================================================
  ! ----------------------------------------------------------------------
  subroutine integral_r(dim,wx,fx,int)
    ! ----------------------------------------------------------------------
    ! Performs integral of function fx using the quadrature rule

    implicit none

    integer*4 :: dim
    real*8 , dimension(1:dim) :: wx
    real*8 , dimension(1:dim) :: fx
    real*8, INTENT(OUT)       :: int

    int=0.d0

    int=sum(wx(:)*fx(:))

  end subroutine integral_r

  !======================================================================
  !======================================================================
  ! ----------------------------------------------------------------------
  subroutine integral2d_r(dimx,dimy,wx,wy,fxy,int)
    ! ----------------------------------------------------------------------
    ! Performs integral of function fx using the quadrature rule

    implicit none

    integer*4 :: dimx, dimy, i, j
    real*8 , dimension(1:dimy) :: wy
    real*8 , dimension(1:dimx) :: wx
    real*8 , dimension(1:dimx,1:dimy) :: fxy
    real*8, INTENT(OUT)       :: int 
    real*8                    :: temp 

    int=0.d0

    !   do i=1, dimx
    !     do j=1, dimy
    !       int=int+fxy(i,j)*wx(i)*wy(j)
    !     end do
    !   end do

    do i=1, dimx
       temp=0.d0
       do j=1, dimy
          temp=temp+fxy(i,j)*wy(j)
       end do
       int=int+temp*wx(i)
    end do

    return

  end subroutine integral2d_r

  ! ----------------------------------------------------------------------
  subroutine integral3d_r(dimx,dimy,dimz,wx,wy,wz,fxyz,int)
    ! ----------------------------------------------------------------------
    ! Performs integral of function fx using the quadrature rule

    implicit none

    integer*4 :: dimx, dimy, dimz,i, j,k
    real*8 , dimension(1:dimz) :: wz
    real*8 , dimension(1:dimy) :: wy
    real*8 , dimension(1:dimx) :: wx
    real*8 , dimension(1:dimx,1:dimy,1:dimz) :: fxyz
    real*8, INTENT(OUT)       :: int  
    real*8                    :: temp, temp1 

    int=0.d0

    do i=1, dimx
       temp=0.d0
       do j=1, dimy
          temp1=0.d0
          do k=1, dimz
             temp1=temp1+fxyz(i,j,k)*wz(k)
          end do
          temp=temp+temp1*wy(j)
       end do
       int=int+temp*wx(i)
    end do

    return

  end subroutine integral3d_r



  !======================================================================
  !======================================================================

  ! ----------------------------------------------------------------------
  subroutine integral_c(dim,wx,fx,int)
    ! ----------------------------------------------------------------------
    ! Performs integral of function fx using the quadrature rule

    implicit none

    integer*4 :: dim
    real*8 , dimension(1:dim) :: wx
    complex*16 , dimension(1:dim) :: fx
    complex*16, INTENT(OUT)       :: int

    int=dcmplx(0.d0,0.d0)

    int=sum(wx(:)*fx(:))

  end subroutine integral_c



  !======================================================================
  !======================================================================

  ! ----------------------------------------------------------------------
  subroutine integral2d_c(dimx,dimy,wx,wy,fxy,int)
    ! ----------------------------------------------------------------------
    ! Performs integral of function fx using the quadrature rule

    implicit none

    integer*4 :: dimx, dimy, i, j
    real*8 , dimension(1:dimy) :: wy
    real*8 , dimension(1:dimx) :: wx
    complex*16, dimension(1:dimx,1:dimy) :: fxy
    complex*16, INTENT(OUT)       :: int 
    complex*16                    :: temp 

    int=dcmplx(0.d0,0.d0)

    do i=1, dimx
       temp=dcmplx(0.d0,0.d0)
       do j=1, dimy
          temp=temp+fxy(i,j)*wy(j)
       end do
       int=int+temp*wx(i)
    end do

    return

  end subroutine integral2d_c

  ! ----------------------------------------------------------------------
  subroutine integral3d_c(dimx,dimy,dimz,wx,wy,wz,fxyz,int)
    ! ----------------------------------------------------------------------
    ! Performs integral of function fx using the quadrature rule

    implicit none

    integer*4 :: dimx, dimy, dimz,i, j,k
    real*8 , dimension(1:dimz) :: wz
    real*8 , dimension(1:dimy) :: wy
    real*8 , dimension(1:dimx) :: wx
    complex*16 , dimension(1:dimx,1:dimy,1:dimz) :: fxyz
    complex*16 , INTENT(OUT)       :: int
    complex*16                     :: temp, temp1 

    int=dcmplx(0.d0,0.d0)

    do i=1, dimx
       temp=dcmplx(0.d0,0.d0)
       do j=1, dimy
          temp1=dcmplx(0.d0,0.d0)
          do k=1, dimz
             temp1=temp1+fxyz(i,j,k)*wz(k)
          end do
          temp=temp+temp1*wy(j)
       end do
       int=int+temp*wx(i)
    end do

    return

  end subroutine integral3d_c

end module integrals
