MODULE singleel
  USE nrtype
  IMPLICIT NONE
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: singleevals
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: singleevecs
  INTEGER(I4B), PARAMETER :: ndim = 1
  INTEGER(I4B) :: CoulLmax
  REAL(DP), PARAMETER :: Z = 2.d0
  LOGICAL(LGT) :: singleel_solved = .false.
contains
  ! ------------------------------------------------------------------------
  subroutine dvr_setup
    ! ------------------------------------------------------------------------ 
    USE pfedvrmod
    USE dvrmod
    implicit NONE
    INTEGER(I4B) :: i, j, k, m
    !--------------------------------------------------------------------------
    ! Set up dvr points, weights and kinetic energy matrix  for each region
    !--------------------------------------------------------------------------

    DO j=1, ndim
       do i=1, num_reg(j)

          CALL DVR_reg( num_fun(j,i),bounds(j,i,:),mat_reg(j,i)%pt, &
               mat_reg(j,i)%wt,mat_reg(j,i)%df,mat_reg(j,i)%ddf, &
               mat_reg(j,i)%ke_mat )

          !       write(*,*)'j=',j,' i=',i,'  Seting up Kinetic-energy matrice......'

       end do
    END DO

    !--------------------------------------------------------------------------------
    !------------------following Nicolai's convention to get Kinetic Matrice---------
    !--------------------------------------------------------------------------------

    DO j=1, ndim
       DO i=1, num_reg(j)
          do k=1,num_fun(j,i)  

             if(mat_reg(j,i)%wt(k).le.0.d0)then
                write(*,*)'i=',i,'  k=',k,' mat_reg(j,i)%wt(k)=',mat_reg(j,i)%wt(k)
                pause
             endif

             if((i.gt.1).and.(k.eq.1)) then
                mat_reg(j,i)%fac1(k)=1.d0/dsqrt(mat_reg(j,i)%wt(1)+mat_reg(j,i-1)%wt(num_fun(j,i-1)))            
             elseif((i.lt.num_reg(j)).and.(k.eq.num_fun(j,i))) then
                mat_reg(j,i)%fac1(k)=1.d0/dsqrt(mat_reg(j,i)%wt(num_fun(j,i))+mat_reg(j,i+1)%wt(1))
             else
                mat_reg(j,i)%fac1(k)=1.d0/dsqrt(mat_reg(j,i)%wt(k))
             end if
             do m=1,num_fun(j,i)
                if((i.gt.1).and.(m.eq.1)) then
                   mat_reg(j,i)%fac2(m)=1.d0/dsqrt(mat_reg(j,i)%wt(1)+mat_reg(j,i-1)%wt(num_fun(j,i-1)))            
                elseif((i.lt.num_reg(j)).and.(m.eq.num_fun(j,i))) then
                   mat_reg(j,i)%fac2(m)=1.d0/dsqrt(mat_reg(j,i)%wt(num_fun(j,i))+mat_reg(j,i+1)%wt(1))
                else
                   mat_reg(j,i)%fac2(m)=1.d0/dsqrt(mat_reg(j,i)%wt(m))
                end if

                mat_reg(j,i)%ke_mat(k,m)=mat_reg(j,i)%ke_mat(k,m)*mat_reg(j,i)%fac1(k)*mat_reg(j,i)%fac2(m) 

                !            write(*,*)'j=',j,' i=',i,'  k=',k,'  m=',m

             end do
          end do
       END DO

    END DO
    !---------------------------------------------------  
  end subroutine dvr_setup

  !======================================================================
  !======================================================================
  !    setup grids for the whole space
  !======================================================================
  !    By Suxing Hu at LANL, references from Barry and Nicolai's DVR codes
  !
  !    Written on:  06/01/2005
  !    Revised on:  mm/dd/yyyy
  !======================================================================
  ! ----------------------------------------------------------------------
  subroutine setup_grid_mapping(nindex,grid,factor,r0BC)
    ! ----------------------------------------------------------------------
    !  Construct full spatial grid from the individual elements for each dimension
    USE pfedvrmod
    Implicit NONE
    REAL(DP) , dimension(:,:) :: grid
    REAL(DP) , dimension(:,:) :: factor 
    INTEGER(I4B), dimension(:,:,:) :: nindex
    !-------------Total number of grid points for each dimension--------------- 

    INTEGER(I4B) :: i, j, k, m, start, start2, end
    LOGICAL(LGT) :: r0BC

    !=========================================================================
    !-----------dim_reg_bdp(j,i,1) is the starting point # for the 
    !-----------jth-dimension and ith-element------------------------------
    !   
    !-----------dim_reg_bdp(j,i,2) is the ending point # for the 
    !-----------jth-dimension and ith-element------------------------------
    !
    !  Storing these numbers are important for correctly operating on 
    !  the wavefunction; Namely, these numbers will tell what points should
    !  be operated by the jth-dimensional ith element!
    !
    !=========================================================================

    DO j=1, ndim

       start=0
       do i=1,num_reg(j)
          if(i.eq.1) then
             if(r0BC) then
                start=-1
                start2=2
                end=num_fun(j,i)-1
             else
                start2=1
                end=num_fun(j,i)-1     
             end if
          else
             start2=1
             end=num_fun(j,i)-1
          end if

          !****************************************************

          do k=start2,num_fun(j,i)

             m=start+k

             nindex(j,i,k)=m    

             grid(j,m)=mat_reg(j,i)%pt(k) 

             factor(j,m)= mat_reg(j,i)%fac1(k)

          end do
          start=start+end


          !****************************************************

       end do

    END DO

  end subroutine setup_grid_mapping
END MODULE singleel
