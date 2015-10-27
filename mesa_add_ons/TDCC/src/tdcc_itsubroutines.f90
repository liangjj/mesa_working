module tdcc_itsubroutines
contains
  !====================================================================================
  !	SUBROUTINE  timfield
  !====================================================================================
  subroutine timfield( t, fldx)

    USE globalmod
    USE pfedvrmod

    IMPLICIT NONE

    REAL*8, INTENT(IN) :: t

    REAL*8  :: fldx, fldy 


    IF(t.le.durationx)THEN
       fldx=efieldx*dsin(pi*t/durationx)**2*dcos(omegax*t+phasex)
    ELSE
       fldx=0.d0
    ENDIF

    !	IF(t.le.durationy)THEN
    !	 fldy=efieldy*dsin(pi*t/durationy)**2*dcos(omegay*t+phasey)
    !	ELSE
    !	 fldy=0.d0
    !	ENDIF

    return
  end subroutine timfield

  !======================================================================
  !======================================================================
  !   Setup DVR Kinetic-energy matrice for each dimension & region
  !======================================================================
  !    By Suxing Hu at LANL, references from Barry and Nicolai's DVR codes
  !
  !    Written on:  5/31/2005
  !======================================================================

  ! ------------------------------------------------------------------------
  subroutine dvr_setup
    ! ------------------------------------------------------------------------ 

    USE globalmod
    USE pfedvrmod
    USE dvrmod

    implicit NONE

    INTEGER(I4B) :: idim, ireg, k, m
    !--------------------------------------------------------------------------
    ! Set up dvr points, weights and kinetic energy matrix  for each region
    !--------------------------------------------------------------------------

    DO idim=1, ndim
       do ireg=1, num_reg(idim)

          CALL DVR_reg( num_fun(idim,ireg),bounds(idim,ireg,:),mat_reg(idim,ireg)%pt, &
               mat_reg(idim,ireg)%wt,mat_reg(idim,ireg)%df,mat_reg(idim,ireg)%ddf, &
               mat_reg(idim,ireg)%ke_mat )

          !       write(*,*)'idim=',idim,' i=',i,'  Seting up Kinetic-energy matrice......'

       end do
    END DO

    !--------------------------------------------------------------------------------
    !------------------following Nicolai's convention to get Kinetic Matrice---------
    !--------------------------------------------------------------------------------

    DO idim=1, ndim
       DO ireg=1, num_reg(idim)
          do k=1,num_fun(idim,ireg)  

             if(mat_reg(idim,ireg)%wt(k).le.0.d0)then
                write(*,*)'ireg=',ireg,'  k=',k,' mat_reg(idim,ireg)%wt(k)=',mat_reg(idim,ireg)%wt(k)
                pause
             endif

             if((ireg.gt.1).and.(k.eq.1)) then
                mat_reg(idim,ireg)%fac1(k)=1.d0/dsqrt(mat_reg(idim,ireg)%wt(1)+mat_reg(idim,ireg-1)%wt(num_fun(idim,ireg-1)))            
             elseif((ireg.lt.num_reg(idim)).and.(k.eq.num_fun(idim,ireg))) then
                mat_reg(idim,ireg)%fac1(k)=1.d0/dsqrt(mat_reg(idim,ireg)%wt(num_fun(idim,ireg))+mat_reg(idim,ireg+1)%wt(1))
             else
                mat_reg(idim,ireg)%fac1(k)=1.d0/dsqrt(mat_reg(idim,ireg)%wt(k))
             end if
             do m=1,num_fun(idim,ireg)
                if((ireg.gt.1).and.(m.eq.1)) then
                   mat_reg(idim,ireg)%fac2(m)=1.d0/dsqrt(mat_reg(idim,ireg)%wt(1)+mat_reg(idim,ireg-1)%wt(num_fun(idim,ireg-1)))            
                elseif((ireg.lt.num_reg(idim)).and.(m.eq.num_fun(idim,ireg))) then
                   mat_reg(idim,ireg)%fac2(m)=1.d0/dsqrt(mat_reg(idim,ireg)%wt(num_fun(idim,ireg))+mat_reg(idim,ireg+1)%wt(1))
                else
                   mat_reg(idim,ireg)%fac2(m)=1.d0/dsqrt(mat_reg(idim,ireg)%wt(m))
                end if

                mat_reg(idim,ireg)%ke_mat(k,m)=mat_reg(idim,ireg)%ke_mat(k,m)*mat_reg(idim,ireg)%fac1(k)*mat_reg(idim,ireg)%fac2(m) 

                !            write(*,*)'idim=',idim,' ireg=',ireg,'  k=',k,'  m=',m

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
  subroutine setup_grid_mapping
    ! ----------------------------------------------------------------------
    !  Construct full spatial grid from the individual elements for each dimension

    USE globalmod
    USE pfedvrmod
    Implicit NONE
    
    INTEGER(I4B) :: idim, ireg, k, m
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

    DO idim=1, ndim

       start=0
       do ireg=1,num_reg(idim)
          if(ireg.eq.1 .and. r0BC) then
             start=-1
             start2=2
             end=num_fun(idim,ireg)-1
          else
             start2=1
             end=num_fun(idim,ireg)-1
          end if

          !****************************************************

          do k=start2,num_fun(idim,ireg)
             m=start+k
             nindex(idim,ireg,k)=m
             grid(idim,m)=mat_reg(idim,ireg)%pt(k) 
             point_order(idim,m)=k
             factor(idim,m)= mat_reg(idim,ireg)%fac1(k)
          end do
          start=start+end
          !****************************************************

       end do

    END DO

  end subroutine setup_grid_mapping

  !======================================================================
  !======================================================================
  !    Diagonalizing real symmetric matrix A(dim,dim) to get eigevec & eigeval
  !======================================================================
  !    By Suxing Hu at LANL, references from Barry and Nicolai's DVR codes
  !
  !    Written on:  5/31/2005
  !
  !     All Subroutines and Functions from LAPACK have been adopted and 
  !     transformed into F90-format, in this code!
  !======================================================================
  ! ----------------------------------------------------------------------
  subroutine diagonalize(dim,A,eig_vec,eig_val)
    ! ----------------------------------------------------------------------
    !  Diagonalize real symmetric matrix A of dimension dim
    !  The eigenvectors and eigenvalues are returned in
    !  eig_vec and eig_val, repsectively
    ! ----------------------------------------------------------------------

    implicit none

    integer, INTENT(IN) :: dim
    real*8, dimension(1:dim,1:dim), INTENT(IN) :: A
    real*8, dimension(1:dim,1:dim), INTENT(OUT) :: eig_vec
    real*8, dimension(1:dim), INTENT(OUT) :: eig_val

    integer :: info
    integer :: lda 
    integer :: lwork
    real*8, dimension(1:dim) :: w
    real*8, dimension(:), ALLOCATABLE :: work   

    lda=dim
    lwork=256*dim
    ALLOCATE(work(1:lwork))

    !----------Diagonalize to find eigenstates and eigenvectors--------------

    call DSYEV('V','U',dim,A,LDA,W,WORK,LWORK,INFO)

    !-------------DVR coefficeients-------------------------------

    eig_vec(:,:)=A(:,:)

    !--------------------Eigenvalues------------------------------

    IF(INFO==0)THEN
       eig_val(:)=w(:)
    ELSE
       write(*,*)'Diagonalizing fails: INFO=',INFO
       STOP
    ENDIF


    DEALLOCATE(work)


  end subroutine diagonalize

end module tdcc_itsubroutines
