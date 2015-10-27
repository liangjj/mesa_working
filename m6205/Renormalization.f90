!***********************************************************************
!Renormalization
!**begin prologue     Renormalization
!**date written       090119   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, Iterative, Lanczos, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             Time_Propagation
!**purpose            Compute renormalization factors required for the
!***                  renormalization of the sector polynomials and kinetic
!***                  matrix elements in a FEDVR basis.
!***references
!***modules needed    See USE statements below
!***comments          
!***                  
!***                  
!***end prologue      Renormalization
!***********************************************************************
!***********************************************************************
                           MODULE Renormalization
                           USE Data
                           USE Grid_Defined_Types
!***********************************************************************
!***********************************************************************
!                          Explicit Interfaces
!***********************************************************************
!
!
                            INTERFACE Normalized_Matrices                       
                       MODULE PROCEDURE Normalized_Even_Matrices,   &
                                        Normalized_Odd_Matrices  
                            END INTERFACE Normalized_Matrices
!
!***********************************************************************
!***********************************************************************
                              CONTAINS
!***********************************************************************
!***********************************************************************
!deck Normalization.f
!***begin prologue     Normalization
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           lobatto functions
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Normalization Factors
!***references
!***routines called    iosys, util and mdutil
!***end prologue       Normalization

  SUBROUTINE Normalization(reg,n_reg)
  IMPLICIT NONE
  TYPE(regional),  DIMENSION(:)       :: reg
  INTEGER                             :: n_reg
  INTEGER                             :: i
  CHARACTER (LEN=3)                   :: itoc
!           
  IF ( n_reg == 1) THEN
!
!                  Only one region.  No endpoint corrections required.
!
       i = 1
       Call Norm ( reg(i)%wt, reg(i)%normalization, reg(i)%n_pts)

  ELSE
!
!                  First region of multi-region domain.  Correction  at  
!                  right endpoint needed from first function in second region.
!
       i = 1
       Call Norm ( reg(i)%wt, reg(i)%normalization, reg(i)%n_pts,               &      
                   wt_right_end=reg(i+1)%wt(1))
!       
       DO i = 2, n_reg - 1
!
!                  General case.  Correction at both the left and right
!                  endpoints needed.  The correction at the left enpoint
!                  requires the last weight from the previous region while
!                  the right endpoint correction requires the first weight
!                  from the next region.
!
          Call Norm ( reg(i)%wt, reg(i)%normalization, reg(i)%n_pts,            &      
                      wt_left_end=reg(i-1)%wt(reg(i-1)%n_pts),                  &
                      wt_right_end=reg(i+1)%wt(1))                     
!
       END DO
!
!                  Last region.  Correct the left end point using the last
!                  weight from the previous region.
!
!
       i = n_reg
       Call Norm ( reg(i)%wt, reg(i)%normalization, reg(i)%n_pts,            &
                   wt_left_end=reg(i-1)%wt(reg(i-1)%n_pts) )
  END IF
  IF (prnt(3) == .true.) THEN
      DO i = 1, n_reg
         write(iout,*)
         call Print_Matrix(type_real_vector,reg(i)%normalization,title='Normalization region-'//itoc(i))
      END DO
  END IF
END SUBROUTINE Normalization
!**********************************************************************
!**********************************************************************
!deck Norm.f
!***begin prologue     Norm
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           lobatto functions
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Modify the weights at the ends of the interval to
!***                   reflect that there are bridge functions present. 
!***                   
!***

!***references

!***routines called    iosys, util and mdutil
!***end prologue       Norm

  SUBROUTINE Norm(wt,normalization,n,wt_left_end,wt_right_end)
  IMPLICIT NONE
  REAL(idp), DIMENSION (:)             :: wt
  REAL(idp), DIMENSION (:)             :: normalization
  INTEGER                              :: n
  REAL(idp), OPTIONAL                  :: wt_right_end
  REAL(idp), OPTIONAL                  :: wt_left_end
!
!
  IF ( .not.present(wt_left_end).and..not.present(wt_right_end) ) THEN
!
       normalization(:) = Sqrt ( 1.d0 / wt(:) )
!
  ELSE IF ( .not.present(wt_left_end).and.present(wt_right_end) ) THEN
!
!      Modify the last weight and then get the inverse square roots.
!
       normalization(1:n-1) = 1.d0 / sqrt ( wt(:) )     
       normalization(n) = 1.d0 / sqrt ( wt(n) + wt_right_end )     
!
  ELSE IF ( present(wt_left_end).and..not.present(wt_right_end) ) THEN
!
!      
       normalization(1) = 1.d0 / sqrt ( wt_left_end + wt(1) )
       normalization(2:n) = 1.d0 / sqrt ( wt(2:n) )
!
!
  ELSE IF ( present(wt_left_end).and.present(wt_right_end) ) THEN
!      
!      Modify the last weight and then get the inverse square roots.
!
       normalization(1) = 1.d0 / Sqrt ( wt_left_end + wt(1) )
       normalization(n)  = 1.d0 / Sqrt ( wt_right_end + wt(n) )
       normalization(2:n-1)  =  1.d0 / sqrt ( wt(2:n-1) )
  END IF
END SUBROUTINE Norm
!***********************************************************************
!***********************************************************************
!deck Normalized_Functions
!***begin prologue     Normalized_Functions                                                                          
!***date written       960718   (yymmdd                                                                                  
!***revision date               (yymmdd)                                                                                 
!***keywords                                                                                                             
!***author             schneider, b. i.(nsf)                                                                             
!***source                                                                                                               
!***purpose            Take the basic interpolation polynomials in each sector                                           
!***                   and normalize them in their own sector.  This does not properly
!***                   normalize the bridge functions.                                                                              
!***references                                                                                                           
!***routines called    iosys, util and mdutil
!***end prologue
  SUBROUTINE Normalized_Functions(reg,n_reg)
  IMPLICIT NONE
  TYPE(regional), DIMENSION(:)   :: reg
  INTEGER                        :: n_reg
  INTEGER                        :: i
  INTEGER                        :: j
  CHARACTER (LEN=3)              :: itoc
!                                                                                                               
!                                                                                                                
!                  Normalize the functions in their own region.  
  DO i = 1, n_reg
     ALLOCATE( reg(i)%bf(1:reg(i)%n_pts,1:reg(i)%n_pts), reg(i)%dbf(1:reg(i)%n_pts,1:reg(i)%n_pts),    &
               reg(i)%ddbf(1:reg(i)%n_pts,1:reg(i)%n_pts)) 
     DO j = 1, reg(i)%n_pts
        reg(i)%bf(:,j)   = reg(i)%p(:,j) / sqrt( reg(i)%wt(j) )
        reg(i)%dbf(:,j)  = reg(i)%dp(:,j) / sqrt( reg(i)%wt(j) )
        reg(i)%ddbf(:,j) = reg(i)%ddp(:,j) / sqrt( reg(i)%wt(j) )
     END DO
  END DO
!
!                                                                       
  IF (prnt(3) == .true. ) THEN
      Write(iout,1)
      DO i = 1, n_reg
         write(iout,*)
         call Print_Matrix(type_real_matrix,reg(i)%bf,reg(i)%n_pts,reg(i)%n_pts,                      &
                           title='Normalized Polynomials Region-'//itoc(i))
         call Print_Matrix(type_real_matrix,reg(i)%dbf,reg(i)%n_pts,reg(i)%n_pts,                     &
                           title='First Derivative of Normalized Polynomials Region-'//itoc(i))
         call Print_Matrix(type_real_matrix,reg(i)%ddbf,reg(i)%n_pts,reg(i)%n_pts,                    &
                           title='Second Derivative of Normalized Polynomials Region-'//itoc(i))
      END DO
  END IF
1    FORMAT(/,1x,'Normalized DVR Functions and Derivatives')  
END SUBROUTINE Normalized_Functions
!***********************************************************************
!***********************************************************************
!deck Normalized_Even_Matrices
!***begin prologue     Normalized_Even_Matrices
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           lobatto functions
!***author             schneider, b. i.(nsf)
!***source
!***purpose            normalize the kinetic energy matrix elements so that they
!***                   are correct across the joining regions.

!***references

!***routines called    iosys, util and mdutil
!***end prologue       Normalized_Even_Matrices
  SUBROUTINE Normalized_Even_Matrices(even,reg,n_reg)
  IMPLICIT NONE
  TYPE(REGIONAL), DIMENSION(:)         :: reg
  TYPE(EVEN_MATRICES)                  :: even
  INTEGER                              :: n_reg
  INTEGER                              :: i
  CHARACTER(LEN=3)                     :: itoc
!
!
  DO i = 1, n_reg
     ALLOCATE ( reg(i)%even%ham(1:reg(i)%n_pts,1:reg(i)%n_pts ) )
  END DO
!
  IF ( n_reg == 1 ) THEN
!
!      One Region only.  No change needed except normalization.
!
       i = 1
!
!
       Call Normalized_Kinetic_Energy (reg(i)%even%ham,              &
                   reg(i)%even%tr,                                   &
                   reg(i)%normalization,                             &
                   reg(i)%n_pts)
  ELSE
       i = 1
!
!      Region one is a special case.  Only the last right hand element needs to be changed
!      and the elements renormalized.
!
!
       Call Normalized_Kinetic_Energy (reg(i)%even%ham,              &
                   reg(i)%even%tr,                                   &
                   reg(i)%normalization,                             &
                   reg(i)%n_pts,                                     &
                   right_end_mat_el=reg(i+1)%even%tr(1,1) )
!
!      Now do the general case.  Both the left and right elemenents are changed.
! 
       DO i = 2 , n_reg - 1
!
!
          Call Normalized_Kinetic_Energy (reg(i)%even%ham,           &
                      reg(i)%even%tr,                                &
                      reg(i)%normalization,                          &
                      reg(i)%n_pts,                                  &
                      left_end_mat_el=reg(i-1)%even%tr(reg(i-1)%n_pts,reg(i-1)%n_pts), &
                      right_end_mat_el=reg(i+1)%even%tr(1,1) )
!
       END DO
!
       i = n_reg
!
!      Special case of last region.  Only the left hand elements needs to be changed.
!
!
!
          Call Normalized_Kinetic_Energy (reg(i)%even%ham,           &
                      reg(i)%even%tr,                                &
                      reg(i)%normalization,                          &
                      reg(i)%n_pts,                                  &
                      left_end_mat_el=reg(i-1)%even%tr(reg(i-1)%n_pts,reg(i-1)%n_pts))
!
  END IF
!
  DO i = 1, n_reg
     DEALLOCATE ( reg(i)%even%tr )
  END DO
!
  IF (prnt(4) == .true. ) THEN
      DO i = 1 , n_reg
         write(iout,*)
         call Print_Matrix(type_real_matrix,reg(i)%even%ham,reg(i)%n_pts,reg(i)%n_pts,                           &
                           title='Normalized Even Sector_Matrix Region-'//itoc(i))
      END DO
  END IF
END SUBROUTINE Normalized_Even_Matrices
!***********************************************************************
!***********************************************************************
!deck Normalized_Odd_Matrices
!***begin prologue     Normalized_Odd_Matrices
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           lobatto functions
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Renormalization
!***                   
!***

!***references

!***routines called    iosys, util and mdutil
!***end prologue       Odd_Renormalization

  SUBROUTINE Normalized_Odd_Matrices(odd,reg,n_reg)
  IMPLICIT NONE
  TYPE(odd_matrices)                   :: odd
  TYPE(regional),  DIMENSION(:)        :: reg
  INTEGER                              :: n_reg
  INTEGER                              :: i
  CHARACTER(LEN=3)                     :: itoc
!
!
  DO i = 1, n_reg
     ALLOCATE ( reg(i)%odd%ham(1:reg(i)%n_pts,1:reg(i)%n_pts ) )
  END DO
!
  IF ( n_reg == 1 ) THEN
!
!      One Region only.  No change needed except normalization.
!
       i = 1
!
       Call Normalized_Kinetic_Energy (reg(i)%odd%ham,              &
                   reg(i)%odd%tr,                                   &
                   reg(i)%normalization,                            &
                   reg(i)%n_pts)
  ELSE
       i = 1
!
!      Region one is a special case.  Only the last right hand element needs to be changed
!      and the elements renormalized.
!
       Call Normalized_Kinetic_Energy (reg(i)%odd%ham,              &
                   reg(i)%odd%tr,                                   &
                   reg(i)%normalization,                            &
                   reg(i)%n_pts,                                    &
                   right_end_mat_el=reg(i+1)%odd%tr(1,1) )
!
!      Now do the general case.  Both the left and right elemenents are changed.
! 
       DO i = 2 , n_reg - 1
!
          Call Normalized_Kinetic_Energy (reg(i)%odd%ham,           &
                      reg(i)%odd%tr,                                &
                      reg(i)%normalization,                         &
                      reg(i)%n_pts,                                 &
                      left_end_mat_el=reg(i-1)%odd%tr(reg(i-1)%n_pts,reg(i-1)%n_pts), &
                      right_end_mat_el=reg(i+1)%odd%tr(1,1) )
       END DO
       i = n_reg
!
!      Special case of last region.  Only the left hand elements needs to be changed.
!
!
          Call Normalized_Kinetic_Energy (reg(i)%odd%ham,           &
                      reg(i)%odd%tr,                                &
                      reg(i)%normalization,                         &
                      reg(i)%n_pts,                                 &
                      left_end_mat_el=reg(i-1)%odd%tr(reg(i-1)%n_pts,reg(i-1)%n_pts) )
  END IF
  IF (prnt(4) == .true. ) THEN
      DO i = 1 , n_reg
         write(iout,*)
         call Print_Matrix(type_real_matrix,reg(i)%odd%ham,reg(i)%n_pts,reg(i)%n_pts,                           &
                           title='Normalized Odd Sector Matrix Region-'//itoc(i))
      END DO
  END IF
!
  DO i = 1, n_reg
     DEALLOCATE(reg(i)%odd%tr)
  END DO
END SUBROUTINE Normalized_Odd_Matrices
!***********************************************************************
!***********************************************************************
!deck Normalized_Kinetic_Energy.f
!***begin prologue     Normalized_Kinetic_Energy
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           lobatto functions
!***author             schneider, b. i.(nsf)
!***source
!***purpose            To do the renormalization of the matrix elements
!***

!***references

!***routines called    iosys, util and mdutil
!***end prologue       Normalized_Kinetic_Energy

  SUBROUTINE Normalized_Kinetic_Energy( ham, tr, norm, n , left_end_mat_el, right_end_mat_el )
  IMPLICIT NONE
  REAL(idp),  DIMENSION(:,:)      :: ham
  REAL(idp),  DIMENSION(:,:)      :: tr
  REAL(idp),  DIMENSION (:)       :: norm
  REAL(idp),  OPTIONAL            :: left_end_mat_el
  REAL(idp),  OPTIONAL            :: right_end_mat_el
  INTEGER                         :: n
  INTEGER                         :: i
  INTEGER                         :: j
!
!
!
  IF ( .not.present(left_end_mat_el).and..not.present(right_end_mat_el) ) THEN
!
!      No bridge functions, so its straightforward.
!
       DO i = 1, n
          ham(i,1:i) = norm(i) * tr(i,1:i) * norm(1:i)
          ham(1:i,i) = ham(i,1:i)
       END DO
!
!
  ELSE IF ( .not.present(left_end_mat_el).and.present(right_end_mat_el) ) THEN
!
!      There is a bridge function at the right hand boundary that requires the
!      the matrix element involving the first function in the next sector.  This
!      is in the right_end_mat_el.
!      
       DO i = 1, n - 1
          ham(i,1:i) = norm(i) * tr(i,1:i) * norm(1:i)
          ham(1:i,i) = ham(i,1:i)
       END DO
       ham(n,1:n-1) = norm(n) * tr(n,1:n-1) * norm(1:n-1)
       ham(1:n-1,n) = ham(n,1:n-1)
       ham(n,n) = norm(n) * ( tr(n,n) + right_end_mat_el ) * norm(n)
!         
  ELSE IF ( present(left_end_mat_el).and..not.present(right_end_mat_el) ) THEN
!
!      This is the last sector.  Here we require information involving the bridge function
!      at the left hand boundary.
!
       ham(1,1) = norm(1) * ( tr(1,1) + left_end_mat_el ) * norm(1)       
       DO i = 2, n
          ham(i,1:i) = norm(i) * tr(i,1:i) * norm(1:i)
          ham(1:i,i) = ham(i,1:i)
       END DO
!
  ELSE IF ( present(left_end_mat_el).and.present(right_end_mat_el) ) THEN
!
!      This is the general case when there are contributions at both ends of the sector.
!      at the left hand boundary.
!       
       ham(1,1) = norm(1) * ( tr(1,1) + left_end_mat_el ) * norm(1)       
       DO i = 2, n - 1
          ham(i,1:i) = norm(i) * tr(i,1:i) * norm(1:i)
          ham(1:i,i) = ham(i,1:i)
       END DO
       ham(n,1:n-1) = norm(n) * tr(n,1:n-1) * norm(1:n-1)
       ham(1:n-1,n) = ham(n,1:n-1)
       ham(n,n) = norm(n) * ( tr(n,n) + right_end_mat_el ) * norm(n)
  END IF
END SUBROUTINE Normalized_Kinetic_Energy
!***********************************************************************
!***********************************************************************
           END MODULE Renormalization
!***********************************************************************
!***********************************************************************
