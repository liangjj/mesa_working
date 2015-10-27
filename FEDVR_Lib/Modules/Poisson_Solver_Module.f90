!***********************************************************************
! DVR_Solver_Module
!**begin prologue     DVR_Solver_Module
!**date written       090119   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, Iterative, Lanczos, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             
!**purpose            Read in FEDVR Data
!***references
!***modules needed    See USE statements below
!***comments          
!***                  
!***                  
!***                  
!***                  
!***end prologue      DVR_Solver_Module
!***********************************************************************
!***********************************************************************
                           MODULE DVR_Solver_Module
                           USE Data_Module
                           USE FEDVR_Shared
                           USE FEDVR_Derived_Types
!***********************************************************************
!***********************************************************************
                              CONTAINS
!***********************************************************************
!***********************************************************************
!deck Poisson.f
!***begin prologue     Poisson
!***date written       030131   (yymmdd)
!***revision date               (yymmdd)
!***keywords           two electron radial integrals
!***author             schneider, b. i.(nsf)
!***source
!***purpose            solve poisson equation using dvr functions

!***references

!***routines called    iosys, util and mdutil
!***end prologue       Poisson
  SUBROUTINE poisson(ke,ang_pot,f,pt,wt,nphy)
  USE dvr_global
  IMPLICIT NONE
  INTEGER                                :: nphy
  REAL*8,   DIMENSION(nphy)              :: pt, wt, ang_pot
  REAL*8,   DIMENSION(nphy,nphy)         :: ke, f
  REAL*8,   DIMENSION(:), ALLOCATABLE    :: rho, v
  REAL*8,   DIMENSION(:,:), ALLOCATABLE  :: t
  REAL*8,   DIMENSION(:,:), ALLOCATABLE  :: work
  REAL*8                                 :: scale, ptlst, a 
  INTEGER,  DIMENSION(:), ALLOCATABLE    :: ipvt
  INTEGER                                :: l, l_barrier, i, info
  INTEGER                                :: l_fac
  CHARACTER(LEN=80)                      :: title
!
  scale = one / dscale
  l_fac=2*val+1
  count = 0
  DO i = 1, dvr_mat(0)%n_final
     DO j = 1, i
        count = count + 1
        dvr_mat(0)%lower(count) = dvr_mat(val)%ham(i,j)
     END DO
  END DO
  dvr_mat(0)%lower(:) = scale * dvr_mat(0)%lower(:)
  dvr_mat(0)%number_of_right_hand_sides = int_one
  IF (dentyp == 'exponential') THEN
      dvr_mat(0)%rhs(:,1) = exp( -grid%grid_points(:))
  ELSE IF (dentyp == 'harmonic') then
      dvr_mat(0)%rhs(:,1) = exp( -grid%grid_points(:) * grid_grid_points(:))
  ELSE IF (dentyp == 'one') then
      dvr_mat(0)%rhs(:,1) = one
  ELSE IF (dentyp == 'linear') then
      dvr_mat(0)%rhs(:,1) = grid%grid_points(:)
  END IF
  dvr_mat(0)%rhs(:,1) = - l_fac * dvr_mat(0)%rhs(:,1)
  IF(prn(9)) then
     title='density'
     call prntfmn(title,dvr_mat(0)%rhs(:,1),physical_points,1,physical_points,1,iout,'e')
  END IF
!
!    solve linear equations
!
  Call DSPSV('u',dvr_mat(0)%n_final,dvr_mat(0)%number_of_right_hand_sides,     &
              dvr_mat(0)%lower,dvr_mat(0)%ipvt,dvr_mat(0)%rhs,                 &
              dvr_mat(0)%n_final,info)
1 Format(/1x,'Solve Linear Systems Angular Momentum = ',i3)
!
     l_fac=2*l+1
     work(:,1)=pt**(l+1)
     ptlst=1.d0/(pt(nphy)**l_fac)
     a=0.d0
     do i=1,nphy
        a = a + rho(i)*wt(i)*work(i,1)
     end do
     a=a/ptlst
     v(nphy)= a*work(nphy,1)
     do i=1,nphy-1
        v(i) = rho(i)*f(i,i) + a*work(i,1)
     end do
     v = v / pt 
     title='potential'
     call prntrm(title,v,nphy,1,nphy,1,iout)
  end do
  DEALLOCATE(t,work,ipvt)
  DEALLOCATE(rho,v)
END SUBROUTINE poisson
!***********************************************************************
!***********************************************************************
!deck Right_Hand_Side
!***begin prologue     Right_Hand_Side
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Form the right hand sides of the Poisson equation
!***                   for some special test cases.
!***                   
!***                   
!***                   
!***references
!***routines called    iosys, util and mdutil
!***end prologue       
!
  SUBROUTINE Right_Hand_Side(grid,dvr_mat)
  IMPLICIT NONE
  TYPE(coordinates)                   :: grid  
  TYPE (dvr_matrices), DIMENSION(0:)  :: dvr_mat
  INTEGER                             :: i
  IF (dvr_mat(0)%type_inhomo == 'one') THEN
     DO i = 1, dvr_mat(0)%number_of_right_hand_sides
        dvr_mat(0)%rhs(:,i) = one * sqrt(grid%grid_weights(:))
     END DO
  ELSE IF(dvr_mat(0)%type_inhomo == 'linear') THEN
     DO i = 1, dvr_mat(0)%number_of_right_hand_sides
        dvr_mat(0)%rhs(:,i) = sqrt(grid%grid_weights(:)) * grid%grid_points(:)
     END DO
  ELSE IF(dvr_mat(0)%type_inhomo == 'exponential') THEN
     DO i = 1, dvr_mat(0)%number_of_right_hand_sides
        dvr_mat(0)%rhs(:,i) = exp(-grid%grid_points(:)) * sqrt(grid%grid_weights(:)) 
     END DO
  ELSE IF(dvr_mat(0)%type_inhomo == 'coulomb') THEN
     DO i = 1, dvr_mat(0)%number_of_right_hand_sides
        dvr_mat(0)%rhs(:,i) = one * sqrt(grid%grid_weights(:))  / grid%grid_points(:i) 
     END DO
  END IF
END SUBROUTINE Right_Hand_Side
!***********************************************************************
!***********************************************************************
!deck Solver.f
!***begin prologue     Solver
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            Solve the Poisson equation using packed solvers.
!***
!***references
!***routines called
!***end prologue       Solver
  SUBROUTINE Solver(grid, dvr_mat,val)
  IMPLICIT NONE
  TYPE(coordinates)                        :: grid  
  TYPE (dvr_matrices), DIMENSION(0:)       :: dvr_mat
  REAL(idp)                                :: scale  
  INTEGER                                  :: info
  INTEGER                                  :: count
  INTEGER                                  :: i
  INTEGER                                  :: j
  INTEGER                                  :: val
  CHARACTER(LEN=3)                         :: itoc
  CHARACTER(LEN=80)                        :: title
  scale = one / dscale
  count = 0
  DO i = 1 , physical_points
     DO j = 1, i
        count = count + 1
        dvr_mat(0)%lower(count) = scale * dvr_mat(val)%ham(i,j)
     END DO
  END DO
  dvr_mat(0)%computed_solution(:,:) = dvr_mat(0)%rhs(:,:)
  Call DSPSV('u',physical_points,dvr_mat(0)%number_of_right_hand_sides,     &
             dvr_mat(0)%lower,dvr_mat(0)%ipvt,dvr_mat(0)%computed_solution, &
             physical_points,info)

1 Format(/1x,'Solve Linear Systems Angular Momentum = ',i3)
END SUBROUTINE Solver
!***********************************************************************
!***********************************************************************
!deck Compare.f
!***begin prologue     Compare
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            Compare computed and exact solutions to the Poisson
!***                   equation for some specific cases.
!***
!***references
!***routines called
!***end prologue       Compare
  SUBROUTINE Compare(grid, dvr_mat,val)
  IMPLICIT NONE
  TYPE(coordinates)                        :: grid  
  TYPE (dvr_matrices), DIMENSION(0:)       :: dvr_mat
  REAL(idp)                                :: scale  
  INTEGER                                  :: info
  INTEGER                                  :: count
  INTEGER                                  :: i
  INTEGER                                  :: j
  INTEGER                                  :: val
  INTEGER                                  :: len
  INTEGER                                  :: lenth
  CHARACTER(LEN=3)                         :: itoc
  CHARACTER(LEN=80)                        :: title
  len=lenth(grid%label)
  DO i = 1, dvr_mat(0)%number_of_right_hand_sides
     dvr_mat(0)%rhs(:,i) = dvr_mat(0)%computed_solution(:,i) / sqrt( grid%grid_weights(:) )
  END DO
  IF (keyword == 'cartesian') THEN
!
!           Some tests for cases on (0,1) whyere the solution vanishes at the two ends.
!
      IF (dvr_mat(0)%type_inhomo == 'one') THEN
!
!         Solution is x *  ( x - 1 ) / 2
!
          DO i = 1, dvr_mat(0)%number_of_right_hand_sides
             dvr_mat(0)%exact_solution(:,i) =  grid%grid_points(:) *                  &
                                             ( grid%grid_points(:) - one ) * half
          END DO
!
      ELSE IF(dvr_mat(0)%type_inhomo == 'linear') THEN
!
!         Solution is x *  ( x*x - 1 ) / 6
!
          DO i = 1, dvr_mat(0)%number_of_right_hand_sides
             dvr_mat(0)%exact_solution(:,i) =  grid%grid_points(:) *                  &
                                             ( grid%grid_points(:) *                  &
                                               grid%grid_points(:) - one ) * sixth
          END DO
      ELSE IF(dvr_mat(0)%type_inhomo == 'exponential') THEN
!
!         Solution is x - 1 + exp(-x) - x * exp(-1)
!
          DO i = 1, dvr_mat(0)%number_of_right_hand_sides
             dvr_mat(0)%exact_solution(:,i) =  grid%grid_points(:) - one              &
                                                       +                              &
                                               exp(-grid%grid_points(:))              &
                                                       -                              &
                                               grid%grid_points(:) * exp (- one )                    
          END DO
      ELSE IF(dvr_mat(0)%type_inhomo == 'coulomb') THEN      
!
!         Solution is x * ( ln (x) - ln(1) )
!
          DO i = 1, dvr_mat(0)%number_of_right_hand_sides
             dvr_mat(0)%exact_solution(:,i) =  grid%grid_points(:) * (                 &
                                               log ( grid%grid_points(:) ) - log(one) )
          END DO
      END IF
!
  ELSE IF(keyword =='spherical') THEN
!
!         Some tests for cases where solution is zero at r = R
!
      IF(grid%label(1:len) == 'r') THEN
!
         IF (dvr_mat(0)%type_inhomo == 'one') THEN
!
!            Solution is (r * r - R * R) / 6
!
            DO i = 1, dvr_mat(0)%number_of_right_hand_sides
               dvr_mat(0)%exact_solution(:,i) =  ( grid%grid_points(:) *                &
                                                   grid%grid_points(:)                  &
                                                 - R_max * R_max )* sixth
         
            END DO
!
         ELSE IF(dvr_mat(0)%type_inhomo == 'linear') THEN
!
!            Solution is (r * r * r - R * R * R) / 12
!
            DO i = 1, dvr_mat(0)%number_of_right_hand_sides
               dvr_mat(0)%exact_solution(:,i) =  ( grid%grid_points(:) *                &
                                                   grid%grid_points(:) *                &
                                                   grid%grid_points(:)                  &
                                                 - R_max * R_max *R_max )* sixth * half
             END DO
         ELSE IF(dvr_mat(0)%type_inhomo == 'coulomb') THEN      
!
!            Solution is (r - R) / 2
!
            DO i = 1, dvr_mat(0)%number_of_right_hand_sides
               dvr_mat(0)%exact_solution(:,i) =  ( grid%grid_points(:) - R_max )* half
            END DO
         ELSE
!
            Call lnkerr('not valid')
!
         END IF
      ELSE
!
         Call lnkerr('not valid')
!
      END IF
!
  ELSE
          Call lnkerr('not valid')
  END IF
  title='Numerical Solution'
  Call Prntfm(title,dvr_mat(0)%rhs,physical_points,                                              &
                    dvr_mat(0)%number_of_right_hand_sides, physical_points,                      &
                    dvr_mat(0)%number_of_right_hand_sides,iout,'e')
  title='Exact Solution'
  Call Prntfm(title,dvr_mat(0)%exact_solution,physical_points,                                   &
                    dvr_mat(0)%number_of_right_hand_sides, physical_points,                      &
                    dvr_mat(0)%number_of_right_hand_sides,iout,'e')
END SUBROUTINE Compare
!***********************************************************************
!***********************************************************************
           END MODULE DVR_Solver_Module
!***********************************************************************
!***********************************************************************
