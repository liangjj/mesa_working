!======================================================================
  SUBROUTINE test_coulomc 
!======================================================================
!                                                                  
!   This program sets up the Coulomb matrix and determines
!   all eigenvalues and eigenvectors of the asymmetric matrix
!   using the LAPACK90 routine, LA_GEGV.
!	                                                         
!   SUBROUTINE contained:
!       getab solveab
!
!   Calling sequence:
!         test_coulomc
!         ------------
!         //        \\
!       getab     solveab
!         |       -------
!       coulom    /   |  \
!                /    |   \    
!            la_gegv sort bvalu2
!                           ||
!                         interv
!
!----------------------------------------------------------------------
!
    USE spline_param

    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(ns,ns) :: a, b

    ! .. builds a and b full matrices for the coulomb problem
    CALL getab

    ! .. solve A x = E B x
    CALL solveab

    CONTAINS

    !==================================================================
      SUBROUTINE getab
    !==================================================================
    !   builds full a and b matrices for the coulomb problem 
    !   with the correct boundary condition at the origin
    !------------------------------------------------------------------
    ! 
        IMPLICIT NONE
        CHARACTER(LEN=1) :: iprint
        INTEGER :: i, j

        ! .. builds full a and b matrices 
        CALL coulom(a, b)
    
        ! .. apply boundary condition at the origin

        a(1,2:ks) = 0.d0
        a(2:ks,1) = 0.d0
        a(1,1)=-1.0d10

        b(1,2:ks) = 0.d0
        b(2:ks,1) = 0.d0
        b(1,1) = 1.0d0

        ! .. prints the matrics a and b

        PRINT *, 'Do you want the matrix printed?'
        PRINT *, '          Y or y  ---  yes'
        PRINT *, '          Others  ---  no'

        READ(5,*) iprint
        if (iprint == 'y' .or. iprint == 'Y') then
          PRINT *, 'A matrix'
          do i = 1,ns
            PRINT '(i5/(6f12.7))', i,(a(i,j),j=1,ns)
          end do

          PRINT *, 'B matrix'
          do i = 1,ns
            PRINT '(i5/(6f12.7))', i,(b(i,j),j=1,ns)
          end do
        end if
        PRINT *
      END SUBROUTINE getab

    !==================================================================
      SUBROUTINE solveab
    !==================================================================
    !   Solve A x = E B x
    !------------------------------------------------------------------
    !
        USE spline_grid
        USE F90_LAPACK, ONLY: LA_GEGV

        IMPLICIT NONE
        REAL(KIND=8), DIMENSION(ns,ns) :: aa, bb
        REAL(KIND=8), DIMENSION(ns) :: eigvalr, eigvali, eigvald
        INTEGER, DIMENSION(ns) :: num
        REAL(KIND=8), DIMENSION(ns,ns) :: vl, eigvec
        INTEGER :: INFO, i, j, m, ii

        REAL(KIND=8) :: x, y, bvalu2

        ! .. backup a and b
        aa = a
        bb = b

        ! .. evaluates the eigenvalues and eigenvectors
        CALL la_gegv(a, b, eigvalr, eigvali, eigvald, vl, eigvec, INFO )

        PRINT *, 'Eigenvalues'
        if (INFO /= 0) then
          PRINT *, 'Error in Eigenvalue routine'
        else
          do i=1,ns
            eigvalr(i) = eigvalr(i)/eigvald(i)
            eigvali(i) = eigvali(i)/eigvald(i)
          end do

          CALL sort(ns,eigvalr,eigvali,num)
          OPEN(UNIT=1, FILE='plot.dat', STATUS='unknown')

          PRINT *, 'Sorted eigenvalues'
          WRITE(6,'(A,9X,A,12X,A,8X,A,8X,A)')    &
                'New','Real part','Im. part','Old','n*'

          do i = 2,ns
            if(eigvalr(i) >= 0) then
              WRITE(6,'(I3,1X,F20.10,F20.10,I7,F20.10)')    &
                i, eigvalr(i),eigvali(i),num(i),1/sqrt(eigvalr(i))
            else
              WRITE(6,'(I3,1X,F20.10,F20.10,I7)')    &
                  i, eigvalr(i),eigvali(i),num(i)
            end if
          end do

          PRINT *, 'Positive eigenvectors'
          do ii = 2,ns
            if (eigvalr(ii) >= 0 .and. eigvali(ii) == 0) then
              WRITE(1,*) ' Eigenvector:   ',ii,eigvalr(ii)
              j = num(ii)
              x = 0
              y = eigvec(1,j)
              WRITE(1,'(F10.4, F12.6)') x,y

              ! .. calculates the wavefunctions at the gauss points
              do i = 1,nv
                do m = 1,ks
                  x = gr(i,m)
                  y = bvalu2(t,eigvec(1,j),ns,ks,x,0)
                  WRITE(1,'(F10.4, F12.6)') x,y
                end do
              end do
              PRINT *
            end if
          end do
        end if
      END SUBROUTINE solveab

  END SUBROUTINE test_coulomc 

!====================================================================
  SUBROUTINE sort(ns,eigvalr,eigvali,num)
!====================================================================
!
!   Sorts eigenvalues in decreasing order.
!   
!------------------------------------------------------------
!
!   on entry
!   --------
!       eigvalr  the real part of eigenvalues.
!       eigvali  the imaginary part of eigenvalues.
!
!   on exit
!   -------
!       eigvalr  the sorted list of real part of eigenvalues.
!       eigvali  the sorted list of imaginary part of eigenvalues.
!       num     the number of the eigenvalue in the unsorted list.
!
!------------------------------------------------------------
!
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(ns), INTENT(INOUT) :: eigvalr, eigvali
    INTEGER, DIMENSION(ns), INTENT(OUT) :: num
    INTEGER, INTENT(IN) :: ns

    ! .. local variables
    REAL(KIND=8) :: tempr, tempi
    INTEGER :: ii, jj, tempn

    do jj =1,ns
      num(jj) = jj
    end do

    do jj = 1,ns-1
      do ii=1,ns-1
        if (eigvalr(ii) > eigvalr(ii+1)) then
	  tempn = num(ii)
	  tempr = eigvalr(ii)
	  tempi = eigvali(ii)
	  num(ii) = num(ii+1)
	  eigvalr(ii) = eigvalr(ii+1)
	  eigvali(ii) = eigvali(ii+1)
	  num(ii+1) = tempn
	  eigvalr(ii+1) = tempr
	  eigvali(ii+1) = tempi
        end if
      end do
    end do

  END SUBROUTINE sort
