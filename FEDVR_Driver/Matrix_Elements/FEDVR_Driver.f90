!deck FEDVR_DRIVER
!**begin prologue     FEDVR_DRIVER
!**date written       060711   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           
!**author             schneider, barry (nsf)
!**source
!**purpose            test FEDVR_DRIVER code
!**description        calls routines to input data, construct dvr matrices
!**                   
!**references
!**routines called        
!                       ----                    --------
!
!**end prologue       FEDVR_DRIVER
  PROGRAM FEDVR_DRIVER
  USE DVR_Module
  USE Read_DVR_Module
  USE Matrix_Scale_and_Assemble
  USE Matrix_Diagonalization
  USE Poisson_Module
  IMPLICIT NONE
  INTEGER                          :: input
  INTEGER                          :: output
  INTEGER                          :: i
  INTEGER                          :: sets
  INTEGER                          :: intkey
  INTEGER                          :: no_sets
  LOGICAL                          :: test_key
  LOGICAL                          :: dollar

  COMMON /io/ input, output
!
!  Get the input and output file numbers which appear in the 
!  Module input_output.f90 and put them into input and output which appears
!  in the ONLY common block in the code.  This is needed to pass into library
!  routines and for no other purpose.
!
  input = inp
  output = iout
!
!  Open the input and output files
!
  OPEN(input,file='FEDVR.inp',status='old')
  OPEN(output,file='FEDVR.out',status='unknown')
!
  write(iout,1)
  write(iout,2)
  write(iout,1)
  IF ( dollar('$begin_calculation',card,cpass,inp) )THEN
     no_sets = intkey(card,'Number_of_Data_Sets',1,' ') ! number of data sets
     Write(iout,*) '** Number of Data Sets ** = ',no_sets
  ELSE
     Call lnkerr('No initialization of data')
  END IF
Loop_over_Data_Sets :  &
  DO sets = 1, no_sets
     Write(iout,*) '                     ** Processing Data Set ** = ',sets
!
     Call Set_General_Keywords(sets)  ! here we just read in some general keyword data
                                      ! see module Read_DVR_Module
!
     Loop_Over_Coordinates : &
     DO i = 1, spdim
        write(iout,3) keyword
        len=lenth(reg_grid(i)%label)
        write(iout,4) reg_grid(i)%label(1:len)
!-------------------------------------------------------------------------------------------------------------------------
        Call Read_Data(reg_grid(i)) !  Read input data for 
                                    !  each type of coordinate
                                    !  or DVR. Contained in module Read_DVR_Module
!-------------------------------------------------------------------------------------------------------------------------
        Call Read_Potential_Data    ! read data for one-body 
                                    ! potential. Contained in module Read_DVR_Module
!-------------------------------------------------------------------------------------------------------------------------
        Call Lobatto_Functions (reg_grid(i))  ! see module DVR_Polynomials_Module, calls
                                              ! either a lobatto or fourier routine
!-------------------------------------------------------------------------------------------------------------------------
        Call Coordinate_Factors(reg_grid(i))  ! also in DVR_Module.  routine computes
                                              ! the grid and a number of functions
                                              ! of the grid that are used in the
                                              ! construction of the matrix elements
!-------------------------------------------------------------------------------------------------------------------------
        Call KE_DVR_Matrices(reg_grid(i))    ! Here is where we calculate the
                                             ! the raw sector functions and the kinetic energy matrices.
                                             ! For some types of DVR basis sets there
                                             ! is only one sector.  The result is 
                                             ! dependent on both the type of DVR and
                                             ! the coordinate system.  The reg_mat arrays contain
                                             ! the results and these arrays have diffrent names.
!-------------------------------------------------------------------------------------------------------------------------
        Call PE_DVR_Matrix(reg_grid(i))       ! Calculates the sector potential energy matrix elements.
!-------------------------------------------------------------------------------------------------------------------------
        Call Final_KE_DVR_Matrices(reg_grid(i))
        Call H_0(reg_grid(i))
        Call Form_Matrix(reg_grid(i))
        IF(diag) THEN
           Call Diagonalize_Global_Matrices(reg_grid(i))
        END IF
        IF (poisson) THEN
            Call Poisson_Equation(reg_grid(i))
        END IF
        IF(diag.or.poisson) THEN
            DEALLOCATE(dvr_mat)           
        END IF
        Call IOsys('rewind all on '//FEDVR_File//' read-and-write', &
                    0,0,0,' ')
        Call IOsys('close '//FEDVR_File,0,0,0,' ')
        Write(iout,*) 'Computation finished for set = ',sets, &
                      ' and dimension = ',i
     END DO Loop_Over_Coordinates
        DEALLOCATE(reg_grid)
  END DO Loop_over_Data_Sets
  CLOSE(input)
  CLOSE(output)
!
1 Format('           **************************************************************' &
         '****************')
2 Format(15x,'Calculation of FEDVR Points, Weights, Polynomials and '                &
             'Matrix Elements')
3 FORMAT(/,25x,'Coordinate System = ',a32)
4 FORMAT(/,25x,'Coordinate = ',a32)
  stop
END PROGRAM FEDVR_DRIVER
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
