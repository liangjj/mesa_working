!
MODULE Matrix_Defined_Types
!***begin prologue     Matrix_Defined_Types
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source             
!***purpose            matrix defined type
!***description        
!***                   
!
!***references

!***routines called    
!***end prologue       Matrix_Defined_Types
  USE accuracy
  IMPLICIT NONE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!---------------------------------------------------------------------
!                    Derived types and Allocated variables for 
!                             the dvr quantities
!---------------------------------------------------------------------
!
!
  TYPE Real_Matrix
       CHARACTER(LEN=16)     :: type
  END TYPE Real_Matrix

  TYPE Complex_Matrix
       CHARACTER(LEN=16)     :: type
  END TYPE Complex_Matrix

  TYPE Real_Vector
       CHARACTER(LEN=16)     :: type
  END TYPE Real_Vector

  TYPE Complex_Vector
       CHARACTER(LEN=16)     :: type
  END TYPE Complex_Vector

  TYPE Real_Triangular_Array
       CHARACTER(LEN=16)     :: type
  END TYPE Real_Triangular_Array

  TYPE Complex_Triangular_Array
       CHARACTER(LEN=16)     :: type
  END TYPE Complex_Triangular_Array

  TYPE(Real_Matrix)              :: type_real_matrix
  TYPE(Complex_Matrix)           :: type_complex_array
  TYPE(Real_Vector)              :: type_real_vector
  TYPE(Complex_Vector)           :: type_complex_vector
  TYPE(Real_Triangular_Array)    :: type_real_triangular_array
  TYPE(Complex_Triangular_Array) :: type_complex_triangular_array
!
!
END MODULE Matrix_Defined_Types
