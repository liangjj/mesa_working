!
MODULE Time_Integrator_Derived_Types
!***begin prologue     Time_Integrator_Derived_Types
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           time, propagation, exponential propagator
!***author             schneider, b. i.(nsf)
!***source             
!***purpose            global shared variables for exponential propagator
!***description        this routine defines the global variables
!***                   and data needed for the exponential propagator code
!
!***references

!***routines called    
!***end prologue       Time_Integrator_Derived_Types
!  USE accuracy
  IMPLICIT NONE
  INTEGER, PARAMETER                                :: idp = SELECTED_REAL_KIND(15,307)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!---------------------------------------------------------------------
!                    Derived types and Allocated variables for 
!                             the Time_Integrator code
!---------------------------------------------------------------------
!
  TYPE Second_Order    
       CHARACTER(LEN=4)                             :: Type
  END TYPE Second_Order

  TYPE Third_Order    
       CHARACTER(LEN=4)                             :: Type
  END TYPE Third_Order

  TYPE Fourth_Order   
       CHARACTER(LEN=4)                             :: Type
  END TYPE Fourth_Order

  TYPE EPIRK2    
       CHARACTER(LEN=4)                             :: Type
  END TYPE EPIRK2

  TYPE EPIRK3A    
       CHARACTER(LEN=4)                             :: Type
  END TYPE EPIRK3A

  TYPE EPIRK3B    
       CHARACTER(LEN=4)                             :: Type
  END TYPE EPIRK3B

  TYPE EPIRK4A    
       CHARACTER(LEN=4)                             :: Type
  END TYPE EPIRK4A

  TYPE EPIRK4B    
       CHARACTER(LEN=4)                             :: Type
  END TYPE EPIRK4B

  TYPE EPIRK4C    
       CHARACTER(LEN=4)                             :: Type
  END TYPE EPIRK4C

  TYPE EPIRK4D    
       CHARACTER(LEN=4)                             :: Type
  END TYPE EPIRK4D


  TYPE Runge_Kutta
       
       TYPE(EPIRK2)                                 :: RK2
       TYPE(EPIRK3A)                                :: RK3A
       TYPE(EPIRK3B)                                :: RK3B
       TYPE(EPIRK4A)                                :: RK4A
       TYPE(EPIRK4B)                                :: RK4B
       TYPE(EPIRK4C)                                :: RK4C
       TYPE(EPIRK4D)                                :: RK4D

  END TYPE Runge_Kutta

  TYPE Quadrature

       TYPE(Second_Order)                           :: SO
       TYPE(Third_Order)                            :: TO
       TYPE(Fourth_Order)                           :: FO

  END TYPE Quadrature

  TYPE Time_Integrator

       TYPE(Runge_Kutta)                           :: RK
       TYPE(Quadrature)                            :: Q

  END TYPE Time_Integrator


  TYPE (Time_Integrator)                           :: Method

!
!
END MODULE Time_Integrator_Derived_Types
