!***********************************************************************
! Runge_Kutta_Module
!**begin prologue     Runge_Kutta_Module
!**date written       082805   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, integrator, Runge_Kutta
!**
!**author             schneider, b. i.(nsf)
!**source             Time_Propagation
!**purpose            To perform the operation of the exponentiation of a matrix on a vector e
!***                  when the matrix exists in diagonal form.
!***description       
!***                  
!***                  
!***                  
!***                  
!***                  
!***                  
!***                  
!***                  
!***                  
!***                  
!***                  
!***                  
!***                  
!**references
!**modules needed     See USE statements below
!**comments           
!**                   
!**                   
!**                   
!**                   
!**end prologue       Runge_Kutta_Module
!***********************************************************************
!***********************************************************************
                           MODULE Runge_Kutta_Module
                           USE Time_Integrator_Derived_Types
!***********************************************************************
!***********************************************************************
!                          Explicit Interfaces
!***********************************************************************
!
                           INTERFACE RK
             MODULE PROCEDURE Runge_Kutta_2,                                   &
                              Runge_Kutta_3A,                                  &
                              Runge_Kutta_3B,                                  &
                              Runge_Kutta_4A,                                  &
                              Runge_Kutta_4B,                                  &
                              Runge_Kutta_4C,                                  &
                              Runge_Kutta_4D

                       END INTERFACE RK
!
!***********************************************************************
!***********************************************************************
                              CONTAINS
!***********************************************************************
!***********************************************************************
!**begin prologue     Runge_Kutta_2
!**date written       960718   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time integrator
!**
!**author             schneider, b. i.(nsf)
!**source             Time_Propagation
!**purpose      
!***            
!**references
!**routines called    iosys, util and mdutil
!**end prologue       Runge_Kutta_2
!***********************************************************************
  SUBROUTINE Runge_Kutta_2(Method,RK,RK2)
  IMPLICIT NONE
  TYPE(Time_Integrator)                  :: Method
  TYPE(Runge_Kutta)                      :: RK
  TYPE(EPIRK2)                           :: RK2
!
END SUBROUTINE Runge_Kutta_2
!***********************************************************************
!***********************************************************************
!**begin prologue     Runge_Kutta_3A
!**date written       960718   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, Lanczos, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             Time_Propagation
!**purpose      
!***            
!**references
!**routines called    iosys, util and mdutil
!**end prologue       Runge_Kutta_3A
!***********************************************************************
  SUBROUTINE Runge_Kutta_3A(Method,RK,RK3A)
  IMPLICIT NONE
  TYPE(Time_Integrator)                  :: Method
  TYPE(Runge_Kutta)                      :: RK
  TYPE(EPIRK3A)                          :: RK3A
!
END SUBROUTINE Runge_Kutta_3A
!***********************************************************************
!***********************************************************************
!**begin prologue     Runge_Kutta_3B
!**date written       960718   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, Lanczos, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             Time_Propagation
!**purpose      
!***            
!**references
!**routines called    iosys, util and mdutil
!**end prologue       Runge_Kutta_3B
!***********************************************************************
  SUBROUTINE Runge_Kutta_3B(Method,RK,RK3B)
  IMPLICIT NONE
  TYPE(Time_Integrator)                  :: Method
  TYPE(Runge_Kutta)                      :: RK
  TYPE(EPIRK3B)                          :: RK3B
!
END SUBROUTINE Runge_Kutta_3B
!***********************************************************************
!***********************************************************************
!**begin prologue     Runge_Kutta_4A
!**date written       960718   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, Lanczos, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             Time_Propagation
!**purpose      
!***            
!**references
!**routines called    iosys, util and mdutil
!**end prologue       Runge_Kutta_4A
!***********************************************************************
 SUBROUTINE Runge_Kutta_4A(Method,RK,RK4A)
  IMPLICIT NONE
  TYPE(Time_Integrator)                  :: Method
  TYPE(Runge_Kutta)                      :: RK
  TYPE(EPIRK4A)                          :: RK4A
!
END SUBROUTINE Runge_Kutta_4A
!***********************************************************************
!***********************************************************************
!**begin prologue     Runge_Kutta_4B
!**date written       960718   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, Lanczos, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             Time_Propagation
!**purpose      
!***            
!**references
!**routines called    iosys, util and mdutil
!**end prologue       Runge_Kutta_4B
!***********************************************************************
  SUBROUTINE Runge_Kutta_4B(Method,RK,RK4B)
  IMPLICIT NONE
  TYPE(Time_Integrator)                  :: Method
  TYPE(Runge_Kutta)                      :: RK
  TYPE(EPIRK4B)                          :: RK4B
!
END SUBROUTINE Runge_Kutta_4B
!***********************************************************************
!***********************************************************************
!**begin prologue     Runge_Kutta_4C
!**date written       960718   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, Lanczos, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             Time_Propagation
!**purpose      
!***            
!**references
!**routines called    iosys, util and mdutil
!**end prologue       Runge_Kutta_4C
!***********************************************************************
  SUBROUTINE Runge_Kutta_4C(Method,RK,RK4C)
  IMPLICIT NONE
  TYPE(Time_Integrator)                  :: Method
  TYPE(Runge_Kutta)                      :: RK
  TYPE(EPIRK4C)                          :: RK4C
!
END SUBROUTINE Runge_Kutta_4C
!***********************************************************************
!***********************************************************************
!**begin prologue     Runge_Kutta_4D
!**date written       960718   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, Lanczos, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             Time_Propagation
!**purpose      
!***            
!**references
!**routines called    iosys, util and mdutil
!**end prologue       Runge_Kutta_4D
!***********************************************************************
  SUBROUTINE Runge_Kutta_4D(Method,RK,RK4D)
  IMPLICIT NONE
  TYPE(Time_Integrator)                  :: Method
  TYPE(Runge_Kutta)                      :: RK
  TYPE(EPIRK4D)                          :: RK4D
!
END SUBROUTINE Runge_Kutta_4D
!***********************************************************************
!***********************************************************************
           END MODULE Runge_Kutta_Module
!***********************************************************************
!***********************************************************************
