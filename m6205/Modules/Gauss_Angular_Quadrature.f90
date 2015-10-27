!***********************************************************************                
                           MODULE Gauss_Angular_Quadrature
!deck Gauss_Angular_Quadrature
!***begin prologue     Gauss_Angular_Quadrature
!***date written       20140706  (yyyymmdd)                                                
!***revision date                (yyyymmdd)                                                
!***keywords           Gauss_Angular_Quadrature
!***author             schneider, b. i.(nist)                                            
!***source                                                                        
!***purpose            Compute the points and weights for an angular quadrature over a sphere
!***                   This module handles either product rules in theta and phi or composite
!***                   rules using the Lebedev quadrature.                     
!***description        The product rules use as default a Gauss-Legendre quadrature from (-1,1) in Cos theta
!***                   and Simpsons rule in phi.  This accounts for the weight functions in the integrand.      
!***                   The Lebedev rules can only handle specific orders as the rules depend on the mapping      
!***                   to specific symmetry groups.  The reader must look at the routines to determine which
!***                   rules are qavailable althogh the invocation of the routine will let you know if you have
!***                   chosen rule that is unvailable.  The module is wrtten with a condsiderable amount of
!***                   overloading so new rules may be added with a minimal of effort.
!***references                                                                          
!***routines called                                                                     
!***end prologue       Available_Table                                 
                           USE input_output
                           USE accuracy
                           USE Data
                           USE Grid_Defined_Types
                           USE Gauss_Quadrature
                           IMPLICIT NONE
!**********************************************************************

!             x= cos(theta) * sin(phi)   phi goes 0 to pi and theta -pi to pi
!             x= sin(theta) * sin(phi)
!             z= cos(phi)
!             Reversing the definitions of phi and theta make no difference since
!             its just a label.  Be careful of the ranges of course.
!**********************************************************************
!          
!***********************************************************************
!***********************************************************************
!                          Explicit Interfaces
!***********************************************************************
!
                            INTERFACE Gauss_Angular_Grid                       
                       MODULE PROCEDURE Theta_Angular_Quadrature,          &
                                        Phi_Angular_Quadrature
                            END INTERFACE Gauss_Angular_Grid                       
!
!
!
!
!************************************************************************                
!!***********************************************************************               
                           Contains
!***********************************************************************                
!***********************************************************************                
  SUBROUTINE Theta_Angular_Quadrature(thet_ang, irange)
  IMPLICIT NONE
  TYPE(THETA)                           :: thet_ang
  INTEGER                               :: i
  LOGICAL, OPTIONAL                     :: irange
!
  IF ( PRESENT(irange) == .false. ) THEN
       thet_ang%edge(2) = -1.0d0
       thet_ang%edge(1) =  1.0d0
  END IF
!
  ALLOCATE( thet_ang%q(1:thet_ang%n_pts),thet_ang%sin_thet(1:thet_ang%n_pts),             &
            thet_ang%wt(1:thet_ang%n_pts) )
  CALL gauss(thet_ang%q,thet_ang%wt,edge=thet_ang%edge,                                   &
             type_quadrature=thet_ang%type_quadrature, fixed_point=thet_ang%fixed_point,  &
             n=thet_ang%n_pts,print=.true.)
  DO i=1,thet_ang%n_pts
     thet_ang%sin_thet(i)=SQRT((1.d0-thet_ang%q(i)*thet_ang%q(i)))
  END DO
!*****************************************************************************80
  END SUBROUTINE Theta_Angular_Quadrature
!*****************************************************************************80
!*****************************************************************************80
  SUBROUTINE Phi_Angular_Quadrature( phi_ang, irange)
  IMPLICIT NONE
  TYPE(PHI)                             :: phi_ang
  INTEGER                               :: i
  LOGICAL, OPTIONAL                     :: irange
!
  IF ( PRESENT(irange) == .false. ) THEN
       phi_ang%edge(2) = two_pi
       phi_ang%edge(1) = 0.d0
  END IF
  ALLOCATE(phi_ang%q(1:phi_ang%n_pts),phi_ang%sin_phi(1:phi_ang%n_pts),                 &
           phi_ang%cos_phi(1:phi_ang%n_pts), phi_ang%wt(1:phi_ang%n_pts) )
  CALL gauss(phi_ang%q,phi_ang%wt,edge=phi_ang%edge,                                    &
             type_quadrature=phi_ang%type_quadrature, fixed_point=phi_ang%fixed_point,  &
             n=phi_ang%n_pts,print=.true.)
  DO  i=1,phi_ang%n_pts
    phi_ang%sin_phi(i)=SIN(phi_ang%q(i))
    phi_ang%cos_phi(i)=COS(phi_ang%q(i))
  END DO
!*****************************************************************************80
  END SUBROUTINE Phi_Angular_Quadrature
!*****************************************************************************80
!*****************************************************************************80
  END MODULE Gauss_Angular_Quadrature
!*****************************************************************************80
