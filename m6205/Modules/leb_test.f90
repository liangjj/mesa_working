!***********************************************************************                
                           MODULE Lebedev
                           USE accuracy
                           USE Data
                           USE Grid_Defined_Types
                           IMPLICIT NONE
!**********************************************************************
  INTEGER, PARAMETER                       :: order_max = 6810
  INTEGER, PARAMETER                       :: rule_max  = 65
  INTEGER                                  :: order
  INTEGER, PARAMETER, DIMENSION(rule_max)  :: table= &
                                                      [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
                                                        1, 1, 1, 1, 1, 0, 1, 0, 0, 1, &
                                                        0, 0, 1, 0, 0, 1, 0, 0, 1, 0, &
                                                        0, 1, 0, 0, 1, 0, 0, 1, 0, 0, &
                                                        1, 0, 0, 1, 0, 0, 1, 0, 0, 1, &
                                                        0, 0, 1, 0, 0, 1, 0, 0, 1, 0, &
                                                        0, 1, 0, 0, 1 ]
!
!          
!***********************************************************************
!***********************************************************************
!                          Explicit Interfaces
!***********************************************************************
!
                            INTERFACE Angular_Grid                       
                       MODULE PROCEDURE Ang_6,                             &
                                        Ang_14
                            END INTERFACE Angular_Grid
!
!
!
!
!***********************************************************************                
!!***********************************************************************               
                           Contains
!***********************************************************************                
!***********************************************************************                
!deck Available_Table
!***begin prologue     Available_Table
!***date written       021231   (yymmdd)                                                
!***revision date               (yymmdd)                                                
!***keywords           lebedev
!***author             schneider, b. i.(nsf)                                            
!***source                                                                        
!***purpose            lebedev rules available
!***                                                                                    
!***description        See below                                                        
!***references                                                                          
!***routines called                                                                     
!***end prologue       Available_Table                                                     
!**********************************************************************
  FUNCTION AVAILABLE_TABLE ( rule )
!
!! Gives the availablity of a given Lebedev rule.
!
!  Modified:
!
!    12 September 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Vyacheslav Lebedev, Dmitri Laikov,
!    A quadrature formula for the sphere of the 131st
!    algebraic order of accuracy,
!    Russian Academy of Sciences Doklady Mathematics,
!    Volume 59, Number 3, 1999, pages 477-481.
!
!
!    AVAILABLE_TABLE, the availability of the rule.
!    * -1, there is no such rule;
!    *  0, there is such a rule, but it is not available in this library.
!    *  1, the rule is available in this library.
!
  IMPLICIT NONE
  INTEGER   :: rule
  INTEGER   :: available_table
!
  IF ( rule < 1 ) THEN
    available_table = - 1
  ELSE IF ( rule_max < rule ) then
    available_table = - 1
  ELSE
    available_table = table(rule)
  END IF
!*****************************************************************************80
  END FUNCTION AVAILABLE_TABLE
!*****************************************************************************80
!*****************************************************************************80
  SUBROUTINE GEN_OH ( ang_grd )
  IMPLICIT NONE
  TYPE(PT_WT)         :: ang_grd
!*****************************************************************************80
!
!! GEN_OH generates points under OH symmetry.
!
!  Discussion:
!
!    Given a point on a sphere, specified by A and B, this routine generates
!    all the equivalent points under OH symmetry, making grid points with
!    weight V.
!
!    The variable NUM is increased by the number of different points
!    generated.
!
!    Depending on CODE, there are from 6 to 48 different but equivalent
!    points that are generated:
!
!      CODE=1:   (0,0,1) etc                                (  6 points)
!      CODE=2:   (0,A,A) etc, A=1/sqrt(2)                   ( 12 points)
!      CODE=3:   (A,A,A) etc, A=1/sqrt(3)                   (  8 points)
!      CODE=4:   (A,A,B) etc, B=sqrt(1-2 A^2)               ( 24 points)
!      CODE=5:   (A,B,0) etc, B=sqrt(1-A^2), A input        ( 24 points)
!      CODE=6:   (A,B,C) etc, C=sqrt(1-A^2-B^2), A, B input ( 48 points)
!
!  Modified:
!
!    11 September 2010
!
!  Author:
!
!    Dmitri Laikov
!
!  Reference:
!
!    Vyacheslav Lebedev, Dmitri Laikov,
!    A quadrature formula for the sphere of the 131st
!    algebraic order of accuracy,
!    Russian Academy of Sciences Doklady Mathematics,
!    Volume 59, Number 3, 1999, pages 477-481.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) CODE, selects the symmetry group.
!
!    Input/output, integer ( kind = 4 ) NUM, presumably a counter for the 
!    total number of points.  It is incremented by the number of points 
!    generated on this call.
!
!    Input, real ( kind = 8 ) A, B, information that may be needed to
!    generate the coordinates of the points (for code = 5 or 6 only).
!
!    Input, real ( kind = 8 ) V, the weight to be assigned the points.
!
!    Output, real ( kind = 8 ) X(NUM), Y(NUM), Z(NUM), W(NUM), the coordinates
!    and weights of the symmetric points generated on this call.
!
!
  IF ( ang_grd%code == 1 ) THEN
       ALLOCATE(ang_grd%pt(1:3,1:6),ang_grd%w(1:6))
       ang_grd%w(1:6)      = ang_grd%wt
       ang_grd%a           = 1.0d0
       ang_grd%pt(1:3,1:6) = 0.d0
       ang_grd%pt(1,1)     =  ang_grd%a
       ang_grd%pt(1,2)     = -ang_grd%a
!--------------------------------------
       ang_grd%pt(2,3)     =  ang_grd%a
       ang_grd%pt(2,4)     = -ang_grd%a
!--------------------------------------
       ang_grd%pt(3,5)     =  ang_grd%a
       ang_grd%pt(3,6)     = -ang_grd%a
!--------------------------------------
       ang_grd%num         = ang_grd%num + 6
  ELSE IF ( ang_grd%code == 2 ) then
       ALLOCATE(ang_grd%pt(1:3,1:12),ang_grd%w(1:12))
       ang_grd%a             = sqrt ( 0.5d0 )
       ang_grd%w(1:12)       = ang_grd%wt
       ang_grd%pt(1:3,1:12)  = 0.d0
!--------------------------------------
       ang_grd%pt(2,1)       =  ang_grd%a
       ang_grd%pt(3,1)       =  ang_grd%a
!--------------------------------------
       ang_grd%pt(2,2)       = -ang_grd%a
       ang_grd%pt(3,2)       =  ang_grd%a
!--------------------------------------
       ang_grd%pt(2,3)       =  ang_grd%a
       ang_grd%pt(3,3)       = -ang_grd%a
!--------------------------------------
       ang_grd%pt(2,4)       = -ang_grd%a
       ang_grd%pt(3,4)       = -ang_grd%a
!--------------------------------------
       ang_grd%pt(1,5)       =  ang_grd%a
       ang_grd%pt(3,5)       =  ang_grd%a
!--------------------------------------
       ang_grd%pt(1,6)       = -ang_grd%a
       ang_grd%pt(3,6)       =  ang_grd%a
!--------------------------------------
       ang_grd%pt(1,7)       =  ang_grd%a
       ang_grd%pt(3,7)       = -ang_grd%a
!--------------------------------------
       ang_grd%pt(1,8)       = -ang_grd%a
       ang_grd%pt(3,8)       = -ang_grd%a
!--------------------------------------
       ang_grd%pt(1,9)       =  ang_grd%a
       ang_grd%pt(2,9)       =  ang_grd%a
!--------------------------------------
       ang_grd%pt(1,10)      = -ang_grd%a
       ang_grd%pt(2,10)      =  ang_grd%a
!--------------------------------------
       ang_grd%pt(1,11)      =  ang_grd%a
       ang_grd%pt(2,11)      = -ang_grd%a
!--------------------------------------
       ang_grd%pt(1,12)      = -ang_grd%a
       ang_grd%pt(2,12)      = -ang_grd%a
!--------------------------------------
       ang_grd%num           = ang_grd%num + 12
  ELSE IF ( ang_grd%code == 3 ) THEN
       ALLOCATE(ang_grd%pt(1:3,1:8),ang_grd%w(1:8))
       ang_grd%a             = sqrt ( 1.0D+000 / 3.0D+000 )
       ang_grd%pt(1:3,1:8)   = 0.d0
       ang_grd%pt(1:3,1)     =  ang_grd%a
!--------------------------------------
       ang_grd%pt(1,2)       = -ang_grd%a
       ang_grd%pt(2:3,2)     =  ang_grd%a
!--------------------------------------
       ang_grd%pt(1,3)       =  ang_grd%a
       ang_grd%pt(2,3)       = -ang_grd%a
       ang_grd%pt(3,3)       =  ang_grd%a
!--------------------------------------
       ang_grd%pt(1:2,4)     = -ang_grd%a
       ang_grd%pt(3,4)       =  ang_grd%a
!--------------------------------------
       ang_grd%pt(1:2,5)     =  ang_grd%a
       ang_grd%pt(3,5)       = -ang_grd%a
!--------------------------------------
       ang_grd%pt(1,6)       = -ang_grd%a
       ang_grd%pt(2,6)       =  ang_grd%a
       ang_grd%pt(3,6)       = -ang_grd%a
!--------------------------------------
       ang_grd%pt(1,7)       =  ang_grd%a
       ang_grd%pt(2:3,7)     = -ang_grd%a
!--------------------------------------
       ang_grd%pt(1:3,8)     = -ang_grd%a
!--------------------------------------
       ang_grd%num           = ang_grd%num  + 8
  ELSE IF ( ang_grd%code == 4 ) THEN
       ALLOCATE(ang_grd%pt(1:3,1:24),ang_grd%w(1:24))
       ang_grd%pt(1:3,1:24)  = 0.d0
       ang_grd%w(1:24)       = ang_grd%wt
       ang_grd%b             = sqrt ( 1.0D+000 - 2.0D+000 * ang_grd%a * ang_grd%a )
       ang_grd%pt(1,1)     =  ang_grd%a
       ang_grd%pt(2,1)     =  ang_grd%a
       ang_grd%pt(3,1)     =  ang_grd%b
!--------------------------------------
       ang_grd%pt(1,2)       = -ang_grd%a
       ang_grd%pt(2,2)       =  ang_grd%a
       ang_grd%pt(3,2)       =  ang_grd%b
!--------------------------------------
       ang_grd%pt(1,3)       =  ang_grd%a
       ang_grd%pt(2,3)       = -ang_grd%a
       ang_grd%pt(3,3)       =  ang_grd%b
!--------------------------------------
       ang_grd%pt(1,4)       = -ang_grd%a
       ang_grd%pt(2,4)       = -ang_grd%a
       ang_grd%pt(3,4)       =  ang_grd%b
!--------------------------------------
       ang_grd%pt(1,5)       =  ang_grd%a
       ang_grd%pt(2,5)       =  ang_grd%a
       ang_grd%pt(3,5)       = -ang_grd%b
!--------------------------------------
       ang_grd%pt(1,6)       = -ang_grd%a
       ang_grd%pt(2,6)       =  ang_grd%a
       ang_grd%pt(3,6)       = -ang_grd%b
!--------------------------------------
       ang_grd%pt(1,7)       =  ang_grd%a
       ang_grd%pt(2,7)       = -ang_grd%a
       ang_grd%pt(3,7)       = -ang_grd%b
!--------------------------------------
       ang_grd%pt(1,8)       = -ang_grd%a
       ang_grd%pt(2,8)       = -ang_grd%a
       ang_grd%pt(3,8)       = -ang_grd%b
!--------------------------------------
       ang_grd%pt(1,9)       =  ang_grd%a
       ang_grd%pt(2,9)       =  ang_grd%b
       ang_grd%pt(3,9)       =  ang_grd%a
!--------------------------------------
       ang_grd%pt(1,10)      = -ang_grd%a
       ang_grd%pt(2,10)      =  ang_grd%b
       ang_grd%pt(3,10)      =  ang_grd%a
!--------------------------------------
       ang_grd%pt(1,11)      =  ang_grd%a
       ang_grd%pt(2,11)      = -ang_grd%b
       ang_grd%pt(3,11)      =  ang_grd%a
!--------------------------------------
       ang_grd%pt(1,12)      = -ang_grd%a
       ang_grd%pt(2,12)      = -ang_grd%b
       ang_grd%pt(3,12)      =  ang_grd%a
!--------------------------------------
       ang_grd%pt(1,13)      =  ang_grd%a
       ang_grd%pt(2,13)      =  ang_grd%b
       ang_grd%pt(3,13)      = -ang_grd%a
!--------------------------------------
       ang_grd%pt(1,14)      = -ang_grd%a
       ang_grd%pt(2,14)      =  ang_grd%b
       ang_grd%pt(3,14)      = -ang_grd%a
!--------------------------------------
       ang_grd%pt(1,15)      =  ang_grd%a
       ang_grd%pt(2,15)      = -ang_grd%b
       ang_grd%pt(3,15)      = -ang_grd%a
!--------------------------------------
       ang_grd%pt(1,16)      = -ang_grd%a
       ang_grd%pt(2,16)      = -ang_grd%b
       ang_grd%pt(3,16)      = -ang_grd%a
!--------------------------------------
       ang_grd%pt(1,17)      =  ang_grd%b
       ang_grd%pt(2,17)      =  ang_grd%a
       ang_grd%pt(3,17)      =  ang_grd%a
!--------------------------------------
       ang_grd%pt(1,18)      = -ang_grd%b
       ang_grd%pt(2,18)      =  ang_grd%a
       ang_grd%pt(3,18)      =  ang_grd%a
!--------------------------------------
       ang_grd%pt(1,19)      =  ang_grd%b
       ang_grd%pt(2,19)      = -ang_grd%a
       ang_grd%pt(3,19)      =  ang_grd%a
!--------------------------------------
       ang_grd%pt(1,20)      = -ang_grd%b
       ang_grd%pt(2,20)      = -ang_grd%a
       ang_grd%pt(3,20)      =  ang_grd%a
!--------------------------------------
       ang_grd%pt(1,21)      =  ang_grd%b
       ang_grd%pt(2,21)      =  ang_grd%a
       ang_grd%pt(3,21)      = -ang_grd%a
!--------------------------------------
       ang_grd%pt(1,22)      = -ang_grd%b
       ang_grd%pt(2,22)      =  ang_grd%a
       ang_grd%pt(3,22)      = -ang_grd%a
!--------------------------------------
       ang_grd%pt(1,23)      =  ang_grd%b
       ang_grd%pt(2,23)      = -ang_grd%a
       ang_grd%pt(3,23)      = -ang_grd%a
!--------------------------------------
       ang_grd%pt(1,24)      = -ang_grd%b
       ang_grd%pt(2,24)      = -ang_grd%a
       ang_grd%pt(3,24)      = -ang_grd%a
!--------------------------------------
       ang_grd%num           = ang_grd%num + 24
  ELSE IF ( ang_grd%code == 5 ) THEN
    ALLOCATE(ang_grd%pt(1:3,1:24),ang_grd%w(1:24))
    ang_grd%b = sqrt ( 1.0D+000 - ang_grd%a * ang_grd%a )
    ang_grd%pt(1:3,1:24)  = 0.d0
    ang_grd%w(1:24)       = ang_grd%wt
    ang_grd%pt(1,1)       =  ang_grd%a
    ang_grd%pt(2,1)       =  ang_grd%b
!--------------------------------------
    ang_grd%pt(1,2)       = -ang_grd%a
    ang_grd%pt(2,2)       =  ang_grd%b
!--------------------------------------
    ang_grd%pt(1,3)       =  ang_grd%a
    ang_grd%pt(2,3)       = -ang_grd%b
!--------------------------------------
    ang_grd%pt(1,4)       = -ang_grd%a
    ang_grd%pt(2,4)       = -ang_grd%b
!--------------------------------------
    ang_grd%pt(1,5)       =  ang_grd%b
    ang_grd%pt(2,5)       =  ang_grd%a
!--------------------------------------
    ang_grd%pt(1,6)       = -ang_grd%b
    ang_grd%pt(2,6)       =  ang_grd%a
!--------------------------------------
    ang_grd%pt(1,7)       =  ang_grd%b
    ang_grd%pt(2,7)       = -ang_grd%a
!--------------------------------------
    ang_grd%pt(1,8)       = -ang_grd%b
    ang_grd%pt(2,8)       = -ang_grd%a
!--------------------------------------
    ang_grd%pt(1,9)       =  ang_grd%a
    ang_grd%pt(3,9)       =  ang_grd%b
!--------------------------------------
    ang_grd%pt(1,10)      = -ang_grd%a
    ang_grd%pt(3,10)      =  ang_grd%b
!--------------------------------------
    ang_grd%pt(1,11)      =  ang_grd%a
    ang_grd%pt(3,11)      = -ang_grd%b
!--------------------------------------
    ang_grd%pt(1,12)      = -ang_grd%a
    ang_grd%pt(3,12)      = -ang_grd%b
!--------------------------------------
    ang_grd%pt(1,13)      =  ang_grd%b
    ang_grd%pt(3,13)      =  ang_grd%a
!--------------------------------------
    ang_grd%pt(1,14)      = -ang_grd%b
    ang_grd%pt(3,14)      =  ang_grd%a
!--------------------------------------
    ang_grd%pt(1,15)      =  ang_grd%b
    ang_grd%pt(3,15)      = -ang_grd%a
!--------------------------------------
    ang_grd%pt(1,16)      = -ang_grd%b
    ang_grd%pt(3,16)      = -ang_grd%a
!--------------------------------------
    ang_grd%pt(2,17)      =  ang_grd%a
    ang_grd%pt(3,17)      =  ang_grd%b
!--------------------------------------
    ang_grd%pt(2,18)      = -ang_grd%a
    ang_grd%pt(3,18)      =  ang_grd%b
!--------------------------------------
    ang_grd%pt(2,19)      =  ang_grd%a
    ang_grd%pt(3,19)      = -ang_grd%b
!--------------------------------------
    ang_grd%pt(2,20)      = -ang_grd%a
    ang_grd%pt(3,20)      = -ang_grd%b
!--------------------------------------
    ang_grd%pt(2,21)      =  ang_grd%b
    ang_grd%pt(3,21)      =  ang_grd%a
!--------------------------------------
    ang_grd%pt(2,22)      = -ang_grd%b
    ang_grd%pt(3,22)      =  ang_grd%a
!--------------------------------------
    ang_grd%pt(2,23)      =  ang_grd%b
    ang_grd%pt(3,23)      = -ang_grd%a
!--------------------------------------
    ang_grd%pt(2,24)      = -ang_grd%b
    ang_grd%pt(3,24)      = -ang_grd%a
!--------------------------------------
    ang_grd%num           = ang_grd%num + 24 
!--------------------------------------
  ELSE IF ( ang_grd%code == 6 ) THEN
    ALLOCATE(ang_grd%pt(1:3,1:48),ang_grd%w(1:48))
    ang_grd%c = sqrt ( 1.0D+000 - ang_grd%a * ang_grd%a - ang_grd%b * ang_grd%b )
    ang_grd%pt(1:3,1:48)  = 0.d0
    ang_grd%w(1:48)       = ang_grd%wt
!--------------------------------------
    ang_grd%pt(1,1)       =  ang_grd%a
    ang_grd%pt(2,1)       =  ang_grd%b
    ang_grd%pt(3,1)       =  ang_grd%c
!--------------------------------------
    ang_grd%pt(1,2)       = -ang_grd%a
    ang_grd%pt(2,2)       =  ang_grd%b
    ang_grd%pt(3,2)       =  ang_grd%c
!--------------------------------------
    ang_grd%pt(1,3)       =  ang_grd%a
    ang_grd%pt(2,3)       = -ang_grd%b
    ang_grd%pt(3,3)       =  ang_grd%c
!--------------------------------------
    ang_grd%pt(1,4)       = -ang_grd%a
    ang_grd%pt(2,4)       = -ang_grd%b
    ang_grd%pt(3,4)       =  ang_grd%c
!--------------------------------------
    ang_grd%pt(1,5)       =  ang_grd%a
    ang_grd%pt(2,5)       =  ang_grd%b
    ang_grd%pt(3,5)       = -ang_grd%c
!--------------------------------------
    ang_grd%pt(1,6)       = -ang_grd%a
    ang_grd%pt(2,6)       =  ang_grd%b
    ang_grd%pt(3,6)       = -ang_grd%c
!--------------------------------------
    ang_grd%pt(1,7)       =  ang_grd%a
    ang_grd%pt(2,7)       = -ang_grd%b
    ang_grd%pt(3,7)       = -ang_grd%c
!--------------------------------------
    ang_grd%pt(1,8)       = -ang_grd%a
    ang_grd%pt(2,8)       = -ang_grd%b
    ang_grd%pt(3,8)       = -ang_grd%c
!--------------------------------------
    ang_grd%pt(1,9)       =  ang_grd%a
    ang_grd%pt(2,9)       =  ang_grd%c
    ang_grd%pt(3,9)       =  ang_grd%b
!--------------------------------------
    ang_grd%pt(1,10)      = -ang_grd%a
    ang_grd%pt(2,10)      =  ang_grd%c
    ang_grd%pt(3,10)      =  ang_grd%b
!--------------------------------------
    ang_grd%pt(1,11)      =  ang_grd%a
    ang_grd%pt(2,11)      = -ang_grd%c
    ang_grd%pt(3,11)      =  ang_grd%b
!--------------------------------------
    ang_grd%pt(1,12)      = -ang_grd%a
    ang_grd%pt(2,12)      = -ang_grd%c
    ang_grd%pt(3,12)      =  ang_grd%b
!--------------------------------------
    ang_grd%pt(1,13)      =  ang_grd%a
    ang_grd%pt(2,13)      =  ang_grd%c
    ang_grd%pt(3,13)      = -ang_grd%b
!--------------------------------------
    ang_grd%pt(1,14)      = -ang_grd%a
    ang_grd%pt(2,14)      =  ang_grd%c
    ang_grd%pt(3,14)      = -ang_grd%b
!--------------------------------------
    ang_grd%pt(1,15)      =  ang_grd%a
    ang_grd%pt(2,15)      = -ang_grd%c
    ang_grd%pt(3,15)      = -ang_grd%b
!--------------------------------------
    ang_grd%pt(1,16)      = -ang_grd%a
    ang_grd%pt(2,16)      = -ang_grd%c
    ang_grd%pt(3,16)      = -ang_grd%b
!--------------------------------------
    ang_grd%pt(1,17)      =  ang_grd%b
    ang_grd%pt(2,17)      =  ang_grd%a
    ang_grd%pt(3,17)      =  ang_grd%c
!--------------------------------------
    ang_grd%pt(1,18)      = -ang_grd%b
    ang_grd%pt(2,18)      =  ang_grd%a
    ang_grd%pt(3,18)      =  ang_grd%c
!--------------------------------------
    ang_grd%pt(1,19)      =  ang_grd%b
    ang_grd%pt(2,19)      = -ang_grd%a
    ang_grd%pt(3,19)      =  ang_grd%c
!--------------------------------------
    ang_grd%pt(1,20)      = -ang_grd%b
    ang_grd%pt(2,20)      = -ang_grd%a
    ang_grd%pt(3,20)      =  ang_grd%c
!--------------------------------------
    ang_grd%pt(1,21)      =  ang_grd%b
    ang_grd%pt(2,21)      =  ang_grd%a
    ang_grd%pt(3,21)      = -ang_grd%c
!--------------------------------------
    ang_grd%pt(1,22)      = -ang_grd%b
    ang_grd%pt(2,22)      =  ang_grd%a
    ang_grd%pt(3,22)      = -ang_grd%c
!--------------------------------------
    ang_grd%pt(1,23)      =  ang_grd%b
    ang_grd%pt(2,23)      = -ang_grd%a
    ang_grd%pt(3,23)      = -ang_grd%c
!--------------------------------------
    ang_grd%pt(1,24)      = -ang_grd%b
    ang_grd%pt(2,24)      = -ang_grd%a
    ang_grd%pt(3,24)      = -ang_grd%c
!--------------------------------------
    ang_grd%pt(1,25)      =  ang_grd%b
    ang_grd%pt(2,25)      =  ang_grd%c
    ang_grd%pt(3,25)      =  ang_grd%a
!--------------------------------------
    ang_grd%pt(1,26)      = -ang_grd%b
    ang_grd%pt(2,26)      =  ang_grd%c
    ang_grd%pt(3,26)      =  ang_grd%a
!--------------------------------------
    ang_grd%pt(1,27)      =  ang_grd%b
    ang_grd%pt(2,27)      = -ang_grd%c
    ang_grd%pt(3,27)      =  ang_grd%a
!--------------------------------------
    ang_grd%pt(1,28)      = -ang_grd%b
    ang_grd%pt(2,28)      = -ang_grd%c
    ang_grd%pt(3,28)      =  ang_grd%a
!--------------------------------------
    ang_grd%pt(1,29)      =  ang_grd%b
    ang_grd%pt(2,29)      =  ang_grd%c
    ang_grd%pt(3,29)      = -ang_grd%a
!--------------------------------------
    ang_grd%pt(1,30)      = -ang_grd%b
    ang_grd%pt(2,30)      =  ang_grd%c
    ang_grd%pt(3,30)      = -ang_grd%a
!--------------------------------------
    ang_grd%pt(1,31)      =  ang_grd%b
    ang_grd%pt(2,31)      = -ang_grd%c
    ang_grd%pt(3,31)      = -ang_grd%a
!--------------------------------------
    ang_grd%pt(1,32)      = -ang_grd%b
    ang_grd%pt(2,32)      = -ang_grd%c
    ang_grd%pt(3,32)      = -ang_grd%a
!--------------------------------------
    ang_grd%pt(1,33)      =  ang_grd%c
    ang_grd%pt(2,33)      =  ang_grd%a
    ang_grd%pt(3,33)      =  ang_grd%b
!--------------------------------------
    ang_grd%pt(1,34)      = -ang_grd%c
    ang_grd%pt(2,34)      =  ang_grd%a
    ang_grd%pt(3,34)      =  ang_grd%b
!--------------------------------------
    ang_grd%pt(1,35)      =  ang_grd%c
    ang_grd%pt(2,35)      = -ang_grd%a
    ang_grd%pt(3,35)      =  ang_grd%b
!--------------------------------------
    ang_grd%pt(1,36)      = -ang_grd%c
    ang_grd%pt(2,36)      = -ang_grd%a
    ang_grd%pt(3,36)      =  ang_grd%b
!--------------------------------------
    ang_grd%pt(1,37)      =  ang_grd%c
    ang_grd%pt(2,37)      =  ang_grd%a
    ang_grd%pt(3,37)      = -ang_grd%b
!--------------------------------------
    ang_grd%pt(1,38)      = -ang_grd%c
    ang_grd%pt(2,38)      =  ang_grd%a
    ang_grd%pt(3,38)      = -ang_grd%b
!--------------------------------------
    ang_grd%pt(1,39)      =  ang_grd%c
    ang_grd%pt(2,39)      = -ang_grd%a
    ang_grd%pt(3,39)      = -ang_grd%b
!--------------------------------------
    ang_grd%pt(1,40)      = -ang_grd%c
    ang_grd%pt(2,40)      = -ang_grd%a
    ang_grd%pt(3,40)      = -ang_grd%b
!--------------------------------------
    ang_grd%pt(1,41)      =  ang_grd%c
    ang_grd%pt(2,41)      =  ang_grd%b
    ang_grd%pt(3,41)      =  ang_grd%a
!--------------------------------------
    ang_grd%pt(1,42)      = -ang_grd%c
    ang_grd%pt(2,42)      =  ang_grd%b
    ang_grd%pt(3,42)      =  ang_grd%a
!--------------------------------------
    ang_grd%pt(1,43)      =  ang_grd%c
    ang_grd%pt(2,43)      = -ang_grd%b
    ang_grd%pt(3,43)      =  ang_grd%a
!--------------------------------------
    ang_grd%pt(1,44)      = -ang_grd%c
    ang_grd%pt(2,44)      = -ang_grd%b
    ang_grd%pt(3,44)      =  ang_grd%a
!--------------------------------------
    ang_grd%pt(1,45)      =  ang_grd%c
    ang_grd%pt(2,45)      =  ang_grd%b
    ang_grd%pt(3,45)      = -ang_grd%a
!--------------------------------------
    ang_grd%pt(1,46)      = -ang_grd%c
    ang_grd%pt(2,46)      =  ang_grd%b
    ang_grd%pt(3,46)      = -ang_grd%a
!--------------------------------------
    ang_grd%pt(1,47)      =  ang_grd%c
    ang_grd%pt(2,47)      = -ang_grd%b
    ang_grd%pt(3,47)      = -ang_grd%a
!--------------------------------------
    ang_grd%pt(1,48)      = -ang_grd%c
    ang_grd%pt(2,48)      = -ang_grd%b
    ang_grd%pt(3,48)      = -ang_grd%a
!------------------------------
    ang_grd%num      = ang_grd%num + 48

  ELSE

    Call lnkerr('GEN_OH - Fatal error')
 
  END IF
!*****************************************************************************80
!*****************************************************************************80
  END SUBROUTINE GEN_OH 
!*****************************************************************************80
!*****************************************************************************80

  Subroutine Generate_Points_Weights(leb_rule)

!*****************************************************************************80
  IMPLICIT NONE
  TYPE(RULE)                 :: leb_rule

  if ( order == 6 ) then
      Call Ang_Grid(leb_rule%rule_06)
  else if ( order == 14 ) then
      Call Ang_Grid(leb_rule%rule_14)
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LD_BANG_GRD%PT_ORDER - Fatal error!'
    write ( *, '(a)' ) '  Unexpected value of ORDER.'
    stop
  end if

  return
!*****************************************************************************80
  END Subroutine Generate_Points_Weights
!*****************************************************************************80

  Subroutine Ang_6 (rule_06 )
  IMPLICIT NONE
  TYPE(ld06)          ::rule_06 
!*****************************************************************************80
  RETURN
  END Subroutine Ang_6
!
  Subroutine Ang_14 (rule_14 )
  implicit none
  TYPE(ld14)          ::rule_14 
  return
  END Subroutine Ang_14
function order_table ( rule )
  implicit none
  integer ( kind = 4 ), parameter :: rule_max = 65

  integer ( kind = 4 ) order_table
  integer ( kind = 4 ) rule
  integer ( kind = 4 ), save :: table(rule_max) = (/ &
       6,   14,   26,   38,   50,   74,   86,  110,  146,  170, &
     194,  230,  266,  302,  350,  386,  434,  482,  530,  590, &
     650,  698,  770,  830,  890,  974, 1046, 1118, 1202, 1274, &
    1358, 1454, 1538, 1622, 1730, 1814, 1910, 2030, 2126, 2222, &
    2354, 2450, 2558, 2702, 2810, 2930, 3074, 3182, 3314, 3470, &
    3590, 3722, 3890, 4010, 4154, 4334, 4466, 4610, 4802, 4934, &
    5090, 5294, 5438, 5606, 5810 /)

  if ( rule < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ORDER_TABLE - Fatal error!'
    write ( *, '(a)' ) '  RULE < 1.'
    stop
  else if ( rule_max < rule ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ORDER_TABLE - Fatal error!'
    write ( *, '(a)' ) '  RULE_MAX < RULE.'
    stop
  end if

  order_table = table(rule)

  END function order_table 
function precision_table ( rule )
!
  implicit none

  integer ( kind = 4 ), parameter :: rule_max = 65

  integer ( kind = 4 ) precision_table
  integer ( kind = 4 ) rule
  integer ( kind = 4 ), save :: table(rule_max) = (/ &
      3,   5,   7,   9,  11,  13,  15,  17,  19,  21, &
     23,  25,  27,  29,  31,  33,  35,  37,  39,  41, &
     43,  45,  47,  49,  51,  53,  55,  57,  59,  61, &
     63,  65,  67,  69,  71,  73,  75,  77,  79,  81, &
     83,  85,  87,  89,  91,  93,  95,  97,  99, 101, &
    103, 105, 107, 109, 111, 113, 115, 117, 119, 121, &
    123, 125, 127, 129, 131 /)

  if ( rule < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PRECISION_TABLE - Fatal error!'
    write ( *, '(a)' ) '  RULE < 1.'
    stop
  else if ( rule_max < rule ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PRECISION_TABLE - Fatal error!'
    write ( *, '(a)' ) '  RULE_MAX < RULE.'
    stop
  end if

  precision_table = table(rule)

  END function precision_table 
  END MODULE Lebedev
