!***********************************************************************
! Channel_Module
!**begin prologue     Channel_Module
!**date written       082805   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, Arnoldi, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             CC_Prop
!**purpose            Channel variables and information
!***description       
!***                  
!**references
!**modules needed     See USE statements below
!**                   
!**end prologue       Channel_Module
!***********************************************************************
!***********************************************************************
                           MODULE Channel_Module
                           IMPLICIT NONE
!***********************************************************************
!***********************************************************************

      INTEGER                                               :: Max_L, Max_M
      INTEGER                                               :: nc
      LOGICAL                                               :: b_spline
      TYPE channel_labels
        INTEGER, DIMENSION(:), ALLOCATABLE                  :: labels  
      END TYPE channel_labels
      TYPE(channel_labels), DIMENSION(:),  ALLOCATABLE      :: channel
END MODULE Channel_Module

