!deck dvr_driver
!**begin prologue     dvr_driver
!**date written       060711   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           dvr!**author             schneider, barry (nsf)
!**source
!**purpose            main dvr driver
!**description        calls routines to input data, construct dvr matrices
!**                   
!**references
!**routines called      Name                    Location    
!                       ----                    --------
!                     Input_Prop             Prop_Sources
!                     Space_Prop             Prop_Sources
!                     Regional_Matrices      Prop_Modules:regional_module
!
!**modules used         Name                    Location               
!                       ----                    --------
!                    dvrprop_global          Modules   
!                    Regional_Matrices       Prop_Modules
!
!**end prologue       dvr_driver
  PROGRAM dvr_driver
  USE dvrprop_global
  USE regional_module
  IMPLICIT NONE
  REAL*4                                   :: secnds
  REAL*4, DIMENSION(10)                    :: del_t
  INTEGER                                  :: iostat 
!
! Read in the Input
!
  time(1)=secnds(0.0)
  CALL Input_DVR
  time(2)=secnds(0.0)
  del_t(1) = time(2) - time(1)
!
! Calculate the Global FEDVR Matrices
!
  CALL Space_DVR
  time(3)=secnds(0.0)
  del_t(2) = time(3) - time(2)
!
! Compute the Regional Matrices from the FEDVR Global Matrices
!
  CALL Regional_Matrices
  time(4)=secnds(0.0)
  del_t(3) = time(4)-time(3)
  Write(iout,1) del_t(1:3)
!
1 FORMAT('***********************************************'                 &
            '*************************'                                    &
            /,10X,'Time to input basic data               = ',f15.8,/,     &
            /,10X,'Time to compute global FEDVR spatial ',                 &
            /,10x,'matrices                               = ',f15.8,/,     &
            /,10x,'Time to compute the regional spatial '                  &
            /,10x,'matrices                               = ',f15.8)
  stop
END PROGRAM dvr_driver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
