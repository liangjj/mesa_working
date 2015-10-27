MODULE Prolate_Functions
!***begin prologue     Prolate_Functions
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)                                                                                            
!***keywords           prolate functions                                                              
!***author             schneider, b. i.(nsf)                                                                                        
!***source                                                                                                                    
!***purpose            
!***description        
!***                   
!***                   
!***                   
!***                   
!***                   
!***                                                                                                  
!                                                                                                                                   
!***references                                                                                                                      
!***routines called                                                                                                                 
!***end prologue       Prolate_Functions                                                                               
!
  USE accuracy
  USE Data_Module
  IMPLICIT NONE
  INTEGER                                 :: lorder
  INTEGER                                 :: morder
  INTEGER                                 :: mabs
  INTEGER                                 :: lsum
  INTEGER                                 :: msum
  INTEGER                                 :: meo
  REAL(idp)                               :: a
  REAL(idp)                               :: R
  REAL(idp)                               :: radius_moeq
  REAL(idp)                               :: xi1
  REAL(idp)                               :: xi2
  REAL(idp)                               :: eta1
  REAL(idp)                               :: eta2
  REAL(idp)                               :: varphi1
  REAL(idp)                               :: varphi2
  REAL(idp)                               :: x1
  REAL(idp)                               :: y1
  REAL(idp)                               :: z1
  REAL(idp)                               :: rho1
  REAL(idp)                               :: x2
  REAL(idp)                               :: y2
  REAL(idp)                               :: z2
  REAL(idp)                               :: rho2
  REAL(idp)                               :: xi_small
  REAL(idp)                               :: xi_large
  REAL(idp)                               :: facm
  REAL(idp)                               :: vardm
  REAL(idp)                               :: dl21
  REAL(idp)                               :: temp
  REAL(idp)                               :: csum_real
  REAL(idp)                               :: csum_imag
  REAL(idp)                               :: varphi_diff
  REAL(idp)                               :: ctemp_real
  REAL(idp)                               :: ctemp_imag
  REAL(idp)                               :: dx
  REAL(idp)                               :: dy
  REAL(idp)                               :: dz
  REAL(idp)                               :: rsqr
  REAL(idp)                               :: r1
  REAL(idp)                               :: r2
  REAL(idp)                               :: r_12
  REAL(idp)                               :: r_12_invs
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE Prolate_Functions
