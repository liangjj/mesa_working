!***********************************************************************
! CC_Prop_Potential_Module
!**begin prologue     CC_Prop_Potential_Module
!**date written       082805   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             Time_Propagation
!**purpose            Contains all of the major subroutines to generate
!***                  a coupling matrix.  Explicit interfaces are used to allow
!***                  a transparent use of generic subroutines which work
!***                  for both real and complex vectors.  This feature
!***                  permits a single code to be used for both real and
!***                  imaginary time propagation.
!***description       
!***                  
!***                  
!***                  
!***                  
!**references
!**modules needed     See USE statements below
!**end prologue       CC_Prop_module
!***********************************************************************
                     MODULE CC_Prop_Potential_Module
                     USE Iterative_Global
                     USE dvrprop_global
                     USE dvr_global
                     USE dvr_shared
                     USE regional_module
                     USE spatial_wavefunction_module
                     USE initial_state_module
                     USE moment_module
                     USE auto_correlation_module
                     USE non_linear_potential_module
                     USE matrix_vector_multiply_module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                            CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_couple.f
!**begin prologue     v_couple
!**date written       960723   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           time development
!**author             schneider, barry (nsf)
!**source
!**purpose            time-dependent potential
!**
!**description        input data for coupling potential, including the
!**                   mean field interaction for a contact potential.
!**references
!**routines called
!**end prologue       v_couple
  SUBROUTINE v_couple
  IMPLICIT NONE
  INTEGER                      :: natom, intkey
  REAL*8                       :: fpkey, sctlen
  CHARACTER (LEN=80)           :: chrkey
  CHARACTER (LEN=2)            :: atom
  LOGICAL                      :: dollar, fretyp, logkey
  IF(dollar('$v_couple',card,cpass,inp) ) THEN
     IF(system == 'spherical') THEN
        WRITE(iout,1)
        WRITE(iout,3)
        CALL intarr(card,'channel_l_values',l_chan,nchan,' ')
        CALL intarr(card,'channel_m_values',m_chan,nchan,' ')
        DO i=1,nchan
           WRITE(iout,4) l_chan(i), m_chan(i)
        END DO
     ELSE IF(system == 'cylindrical') THEN
        WRITE(iout,2)
        WRITE(iout,5)
        CALL intarr(card,'channel_m_values',m_chan,nchan,' ')
        DO i=1,nchan 
           WRITE(iout,6) m_chan(i)
        END DO
     END IF
  END IF
1    FORMAT(/,1X,'Spherical Coordinates')
2    FORMAT(/,1X,'Cylindrical Coordinates')
3    FORMAT(/,10x,'Channel L Value',10x,'Channel M Value')
4    FORMAT(16x,i3,23x,i3)
5    FORMAT(/,10x,'Channel M Value')
6    FORMAT(16x,i3)
END SUBROUTINE v_couple
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
END MODULE CC_Prop_Potential_Module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
