!deck prop_main
!**begin prologue     prop_main
!**date written       960723   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           time development
!**author             schneider, barry (nsf)
!**source
!**purpose            main propagator program
!**description        calls routines to input data, construct FEDVR propagators
!**                   and then propagate the initial wave function.
!**                   
!**references
!**routines called    drum, Input_Prop, Space_Prop, Regional_Matrices,
!**                   v_couple, Arnoldi_Propagator, Time_Prop, Split_Operator
!**modules            dvrprop_globalt
!**end prologue       prop_main
  PROGRAM prop_main
  USE dvrprop_global
  USE arnoldi_module
  USE split_operator_module
  IMPLICIT NONE
  REAL*4                                   :: secnds
  REAL*4, DIMENSION(10)                    :: del_t
  INTEGER                                  :: iostat
  INTEGER                                  :: i, bigs, bigv 
  INTEGER, DIMENSION(2)                    :: words
  CHARACTER(LEN=4096)                      :: ops
!
  CALL Drum
  CALL IOsys('read character options from rwf',-1,0,0,ops)
  OPEN (UNIT=iplot(1),FILE='misc',                                      &
        ACCESS='sequential',FORM='formatted',                           &
        IOSTAT=IOSTAT,STATUS='unknown')
  IF(IOSTAT /= 0) THEN
     CALL lnkerr('error in file handling')
  END IF
!
! Read in the Input
!
  time(1)=secnds(0.0)
  CALL Input_Prop
  time(2)=secnds(0.0)
  del_t(1) = time(2) - time(1)
!
! Calculate the Global FEDVR Matrices
!
  CALL Space_Prop
  time(3)=secnds(0.0)
  del_t(2) = time(3) - time(2)
!
! Compute the Regional Matrices from the FEDVR Global Matrices
!
  CALL Regional_Matrices
  time(4)=secnds(0.0)
  del_t(3) = time(4)-time(3)
!
! Get potential parameters.  The actual potential itself
! is calculated in the propagation step.
!
  CALL v_couple
  time(5)=secnds(0.0)
  del_t(4) = time(5)-time(4)
  CALL setup
  time(6)=secnds(0.0)
  del_t(5) = time(6)-time(5)
  WRITE(iout,1) del_t(1:5)
  WRITE(iout,4)
!
! ALLOCATE the Memory for the Matrices used to Store the Wavefunction
! and for the Potential
!
  IF ( type_calculation == 'real time' ) THEN
       ALLOCATE( v_tot(n3d), psi_z(n3d), v_scr_z(n3d) )
  ELSE
       ALLOCATE( v_tot(n3d), psi_d(n3d), v_scr_d(n3d) )
  END IF
  IF(algorithm == 'arnoldi') THEN
!                                                                        
!                                                                         
!    Read in any data for the arnoldi iterations.                       
!                                                                         
     CALL arnoldi_data
     time(7)=secnds(0.0)
     del_t(6) = time(7)-time(6)
!                                                                         
     bigv=max(maxvec,ntrial)
     bigs=max(5,n3d)

!    Allocate memory for the arnoldi arrays.              
!                                                                         
     IF ( type_calculation == 'real time' ) THEN
          ALLOCATE(vec_z(n3d,bigv),hvec_z(n3d,maxvec),                         &
                   b_z(maxvec,maxvec),bwrk_z(maxvec,maxvec),                   &
                   eig(maxvec),work_z(bigs*maxvec),rwork(5*maxvec),            &
                   u_z(n3d,maxvec),vscr_z(n3d))
          CALL Arnoldi_Propagation(psi_z,v_scr_z)
          DEALLOCATE(vec_z,hvec_z,b_z,bwrk_z,eig,work_z,rwork,u_z,vscr_z)
          DEALLOCATE(v_tot, psi_z, v_scr_z)
     ELSE
          ALLOCATE(vec_d(n3d,bigv),hvec_d(n3d,maxvec),                         &
                   b_d(maxvec,maxvec),bwrk_d(maxvec,maxvec),                   &
                   eig(maxvec),work_d(bigs*maxvec),rwork(5*maxvec),            &
                   u_d(n3d,maxvec),vscr_d(n3d))
          CALL Arnoldi_Propagation(psi_d,v_scr_d)
          DEALLOCATE(vec_d,hvec_d,b_d,bwrk_d,eig,work_d,rwork,u_d,vscr_d)
          DEALLOCATE(v_tot, psi_d, v_scr_d)
     END IF
     time(8)=secnds(0.0)
     del_t(7) = time(8)-time(7)
     WRITE(iout,2) del_t(6:7)
     WRITE(iout,4)
  ELSE IF(algorithm == 'split_operator') THEN
     IF ( type_calculation == 'real time' ) THEN
          ALLOCATE(exp_z(maxmem_reg,maxmem_reg))
          CALL SO_Sector_Time_Prop(exp_z)
          DEALLOCATE(exp_z)
          time(7) = secnds(0.0)
          del_t(6) = time(7)-time(6)
          ALLOCATE(exp_diag_z(n3d))
          CALL SO_Propagation(psi_z,v_scr_z)
          DEALLOCATE(exp_diag_z)
          DEALLOCATE(v_tot, psi_z, v_scr_z)
          time(8)=secnds(0.0)
          del_t(7) = time(8)-time(7)
     ELSE
          ALLOCATE(exp_d(maxmem_reg,maxmem_reg))
          CALL SO_Sector_Time_Prop(exp_d)
          DEALLOCATE(exp_d)
          time(7) = secnds(0.0)
          del_t(6) = time(7)-time(6)
          ALLOCATE(exp_diag_d(n3d))
          CALL SO_Propagation(psi_d,v_scr_d)
          DEALLOCATE(exp_diag_d)
          time(8)=secnds(0.0)
          del_t(7) = time(8)-time(7)
          DEALLOCATE(v_tot, psi_d, v_scr_d)
     END IF
     WRITE(iout,3) del_t(6:7)
     WRITE(iout,4)
  ELSE
    CALL lnkerr('algorithm error')
  END IF

1 FORMAT('***********************************************'                 &
            '*************************'                                    &
            /,10X,'Time to input basic data               = ',f15.8,/,     &
            /,10X,'Time to compute global FEDVR spatial ',                 &
            /,10x,'matrices                               = ',f15.8,/,     &
            /,10x,'Time to compute the regional spatial '                  &
            /,10x,'matrices                               = ',f15.8,/,     &
            /,10x,'Time to enter potential parameters     = ',f15.8,/,     &
            /,10x,'Time to compute the time grid and plot '                &
            /,10x,'files                                  = ',f15.8)
2 FORMAT('***********************************************'                 &
            '*************************'                                    &
            /,10x,'Time to input Arnoldi data             = ',f15.8,/,     &
            /,10x,'Time for the Arnoldi propagation       = ',f15.8)
3 FORMAT('***********************************************'                 &
            '*************************'                                    &
            /,10x,'Time to compute split operator sector '                 &
            /,10x,'time propagators                        = ',f15.8,/,    & 
            /,10x,'Time for the Split Operator Propagation = ',f15.8)
4 FORMAT('***********************************************'                 &
         '*************************')
  call chainx(0)
  stop
END PROGRAM prop_main
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
