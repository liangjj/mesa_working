!deck prop_main_it
!**begin prologue     prop_main_it
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
!**modules            dvrprop_global_it
!**end prologue       prop_main_it
  PROGRAM prop_main_it
  USE dvrprop_global_it
  IMPLICIT NONE
  REAL*4                                   :: secnds
  REAL*4, DIMENSION(10)                    :: del_t
  INTEGER                                  :: iostat
  CHARACTER(LEN=4096)                      :: ops
!
  CALL Drum
  CALL IOsys('read character options from rwf',-1,0,0,ops)
  OPEN (UNIT=iplot(1),FILE='misc',             &
        ACCESS='sequential',FORM='formatted',  &
        IOSTAT=IOSTAT,STATUS='unknown')
  IF(IOSTAT /= 0) THEN
     CALL lnkerr('error in file handling')
  END IF
  time(1)=secnds(0.0)
  CALL Input_Prop
  time(2)=secnds(0.0)
  del_t(1)=time(2)-time(1)
  CALL Space_Prop
  time(3)=secnds(0.0)
  del_t(2)=time(3)-time(2)
  CALL Regional_Matrices
  time(4)=secnds(0.0)
  del_t(3)=time(4)-time(3)
!
! Get potential parameters.  The actual potential itself
! is calculated in the propagation step.
!
  CALL v_couple
  time(5)=secnds(0.0)
  del_t(4)=time(5)-time(4)
  WRITE(iout,1) del_t(1:4)
  WRITE(iout,4)
  IF(algorithm == 'arnoldi') THEN
     CALL Arnoldi_Propagator
     time(5)=secnds(0.0)
     del_t(4)=time(5)-time(4)
     WRITE(iout,2) del_t(4)
     WRITE(iout,4)
  ELSE IF(algorithm == 'split_operator') THEN
     CALL Time_Prop
     time(5)=secnds(0.0)
     del_t(4)=time(5)-time(4)
     CALL Split_Operator
     time(6)=secnds(0.0)
     del_t(5)=time(6)-time(5)
     WRITE(iout,3) del_t(4:5)
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
            /,10x,'Time to enter potential parameters     = ',f15.8)
2 FORMAT('***********************************************'                 &
            '*************************'                                    &
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
END PROGRAM prop_main_it

