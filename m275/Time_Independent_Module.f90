!
  MODULE Time_Independent_Module
  USE Data
  USE Matrix_Print
  IMPLICIT NONE
  TYPE(Real_Vector)                        :: type_real_vector
  TYPE(Real_Matrix)                        :: type_real_matrix
!
!***********************************************************************
!***********************************************************************
!                          Explicit Interfaces
!***********************************************************************
!
                                       Contains
!____________________________________________________________________________________________!
!____________________________________________________________________________________________!
!_______________________________________TIME__INDEPENDENT____________________________________!
!______________________(FINDING EIGENVALUES AND EIGENVECTORS USING DSTEVX)___________________!
!____________________________________________________________________________________________!
!
!                                     DSTEVX:
! If Y_N='N', only eigenvalues will be computed; 
! If Y_N='V', eigenvalues and eigenvectors wil be computed.
!
! If Range = 'A', all eigenvalues will be computed.
! If Range = 'V', all eigenvalues in range (VL,VU] will be computed.
! If Range = 'I', all eigenvalues from ILth through IUth will be computed.
! 

  SUBROUTINE H_0 
  IMPLICIT NONE
  REAL(idp)                                :: DLAMCH
  INTEGER                                  :: i 
  INTEGER                                  :: len
!
  ALLOCATE(D(1:M_Size), E(1:M_Size), IBLOCK(1:M_Size), ISPLIT(1:M_Size),      &
           WORK_d(1:5*M_Size), W(1:M_Size), IWORK(1:5*M_Size), IFAIL(1:M_Size) )
  IF (Get_Eigenvectors == .true.) THEN
      Y_N = 'V'
      ALLOCATE(Z(M_Size,Number_of_Eigenvectors))
  END IF
  ABSTOL = 2 * DLAMCH('S')
  AD = - 1.d0 / ( 2.d0 * Step_Size * Step_Size )
  D(1:M_Size) = - 2.d0 * AD
  E(1:M_Size) = AD
  Call Potential
  Call cpu_time(t_i)
  Call DSTEBZ('I', 'E', M_Size, VL, VU, IL, IU, ABSTOL, D, E, M_FOUND, NSPLIT, W,    &
               IBLOCK, ISPLIT, WORK_d, IWORK, INFO)
  Call cpu_time(t_f)
  write(iout,*)
  write(iout,*) '               Elapsed Time for eigenvalues= ',t_f-t_i
  IF (INFO /= 0 ) THEN
      write(iout,*)
      Write(iout,*) '               Eigenvalue Routine Failed'
  ELSE
     write(iout,*)
     Write(iout,*) '               Found ',M_FOUND,' Eigenvalues'
     Write(iout,*)
     Call Print_Matrix(type_real_vector,W(1:M_Found),frmt='fr')
     Call IOsys('write real eigenvalues to tdse',M_FOUND,W,0,' ')           
  END IF
  IF (Get_Eigenvectors == .true.) THEN
      Call cpu_time(t_i)
      call DSTEIN(M_Size, D, E, Number_of_Eigenvectors , W, IBLOCK, ISPLIT, Z, M_Size, &
                                WORK_d, IWORK, IFAIL, INFO)
      IF (Prnt(1) == .true. ) THEN
         title = '                   Eigenvectors:'
         Call Print_Matrix(type_real_matrix,Z,M_Size, Number_of_Eigenvectors)
         Call IOsys('write real eigenvectors to tdse',Number_of_Eigenvectors*M_Size,Z,0,' ')           
      END IF
      Call cpu_time(t_f)
      write(iout,*)
      write(iout,*) '           Elapsed Time for eigenvectors= ',t_f-t_i      
  END IF
  DEALLOCATE( IBLOCK, ISPLIT, WORK_d, IWORK, IFAIL ) 
!********************************************************************************
  END SUBROUTINE H_0
!********************************************************************************
  SUBROUTINE Potential
  IMPLICIT NONE
  CHARACTER (LEN=80)               :: type
  CHARACTER (LEN=80)               :: chrkey
  LOGICAL                          :: dollar
  INTEGER                          :: i
  INTEGER                          :: len
  REAL(idp)                        :: XX
!
  ALLOCATE( V(1:M_Size), X(1:M_Size) )
!
  V(:) = 0.d0
!
  IF ( dollar('$potential',card,cpass,inp) )THEN
       type = chrkey (card,'type_potential','none',' ')
  ELSE
       call lnkerr('Potential Keywords Absent:Quit')
  END If
  XX = left_end
  DO i = 1, M_Size
     XX = XX + Step_Size
     X(i) = XX
  END DO
  Call pakstr(type,len)
  Call Potential_Mat( type(1:len), V, X, M_Size)
  D(:) = D(:) + V(:)
!***********************************************************************
!***********************************************************************
  END SUBROUTINE Potential
!***********************************************************************
!***********************************************************************
!deck Potential_Mat
!***begin prologue     Potential_Mat
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           lobatto functions
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Calculate the local potential matrix for a number of
!***                   fairly standard potentials.
!***                   
!***                   
!***                   
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       Potential_Mat
  SUBROUTINE Potential_Mat(type, pot, q,  n_pts)
  IMPLICIT NONE
  REAL(idp), DIMENSION(:)          :: pot
  REAL(idp), DIMENSION(:)          :: q                       
  INTEGER                          :: n_pts
  INTEGER                          :: i
  CHARACTER (LEN=80)               :: type
  REAL(idp)                        :: fpkey
  LOGICAL                          :: dollar
  LOGICAL                          :: logkey
  INTEGER                          :: intkey
  CHARACTER(LEN=80)                :: chrkey
  CHARACTER(LEN=8)                 :: itoc
!
!
!
  IF(type == 'none') then
!
     call none(pot,                                       &
               q,                                         &
               n_pts,                                     &
               prnt)
!
  ELSE IF(type == 'well') then
!
     depth=fpkey(card,'well_depth',zero,' ')
     call vwell(pot,                                      &
                depth,                                    &
                n_pts,                                    &
                prnt)
!
  ELSE IF(type == 'exponential') then
!
     amp(1)=fpkey(card,'amplitude',-1.d0,' ')
     expnt(1)=fpkey(card,'exponent',expnt,' ')
     call vexp(pot,                                       &
               q,                                         &
               amp,                                       &
               expnt,                                     &
               n_pts,                                     &
               prnt)
!
  ELSE IF(type == 'yukawa') then
!
     amp(1)=fpkey(card,'amplitude',-1.d0,' ')
     expnt(1)=fpkey(card,'exponent',expnt,' ')
     call vyukawa(pot,                                    &
                  q,                                      &
                  amp,                                    &
                  expnt,                                  &
                  n_pts,                                  &
                  prnt)
!
  ELSE IF(type == 'power_exponential') then
!
     amp(1)=fpkey(card,'amplitude',1.d0,' ')
     expnt(1)=fpkey(card,'exponent',expnt,' ')
     n_p=intkey(card,'power',0,' ')
     call v_pow_exp(pot,                                  &
                    q,                                    &
                    amp,                                  &
                    expnt,                                &
                    n_p,                                  &
                    n_pts,                                &
                    prnt)
!
  ELSE IF(type == 'sum_exponential') then
!
     call fparr(card,'amplitudes',amp,2,' ')
     call fparr(card,'exponents',expnt,2,' ')
     call vexp_sum(pot,                                   &
                   q,                                     &
                   amp,                                   &
                   expnt,                                 &
                   n_pts,                                 &
                   prnt)
!
  ELSE IF(type == 'coulomb') then
!
     charge=fpkey(card,'charge',-one,' ')
     call vcoul(pot,                                      &
                q,                                        &
                charge,                                   &
                n_pts,                                    &
                prnt)
!
  ELSE IF(type == 'eberlonium') then
!
     charge=fpkey(card,'charge',-one,' ')
     n_p=intkey(card,'power',0,' ')
     amp(1)=fpkey(card,'a',1.d0,' ')
     amp(2)=fpkey(card,'b',1.d0,' ')
     call v_eberlonium(pot,                               &
                       q,                                 &
                       charge,                            &
                       amp(1),                            &
                       amp(2),                            &
                       n_p,                               &
                       n_pts,                             &
                       prnt)
!
  ELSE IF(type == 'inverse_r4') then
!
     call vir4(pot,                                       &
               q,                                         &
               n_pts,                                     &
               prnt)
!
  ELSE IF(type == 'rounded_well') then
!
     nwell=intkey(card,'nwell',ten,' ')
     awell=fpkey(card,'a_well',14.d0,' ')
     call vrwell(pot,                                     &
                 q,                                       &
                 awell,                                   &
                 nwell,                                   &
                 n_pts,                                   &
                 prnt)
!
  ELSE IF(type == 'harmonic_oscillator') then
!
!        enter the mass and frequency in atomic units.
!         
     mass=fpkey(card,'mass',1.d0,' ')
     omega=fpkey(card,'omega',1.d0,' ')
     factor= mass * omega * omega * half
     hbar=one
     write(iout,1) mass, omega
     call vhmo(pot,                                       &
               q,                                         &
               factor,                                    &
               n_pts,                                     &
               prnt)
!
  ELSE IF(type == 'anharmonic_oscillator') then
!
     call vanhmo(pot,                                     &
                 q,                                       &
                 n_pts,                                   &
                 prnt)
!
  ELSE IF(type == 'expres') then
!
     call fparr(card,'amplitude',amp,2,' ')
     call fparr(card,'exponent',expnt,2,' ')
     shift=fpkey(card,'exponent_shift',zero,' ')
     call vres(pot,                                       &
               q,                                         &
               amp,                                       &
               expnt,                                     &
               shift,                                     &
               n_pts,                                     &
               prnt)
!
  ELSE IF(type == 'periodic') then
!
     n_scale=intkey(card,'n_i',10,' ')         
     e_c=fpkey(card,'e_c',.001d0,' ')         
     awell=n_scale/e_c
     call vperiod(pot,                                    &
                  q,                                      &
                  awell,                                  &
                  n_pts,                                  &
                  prnt)
!
  ELSE
     Call lnkerr('screwed up potential. quit')
  END IF
  IF(prnt(7)) then
     write(iout,2) type
     call Print_Matrix(type_real_vector,pot,title='potential')
  END IF
!
!
 1 Format(/,1x,'oscillator mass      = ',e15.8, &
          /,1x,'oscillator-frequency = ',e15.8)
 2 FORMAT(/,20x,'Potential Type = ',a32)
!
END SUBROUTINE Potential_Mat
!********************************************************************************
!********************************************************************************
  END MODULE Time_Independent_Module
!********************************************************************************
